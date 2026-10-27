# -*- coding: utf-8 -*-
"""
Fase 2m: gargantas simetricas con vector temporal (caso E: A1 = 0) en la familia
    G4 = 1 + b X,      G2 = F + h(X),   h = m2 X + l2 X^2 + l3 X^3,      G3 = 0,
y la energia de la prueba d = N (energia_E.py):
    E[N] = int [ -2 N^3 G4 sig^2 + (2N/g) (2 G4_X N (U/N)' - G2_F U')^2 + J_N^2/(4 m_at) ] dl.
Fondo: E_A, E_C, E_u resueltas para (N'', r'', U''); E_B es la ligadura (se vigila).  En la garganta (l = 0)
r' = N' = U' = 0 y E_B da h(X0) = -2 G4(X0)/r0^2; con r0 = 1 y N0 = 1 queda X0 como parametro de tiro para que el
modo de Yukawa creciente de U desaparezca (U -> 0 lejos).  Condicion de garganta (exacta, reduccion_E):
    r''(0) = -U0^2 (G2_X r0^2 + 2 G4_X) / (4 r0 N0^2 Gcal) > 0,   Gcal = G4 - 2 X G4_X.
Uso: python garganta_E.py   (escribe resultados/garganta_E.json y perfiles_E.npz)
"""
import os
import sys
import json
import pickle
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import chorros_proca as CP                              # noqa: E402
from chorros_proca import NN, rr, UU, VV, gsym          # noqa: E402

d0 = pickle.load(open(os.path.join(AQUI, "resultados", "proca_simbolico.pkl"), "rb"))
en = pickle.load(open(os.path.join(AQUI, "resultados", "energia_E.pkl"), "rb"))
fe = d0["fondo_eq"]
b, m2, l2, l3 = sp.symbols("b m2 l2 l3", real=True)
Xs = sp.Symbol("X")
h = m2 * Xs + l2 * Xs ** 2 + l3 * Xs ** 3
Xb = UU[0] ** 2 / (2 * NN[0] ** 2)
Fb = UU[1] ** 2 / (2 * NN[0] ** 2)
teo = {gsym("G2", 0, 0): Fb + h.subs(Xs, Xb), gsym("G2", 0, 1): 1, gsym("G2", 1, 0): sp.diff(h, Xs).subs(Xs, Xb),
       gsym("G2", 2, 0): sp.diff(h, Xs, 2).subs(Xs, Xb), gsym("G2", 3, 0): sp.diff(h, Xs, 3).subs(Xs, Xb),
       gsym("G4", 0, 0): 1 + b * Xb, gsym("G4", 1, 0): b}
cero = {VV[k]: 0 for k in range(CP.NMAXB)}


def T(e):
    e = e.xreplace(cero).xreplace(teo)
    return e.xreplace({s: 0 for s in e.free_symbols if CP.es_G(s)})


ARGS = (NN[0], NN[1], rr[0], rr[1], UU[0], UU[1], b, m2, l2, l3)
inc = [NN[2], rr[2], UU[2]]
eqs = [sp.expand(T(fe[n])) for n in ("A", "C", "u")]
Mx = sp.Matrix([[e.coeff(u) for u in inc] for e in eqs])
bx = sp.Matrix([sp.expand(e - sum(e.coeff(u) * u for u in inc)) for e in eqs])
f_M = sp.lambdify(ARGS, Mx, "numpy")
f_b = sp.lambdify(ARGS, bx, "numpy")
f_B = sp.lambdify(ARGS, T(fe["B"]), "numpy")
# densidades de la prueba d = N (necesitan segundas derivadas: se pasan N2, r2, U2 y terceras via diferencias)
ARGS2 = ARGS[:6] + (NN[2], rr[2], UU[2], NN[3], rr[3], UU[3]) + ARGS[6:]
Gc = T(en["Gc"])
G4 = T(gsym("G4", 0, 0))
sig = en["sig"]
g = T(gsym("G2", 0, 1) + gsym("G2", 0, 2) * UU[1] ** 2 / NN[0] ** 2)
WN = T(en["WN"])
kap = T(en["kap"])
vol = -2 * NN[0] ** 3 * G4 * sig ** 2
ele = sp.simplify(WN ** 2 / (4 * kap))
mt = T(en["mt"])
JN = T(en["JN"])
cww = T(sp.expand(d0["L2"]).coeff(d0["pj"][("w", (0, 0))], 2).xreplace(
    {s: 0 for s in d0["L2"].free_symbols if s.name[:2] in ("a_", "be", "c_", "d_", "al", "w_")}))
f_dens = sp.lambdify(ARGS2, [vol, ele, JN, mt, Gc, G4, cww], "numpy")


def rhs(l, y, p):
    N, N1, r, r1, U, U1 = y
    a = (N, N1, r, r1, U, U1) + tuple(p)
    x = np.linalg.solve(np.array(f_M(*a), dtype=float), -np.array(f_b(*a), dtype=float).ravel())
    return [N1, x[0], r1, x[1], U1, x[2]]


def h_num(X, p):
    _, m2_, l2_, l3_ = p
    return m2_ * X + l2_ * X ** 2 + l3_ * X ** 3


def integra(X0, p, L=40.0, denso=False):
    b_ = p[0]
    hX0 = h_num(X0, p)
    G40 = 1 + b_ * X0
    if hX0 >= 0 or G40 <= 0:
        return None
    r0 = np.sqrt(-2 * G40 / hX0)
    y0 = [1.0, 0.0, r0, 0.0, np.sqrt(2 * X0), 0.0]

    def explota(l, y, p):
        return abs(y[4]) - 3 * np.sqrt(2 * X0) - 1e-9
    explota.terminal = True

    def cruza(l, y, p):
        return y[4]
    cruza.terminal = True
    sol = solve_ivp(rhs, (0, L), y0, args=(p,), rtol=1e-11, atol=1e-13, events=(explota, cruza), dense_output=denso,
                    max_step=0.05)
    return sol


def tiro(X0, p):
    """+1 si U explota hacia arriba, -1 si cruza el cero (o explota hacia abajo)."""
    s = integra(X0, p)
    if s is None:
        return np.nan
    if s.t_events[1].size:
        return -1.0
    return 1.0 if s.y[4, -1] > 0 else -1.0


if __name__ == "__main__":
    import itertools
    out = []
    # barrido de teorias: b, m2, l2, l3
    rejilla = list(itertools.product([0.0, -0.5, -1.0, -2.0], [1.0], [-1.0, -3.0, -10.0], [0.0, 1.0, 5.0]))
    for p in rejilla:
        Xs_ = np.linspace(0.02, 3.0, 80)
        sgn = [tiro(X, p) for X in Xs_]
        for i in range(len(Xs_) - 1):
            if np.isfinite(sgn[i]) and np.isfinite(sgn[i + 1]) and sgn[i] != sgn[i + 1]:
                out.append(dict(p=p, X0a=float(Xs_[i]), X0b=float(Xs_[i + 1])))
        print(p, [o for o in out if o["p"] == p], flush=True)
    json.dump(out, open(os.path.join(AQUI, "resultados", "garganta_E_barrido.json"), "w"), indent=1)
