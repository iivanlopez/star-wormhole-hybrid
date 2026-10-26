# -*- coding: utf-8 -*-
"""
Fase 2g, parte II: coeficientes de la accion reducida de dos campos (fluido + monopolo) para el calculo numerico.

En el gauge q = x = 0 (reduccion_fluido.py):  L2 = -(N^2 a' + beta_t) W + kap_a a^2 + a M - kap_b beta^2 + L_rest(c, d),
y la ecuacion de c de L_rest no contiene derivadas de c (se comprueba aqui, fuera de la capa de masas):
    E_c^rest = k_c c + e0 d + e1 d' + e_tt d_tt,
de modo que  L_rest(c, d) = L_rest(0, d) + k_c c^2/2 + c (e0 d + e1 d' + e_tt d_tt)  (salvo derivadas totales), con
L_rest(0, N h) = k_d N^2 h_t^2 + 2 N^3 Gcal (h' - sig h)^2 (familia_c0.py).  Este script extrae k_c, e0, e1, e_tt con
el fluido (el fluido aporta -p~ N r^2 a k_c y 2 p~ N r a e0), las soluciones de las ligaduras del fondo con fluido para
G2 y G2_X, y los compila para chorros numericos.  Escribe resultados/coeficientes_fluido.pkl.
Uso: python coeficientes_fluido.py
"""
import os
import sys
import time
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
import chorros as CH                                       # noqa: E402
from chorros import Phi, rr, PP, Ee, eP, Xs              # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


def S(n):
    return sp.Symbol(n, real=True)


fl = pickle.load(open(os.path.join(AQUI, "resultados", "fluido_simbolico.pkl"), "rb"))
rd = pickle.load(open(os.path.join(AQUI, "resultados", "reduccion_fluido.pkl"), "rb"))
Pf = fl["Pf"]
NPF = len(Pf)
Lrest, J = rd["Lrest"], rd["J"]
Xb = {Xs: -PP[1] ** 2 / 2}
N = sp.sqrt(Ee)


def Dbg(e, amps=()):
    out = CH.D(e, amps)
    for k in range(NPF - 1):
        if e.has(Pf[k]):
            out += sp.diff(e, Pf[k]) * (-Pf[k + 1] * Phi[1] / Ee)
    return sp.expand(out)


NMAX = max(i for (i, j) in J["c"])
CAMPOS = ("a", "be", "c", "d")
FILAS = [tuple(J[f][(i, j)] for j in range(NMAX + 1)) for f in CAMPOS for i in range(NMAX + 1)]


def Dl(e):
    return Dbg(e, tuple(FILAS))


def Dt(e):
    out = 0
    for f in CAMPOS:
        for (i, j), s in J[f].items():
            if i < NMAX and e.has(s):
                out += sp.diff(e, s) * J[f][(i + 1, j)]
    return sp.expand(out)


def euler_lagrange(Lg, f):
    tot = 0
    for (i, j), s in J[f].items():
        if i + j > 3 or not Lg.has(s):
            continue
        Q = sp.expand(sp.diff(Lg, s)) * (-1) ** (i + j)
        for _ in range(i):
            Q = Dt(Q)
        for _ in range(j):
            Q = Dl(Q)
        tot += Q
    return sp.expand(tot)


Ec = euler_lagrange(Lrest, "c")
c, d = J["c"], J["d"]
k_c = sp.expand(Ec.coeff(c[(0, 0)]))
e0 = sp.expand(Ec.coeff(d[(0, 0)]))
e1 = sp.expand(Ec.coeff(d[(0, 1)]))
e_tt = sp.expand(Ec.coeff(d[(2, 0)]))
resto = sp.expand(Ec - k_c * c[(0, 0)] - e0 * d[(0, 0)] - e1 * d[(0, 1)] - e_tt * d[(2, 0)])
assert resto == 0, resto.free_symbols
log("E_c^rest = k_c c + e0 d + e1 d' + e_tt d_tt exactamente (fuera de la capa de masas, con fluido)")
# comparacion con la Fase 2f (sin fluido) y aportes del fluido
mr = pickle.load(open(os.path.join(F2F, "resultados", "monopolo_reduccion.pkl"), "rb"))
kc_vac = sp.expand(mr["k_c"].xreplace(Xb))
dk = sp.simplify(k_c - kc_vac + Pf[0] * N * rr[0] ** 2)
assert dk == 0, dk
log("k_c = k_c(Fase 2f) - p~ N r^2")
for nombre, e in (("k_c", k_c), ("e0", e0), ("e1", e1), ("e_tt", e_tt)):
    log(f"  {nombre}: {len(sp.Add.make_args(e))} terminos, simbolos {sorted(s.name for s in e.free_symbols)}")

# ligaduras del fondo con fluido: G2_00 de E_A = 0, G2_01 de E_B = 0
fe = fl["fondo_eq"]
g200, g201 = S("G2_00"), S("G2_01")
EA = sp.expand(fe["A"].xreplace(Xb))
EB = sp.expand(fe["B"].xreplace(Xb))
sA = sp.solve(EA, g200)
sB = sp.solve(EB, g201)
assert len(sA) == 1 and len(sB) == 1
G2_00 = sp.expand(sA[0])
G2_01 = sp.expand(sB[0].subs(g200, G2_00))
log("ligaduras con fluido: G2_00 de E_A = 0, G2_01 de E_B = 0")
pickle.dump(dict(k_c=k_c, e0=e0, e1=e1, e_tt=e_tt, G2_00=G2_00, G2_01=G2_01, gam=rd["gam"], al1=rd["al1"],
                 al0=rd["al0"], Pf=Pf), open(os.path.join(AQUI, "resultados", "coeficientes_fluido.pkl"), "wb"))
log("guardado resultados/coeficientes_fluido.pkl")
