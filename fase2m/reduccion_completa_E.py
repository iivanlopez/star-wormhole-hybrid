# -*- coding: utf-8 -*-
"""
Fase 2m (continuacion): reduccion completa del monopolo con vector temporal (caso E: A1 = 0, G3 = 0).

Partimos de L_eff = -W^2/(4 kap) + L_rest (F ya eliminado, energia_E.py) en el gauge radial c = 0.
Quedan d (par), at = dA0 - U a (par, auxiliar) y wt = dA1 + U beta/N^2 (impar, auxiliar).  Afirmacion:
    L_eff = -(2 Gcal/N) d_t^2 + ct wt^2 + wt S_t + mt at^2 + at J(d) + L_s(d)      (salvo derivadas totales)
con
    S = -(4 r N^2/U) (Gcal d/N)',   Gcal = G4 - 2X G4_X,
    ct = c_ww - (e_ww)'/2,  mt = c_aa - (c_aa1)'/2,  J(d) estatico (sin velocidades),
es decir: at no se acopla a wt ni a velocidades (los terminos at d_tt + at_t d_t y at wt_t, at_t wt, at' wt_t,
at_t wt' se cancelan), y wt solo ve S_t.  Eliminando los auxiliares:
    S_red[d] = int dt dl [ T(d_t) - V(d) ],
    T(d_t) = -(2 Gcal/N) d_t^2 + k2 ((Gcal d_t/N)')^2,     k2 = -4 r^2 N^4/(U^2 ct),
    V(d)   = -L_s(d) + J(d)^2/(4 mt).
Comprobaciones (aritmetica exacta en chorros de la capa de masas, caso E):
  1. las ecuaciones de Euler-Lagrange de L_eff - L_afirmada se anulan (d, at, wt);
  2. (en esfera_E.py, numerico) el modo gauge residual d = r' Z(t) (c = 0 deja Z = Z(t)) esta en el nucleo de T y
     de V:  K r' = 0 y V r' = -E_s(r') + J^dag(J(r')/(2 mt)) = 0.
Consecuencia (exacta): donde ct > 0 (y Gcal, N > 0) los dos terminos de T son negativos; como ct(0) > 0 en toda
garganta (longitudinal_E.py), toda perturbacion con soporte cerca de la garganta tiene energia cinetica negativa:
el monopolo tiene un FANTASMA.
Uso: python reduccion_completa_E.py
"""
import os
import sys
import json
import time
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import chorros_proca as CP                              # noqa: E402
from chorros_proca import NN, rr, UU, VV, gsym          # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


d0 = pickle.load(open(os.path.join(AQUI, "resultados", "proca_simbolico.pkl"), "rb"))
en = pickle.load(open(os.path.join(AQUI, "resultados", "energia_E.pkl"), "rb"))
re = pickle.load(open(os.path.join(AQUI, "resultados", "reduccion_E.pkl"), "rb"))
fe = d0["fondo_eq"]
q, Lr, kap, Wj = re["q"], re["Lr"], en["kap"], en["Wj"]
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
cv = {VV[k]: 0 for k in range(CP.NMAXB)}
c0 = {q[(v, o)]: 0 for v in ("c", "a", "be") for o in ORD}
Leff = sp.expand((-Wj ** 2 / (4 * kap) + Lr).xreplace(cv).xreplace(c0))
G3z = {s: 0 for s in Leff.free_symbols if s.name.startswith("G3")}
Leff = sp.expand(Leff.xreplace(G3z))
D = CP.D
dj = {o: q[("d", o)] for o in ORD}
aj = {o: q[("at", o)] for o in ORD}
wj = {o: q[("wt", o)] for o in ORD}
G4, G4X = gsym("G4", 0, 0), gsym("G4", 1, 0)
Gc = G4 - UU[0] ** 2 / NN[0] ** 2 * G4X
N, r, U = NN[0], rr[0], UU[0]
CHK = {}

# ------------------------------------------------------------------ piezas de la forma afirmada
coef = lambda e, *s: sp.expand(e.coeff(s[0]) if len(s) == 1 else e.coeff(s[0]).coeff(s[1]))   # noqa: E731
c_ww = sp.expand(Leff.coeff(wj[(0, 0)], 2).xreplace({x: 0 for x in q.values()}))
e_ww = coef(Leff, wj[(0, 0)], wj[(0, 1)])
ct = sp.expand(c_ww - D(e_ww) / 2)
c_aa = sp.expand(Leff.coeff(aj[(0, 0)], 2).xreplace({x: 0 for x in q.values()}))
c_aa1 = coef(Leff, aj[(0, 0)], aj[(0, 1)])
mt = sp.expand(c_aa - D(c_aa1) / 2)
k_dd = sp.expand(Leff.coeff(dj[(1, 0)], 2).xreplace({x: 0 for x in q.values()}))
assert sp.simplify(k_dd + 2 * Gc / N) == 0
log("coeficiente de d_t^2 = -2 Gcal/N  (exacto)")
# parte estatica en (d, at): L con at y d estaticos, wt = 0
est = {x: 0 for x in q.values() if x.name.split("_")[1][0] != "0"}
est.update({x: 0 for x in wj.values()})
Lst = sp.expand(Leff.xreplace(est))
Jd = sp.expand(sum(Lst.coeff(aj[(0, 0)]).xreplace({aj[(0, 1)]: 0}).coeff(dj[o]) * dj[o] for o in ORD))
J1 = sp.expand(sum(Lst.coeff(aj[(0, 1)]).xreplace({aj[(0, 0)]: 0}).coeff(dj[o]) * dj[o] for o in ORD))
d3 = sp.Symbol("d_03", real=True)
Jtot = sp.expand(Jd - D(J1, [[dj[(0, 0)], dj[(0, 1)], dj[(0, 2)], d3]]))   # at' J1 -> -at J1'
Ls = sp.expand(Lst.xreplace({aj[o]: 0 for o in ORD}))
# S
S_af = -4 * r * N ** 2 / U * D(Gc * dj[(0, 0)] / N, [[dj[(0, 0)], dj[(0, 1)], dj[(0, 2)], d3]])
S_af = sp.expand(S_af)
assert not S_af.has(d3)
St = sp.expand(S_af.xreplace({dj[(0, 0)]: dj[(1, 0)], dj[(0, 1)]: dj[(1, 1)]}))      # d -> d_t
Laf = (k_dd * dj[(1, 0)] ** 2 + ct * wj[(0, 0)] ** 2 + wj[(0, 0)] * St + mt * aj[(0, 0)] ** 2
       + aj[(0, 0)] * Jtot.xreplace({d3: 0}) + Ls)
assert not Jtot.has(d3)
resto = sp.expand(Leff - Laf)
log(f"L_eff - L_afirmada: {len(sp.Add.make_args(resto))} terminos (deben ser una derivada total)")

# ------------------------------------------------------------------ 1. EL del resto en chorros de la capa de masas
iw = sp.Symbol("iw")
NMAX = 6
NV = ("d", "at", "wt")
amp = {(v, k): sp.Symbol(f"{v}{k}", real=True) for v in NV for k in range(NMAX + 1)}
modo = {q[(v, o)]: iw ** o[0] * amp[(v, o[1])] for v in NV for o in ORD}
AMPS = [[amp[(v, k)] for k in range(NMAX + 1)] for v in NV]


def EL(Lx):
    out = {}
    for v in NV:
        tot = 0
        for o in ORD:
            Q = sp.diff(Lx, q[(v, o)]).xreplace(modo) * (-1) ** (o[0] + o[1]) * iw ** o[0]
            for _ in range(o[1]):
                Q = D(sp.sympify(Q), AMPS)
            tot += Q
        out[v] = sp.expand(tot)
    return out


Eres = EL(resto)
mx = 0
for sem in range(2):
    v = CP.chorro_en_capa(fe, 1200 + sem, orden=5, caso="E")
    for nom, e in Eres.items():
        x = v.ev(e)
        if x != 0:
            mx = max(mx, max(abs(float(c)) for c in sp.Poly(x, *(list(amp.values()) + [iw])).coeffs()))
CHK["resto_derivada_total"] = mx
log(f"1. EL(L_eff - L_afirmada) en 2 chorros de la capa de masas: max {mx:.1e}")
assert mx == 0

# ------------------------------------------------------------------ 2. el modo gauge residual d = r' en el nucleo
# T = kT(d_t), V(d) = -Ls + J^2/(4 mt).  Operadores (EL estaticos, d -> funcion de l):
k2 = -4 * r ** 2 * N ** 4 / (U ** 2 * ct)
ds = [sp.Symbol(f"x{k}", real=True) for k in range(9)]                     # x = d(l) y derivadas


def Dx(e):
    out = D(e)
    for k in range(8):
        out += sp.diff(e, ds[k]) * ds[k + 1]
    return sp.expand(out)


def el_l(dens, orden=2):
    """Euler-Lagrange en l de una densidad en (x, x', x'')."""
    out = sp.diff(dens, ds[0]) - Dx(sp.diff(dens, ds[1]))
    if orden >= 2:
        out += Dx(Dx(sp.diff(dens, ds[2])))
    return sp.expand(out)


a_x = {dj[(0, 0)]: ds[0], dj[(0, 1)]: ds[1], dj[(0, 2)]: ds[2]}
Qx = Gc * ds[0] / N
Tdens = k_dd * ds[0] ** 2 + k2 * Dx(Qx) ** 2
Vdens = -Ls.xreplace(a_x) + Jtot.xreplace(a_x) ** 2 / (4 * mt)
pickle.dump(dict(ct=ct, mt=mt, Jtot=Jtot, Ls=Ls, k2=k2, k_dd=k_dd, Gc=Gc, Tdens=Tdens, Vdens=Vdens, ds=ds),
            open(os.path.join(AQUI, "resultados", "reduccion_completa_E.pkl"), "wb"))
log("guardado resultados/reduccion_completa_E.pkl")
# 2. La comprobacion del modo gauge residual (K r' = V r' = 0) se hace en esfera_E.py sobre las gargantas
#    numericas: en chorros exactos los cocientes por ct y mt derivados dos veces no terminan en un tiempo razonable.
#    (Logicamente es redundante: L_afirmada = L_eff salvo derivadas totales (1) y L2 es invariante gauge
#    (verificar_proca.py); comprueba el algebra de la eliminacion de los auxiliares.)

# ------------------------------------------------------------------ resumen simbolico
CHK["S"] = "S = -(4 r N^2/U) (Gcal d/N)'"
CHK["k2"] = "k2 = -4 r^2 N^4/(U^2 ct)"
CHK["T"] = "T = int [ -(2 Gcal/N) d_t^2 + k2 ((Gcal d_t/N)')^2 ]"
json.dump(CHK, open(os.path.join(AQUI, "resultados", "reduccion_completa_E.json"), "w"), indent=1, default=str)
log("comprobaciones superadas")
