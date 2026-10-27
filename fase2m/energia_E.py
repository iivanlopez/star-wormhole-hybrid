# -*- coding: utf-8 -*-
"""
Fase 2m: energia estatica del monopolo con vector temporal (caso E) tras eliminar F (reduccion_E.py).

  L_eff = -W^2/(4 kap) + L_rest(c, d, at, wt)     (en chorros; W con iw wt -> wt_t).
1. Comprueba que L2 - [kap F^2 - F W + L_rest] es una derivada total (sus ecuaciones de Euler-Lagrange se anulan en
   chorros de la capa de masas).
2. Sector estatico par (wt = 0 por paridad TC): fija el gauge radial con c = 0.  Lee la forma de at en L_eff estatico
   (los terminos at'^2 se cancelan entre -W^2/4kap y L_rest) y elimina at por su ecuacion.
3. Evalua la energia de la familia d = N f:  E[f] = -L_eff(c = 0, d = N f, at = at*(f)) y la compara con
   -2 N^3 Gcal (f' - sig f)^2, Gcal = G4 - 2 X G4_X (X = U^2/2N^2): da la parte nueva del vector.
Uso: python energia_E.py
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


r = pickle.load(open(os.path.join(AQUI, "resultados", "reduccion_E.pkl"), "rb"))
d0 = pickle.load(open(os.path.join(AQUI, "resultados", "proca_simbolico.pkl"), "rb"))
fe = d0["fondo_eq"]
L, Lr, Wamp, kap, q, amp, iw = (r[k] for k in ("L", "Lr", "W", "kap", "q", "amp", "iw"))
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
NV = ("a", "be", "c", "d", "at", "wt")
cv = {VV[k]: 0 for k in range(CP.NMAXB)}
# W en chorros
Wj = 0
for t1 in sp.Add.make_args(sp.expand(Wamp.xreplace(cv))):
    fac = 1
    rest = t1
    if t1.has(iw):
        rest = t1 / iw
        assert rest.has(amp[("wt", 0)])
        Wj += rest.xreplace({amp[("wt", 0)]: q[("wt", (1, 0))]})
        continue
    hit = [(v, k) for v in NV for k in range(3) if t1.has(amp[(v, k)])]
    assert len(hit) == 1, t1
    v, k = hit[0]
    Wj += t1.xreplace({amp[(v, k)]: q[(v, (0, k))]})
Wj = sp.expand(Wj)
Fj = NN[0] ** 2 * q[("a", (0, 1))] + q[("be", (1, 0))]
resto = sp.expand(L - kap * Fj ** 2 + Fj * Wj - Lr)
iwS = sp.Symbol("iw")
NMAX = 5
ampS = {(v, k): sp.Symbol(f"{v}{k}", real=True) for v in NV for k in range(NMAX + 1)}
modo = {q[(v, o)]: iwS ** o[0] * ampS[(v, o[1])] for v in NV for o in ORD}
AMPS = [[ampS[(v, k)] for k in range(NMAX + 1)] for v in NV]


def EL(Lx, vars_=NV):
    out = {}
    for v in vars_:
        tot = 0
        for o in ORD:
            Q = sp.diff(Lx, q[(v, o)]).xreplace(modo) * (-1) ** (o[0] + o[1]) * iwS ** o[0]
            for _ in range(o[1]):
                Q = CP.D(sp.sympify(Q), AMPS)
            tot += Q
        out[v] = sp.expand(tot)
    return out


Eres = EL(resto)
CHK = {}
mx = 0
for sem in range(2):
    v = CP.chorro_en_capa(fe, 900 + sem, orden=4, caso="E")
    for nom, e in Eres.items():
        x = v.ev(e)
        mx = max(mx, max([abs(float(c)) for c in sp.Poly(x, *(list(ampS.values()) + [iwS])).coeffs()] or [0]))
CHK["resto_derivada_total"] = mx
log(f"L2 - (kap F^2 - F W + L_rest): ecuaciones de Euler-Lagrange en 2 chorros de la capa de masas: max {mx:.1e}")
assert mx == 0

# ------------------------------------------------------------------ sector estatico
est = {q[(v, o)]: 0 for v in NV for o in ORD if o[0] > 0}
est.update({q[(v, o)]: 0 for v in ("a", "be", "wt", "c") for o in ORD})
Leff = sp.expand((-Wj ** 2 / (4 * kap) + Lr).xreplace(est).xreplace(cv))
at0, at1, at2 = q[("at", (0, 0))], q[("at", (0, 1))], q[("at", (0, 2))]
d_0, d_1, d_2 = q[("d", (0, 0))], q[("d", (0, 1))], q[("d", (0, 2))]
pol = sp.Poly(Leff, at0, at1, at2, d_0, d_1, d_2)
forma = {m: sp.factor_terms(sp.simplify(c)) for m, c in zip(pol.monoms(), pol.coeffs())}
for m, c in forma.items():
    log(f"  L_eff estatico (c = 0): {m}  ->  {str(c)[:300]}")
CHK["L_eff_estatico"] = {str(m): str(c) for m, c in forma.items()}
pickle.dump(dict(Leff=Leff, Wj=Wj, forma=forma), open(os.path.join(AQUI, "resultados", "energia_E.pkl"), "wb"))
json.dump(CHK, open(os.path.join(AQUI, "resultados", "energia_E.json"), "w"), indent=1)

# ------------------------------------------------------------------ eliminacion de at y energia de la prueba d = N
d_3 = sp.Symbol('d_03', real=True)
c_aa = Leff.coeff(at0, 2).xreplace({at1: 0, at2: 0})
c_aa1 = sp.expand(Leff.coeff(at0).coeff(at1))
assert sp.simplify(Leff.coeff(at1, 2)) == 0 and sp.simplify(Leff.coeff(at2)) == 0
mt = sp.expand(c_aa - CP.D(c_aa1) / 2)                    # at at' = (at^2)'/2
Jat = sp.expand(Leff.coeff(at0).xreplace({at1: 0}))       # lineal en d
J1 = sp.expand(Leff.coeff(at1).xreplace({at0: 0}))        # at' J1  ->  -at J1'
Jtot = sp.expand(Jat - CP.D(J1, [[d_0, d_1, d_2, d_3]]))
Ld = sp.expand(Leff.xreplace({at0: 0, at1: 0, at2: 0}))
fN = {d_0: NN[0], d_1: NN[1], d_2: NN[2], d_3: NN[3]}
JN = sp.factor_terms(sp.simplify(Jtot.xreplace(fN)))
mN = sp.factor_terms(sp.simplify(mt))
log(f"m_at (coeficiente de at^2) = {str(mN)[:400]}")
log(f"J en d = N: {JN}")
Gc = gsym("G4", 0, 0) - UU[0] ** 2 / NN[0] ** 2 * gsym("G4", 1, 0)
sig = rr[1] / rr[0] - NN[1] / NN[0]
WN = sp.expand(Wj.xreplace(est).xreplace(cv).xreplace({at0: 0, at1: 0}).xreplace(fN))
dens = sp.expand(-Ld.xreplace(fN))                         # energia sin la correccion de at
delta = sp.factor_terms(sp.simplify(dens + 2 * NN[0] ** 3 * Gc * sig ** 2 - WN ** 2 / (4 * kap)))
log(f"-L_d(N) = -2 N^3 Gcal sig^2 + W_N^2/(4 kap) + Delta,  Delta = {delta}")
CHK["m_at"] = str(mN)
CHK["J_N"] = str(JN)
CHK["Delta"] = str(delta)
CHK["W_N"] = str(sp.factor_terms(WN))
pickle.dump(dict(Leff=Leff, Wj=Wj, mt=mt, Jtot=Jtot, Ld=Ld, JN=JN, WN=WN, delta=delta, kap=kap, Gc=Gc, sig=sig),
            open(os.path.join(AQUI, "resultados", "energia_E.pkl"), "wb"))
json.dump(CHK, open(os.path.join(AQUI, "resultados", "energia_E.json"), "w"), indent=1)
log("guardado resultados/energia_E.pkl")
