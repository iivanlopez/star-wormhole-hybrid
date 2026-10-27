# -*- coding: utf-8 -*-
"""
Fase 2m: estructura del monopolo con vector temporal (caso E: A1 = 0, G3 = 0; paridad TC).

Variables invariantes bajo las reparametrizaciones temporales (T):
    at = al - U a        (delta A0 menos el arrastre del lapso),     wt = w + U beta / N^2.
En ellas, L2 debe depender de (a, beta) solo a traves de F = N^2 a' + beta_t.  Este script comprueba en chorros exactos de
la capa de masas (caso E) que
    L2 = kap F^2 - F W + L_rest(c, d, at, wt)      (salvo derivadas totales),
con kap = U^2 r^2 (G2_F + 2 F G2_FF) / (2 N^5)   y W sin derivadas temporales,
leyendo E_beta = -i w (2 kap F - W).  Luego:
  - la ligadura Pi = 2 kap F - W = 0 da F = W/(2 kap) y L_eff = -W^2/(4 kap) + L_rest: la energia del campo
    electrico perturbado entra con signo POSITIVO, +W^2/(4 kap);
  - lee W = gam c + al1 d' + al0 d + (terminos en at, wt) y la evalua en la prueba d = N, c = at = wt = 0;
  - clasifica at, wt (auxiliares o dinamicos) por los coeficientes de sus derivadas temporales en L_rest.
Escribe resultados/reduccion_E.pkl.
Uso: python reduccion_E.py
"""
import os
import sys
import json
import time
import pickle
from math import comb
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import chorros_proca as CP                              # noqa: E402
from chorros_proca import NN, rr, UU, VV, gsym          # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


d = pickle.load(open(os.path.join(AQUI, "resultados", "proca_simbolico.pkl"), "rb"))
L2, pj, fe = d["L2"], d["pj"], d["fondo_eq"]
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
LE = sp.expand(L2.xreplace({VV[k]: 0 for k in range(CP.NMAXB)}))
LE = sp.expand(LE.xreplace({s: 0 for s in LE.free_symbols if s.name.startswith("G3")}))
NV = ("a", "be", "c", "d", "at", "wt")
q = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in NV for o in ORD}


def leib(fbg, var, o):
    nt, nl = o
    out, fk = 0, fbg
    for k in range(nl + 1):
        out += comb(nl, k) * fk * q[(var, (nt, nl - k))]
        fk = CP._D0(fk)
    return out


rep = {}
for o in ORD:
    for v in ("a", "be", "c", "d"):
        rep[pj[(v, o)]] = q[(v, o)]
    rep[pj[("al", o)]] = q[("at", o)] + leib(UU[0], "a", o)
    rep[pj[("w", o)]] = q[("wt", o)] - leib(UU[0] / NN[0] ** 2, "be", o)
L = sp.expand(LE.xreplace(rep))
log(f"L2 (caso E) en las variables invariantes: {len(sp.Add.make_args(L))} terminos")

# ecuaciones de Euler-Lagrange para modos e^{i w t}
iw = sp.Symbol("iw")
NMAX = 5
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
                Q = CP.D(sp.sympify(Q), AMPS)
            tot += Q
        out[v] = sp.expand(tot)
    return out


E = EL(L)
log("ecuaciones lineales: " + ", ".join(f"{v}: {len(sp.Add.make_args(e))}" for v, e in E.items()))
kap = UU[0] ** 2 * rr[0] ** 2 * (gsym("G2", 0, 1) * NN[0] ** 2 + gsym("G2", 0, 2) * UU[1] ** 2) / (2 * NN[0] ** 7)
Fm = NN[0] ** 2 * amp[("a", 1)] + iw * amp[("be", 0)]
Ebe = E["be"]
# parte en (a, be) de E_beta frente a -i w 2 kap F
ab = [amp[(v, k)] for v in ("a", "be") for k in range(NMAX + 1)]
parte_ab = sp.expand(sum(Ebe.coeff(s) * s for s in ab))
resto = sp.expand(Ebe - parte_ab)
Wexp = sp.expand(resto / iw)                            # E_be = -iw (2 kap F - W)  =>  resto = iw W
CHK = {}
dif_ab = sp.expand(parte_ab + iw * 2 * kap * Fm)
otros = [amp[(v, k)] for v in ("c", "d", "at", "wt") for k in range(NMAX + 1)]
res = []
for sem in range(3):
    v = CP.chorro_en_capa(fe, 700 + sem, orden=4, caso="E")
    x = v.ev(dif_ab)
    y = v.ev(Wexp)
    # W no debe depender de iw
    res.append((float(sp.Abs(sum(abs(c) for c in sp.Poly(x, *(ab + [iw])).coeffs())) if x != 0 else 0),
                sp.Poly(y, iw).degree() if y.has(iw) else 0))
    # identidad de Bianchi: E_a = -D(N^2 Pi) con Pi = -E_be / iw
    Pi = sp.expand(-Ebe / iw)
    bia = v.ev(sp.expand(E["a"] + CP.D(NN[0] ** 2 * Pi, AMPS)))
    res[-1] = res[-1] + (bia == 0,)
    log(f"chorro E {sem}: |parte (a, beta) de E_beta + 2 i w kap F| = {res[-1][0]:.1e}, grado de W en iw = "
        f"{res[-1][1]}, E_a = -(N^2 Pi)': {res[-1][2]}")
CHK["kap_F2_menos_FW"] = res
assert all(r[0] == 0 and r[1] <= 1 and r[2] for r in res)
W1 = sp.expand(Wexp.coeff(iw))                           # parte de W con una derivada temporal
assert W1.free_symbols & set(otros) <= {amp[("wt", k)] for k in range(NMAX + 1)}, W1.free_symbols
log(f"parte de W con derivada temporal (solo wt, impar TC): {sp.factor_terms(W1)}")
CHK["W_iw"] = str(sp.factor_terms(W1))

# W en la capa de masas: forma y valor en la prueba d = N
Wc = {s: sp.factor_terms(Wexp.coeff(s)) for s in otros if Wexp.coeff(s) != 0}
log("W = " + " + ".join(f"({str(c)[:120]}) {s}" for s, c in Wc.items()))
prueba = {amp[(v, k)]: 0 for v in ("c", "at", "wt") for k in range(NMAX + 1)}
e = NN[0]
for k in range(NMAX + 1):
    prueba[amp[("d", k)]] = e
    e = CP._D0(e)
WN = sp.factor_terms(sp.expand(Wexp.xreplace(prueba)))
log(f"W en la prueba d = N, c = at = wt = 0 (fuera de la capa de masas): {WN}")
# parte cinetica de L_rest (a = be = 0)
Lr = sp.expand(L.xreplace({q[(v, o)]: 0 for v in ("a", "be") for o in ORD}))
cin = {}
for t1 in sp.Add.make_args(Lr):
    sym = [s for s in t1.free_symbols if s in q.values() and s.name.split("_")[1][0] == "1"]
    if len(sym) == 2 or any(t1.coeff(s, 2) != 0 for s in sym):
        key = tuple(sorted(str(s) for s in sym))
        cin[key] = cin.get(key, 0) + t1
log("terminos cineticos de L_rest: " + "; ".join(f"{k}: {sp.factor_terms(v)}" for k, v in cin.items()))
pickle.dump(dict(L=L, E=E, W=Wexp, kap=kap, q=q, amp=amp, iw=iw, Lr=Lr, WN=WN),
            open(os.path.join(AQUI, "resultados", "reduccion_E.pkl"), "wb"))
CHK["W_prueba_N"] = str(WN)
CHK["cineticos_L_rest"] = {str(k): str(sp.factor_terms(v)) for k, v in cin.items()}
json.dump(CHK, open(os.path.join(AQUI, "resultados", "reduccion_E.json"), "w"), indent=1, default=str)
log("guardado resultados/reduccion_E.pkl")
