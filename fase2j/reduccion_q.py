# -*- coding: utf-8 -*-
"""
Fase 2j: el monopolo de la teoria quintica de Horndeski (G2, G3, G4, G5; F4 = F5 = 0) en el gauge unitario, con el mismo
metodo que fase2f/monopolo_reduccion.py y fase2g/familia_c0.py, y la prueba d = N del teorema de la Fase 2g.

Comprobaciones (aritmetica exacta en chorros de la capa de masas, chorros_q.py):
  1. E_beta (unitario) = k_beta beta + i w W,  W = gam c + al1 d' + al0 d;  k_beta = 0 en la capa de masas;
     relacion al0 = -Phi' al1 (identica o en la capa de masas);
  2. los coeficientes de a en E_a y E_c se anulan y E_c depende de (a', beta) solo via F = e^{2Phi} a' + beta_t:
     L2 = -F W + L_rest(c, d);
  3. c es auxiliar: la ecuacion de c de L_rest no contiene derivadas de c;
  4. la prueba d = N f, c = 0 (cumple W = 0 si al0 = -Phi' al1):  -L_rest(0, N f)|_estatico frente a
     -2 N^3 Gtil (f' - sig f)^2 con Gtil = -N al1/(4 r), sig = r'/r - Phi', y la parte cinetica k_d.
Escribe resultados/reduccion_q.pkl y .json.
Uso: python reduccion_q.py
"""
import os
import sys
import json
import time
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_q import Phi, rr, PP, Ee, Xs, D, chorro_en_capa   # noqa: E402

T0 = time.time()
RES = {}


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


def S(n):
    return sp.Symbol(n, real=True)


d = pickle.load(open(os.path.join(AQUI, "resultados", "q_simbolico.pkl"), "rb"))
fe, lin, L2, iw, amp = (d[k] for k in ("fondo_eq", "lineales", "L2", "iw", "amp"))
Xb = {Xs: -PP[1] ** 2 / 2}
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
pj = {(v, o): S(f"{v}_{o[0]}{o[1]}") for v in ("a", "be", "c", "d", "q") for o in ORD}
horn = {s: 0 for e in list(fe.values()) + [L2] for s in e.free_symbols if s.name.startswith(("F4_", "F5_"))}
fe = {k: sp.expand(v.xreplace(horn)) for k, v in fe.items()}
L2 = sp.expand(L2.xreplace(horn))
lin = {k: sp.expand(v.xreplace(horn)) for k, v in lin.items()}
N = sp.sqrt(Ee)

# ------------------------------------------------------------------ 1. ligadura de momento
uni = {amp[("q", k)]: 0 for k in range(6)}
Ebe = sp.expand(lin["be"].subs(uni))
k_beta = Ebe.coeff(amp[("be", 0)])
resto = sp.expand(Ebe - k_beta * amp[("be", 0)])
gam = sp.expand(resto.coeff(amp[("c", 0)]) / iw)
al1 = sp.expand(resto.coeff(amp[("d", 1)]) / iw)
al0 = sp.expand(resto.coeff(amp[("d", 0)]) / iw)
assert sp.expand(resto - iw * (gam * amp[("c", 0)] + al1 * amp[("d", 1)] + al0 * amp[("d", 0)])) == 0
assert not any(x.has(iw) for x in (gam, al1, al0))
rel = sp.simplify((al0 + Phi[1] * al1).xreplace(Xb))
log(f"1. E_beta = k_beta beta + i w (gam c + al1 d' + al0 d): gam {len(sp.Add.make_args(gam))}, al1 "
    f"{len(sp.Add.make_args(al1))}, al0 {len(sp.Add.make_args(al0))} terminos;  al0 + Phi' al1 = {rel}")
RES["al0_mas_Phi_al1"] = str(rel)
log(f"   al1 = {sp.factor(al1.xreplace(Xb))}")
RES["al1"] = str(sp.factor(al1.xreplace(Xb)))
Ec = sp.expand(lin["c"].subs(uni))
Ea = sp.expand(lin["a"].subs(uni))
cA1 = Ec.coeff(amp[("a", 1)])
cB0 = sp.expand(Ec.coeff(amp[("be", 0)]) / iw)
for sem in range(3):
    v = chorro_en_capa(fe, 400 + sem, orden=3)
    assert v.ev(k_beta) == 0
    assert v.ev(Ea.coeff(amp[("a", 0)])) == 0 and v.ev(Ec.coeff(amp[("a", 0)])) == 0
    assert v.ev(cA1 - Ee * cB0) == 0
    if rel != 0:
        log(f"   al0 + Phi' al1 en la capa de masas (chorro {sem}): {v.ev(rel.subs(Xs, -PP[1] ** 2 / 2))}")
log("2. en la capa de masas: k_beta = 0, a solo via F = e^{2Phi} a' + beta_t  ->  L2 = -F W + L_rest(c, d)")

# ------------------------------------------------------------------ 3. c auxiliar
L_rest = sp.expand(L2.xreplace({pj[(v, o)]: 0 for v in ("a", "be", "q") for o in ORD}))
NMAX = 7
dJ = {(i, j): S(f"d{i}{j}") for i in range(NMAX + 1) for j in range(NMAX + 1)}
cJ = {(i, j): S(f"c{i}{j}") for i in range(NMAX + 1) for j in range(NMAX + 1)}
AJ = [[dJ[(i, j)] for j in range(NMAX + 1)] for i in range(NMAX + 1)]
CJ = [[cJ[(i, j)] for j in range(NMAX + 1)] for i in range(NMAX + 1)]


def Dl(e, extra=()):
    return D(e, tuple(AJ) + tuple(CJ) + tuple(extra))


def Dt(e):
    out = 0
    for (i, j), s in dJ.items():
        if i < NMAX and e.has(s):
            out += sp.diff(e, s) * dJ[(i + 1, j)]
    for (i, j), s in cJ.items():
        if i < NMAX and e.has(s):
            out += sp.diff(e, s) * cJ[(i + 1, j)]
    return sp.expand(out)


Lr = L_rest.xreplace({pj[("d", o)]: dJ[o] for o in ORD}).xreplace({pj[("c", o)]: cJ[o] for o in ORD})


def euler_lagrange(L, campo, extra=()):
    tot = 0
    for (i, j) in [(i, j) for i in range(3) for j in range(4) if i + j <= 3]:
        s = campo[(i, j)]
        if not L.has(s):
            continue
        Q = sp.expand(sp.diff(L, s)) * (-1) ** (i + j)
        for _ in range(i):
            Q = Dt(Q)
        for _ in range(j):
            Q = Dl(Q, extra)
        tot += Q
    return sp.expand(tot)


Ec_rest = euler_lagrange(Lr, cJ)
derivs_c = sorted({str(s) for s in Ec_rest.free_symbols if s in set(cJ.values()) and s != cJ[(0, 0)]})
nulos = []
for sem in range(2):
    v = chorro_en_capa(fe, 410 + sem, orden=4)
    nulos.append({s for s in derivs_c if v.ev(Ec_rest.coeff(S(s))) == 0})
aux = all(n == set(derivs_c) for n in nulos)
RES["c_auxiliar"] = dict(derivadas=derivs_c, auxiliar=aux)
log(f"3. E-L de L_rest para c: derivadas de c presentes {derivs_c}; coeficientes nulos en la capa de masas: {aux}")
k_c = sp.expand(Ec_rest.coeff(cJ[(0, 0)]))

# ------------------------------------------------------------------ 4. la prueba d = N f, c = 0
Lstat = sp.expand(Lr.xreplace({cJ[(i, j)]: 0 for (i, j) in cJ}))
cin = sp.expand(Lstat.coeff(dJ[(1, 0)], 2))
Ld0 = sp.expand(Lstat.xreplace({dJ[(i, j)]: 0 for (i, j) in dJ if i > 0}))
f = sp.symbols("f0:5", real=True)


def Df(e):
    out = D(e)
    for k in range(4):
        out += sp.diff(e, f[k]) * f[k + 1]
    return sp.expand(out)


dd = [N * f[0]]
for k in range(3):
    dd.append(Df(dd[-1]))
menosL = sp.expand(-Ld0.xreplace({dJ[(0, k)]: dd[k] for k in range(3)}).xreplace(Xb))
Gtil = sp.expand((-N * al1 / (4 * rr[0])).xreplace(Xb))
sig = rr[1] / rr[0] - Phi[1]
dif_id = sp.simplify(menosL + 2 * N ** 3 * Gtil * (f[1] - sig * f[0]) ** 2)
RES["identidad_dN_exacta"] = str(dif_id)
log(f"4. -L_rest(0, N f) + 2 N^3 Gtil (f' - sig f)^2 = {sp.factor(dif_id) if dif_id != 0 else 0}   (Gtil = -N al1/4r)")
log(f"   Gtil = {sp.factor(Gtil)}")
log(f"   cinetica k_d = {sp.factor(cin.xreplace(Xb))}")
RES["Gtil"] = str(sp.factor(Gtil))
RES["k_d"] = str(sp.factor(cin.xreplace(Xb)))
if dif_id != 0:
    # en la capa de masas
    for sem in range(3):
        v = chorro_en_capa(fe, 420 + sem, orden=4)
        for k in range(5):
            v[f[k]] = v.r()
        log(f"   en la capa de masas (chorro {sem}): {v.ev(dif_id)}")

# ------------------------------------------------------------------ 5. accion reducida S_red[d] = int [K d_t^2 - P d'^2 - U d^2]
NF = 8
g_ = sp.symbols(f"gam0:{NF}", real=True)
a1_ = sp.symbols(f"al1_0:{NF}", real=True)
a0_ = sp.symbols(f"al0_0:{NF}", real=True)
EXTRA = (g_, a1_, a0_)
c_expr = {}
base = sp.expand(-(a1_[0] * dJ[(0, 1)] + a0_[0] * dJ[(0, 0)]) / g_[0])
for i in range(3):
    e = base
    for _ in range(i):
        e = Dt(e)
    for j in range(3 - i):
        c_expr[(i, j)] = e
        e = Dl(e, EXTRA)
Lred = sp.expand(Lr.xreplace({cJ[o]: c_expr[o] for o in ORD}))
Ered = euler_lagrange(Lred, dJ, EXTRA)
coef = {}
for (i, j), s in dJ.items():
    c_ = Ered.coeff(s)
    if c_ != 0:
        coef[(i, j)] = sp.expand(c_)
log("5. chorros de d presentes en la ecuacion reducida: " + ", ".join(f"d{i}{j}" for (i, j) in sorted(coef)))
Kred = sp.expand(-coef[(2, 0)] / 2)
Pred = sp.expand(coef[(0, 2)] / 2)
Ured = sp.expand(-coef[(0, 0)] / 2)
# segundo orden: los coeficientes de d''' ... se anulan en la capa de masas (con las expresiones de gam, al1, al0)
alto = [k for k in coef if k[0] == 0 and k[1] >= 3]
if alto:
    v = chorro_en_capa(fe, 430, orden=6)
    sub = {}
    # derivadas de gam, al1, al0 con el operador D en el chorro
    for simb, expr in ((g_, gam), (a1_, al1), (a0_, al0)):
        e = expr
        for k_ in range(5):
            sub[simb[k_]] = v.ev(e)
            e = D(e)
    for k in alto:
        log(f"   coeficiente de d{k[0]}{k[1]} en la capa de masas: {v.ev(coef[k].xreplace(sub))}")
log(f"   K: {len(sp.Add.make_args(Kred))}, P: {len(sp.Add.make_args(Pred))}, U: {len(sp.Add.make_args(Ured))} terminos")

with open(os.path.join(AQUI, "resultados", "reduccion_q.pkl"), "wb") as fh:
    pickle.dump(dict(gam=gam, al1=al1, al0=al0, k_c=k_c, L_rest=L_rest, Gtil=Gtil, menosL=menosL, dif_id=dif_id,
                     K=Kred, P=Pred, U=Ured, coef=coef, extra=EXTRA), fh)
json.dump(RES, open(os.path.join(AQUI, "resultados", "reduccion_q.json"), "w"), indent=1)
log("guardado resultados/reduccion_q.pkl y .json")
