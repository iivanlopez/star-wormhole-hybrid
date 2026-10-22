# -*- coding: utf-8 -*-
"""
Accion reducida del sector l = 0 de la teoria cuartica mas alla de Horndeski (fondo estatico general), a partir de la
lagrangiana cuadratica L2 de generar_bH.py (resultados/bH_simbolico.pkl).

Gauge unitario q = 0 (phi' != 0).  Quedan a (lapso), beta (g_tl), c (g_ll / 2) y d (delta r); el gauge residual
tau(t, l) actua solo sobre (a, beta): a -> a + tau_t, beta -> beta - e^{2 Phi} tau_l.  En la capa de masas L2 es lineal
en (a, beta) y, por la invariancia gauge, depende de ellos solo a traves de F = e^{2Phi} a' + beta_t:
    L2 = -F W + L_rest(c, d) + (derivadas totales),      W = gam c + al1 d' + al0 d,
donde W es la ligadura de momento (E_beta = i omega W en el gauge unitario) y L_rest = L2(a = beta = q = 0).  F es un
multiplicador: W = 0 da c = -(al1 d' + al0 d)/gam y la accion reducida es
    S_red[d] = int L_rest(c*(d), d) = int [ K d_t^2 - P d_l^2 - U d^2 ]      (forma canonica).
K, P, U se leen de la ecuacion de Euler-Lagrange de S_red:  E(d) = -2 [ K d_tt - (P d_l)_l + U d ].
gam, al1, al0 se tratan como funciones abstractas de l (simbolos gam_k, al1_k, al0_k para sus derivadas) y sus
expresiones en el fondo se exportan aparte.  Escribe resultados/monopolo_reduccion.pkl y monopolo_gen.py.
Comprobaciones (aritmetica exacta en chorros de la capa de masas, chorros.py):
  - el coeficiente de beta en E_beta y el de a en E_a, E_c se anulan;  E_c contiene a' y beta solo como F;
  - L_rest no contiene derivadas de c en su ecuacion de Euler-Lagrange (c auxiliar);
  - E(d) es de segundo orden (los coeficientes de las derivadas de orden 3..6 se anulan) y el de d' es 2 P'.
Uso: python monopolo_reduccion.py
"""
import os
import sys
import time
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros import Phi, rr, PP, Ee, Xs, D, chorro_en_capa, es_G       # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


d = pickle.load(open(os.path.join(AQUI, "resultados", "bH_simbolico.pkl"), "rb"))
fe, lin, L2, iw, amp = (d[k] for k in ("fondo_eq", "lineales", "L2", "iw", "amp"))
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
pj = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in ("a", "be", "c", "d", "q") for o in ORD}

# ------------------------------------------------------------------ 1. ligadura de momento en el gauge unitario
uni = {amp[("q", k)]: 0 for k in range(6)}
Ebe = sp.expand(lin["be"].subs(uni))
k_beta = Ebe.coeff(amp[("be", 0)])
resto = sp.expand(Ebe - k_beta * amp[("be", 0)])
gam = sp.expand(resto.coeff(amp[("c", 0)]) / iw)
al1 = sp.expand(resto.coeff(amp[("d", 1)]) / iw)
al0 = sp.expand(resto.coeff(amp[("d", 0)]) / iw)
assert sp.expand(resto - iw * (gam * amp[("c", 0)] + al1 * amp[("d", 1)] + al0 * amp[("d", 0)])) == 0
assert not any(x.has(iw) for x in (gam, al1, al0))
Ec = sp.expand(lin["c"].subs(uni))
Ea = sp.expand(lin["a"].subs(uni))
cA1 = Ec.coeff(amp[("a", 1)])
cB0 = sp.expand(Ec.coeff(amp[("be", 0)]) / iw)
log(f"E_beta (unitario) = k_beta beta + i w (gam c + al1 d' + al0 d): gam {len(sp.Add.make_args(gam))}, "
    f"al1 {len(sp.Add.make_args(al1))}, al0 {len(sp.Add.make_args(al0))} terminos")
for sem in range(3):
    v = chorro_en_capa(fe, 400 + sem, orden=3)
    assert v.ev(k_beta) == 0
    assert v.ev(Ea.coeff(amp[("a", 0)])) == 0 and v.ev(Ec.coeff(amp[("a", 0)])) == 0
    assert v.ev(cA1 - Ee * cB0) == 0                     # E_c contiene e^{2 Phi} a' + beta_t = F
log("en la capa de masas: k_beta = 0, coeficiente de a en E_a y E_c = 0, E_c depende de (a', beta) solo via F "
    "(3 chorros exactos)")

# ------------------------------------------------------------------ 2. L_rest(c, d) y su ecuacion para c
L_rest = sp.expand(L2.xreplace({pj[(v, o)]: 0 for v in ("a", "be", "q") for o in ORD}))
log(f"L_rest(c, d): {len(sp.Add.make_args(L_rest))} terminos")
NMAX = 7
dJ = {(i, j): sp.Symbol(f"d{i}{j}", real=True) for i in range(NMAX + 1) for j in range(NMAX + 1)}
cJ = {(i, j): sp.Symbol(f"c{i}{j}", real=True) for i in range(NMAX + 1) for j in range(NMAX + 1)}
AJ = [[dJ[(i, j)] for j in range(NMAX + 1)] for i in range(NMAX + 1)]       # filas: orden en t
CJ = [[cJ[(i, j)] for j in range(NMAX + 1)] for i in range(NMAX + 1)]


def Dl(e, extra=()):
    """d/dl: fondo + chorros de d y c (indice l) + funciones abstractas 'extra'."""
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
v = chorro_en_capa(fe, 410, orden=4)
nulos = [s for s in derivs_c if v.ev(Ec_rest.coeff(sp.Symbol(s, real=True))) == 0]
assert set(nulos) == set(derivs_c), (derivs_c, nulos)
k_c = sp.expand(Ec_rest.coeff(cJ[(0, 0)]))
log(f"E-L de L_rest para c: las derivadas de c {derivs_c} tienen coeficiente nulo en la capa de masas (c auxiliar)")

# ------------------------------------------------------------------ 3. sustitucion c = -(al1 d' + al0 d)/gam y E-L
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
log(f"L_red(d): {len(sp.Add.make_args(Lred))} terminos")
Ered = euler_lagrange(Lred, dJ, EXTRA)
log(f"E-L de L_red: {len(sp.Add.make_args(Ered))} terminos")
coef = {}
for (i, j), s in dJ.items():
    c_ = Ered.coeff(s)
    if c_ != 0:
        coef[(i, j)] = sp.expand(c_)
log("chorros de d presentes en E(d): " + ", ".join(f"d{i}{j}" for (i, j) in sorted(coef)))
K = sp.expand(-coef[(2, 0)] / 2)
P = sp.expand(coef[(0, 2)] / 2)
U = sp.expand(-coef[(0, 0)] / 2)
log(f"K: {len(sp.Add.make_args(K))}, P: {len(sp.Add.make_args(P))}, U: {len(sp.Add.make_args(U))} terminos")
with open(os.path.join(AQUI, "resultados", "monopolo_reduccion.pkl"), "wb") as fh:
    pickle.dump(dict(gam=gam, al1=al1, al0=al0, k_c=k_c, K=K, P=P, U=U, coef=coef, L_rest=L_rest,
                     extra=EXTRA), fh)
log("guardado resultados/monopolo_reduccion.pkl")
