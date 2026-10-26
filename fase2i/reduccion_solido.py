# -*- coding: utf-8 -*-
"""
Fase 2i, parte I: estructura del sector l = 0 con el solido (GLPV cuartica + F(X_r, X_t, f)) en el gauge unitario.

A. Gauge unitario (q = 0); el gauge residual tau solo actua sobre (a, beta), porque el solido estatico no cambia con el
   tiempo (y no se mueve).  Se comprueba con chorros exactos en la capa de masas con solido que
       L2 = -(N^2 a' + beta_t)(W + P_m) + L_rest(c, d, y)        (salvo derivadas totales),
       W = gam c + al1 d' + al0 d  (Fase 2f),   P_m = (r^2/N)(rho~ + p~_r) xi,   xi = -y/f0'  (desplazamiento),
   sin terminos a^2, beta^2 ni a beta (a diferencia del fluido, el solido no depende del lapso).  Las ecuaciones de a
   y beta dan la misma ligadura,  W = -(r^2/N)(rho~ + p~_r) xi,  y el monopolo tiene dos campos, (xi, d).
B. Ecuacion de c en L_rest: E_c = k_c c + e0 d + e1 d' + e_tt d_tt + g0 y + g1 y' (c auxiliar).
C. Forma lagrangiana (comovil) de la energia estatica del solido, con xi como variable: el elemento de material
   que en el fondo esta en l pasa a l + xi, y
       L_m = r^2 N exp(eps_r + 2 eps_t + nu) G(eps_r, eps_t)     (G(e_r, e_t) = F(X_r e^{-2 e_r}, X_t e^{-2 e_t}, f)),
       eps_r = ln[(1 + xi') sqrt(1 + 2 c(l + xi))],  eps_t = (1/2) ln[C(l + xi)/r(l)^2],  nu = Phi(l + xi) - Phi(l);
   desarrollada a segundo orden, sin integrar por partes (la integral es sobre el dominio material, fijo: las
   superficies libres no necesitan terminos de superficie).  Sus ecuaciones de Euler-Lagrange en (xi, c, d) coinciden
   con las de la forma euleriana (y = -f0' xi) en chorros con el solido en equilibrio (E_f = 0).
   Coeficientes: G0 = -rho~, G_r = rho~ + p~_r, G_t = 2 (rho~ + p~_t), G_rr = 4 X_r^2 F_rr + 4 X_r F_r,
   G_rt = 4 X_r X_t F_rt, G_tt = 4 X_t^2 F_tt + 4 X_t F_t  (modulos en deformaciones logaritmicas).
Escribe resultados/reduccion_solido.pkl y .json.
Uso: python reduccion_solido.py
"""
import os
import sys
import json
import time
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
sys.path.insert(0, AQUI)
from chorros import Phi, rr, PP, Ee, eP, Xs              # noqa: E402
from simbolos_solido import (ff, Fs, Xrb, Xtb, D, chorro_capa_solido, chorro_equilibrio,  # noqa: E402
                             rho, p_r, p_t)

T0 = time.time()
RES = {}


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


def S(n):
    return sp.Symbol(n, real=True)


d_ = pickle.load(open(os.path.join(AQUI, "resultados", "solido_simbolico.pkl"), "rb"))
L2g, L2m, fe = d_["L2_g"], d_["L2_m"], d_["fondo_eq"]
mr = pickle.load(open(os.path.join(F2F, "resultados", "monopolo_reduccion.pkl"), "rb"))
Xb = {Xs: -PP[1] ** 2 / 2}
N = sp.sqrt(Ee)

# ------------------------------------------------------------------ A. lagrangiana en el gauge q = 0, chorros (t, l)
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
cero = {S(f"q_{o[0]}{o[1]}"): 0 for o in ORD}
L = sp.expand((L2g + L2m.subs(eP, N)).xreplace(cero).xreplace(Xb))
NMAX = 6
CAMPOS = ("a", "be", "c", "d", "y")
J = {f: {(i, j): S(f"{f}{i}{j}") for i in range(NMAX + 1) for j in range(NMAX + 1)} for f in CAMPOS}
L = L.xreplace({S(f"{f}_{o[0]}{o[1]}"): J[f][o] for f in CAMPOS for o in ORD})
FILAS = [tuple(J[f][(i, j)] for j in range(NMAX + 1)) for f in CAMPOS for i in range(NMAX + 1)]


def Dl(e):
    return D(e, tuple(FILAS))


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


EL = {f: euler_lagrange(L, f) for f in CAMPOS}
log("A. ecuaciones de Euler-Lagrange (gauge q = 0): " + ", ".join(f"{f}: {len(sp.Add.make_args(e))}" for f, e in EL.items()))
gam = sp.expand(mr["gam"].xreplace(Xb))
al1 = sp.expand(mr["al1"].xreplace(Xb))
al0 = sp.expand(mr["al0"].xreplace(Xb))
c, d, a, b, y = J["c"], J["d"], J["a"], J["be"], J["y"]
rpr = sp.expand(rho + p_r)                              # rho~ + p~_r = -2 ff1^2 Fs_100
W = gam * c[(0, 0)] + al1 * d[(0, 1)] + al0 * d[(0, 0)]
Pm = -rr[0] ** 2 / N * rpr * y[(0, 0)] / ff[1]         # = (r^2/N)(rho~ + p~_r) xi,  xi = -y/ff1
F_ = Ee * a[(0, 1)] + b[(1, 0)]
Lrest = sp.expand(L.xreplace({s: 0 for f in ("a", "be") for s in J[f].values()}))
Lhyp = sp.expand(-F_ * (W + Pm) + Lrest)
ELh = {f: euler_lagrange(Lhyp, f) for f in CAMPOS}
for sem in range(3):
    v = chorro_capa_solido(fe, 1500 + sem, orden=6)
    for f in CAMPOS:
        val = v.ev(sp.expand(EL[f] - ELh[f]))
        assert val == 0, (f, val)
RES["forma"] = ("L2 = -(N^2 a' + beta_t)(W + P_m) + L_rest(c, d, y), P_m = (r^2/N)(rho~ + p~_r) xi: Euler-Lagrange "
                "identicas en 3 chorros exactos en la capa de masas con solido")
log("   forma comprobada: L2 = -(N^2 a' + beta_t)(W + P_m) + L_rest(c, d, y),  P_m = (r^2/N)(rho~+p~_r) xi  "
    "(E-L exactas, 3 chorros)")

# ------------------------------------------------------------------ B. ecuacion de c en L_rest
Ec = euler_lagrange(Lrest, "c")
k_c = sp.expand(Ec.coeff(c[(0, 0)]))
e0 = sp.expand(Ec.coeff(d[(0, 0)]))
e1 = sp.expand(Ec.coeff(d[(0, 1)]))
e_tt = sp.expand(Ec.coeff(d[(2, 0)]))
g0 = sp.expand(Ec.coeff(y[(0, 0)]))
g1 = sp.expand(Ec.coeff(y[(0, 1)]))
resto = sp.expand(Ec - k_c * c[(0, 0)] - e0 * d[(0, 0)] - e1 * d[(0, 1)] - e_tt * d[(2, 0)] - g0 * y[(0, 0)]
                  - g1 * y[(0, 1)])
assert resto == 0, resto.free_symbols
kc_vac = sp.expand(mr["k_c"].xreplace(Xb))
dk = sp.simplify(k_c - kc_vac)
log(f"B. E_c^rest = k_c c + e0 d + e1 d' + e_tt d_tt + g0 y + g1 y' exactamente;  k_c - k_c(Fase 2f) = {sp.factor(dk)}")
log(f"   g0 = {sp.factor(g0)},  g1 = {sp.factor(g1)}")
RES["c"] = dict(dk_c=str(sp.factor(dk)), g0=str(sp.factor(g0)), g1=str(sp.factor(g1)))

# ------------------------------------------------------------------ C. forma lagrangiana de la energia estatica
sJ = {o: S(f"s{o[0]}{o[1]}") for o in [(0, j) for j in range(NMAX + 1)]}
cJ = {o: c[o] for o in [(0, j) for j in range(NMAX + 1)]}
dJ = {o: d[o] for o in [(0, j) for j in range(NMAX + 1)]}
s0, s1 = sJ[(0, 0)], sJ[(0, 1)]
c0, c1 = cJ[(0, 0)], cJ[(0, 1)]
d0_, d1 = dJ[(0, 0)], dJ[(0, 1)]
r0, r1, r2 = rr[0], rr[1], rr[2]
er1 = s1 + c0
er2 = -s1 ** 2 / 2 + c1 * s0 - c0 ** 2
et1 = (r1 * s0 + d0_) / r0
et2 = ((r1 ** 2 + r0 * r2) * s0 ** 2 + 2 * r1 * d0_ * s0 + 2 * r0 * d1 * s0) / (2 * r0 ** 2) - (r1 * s0 + d0_) ** 2 / r0 ** 2
nu1 = Phi[1] * s0
nu2 = Phi[2] * s0 ** 2 / 2
z1 = er1 + 2 * et1 + nu1
z2 = er2 + 2 * et2 + nu2
G0 = Fs(0, 0, 0)
Gr = -2 * Xrb * Fs(1, 0, 0)
Gt = -2 * Xtb * Fs(0, 1, 0)
Grr = 4 * Xrb ** 2 * Fs(2, 0, 0) + 4 * Xrb * Fs(1, 0, 0)
Grt = 4 * Xrb * Xtb * Fs(1, 1, 0)
Gtt = 4 * Xtb ** 2 * Fs(0, 2, 0) + 4 * Xtb * Fs(0, 1, 0)
Llag = sp.expand(r0 ** 2 * N * (G0 * (z2 + z1 ** 2 / 2) + Gr * er2 + Gt * et2 + (Gr * er1 + Gt * et1) * z1
                                + (Grr * er1 ** 2 + 2 * Grt * er1 * et1 + Gtt * et1 ** 2) / 2))
# forma euleriana estatica del solido (a = beta = 0, sin derivadas temporales), y = -ff1 xi
L2m_l = L2m.subs(eP, N).xreplace(cero)
est = {S(f"{f}_{o[0]}{o[1]}"): 0 for f in ("a", "be") for o in ORD}
est.update({S(f"{f}_{o[0]}{o[1]}"): 0 for f in ("c", "d", "y") for o in ORD if o[0] > 0})
L2m_s = sp.expand(L2m_l.xreplace(est))
yx = {0: -ff[1] * s0, 1: -ff[2] * s0 - ff[1] * s1, 2: -ff[3] * s0 - 2 * ff[2] * s1 - ff[1] * sJ[(0, 2)]}
LE = sp.expand(L2m_s.xreplace({S("y_00"): yx[0], S("y_01"): yx[1], S("y_02"): yx[2],
                               S("c_00"): c0, S("c_01"): c1, S("c_02"): cJ[(0, 2)],
                               S("d_00"): d0_, S("d_01"): d1, S("d_02"): dJ[(0, 2)]}))
FS = [tuple(sJ[(0, j)] for j in range(NMAX + 1)), tuple(cJ[(0, j)] for j in range(NMAX + 1)),
      tuple(dJ[(0, j)] for j in range(NMAX + 1))]


def el_estatica(Lg, campo):
    tot = 0
    for j in range(4):
        s = campo[(0, j)]
        if not Lg.has(s):
            continue
        Q = sp.expand(sp.diff(Lg, s)) * (-1) ** j
        for _ in range(j):
            Q = D(Q, FS)
        tot += Q
    return sp.expand(tot)


fe_f = sp.expand(fe["f"])
for nombre, campo in (("xi", sJ), ("c", cJ), ("d", dJ)):
    dif = sp.expand(el_estatica(LE, campo) - el_estatica(Llag, campo))
    for sem in range(3):
        v = chorro_equilibrio(fe_f, 1700 + sem, orden=6)
        val = v.ev(dif)
        assert val == 0, (nombre, val)
RES["forma_lagrangiana"] = ("E-L de la forma lagrangiana en (xi, c, d) = E-L de la euleriana con y = -f0' xi, en 3 chorros "
                            "con el solido en equilibrio (metrica arbitraria)")
log("C. forma lagrangiana (comovil) de la energia estatica: E-L en (xi, c, d) iguales a las de la euleriana "
    "(3 chorros con E_f = 0, metrica arbitraria)")

with open(os.path.join(AQUI, "resultados", "reduccion_solido.pkl"), "wb") as fh:
    pickle.dump(dict(Lrest=Lrest, J=J, gam=gam, al1=al1, al0=al0, Pm=Pm, k_c=k_c, e0=e0, e1=e1, e_tt=e_tt, g0=g0,
                     g1=g1, Llag=Llag, sJ=sJ), fh)
json.dump(RES, open(os.path.join(AQUI, "resultados", "reduccion_solido.json"), "w"), indent=1)
log("guardado resultados/reduccion_solido.pkl y .json")
