# -*- coding: utf-8 -*-
"""
Fase 2i, parte I: solido elastico en el sector l = 0 de la teoria cuartica mas alla de Horndeski.

En simetria esferica un solido se describe con su coordenada radial material f(t, l) (coordenadas del cuerpo
phi^I = f n^I) y una lagrangiana F(X_r, X_t, f), con
    X_r = g^{mu nu} d_mu f d_nu f      (estiramiento radial al cuadrado),
    X_t = f^2/C                        (estiramiento tangencial al cuadrado, C = r^2),
y la dependencia explicita en f describe un material no homogeneo (capas):
    S = (1/16 pi) int sqrt(-g) [ L_GLPV + F(X_r, X_t, f) ],
    rho~ = -F,   p~_r = F - 2 X_r F_r,   p~_t = F - X_t F_t           (p~ = 16 pi p en unidades G = c = 1).
Un fluido es el caso F = F(n), n = sqrt(X_r) X_t (numero de particulas): entonces p~_r = p~_t.
Fondo estatico: f = f0(l) (la normalizacion de las etiquetas es libre); la ecuacion de f es el equilibrio
hidrostatico anisotropo  p~_r' = -(rho~ + p~_r) Phi' - 2 (p~_r - p~_t) r'/r.  Perturbacion: f = f0 + y(t, l).
Con la metrica  A = e^{2 Phi}(1 + 2a), b = beta, B = 1 + 2c, C = r^2 + 2 r d  (Fase 2f):
    L_m = sqrt(A B + b^2) C F(X_r, X_t, f),   X_r = (-B f_t^2 + 2 b f_t f_l + A f_l^2)/(A B + b^2).
Simbolos nuevos: ff_k = d^k f0/dl^k y Fs_ijk = d^i_{X_r} d^j_{X_t} d^k_f F en el fondo, con
    D_l Fs_ijk = Fs_(i+1)jk X_r' + Fs_i(j+1)k X_t' + Fs_ij(k+1) ff_1,   X_r = ff_1^2,  X_t = ff_0^2/r^2.

Comprobaciones:
  1. ecuaciones del fondo de la materia: E_A = -r^2 rho~/(2N), E_B = N r^2 p~_r/2, E_C = N p~_t, E_b = 0, y
     ff_1 E_f = N r^2 [p~_r' + (rho~ + p~_r) Phi' + 2 (p~_r - p~_t) r'/r]  (equilibrio anisotropo);
  2. identidad de Noether del fondo para la materia sola (con la ecuacion de f): identicamente nula;
  3. limite RG (G4 = 1, sin escalar): E_B y E_A dan G^l_l = p~_r/2 y G^t_t = -rho~/2;
  4. modos gauge de la teoria completa con el solido (zeta: a = Phi' zeta, beta = zeta_t, c = zeta', d = r' zeta,
     q = phi' zeta, y = f0' zeta; tau: a = tau_t, beta = -E tau', y = 0) sobre chorros exactos en la capa de masas
     con solido: las 6 ecuaciones lineales se anulan para todo omega.
Escribe resultados/solido_simbolico.pkl y resultados/generar_solido.json.
Uso: python generar_solido.py
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
import chorros as CH                                       # noqa: E402
from chorros import Phi, rr, PP, Ee, eP, Xs              # noqa: E402
from simbolos_solido import ff, Fs, Xrb, Xtb, D, chorro_capa_solido, rho, p_r, p_t   # noqa: E402

T0 = time.time()
CHK = {}


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


d0 = pickle.load(open(os.path.join(F2F, "resultados", "bH_simbolico.pkl"), "rb"))
fe_g, lin_g, L2_g, amp_g, iw = d0["fondo_eq"], d0["lineales"], d0["L2"], d0["amp"], d0["iw"]

# ------------------------------------------------------------------ 1. lagrangiana del solido y fondo
A_, b_, B_, C_, ft, fl, fv = sp.symbols("A_ b_ B_ C_ ft fl fv", real=True)
U = (A_, b_, B_, C_, ft, fl, fv)
Xr = (-B_ * ft ** 2 + 2 * b_ * ft * fl + A_ * fl ** 2) / (A_ * B_ + b_ ** 2)
Xt = fv ** 2 / C_
GRADO = 2
Fpol = 0
for i in range(GRADO + 1):
    for j in range(GRADO + 1 - i):
        for k in range(GRADO + 1 - i - j):
            Fpol += Fs(i, j, k) * (Xr - Xrb) ** i * (Xt - Xtb) ** j * (fv - ff[0]) ** k / (
                sp.factorial(i) * sp.factorial(j) * sp.factorial(k))
Lm = sp.sqrt(A_ * B_ + b_ ** 2) * C_ * Fpol
BG = {A_: Ee, b_: 0, B_: 1, C_: rr[0] ** 2, ft: 0, fl: ff[1], fv: ff[0]}


def en_bg(e):
    return sp.expand(sp.simplify(e.subs(BG)).subs(sp.sqrt(Ee), eP))


fe_m = {"A": en_bg(sp.diff(Lm, A_)), "b": en_bg(sp.diff(Lm, b_)), "B": en_bg(sp.diff(Lm, B_)),
        "C": en_bg(sp.diff(Lm, C_)), "p": sp.Integer(0)}
fe_m["f"] = sp.expand(en_bg(sp.diff(Lm, fv)) - D(en_bg(sp.diff(Lm, fl))))
aE = {eP: sp.sqrt(Ee)}
assert sp.simplify((fe_m["A"] + rr[0] ** 2 * rho / (2 * eP)).subs(aE)) == 0
assert sp.simplify((fe_m["B"] - eP * rr[0] ** 2 * p_r / 2).subs(aE)) == 0
assert sp.simplify((fe_m["C"] - eP * p_t).subs(aE)) == 0
assert fe_m["b"] == 0
# equilibrio anisotropo: ff1 E_f = N r^2 [p_r' + (rho + p_r) Phi' + 2 (p_r - p_t) r'/r]
tov = eP * rr[0] ** 2 * (D(p_r) + (rho + p_r) * Phi[1] + 2 * (p_r - p_t) * rr[1] / rr[0])
assert sp.simplify(sp.expand(ff[1] * fe_m["f"] - tov).subs(aE)) == 0
CHK["fondo_materia"] = "E_A = -r^2 rho~/(2N), E_B = N r^2 p~_r/2, E_C = N p~_t, E_b = 0; ff1 E_f = equilibrio anisotropo"
log("1. fondo del solido: E_A = -r^2 rho~/(2N), E_B = N r^2 p~_r/2, E_C = N p~_t, E_b = 0;"
    " ff1 E_f = N r^2 [p~_r' + (rho~ + p~_r) Phi' + 2 (p~_r - p~_t) r'/r]")

# 2. Noether de la materia sola: E_A A' + E_C C' + E_f f0' - 2 D_l E_B = 0 (B = 1)
noe = sp.expand((fe_m["A"] * 2 * Phi[1] * Ee + fe_m["C"] * 2 * rr[0] * rr[1] + fe_m["f"] * ff[1]
                 - 2 * D(fe_m["B"])).subs(aE))
assert sp.simplify(noe) == 0
CHK["noether_solido"] = "0 (identica, con la ecuacion de f)"
log("2. identidad de Noether del fondo para el solido: identicamente nula")

fe_t = {k: sp.expand(fe_g.get(k, 0) + fe_m[k].subs(aE)) for k in ("A", "b", "B", "C", "p")}
fe_t["f"] = sp.expand(fe_m["f"].subs(aE))

# 3. limite RG
rg = {}
for e in list(fe_t.values()):
    for s in e.free_symbols:
        if CH.es_G(s):
            rg[s] = 0
rg[sp.Symbol("G4_00", real=True)] = 1
EB = sp.expand(fe_t["B"].subs(rg).subs(Xs, -PP[1] ** 2 / 2))
EA = sp.expand(fe_t["A"].subs(rg).subs(Xs, -PP[1] ** 2 / 2))
N = sp.sqrt(Ee)
tovB = (rr[1] ** 2 - 1) / rr[0] ** 2 + 2 * rr[1] * Phi[1] / rr[0] - p_r / 2
tovA = (2 * rr[0] * rr[2] + rr[1] ** 2 - 1) / rr[0] ** 2 + rho / 2
cB = sp.factor(sp.cancel(sp.together(EB / tovB)))
cA = sp.factor(sp.cancel(sp.together(EA / tovA)))
assert cB.free_symbols <= {rr[0], Ee} and cA.free_symbols <= {rr[0], Ee}, (cB, cA)
CHK["RG"] = dict(EB_sobre=str(cB), EA_sobre=str(cA))
log(f"3. limite RG: E_B = ({cB}) [G^l_l - p~_r/2],  E_A = ({cA}) [G^t_t + rho~/2]")

# ------------------------------------------------------------------ 4. accion cuadratica del solido
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
PERT = ("a", "be", "c", "d", "q", "y")
pj = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in PERT for o in ORD}
dU = {A_: 2 * Ee * pj[("a", (0, 0))], b_: pj[("be", (0, 0))], B_: 2 * pj[("c", (0, 0))],
      C_: 2 * rr[0] * pj[("d", (0, 0))], ft: pj[("y", (1, 0))], fl: pj[("y", (0, 1))], fv: pj[("y", (0, 0))]}
L2_m = 0
for i, u in enumerate(U):
    for j, v in enumerate(U):
        if j < i:
            continue
        h = en_bg(sp.diff(Lm, u, v))
        if h == 0:
            continue
        L2_m += (sp.Rational(1, 2) if i == j else 1) * h * dU[u] * dU[v]
L2_m = sp.expand(L2_m.subs(aE))
log(f"4. accion cuadratica del solido: {len(sp.Add.make_args(L2_m))} terminos")

NMAX = 5
amp = dict(amp_g)
for k in range(NMAX + 1):
    amp[("y", k)] = sp.Symbol(f"y{k}", real=True)
modo = {pj[(v, o)]: iw ** o[0] * amp[(v, o[1])] for v in PERT for o in ORD}
AMPS = [[amp[(v, k)] for k in range(NMAX + 1)] for v in PERT]


def DL(e):
    return D(e, AMPS)


lin_m = {}
for v in PERT:
    tot = 0
    for o in ORD:
        Q = sp.diff(L2_m, pj[(v, o)]).subs(modo) * (-1) ** (o[0] + o[1]) * iw ** o[0]
        for _ in range(o[1]):
            Q = DL(Q)
        tot += Q
    lin_m[v] = sp.expand(tot)
log("   ecuaciones lineales del solido: " + ", ".join(f"{v}: {len(sp.Add.make_args(e))}" for v, e in lin_m.items()))
lin_t = {v: sp.expand(lin_g.get(v, 0) + lin_m[v]) for v in PERT}


# ------------------------------------------------------------------ 5. modos gauge en la capa de masas con solido
zt = sp.symbols("zt0:7", real=True)
ta = sp.symbols("ta0:7", real=True)


def Dg(e):
    return D(e, [list(zt), list(ta)])


campos = {"a": Phi[1] * zt[0] + iw * ta[0], "be": iw * zt[0] - Ee * ta[1], "c": zt[1], "d": rr[1] * zt[0],
          "q": PP[1] * zt[0], "y": ff[1] * zt[0]}
subs_g = {}
for v, e in campos.items():
    ek = e
    for k in range(NMAX + 1):
        subs_g[amp[(v, k)]] = ek
        ek = Dg(ek)
for sem in range(2):
    val = chorro_capa_solido(fe_t, 900 + sem, orden=5)
    for k in range(7):
        val[zt[k]] = val.r()
        val[ta[k]] = val.r()
    val[iw] = sp.Rational(3, 7)
    for ec in PERT:
        r_ = val.ev(lin_t[ec].xreplace(subs_g))
        assert r_ == 0, (ec, r_)
CHK["modos_gauge"] = "exactos: 6 ecuaciones, 2 chorros en la capa de masas con solido, omega racional"
log("5. modos gauge (zeta, tau) con solido: las 6 ecuaciones lineales se anulan exactamente (2 chorros)")

with open(os.path.join(AQUI, "resultados", "solido_simbolico.pkl"), "wb") as fh:
    pickle.dump(dict(fondo_eq=fe_t, lineales=lin_t, L2_m=L2_m, L2_g=L2_g, ff=ff, amp=amp, iw=iw, fe_m=fe_m,
                     lin_m=lin_m, rho=rho, p_r=p_r, p_t=p_t), fh)
json.dump(CHK, open(os.path.join(AQUI, "resultados", "generar_solido.json"), "w"), indent=1)
log("guardado resultados/solido_simbolico.pkl")
