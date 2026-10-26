# -*- coding: utf-8 -*-
"""
Fase 2g, parte II: fluido perfecto en el sector l = 0 de la teoria cuartica mas alla de Horndeski.

El fluido (barotropo, irrotacional en l = 0) se describe con un campo chi y la lagrangiana P(Y), Y = -(1/2)(d chi)^2:
    S = (1/16 pi) int sqrt(-g) [ L_GLPV + P(Y) ],
    presion p~ = P,  densidad rho~ = 2 Y P_Y - P,  u_m = -d_m chi/sqrt(2Y),  c_s^2 = P_Y/(P_Y + 2 Y P_YY)
(p~ = 16 pi p, rho~ = 16 pi rho en unidades G = c = 1).  Fondo estatico: chi = t (la normalizacion de chi se absorbe
en P), Y_b = 1/(2 N^2): p y rho son funciones del lapso y el equilibrio hidrostatico p' = -(rho + p) Phi' se cumple
solo (la ecuacion de chi es trivial en el fondo).  Perturbacion: chi = t + x(t, l).

La parte GLPV (ecuaciones del fondo y accion cuadratica de la Fase 2f, validas fuera de la capa de masas) se suma a
la del fluido, que es algebraica en (A, b, B, C) y en las primeras derivadas de chi:
    L_m = sqrt(A B + b^2) C P(Y),   Y = (B chi_t^2 - 2 b chi_t chi_l - A chi_l^2) / (2 (A B + b^2)).
Simbolos nuevos: Pf0, Pf1, ... = P, P_Y, P_YY, ... en Y_b;  D_l Pf_k = Pf_{k+1} Y_b' = -Pf_{k+1} Phi1/E.

Comprobaciones:
  1. identidad de Noether del fondo para la parte del fluido sola (reparametrizaciones de l): se anula identicamente
     (el equilibrio hidrostatico esta incluido);
  2. limite RG (G4 = 1, sin escalar): E_B y E_A dan las ecuaciones de Tolman-Oppenheimer-Volkoff en longitud propia;
  3. modos gauge de la teoria completa con fluido (Lie de xi = tau d_t + zeta d_l: a = Phi' zeta + tau_t,
     beta = zeta_t - E tau', c = zeta', d = r' zeta, q = phi' zeta, x = tau) sobre chorros exactos en la capa de masas
     con fluido: las 6 ecuaciones lineales se anulan para todo omega.
Escribe resultados/fluido_simbolico.pkl (fondo_eq y lineales totales, L2 del fluido) y resultados/generar_fluido.json.
Uso: python generar_fluido.py
"""
import os
import sys
import json
import time
import random
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
import chorros as CH                                       # noqa: E402
from chorros import Phi, rr, PP, Ee, eP, Xs              # noqa: E402

T0 = time.time()
CHK = {}


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


d0 = pickle.load(open(os.path.join(F2F, "resultados", "bH_simbolico.pkl"), "rb"))
fe_g, lin_g, L2_g, amp_g, iw = d0["fondo_eq"], d0["lineales"], d0["L2"], d0["amp"], d0["iw"]

NPF = 8
Pf = sp.symbols(f"Pf0:{NPF}", real=True)


def D(e, amps=()):
    """Derivada total en l: la de chorros.D mas Pf_k -> -Pf_{k+1} Phi1/E."""
    out = CH.D(e, amps)
    for k in range(NPF - 1):
        if e.has(Pf[k]):
            out += sp.diff(e, Pf[k]) * (-Pf[k + 1] * Phi[1] / Ee)
    return sp.expand(out)


# ------------------------------------------------------------------ 1. lagrangiana del fluido y fondo
A_, b_, B_, C_, ct, cl = sp.symbols("A_ b_ B_ C_ ct cl", real=True)
Yb = 1 / (2 * Ee)
Y = (B_ * ct ** 2 - 2 * b_ * ct * cl - A_ * cl ** 2) / (2 * (A_ * B_ + b_ ** 2))
Pdes = sum(Pf[k] * (Y - Yb) ** k / sp.factorial(k) for k in range(4))
Lm = sp.sqrt(A_ * B_ + b_ ** 2) * C_ * Pdes
BG = {A_: Ee, b_: 0, B_: 1, C_: rr[0] ** 2, ct: 1, cl: 0}
U = (A_, b_, B_, C_, ct, cl)


def en_bg(e):
    return sp.expand(sp.simplify(e.subs(BG)).subs(sp.sqrt(Ee), eP))


fe_m = {"A": en_bg(sp.diff(Lm, A_)), "b": en_bg(sp.diff(Lm, b_)), "B": en_bg(sp.diff(Lm, B_)),
        "C": en_bg(sp.diff(Lm, C_)), "p": sp.Integer(0)}
log("fondo del fluido: " + ", ".join(f"E_{k} = {v}" for k, v in fe_m.items()))
rho = 2 * Yb * Pf[1] - Pf[0]
pres = Pf[0]
aE = {eP: sp.sqrt(Ee)}
assert sp.simplify((fe_m["A"] + rr[0] ** 2 * rho / (2 * eP)).subs(aE)) == 0
assert sp.simplify((fe_m["B"] - eP * rr[0] ** 2 * pres / 2).subs(aE)) == 0
assert sp.simplify((fe_m["C"] - eP * pres).subs(aE)) == 0
assert fe_m["b"] == 0
log("   E_A^m = -r^2 rho~/(2N),  E_B^m = N r^2 p~/2,  E_C^m = N p~  (rho~ = 2 Y P_Y - P, p~ = P)")

# 1a. Noether del fluido: E_A A' + E_C C' - 2 D_l E_B = 0 (B = 1)
noe = sp.expand((fe_m["A"] * 2 * Phi[1] * Ee + fe_m["C"] * 2 * rr[0] * rr[1] - 2 * D(fe_m["B"])).subs(eP, sp.sqrt(Ee)))
assert sp.simplify(noe) == 0
CHK["noether_fluido"] = "0 (identica)"
log("1. identidad de Noether del fondo para el fluido: identicamente nula (equilibrio hidrostatico incluido)")

# ------------------------------------------------------------------ 2. accion cuadratica del fluido
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
PERT = ("a", "be", "c", "d", "q", "x")
pj = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in PERT for o in ORD}
dU = {A_: 2 * Ee * pj[("a", (0, 0))], b_: pj[("be", (0, 0))], B_: 2 * pj[("c", (0, 0))],
      C_: 2 * rr[0] * pj[("d", (0, 0))], ct: pj[("x", (1, 0))], cl: pj[("x", (0, 1))]}
L2_m = 0
for i, u in enumerate(U):
    for j, v in enumerate(U):
        if j < i:
            continue
        h = en_bg(sp.diff(Lm, u, v))
        if h == 0:
            continue
        L2_m += (sp.Rational(1, 2) if i == j else 1) * h * dU[u] * dU[v]
L2_m = sp.expand(L2_m)
log(f"2. accion cuadratica del fluido: {len(sp.Add.make_args(L2_m))} terminos")
log(f"   L2_m = {sp.factor_terms(L2_m)}")

# ecuaciones lineales del fluido (modos e^{i omega t}) con los mismos simbolos de amplitud que la Fase 2f
NMAX = 5
amp = dict(amp_g)
for k in range(NMAX + 1):
    amp[("x", k)] = sp.Symbol(f"x{k}", real=True)
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
log("   ecuaciones lineales del fluido: " + ", ".join(f"{v}: {len(sp.Add.make_args(e))}" for v, e in lin_m.items()))

fe_t = {k: sp.expand(fe_g.get(k, 0) + fe_m[k].subs(eP, sp.sqrt(Ee))) for k in ("A", "b", "B", "C", "p")}
lin_t = {v: sp.expand(lin_g.get(v, 0) + lin_m[v].subs(eP, sp.sqrt(Ee))) for v in PERT}

# ------------------------------------------------------------------ 3. limite RG + fluido: TOV
rg = {}
for e in list(fe_t.values()):
    for s in e.free_symbols:
        if CH.es_G(s):
            rg[s] = 0
rg[sp.Symbol("G4_00", real=True)] = 1
EB = sp.expand(fe_t["B"].subs(rg).subs(Xs, -PP[1] ** 2 / 2))
EA = sp.expand(fe_t["A"].subs(rg).subs(Xs, -PP[1] ** 2 / 2))
N = sp.sqrt(Ee)
p_, rho_ = pres, rho.subs(eP, N)
# G^l_l = (r'^2 - 1)/r^2 + 2 r' Phi'/r = p~/2 ;  G^t_t = (2 r r'' + r'^2 - 1)/r^2 = -rho~/2
tovB = (rr[1] ** 2 - 1) / rr[0] ** 2 + 2 * rr[1] * Phi[1] / rr[0] - p_ / 2
tovA = (2 * rr[0] * rr[2] + rr[1] ** 2 - 1) / rr[0] ** 2 + rho_ / 2
cB = sp.factor(sp.cancel(sp.together(EB / tovB)))
cA = sp.factor(sp.cancel(sp.together(EA / tovA)))
assert cB.free_symbols <= {rr[0], Ee} and cA.free_symbols <= {rr[0], Ee}, (cB, cA)
CHK["TOV"] = dict(EB_sobre_tov=str(cB), EA_sobre_tov=str(cA))
log(f"3. limite RG + fluido: E_B = ({cB}) x [G^l_l - p~/2],  E_A = ({cA}) x [G^t_t + rho~/2]  (TOV)")


# ------------------------------------------------------------------ 4. modos gauge en la capa de masas con fluido
def chorro_capa_fluido(semilla, orden=3):
    v = CH.Valores(semilla)
    v[Ee] = sp.Rational(121, 100)
    v[eP] = sp.Rational(11, 10)
    v[Phi[1]] = v.r()
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[rr[1]] = v.r()
    v[PP[0]] = v.r()
    v[PP[1]] = v.r()
    v[Xs] = -v[PP[1]] ** 2 / 2
    for k in range(NPF):
        v[Pf[k]] = v.r()
    g200 = CH.gsym("G2", 0, 0)
    eB = v.ev(fe_t["B"], salvo=(g200,))
    c1 = eB.coeff(g200)
    v[g200] = -sp.expand(eB - c1 * g200) / c1
    eqs = [fe_t["A"], fe_t["C"], fe_t["p"]]
    for k in range(2, orden + 1):
        inc = [Phi[k], rr[k], PP[k]]
        ev = [v.ev(e) for e in eqs]
        Mk = sp.Matrix([[x.coeff(u) for u in inc] for x in ev])
        bk = sp.Matrix([sp.expand(x - sum(x.coeff(u) * u for u in inc)) for x in ev])
        assert all(not x.free_symbols for x in bk) and Mk.det() != 0
        sol = Mk.LUsolve(-bk)
        for u, s in zip(inc, sol):
            s = sp.nsimplify(s) if not s.is_Rational else s
            assert s.is_Rational, s
            v[u] = s
        if k < orden:
            eqs = [D(e) for e in eqs]
    # comprobacion: E_B (y su derivada) en la capa
    assert v.ev(fe_t["B"]) == 0
    return v


# perturbacion gauge: zeta(l), tau(l) (chorros arbitrarios) con dependencia e^{i omega t}
zt = sp.symbols("zt0:6", real=True)
ta = sp.symbols("ta0:6", real=True)


def Dg(e):
    return D(e, [list(zt), list(ta)])


modos = {}
z0 = [zt[0]]
for k in range(1, 5):
    z0.append(Dg(z0[-1]))
t0 = [ta[0]]
for k in range(1, 5):
    t0.append(Dg(t0[-1]))
# a = Phi' zeta + iw tau,  beta = iw zeta - E tau',  c = zeta',  d = r' zeta,  q = phi' zeta,  x = tau
campos = {"a": Phi[1] * zt[0] + iw * ta[0], "be": iw * zt[0] - Ee * ta[1], "c": zt[1], "d": rr[1] * zt[0],
          "q": PP[1] * zt[0], "x": ta[0]}
subs_g = {}
for v, e in campos.items():
    ek = e
    for k in range(NMAX + 1):
        subs_g[amp[(v, k)]] = ek
        ek = Dg(ek)
res = {}
for sem in range(2):
    val = chorro_capa_fluido(700 + sem, orden=4)
    for k in range(6):
        val[zt[k]] = val.r()
        val[ta[k]] = val.r()
    val[iw] = sp.Rational(3, 7)
    for ec in PERT:
        r_ = val.ev(lin_t[ec].xreplace(subs_g))
        res[(sem, ec)] = str(r_)
        assert r_ == 0, (ec, r_)
CHK["modos_gauge"] = "exactos: 6 ecuaciones, 2 chorros en la capa de masas con fluido, omega racional"
log("4. modos gauge (zeta, tau) con fluido: las 6 ecuaciones lineales se anulan exactamente (2 chorros)")

with open(os.path.join(AQUI, "resultados", "fluido_simbolico.pkl"), "wb") as fh:
    pickle.dump(dict(fondo_eq=fe_t, lineales=lin_t, L2_m=L2_m, L2_g=L2_g, Pf=Pf, amp=amp, iw=iw, fe_m=fe_m,
                     lin_m=lin_m), fh)
json.dump(CHK, open(os.path.join(AQUI, "resultados", "generar_fluido.json"), "w"), indent=1)
log("guardado resultados/fluido_simbolico.pkl")
