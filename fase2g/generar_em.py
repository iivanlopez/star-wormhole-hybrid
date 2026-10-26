# -*- coding: utf-8 -*-
"""
Fase 2g, parte IV: campo electrico radial ("carga sin carga" de Wheeler) en el sector l = 0 de la GLPV cuartica.

En simetria esferica el campo de Maxwell no tiene grados de libertad en l = 0: la ley de Gauss fija el flujo Q y,
eliminando F_tl (rutiano), su accion es la de una "materia" algebraica en la metrica:
    L_em = -q2 sqrt(A B + b^2)/C          (q2 = Q^2 en la normalizacion de la accion, rho~ = -p~_r = p~_t = q2/r^4).
Se suma a la parte GLPV de la Fase 2f (valida fuera de la capa de masas), se comprueban la identidad de Noether del
fondo y los modos gauge en la capa de masas con carga, y se analiza la estructura del monopolo en el gauge unitario:
    coeficientes de a^2, beta^2 y de a (c, d) en la accion; forma de la ligadura de momento con carga.
Escribe resultados/em_simbolico.pkl.
Uso: python generar_em.py
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
import chorros as CH                                       # noqa: E402
from chorros import Phi, rr, PP, Ee, eP, Xs              # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


d0 = pickle.load(open(os.path.join(F2F, "resultados", "bH_simbolico.pkl"), "rb"))
fe_g, lin_g, L2_g, amp, iw = d0["fondo_eq"], d0["lineales"], d0["L2"], d0["amp"], d0["iw"]
q2 = sp.Symbol("q2", positive=True)
A_, b_, B_, C_ = sp.symbols("A_ b_ B_ C_", real=True)
Lm = -q2 * sp.sqrt(A_ * B_ + b_ ** 2) / C_
BG = {A_: Ee, b_: 0, B_: 1, C_: rr[0] ** 2}
U = (A_, b_, B_, C_)


def en_bg(e):
    return sp.expand(sp.simplify(e.subs(BG)).subs(sp.sqrt(Ee), eP))


fe_m = {"A": en_bg(sp.diff(Lm, A_)), "b": en_bg(sp.diff(Lm, b_)), "B": en_bg(sp.diff(Lm, B_)),
        "C": en_bg(sp.diff(Lm, C_)), "p": sp.Integer(0)}
log("fondo de la carga: " + ", ".join(f"E_{k} = {v}" for k, v in fe_m.items()))
aE = {eP: sp.sqrt(Ee)}
noe = sp.expand((fe_m["A"] * 2 * Phi[1] * Ee + fe_m["C"] * 2 * rr[0] * rr[1] - 2 * CH.D(fe_m["B"])).subs(aE))
assert sp.simplify(noe) == 0
log("identidad de Noether de la carga sola: nula (el tensor de Maxwell se conserva)")

ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
PERT = ("a", "be", "c", "d", "q")
pj = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in PERT for o in ORD}
dU = {A_: 2 * Ee * pj[("a", (0, 0))], b_: pj[("be", (0, 0))], B_: 2 * pj[("c", (0, 0))], C_: 2 * rr[0] * pj[("d", (0, 0))]}
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
log(f"accion cuadratica de la carga: {sp.factor_terms(L2_m)}")

NMAX = 5
modo = {pj[(v, o)]: iw ** o[0] * amp[(v, o[1])] for v in PERT for o in ORD}
AMPS = [[amp[(v, k)] for k in range(NMAX + 1)] for v in PERT]
lin_m = {}
for v in PERT:
    tot = 0
    for o in ORD:
        Q = sp.diff(L2_m, pj[(v, o)]).subs(modo) * (-1) ** (o[0] + o[1]) * iw ** o[0]
        for _ in range(o[1]):
            Q = CH.D(Q, AMPS)
        tot += Q
    lin_m[v] = sp.expand(tot)
fe_t = {k: sp.expand(fe_g.get(k, 0) + fe_m[k].subs(aE)) for k in ("A", "b", "B", "C", "p")}
lin_t = {v: sp.expand(lin_g.get(v, 0) + lin_m[v]) for v in PERT}


def capa(semilla, orden=4):
    v = CH.Valores(semilla)
    v[Ee] = sp.Rational(121, 100)
    v[eP] = sp.Rational(11, 10)
    v[Phi[1]] = v.r()
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[rr[1]] = v.r()
    v[PP[0]] = v.r()
    v[PP[1]] = v.r()
    v[Xs] = -v[PP[1]] ** 2 / 2
    v[q2] = sp.Rational(v.rnd.randint(100, 900), 1000)
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
        sol = Mk.LUsolve(-bk)
        for u, s in zip(inc, sol):
            assert s.is_Rational
            v[u] = s
        if k < orden:
            eqs = [CH.D(e) for e in eqs]
    return v


# modos gauge con carga (sin campo extra: x no existe)
zt = sp.symbols("zt0:6", real=True)
ta = sp.symbols("ta0:6", real=True)
campos = {"a": Phi[1] * zt[0] + iw * ta[0], "be": iw * zt[0] - Ee * ta[1], "c": zt[1], "d": rr[1] * zt[0],
          "q": PP[1] * zt[0]}
subs_g = {}
for vv, e in campos.items():
    ek = e
    for k in range(NMAX + 1):
        subs_g[amp[(vv, k)]] = ek
        ek = CH.D(ek, [list(zt), list(ta)])
for sem in range(2):
    val = capa(800 + sem, orden=4)
    for k in range(6):
        val[zt[k]] = val.r()
        val[ta[k]] = val.r()
    val[iw] = sp.Rational(2, 9)
    for ec in PERT:
        r_ = val.ev(lin_t[ec].xreplace(subs_g))
        assert r_ == 0, (ec, r_)
log("modos gauge (zeta, tau) con carga: las 5 ecuaciones lineales se anulan exactamente (2 chorros)")

# estructura en el gauge unitario: E_beta y E_a
gauge = {amp[("q", k)]: 0 for k in range(6)}
mr = pickle.load(open(os.path.join(F2F, "resultados", "monopolo_reduccion.pkl"), "rb"))
Xb = {Xs: -PP[1] ** 2 / 2}
gam, al1, al0 = (sp.expand(mr[k].xreplace(Xb)) for k in ("gam", "al1", "al0"))
res = {}
for sem in range(2):
    v = capa(900 + sem, orden=5)
    for nombre in ("be", "a", "c", "d"):
        E = sp.expand(lin_t[nombre].xreplace(gauge))
        pol = sp.Poly(E, *[amp[(f, k)] for f in ("a", "be", "c", "d") for k in range(6)])
        txt = []
        for mono, cf in zip(pol.monoms(), pol.coeffs()):
            nm = [str(amp[(f, k)]) for f in ("a", "be", "c", "d") for k in range(6)]
            txt.append("*".join(n for n, ee in zip(nm, mono) if ee) + f": {v.ev(cf)}")
        res[(sem, nombre)] = txt
    gv, a1v, a0v = (v.ev(e) for e in (gam, al1, al0))
    log(f"chorro {sem}: gam = {float(gv):.6f}, al1 = {float(a1v):.6f}, al0 = {float(a0v):.6f}, "
        f"q2/r^2 N = {float(v[q2] / v[rr[0]] ** 2 * sp.sqrt(v[Ee])):.6f}")
    for nombre in ("be", "a"):
        log(f"   E_{nombre}: " + "; ".join(res[(sem, nombre)]))
pickle.dump(dict(fondo_eq=fe_t, lineales=lin_t, L2_m=L2_m, q2=q2), open(os.path.join(AQUI, "resultados", "em_simbolico.pkl"), "wb"))
log("guardado resultados/em_simbolico.pkl")
