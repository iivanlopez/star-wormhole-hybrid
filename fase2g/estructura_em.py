# -*- coding: utf-8 -*-
"""
Fase 2g, parte IV: estructura del monopolo con carga (generar_em.py) en el gauge unitario, sobre chorros exactos en la
capa de masas con carga.  Comprueba si la ecuacion de a es exactamente (N^2 W)' (la carga no aporta terminos en a ni
en beta, al ser rho + p_r = 0) y lee los coeficientes de a y beta en E_c, E_d.
Uso: python estructura_em.py
"""
import os
import sys
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
import chorros as CH                                       # noqa: E402
from chorros import Phi, rr, PP, Ee, eP, Xs              # noqa: E402

d0 = pickle.load(open(os.path.join(F2F, "resultados", "bH_simbolico.pkl"), "rb"))
amp, iw = d0["amp"], d0["iw"]
d = pickle.load(open(os.path.join(AQUI, "resultados", "em_simbolico.pkl"), "rb"))
fe_t, lin_t, q2 = d["fondo_eq"], d["lineales"], d["q2"]
mr = pickle.load(open(os.path.join(F2F, "resultados", "monopolo_reduccion.pkl"), "rb"))
Xb = {Xs: -PP[1] ** 2 / 2}
gam, al1, al0 = (sp.expand(mr[k].xreplace(Xb)) for k in ("gam", "al1", "al0"))
gauge = {amp[("q", k)]: 0 for k in range(6)}
cs = [amp[("c", k)] for k in range(6)]
ds = [amp[("d", k)] for k in range(6)]
Wexp = gam * cs[0] + al1 * ds[1] + al0 * ds[0]
DW = CH.D(Ee * Wexp, [cs, ds])
Ea = sp.expand(lin_t["a"].xreplace(gauge))
Ec = sp.expand(lin_t["c"].xreplace(gauge))
Ed = sp.expand(lin_t["d"].xreplace(gauge))


def capa(semilla, orden=5):
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
            v[u] = s
        if k < orden:
            eqs = [CH.D(e) for e in eqs]
    return v


for sem in range(2):
    v = capa(950 + sem)
    dif = sp.expand(v.ev(Ea) - v.ev(DW))
    print(f"chorro {sem}: E_a - (N^2 W)' = {dif}")
    for nom, E in (("E_c", Ec), ("E_d", Ed)):
        pol = sp.Poly(E, *[amp[(f, k)] for f in ("a", "be") for k in range(6)])
        print("  ", nom, [(m, float(v.ev(c))) for m, c in zip(pol.monoms(), pol.coeffs()) if any(m)
                          and not v.ev(c).free_symbols])
