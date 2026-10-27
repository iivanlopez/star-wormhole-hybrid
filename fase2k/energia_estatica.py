# -*- coding: utf-8 -*-
"""
Fase 2k: energia de las perturbaciones estaticas en la forma ADM cerrada (adm_estatico.py)
    Lest = Nc r^2 A2 - A3 (r^2 v)' + (2 A4/Nc) r' (v^2 r)' + Nc B4 (2 - 4 r r'' - 2 r'^2),   ' = d/drho.
Perturbaciones estaticas (gauge de longitud propia del corte, h = 0): dNc = nu, dv = w, dr = d.  Energia
V = -int L2 (L2 = parte cuadratica), sobre configuraciones que cumplen las ecuaciones de los auxiliares (nu, w).

Familia dX = 0 (nu = 0; X = qt^2/2Nc^2 no cambia, como la prueba d = N de las Fases 2g-2j):
  w sale de la ligadura de momento dE_v = 0 y d debe cumplir la del lapso dE_Nc = 0.  Se comprueba:
  (1) caso v = 0 (phi = qt, sin desplazamiento): identidad  -L2(0, 0, d) = -2 beta r^2 [(d/r)']^2 + derivada total
      (beta = Nc B4) usando solo E_r;  y dE_Nc(0, 0, d) = -4 beta_N D[d],  D[d] = (r d')' + (1 - r'^2 - r r'') d/r,
      usando E_Nc.  Ninguna funcion libre (A2, A4, ni derivadas) aparece.
  (2) caso general (v != 0): forma de -L2(0, w, d) con w de la ligadura, en terminos de beta r^2[(d/r)']^2 y un resto.
Uso: python energia_estatica.py
"""
import os
import json
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
rho = sp.Symbol("rho", real=True)
eps = sp.Symbol("eps")
F = lambda n: sp.Function(n)(rho)                       # noqa: E731
r, v, p, a3, Q, be = F("r"), F("v"), F("p"), F("a3"), F("Q"), F("be")      # p = Nc A2, a3 = A3, Q = 2A4/Nc, be = Nc B4
pN, a3N, QN, beN = F("pN"), F("a3N"), F("QN"), F("beN")                    # derivadas respecto de Nc
w, d = F("w"), F("d")
D = lambda e, k=1: sp.diff(e, rho, k)                    # noqa: E731
RES = {}


def lag(vv, rr):
    """Lest a lapso fijo (nu = 0): las funciones de Nc son funciones dadas de rho."""
    return p * rr ** 2 - a3 * D(rr ** 2 * vv) + Q * D(rr) * D(vv ** 2 * rr) + be * (2 - 4 * rr * D(rr, 2) - 2 * D(rr) ** 2)


def EL(L, f):
    return sp.expand(sp.euler_equations(L, [f], [rho])[0].lhs)


L0 = lag(v, r)
Er = EL(L0, r)
Ev = EL(L0, v)
# ecuacion del lapso (derivada parcial en Nc con Nc' a traves de a3' = a3N Nc': se usa la forma original sin integrar)
ENc = r ** 2 * pN - a3N * D(r ** 2 * v) + QN * D(r) * D(v ** 2 * r) + beN * (2 - 4 * r * D(r, 2) - 2 * D(r) ** 2)
L2 = sp.expand(sp.diff(lag(v + eps * w, r + eps * d), eps, 2).subs(eps, 0) / 2)
dEv = sp.expand(sp.diff(EL(lag(v + eps * w, r + eps * d), v).subs({}), eps).subs(eps, 0)) if False else None
# ligaduras linealizadas a nu = 0
Ev_eps = sp.euler_equations(lag(v + eps * w, r + eps * d), [w], [rho])[0].lhs   # dL/dw = eps E_v(v + eps w, ...)
dEv = sp.expand(sp.diff(Ev_eps, eps, 2).subs(eps, 0) / 2)                        # parte lineal en (w, d)
dENc = sp.expand(sp.diff(ENc.subs({v: v + eps * w, r: r + eps * d}).doit(), eps).subs(eps, 0))


def es_derivada_total(e, campos):
    return all(sp.simplify(EL(e, f)) == 0 for f in campos)


# ------------------------------------------------------------------ (1) v = 0
sv0 = {v: 0}
L2v0 = sp.expand(L2.subs(sv0).doit().subs(w, 0))
ident = -L2v0 + 2 * be * r ** 2 * D(d / r) ** 2
# usar E_r (v = 0) para eliminar p:  E_r = 0  ->  p = ...
Er0 = sp.expand(Er.subs(sv0).doit())
p_sol = sp.solve(Er0, p)[0]
ident_s = sp.expand(ident.subs(p, p_sol))
td = es_derivada_total(ident_s, [d])
RES["v0_identidad_derivada_total"] = td
print(f"(1) v = 0:  -L2(0,0,d) + 2 beta r^2 [(d/r)']^2 es derivada total (con E_r): {td}", flush=True)
assert td
ENc0 = sp.expand(ENc.subs(sv0).doit())
pN_sol = sp.solve(ENc0, pN)[0]
DD = D(r * D(d)) + (1 - D(r) ** 2 - r * D(r, 2)) * d / r
lig = sp.simplify(dENc.subs(sv0).doit().subs(w, 0).subs(pN, pN_sol) + 4 * beN * DD)
RES["v0_ligadura_lapso"] = str(lig)
print(f"    dE_Nc(0,0,d) + 4 beta_N D[d] = {lig}  (con E_Nc)", flush=True)
assert lig == 0

# ------------------------------------------------------------------ (2) v != 0: estructura
kw = sp.expand(dEv.coeff(w))
print(f"(2) ligadura de momento linealizada: coeficiente de w = {sp.factor(kw)}", flush=True)
Ev_s = sp.solve(Ev, D(Q, 1))
print(f"    E_v de fondo: {sp.factor(Ev)}", flush=True)
RES["kw"] = str(sp.factor(kw))
RES["Ev"] = str(sp.factor(Ev))
# -L2(0, 0, d) con v != 0 frente a -2 beta r^2 [(d/r)']^2 usando E_r
p_solg = sp.solve(Er, p)[0]
rest = sp.expand((-L2.subs(w, 0) + 2 * be * r ** 2 * D(d / r) ** 2).subs(p, p_solg))
# forma canonica del resto: A d'^2 + C d^2 salvo derivadas totales (se identifica con dos pruebas de coeficientes)
A_ = sp.Symbol("A_")
print("    resto de -L2(0,0,d) (con E_r):", sp.factor(sp.simplify(rest)), flush=True)
RES["resto_v"] = str(sp.factor(sp.simplify(rest)))
json.dump(RES, open(os.path.join(AQUI, "resultados", "energia_estatica.json"), "w"), indent=1)
