# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): identidades integrales de la ecuacion de beta con la fuente cubica (nogo_beta.py)
    (beta'/r)' + K beta/r^3 = -(Q r^4 u')'/(4 r^3),     u = v^2/r^2,  K = 1 + r r'' - r'^2 .
(I1) sin dividir:  K beta/r^3 = (3/4) u (Q r')' - [ beta'/r + Q r u'/4 + (3/4) Q r' u ]'      (sobre la ecuacion)
     =>  int K beta/r^3 = (3/4) int (v/r)^2 (Q r')'            (fondo asintoticamente plano, flujo que decae).
(I2) dividiendo por beta > 0, con y = beta'/beta:
     (1/r)(y + Q r^2 u'/(8 beta))^2 + (1 + r'^2)/r^3 = -(3/4) Q u' r'/beta + Q^2 r^3 u'^2/(64 beta^2)
                                                - [beta'/(r beta) + r'/r^2 + Q r u'/(4 beta)]'     (sobre la ecuacion).
Comprobacion: se sustituye beta'' de la ecuacion y se verifica que la diferencia es D_rho de una expresion dada.
Uso: python identidad_cubica.py
"""
import os
import json
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
x = sp.Symbol("rho", real=True)
r, v, Q, be = (sp.Function(n)(x) for n in ("r", "v", "Q", "beta"))
D = lambda f, n=1: sp.diff(f, x, n)                  # noqa: E731
u = v ** 2 / r ** 2
K = 1 + r * D(r, 2) - D(r) ** 2
b2 = sp.solve(sp.Eq(D(D(be) / r) + K * be / r ** 3, -D(Q * r ** 4 * D(u)) / (4 * r ** 3)), D(be, 2))[0]


def sobre_ecuacion(e):
    return sp.simplify(sp.expand(e.subs(D(be, 2), b2)))


RES = {}
# (I1)
F1 = D(be) / r + Q * r * D(u) / 4 + sp.Rational(3, 4) * Q * D(r) * u
dif1 = sobre_ecuacion(K * be / r ** 3 - sp.Rational(3, 4) * u * D(Q * D(r)) + D(F1))
RES["I1"] = str(dif1)
print("(I1) K beta/r^3 - (3/4) u (Q r')' + [beta'/r + Q r u'/4 + (3/4) Q r' u]' =", dif1, flush=True)
assert dif1 == 0
# (I2): dividir la ecuacion por beta e integrar por partes
y = D(be) / be
lhs = (1 / r) * (y + Q * r ** 2 * D(u) / (8 * be)) ** 2 + (1 + D(r) ** 2) / r ** 3
rhs = -sp.Rational(3, 4) * Q * D(u) * D(r) / be + Q ** 2 * r ** 3 * D(u) ** 2 / (64 * be ** 2)
F2 = D(be) / (r * be) + D(r) / r ** 2 + Q * r * D(u) / (4 * be)
dif2 = sobre_ecuacion(lhs - rhs + D(F2))
RES["I2"] = str(dif2)
print("(I2) (1/r)(beta'/beta + Q r^2 u'/(8 beta))^2 + (1 + r'^2)/r^3 + (3/4) Q u' r'/beta - Q^2 r^3 u'^2/(64 beta^2) + [beta'/(r beta) "
      "+ r'/r^2 + Q r u'/(4 beta)]' =", dif2, flush=True)
assert dif2 == 0
json.dump(RES, open(os.path.join(AQUI, "resultados", "identidad_cubica.json"), "w"), indent=1)
