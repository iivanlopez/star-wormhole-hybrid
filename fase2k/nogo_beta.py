# -*- coding: utf-8 -*-
"""
Fase 2k: la ecuacion de beta = Nc B4 en el fondo estatico (forma ADM, adm_estatico.py) y la identidad integral.

Lest ~ p r^2 + A3' r^2 v + Q r'(v^2 r)' + 2 beta (1 + r'^2) + 4 beta' r r'     (p = Nc A2, Q = 2 A4/Nc, beta = Nc B4)
Primera integral (rho no aparece) = ligadura de la componente h (gauge de longitud propia):
    H = 2 beta (r'^2 - 1) + 4 beta' r r' + Q r' (v^2 r)' - p r^2 = 0 .
Con E_r y H (eliminando p) y la ligadura de momento E_v (r A3' = 2 v (Q r')') se obtiene
    (beta'/r)' + K beta / r^3 = S_v,      K = 1 + r r'' - r'^2,
con S_v = 0 si v = 0.  Se comprueba simbolicamente y se da S_v.  Para v = 0: dividiendo por beta e integrando,
    int beta'^2/(r beta^2) = -int K/r^3 = -int (1 + r'^2)/r^3 - [r'/r^2]  < 0
si beta > 0 y beta'/(r beta), r'/r^2 -> 0 en los dos extremos: contradiccion.  (Comprobado tambien: dH/drho =
-(E_r r' + E_Nc Nc' + E_v v').)
Uso: python nogo_beta.py
"""
import os
import json
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
rho = sp.Symbol("rho", real=True)
F = lambda n: sp.Function(n)(rho)                       # noqa: E731
r, v, p, a3, Q, be = F("r"), F("v"), F("p"), F("a3"), F("Q"), F("be")
D = lambda e, k=1: sp.diff(e, rho, k)                    # noqa: E731
Lest = p * r ** 2 + D(a3) * r ** 2 * v + Q * D(r) * D(v ** 2 * r) + 2 * be * (1 + D(r) ** 2) + 4 * D(be) * r * D(r)
H = 2 * be * (D(r) ** 2 - 1) + 4 * D(be) * r * D(r) + Q * D(r) * D(v ** 2 * r) - p * r ** 2
# E_r con las funciones de Nc como funciones dadas de rho (la variacion de r no las toca)
Er = sp.expand(sp.euler_equations(Lest, [r], [rho])[0].lhs)
Ev = sp.expand(sp.euler_equations(Lest, [v], [rho])[0].lhs)
p_H = sp.solve(H, p)[0]
ec = sp.expand(Er.subs(p, p_H) * sp.Rational(1, 1))
# forma normalizada: dividir por el coeficiente de beta''
cb = ec.coeff(D(be, 2))
ec_n = sp.expand(ec / cb)
K = 1 + r * D(r, 2) - D(r) ** 2
objetivo = D(D(be) / r) + K * be / r ** 3
# Er/(coef) frente a r * objetivo: el resto solo contiene v
resto = sp.simplify(ec_n - sp.expand(objetivo * r))
resto_v0 = sp.simplify(resto.subs(v, 0).doit())
print("E_r con p de H, normalizada:  r[(beta'/r)' + K beta/r^3] + resto;  resto(v = 0) =", resto_v0, flush=True)
assert resto_v0 == 0
# usar la ligadura de momento para eliminar a3''
a3pp = sp.solve(Ev, D(a3))
resto_m = sp.simplify(resto.subs(D(a3), a3pp[0]).doit()) if a3pp else resto
print("resto con v != 0 (ligadura de momento usada):", sp.factor(resto_m), flush=True)
json.dump({"resto_v0": str(resto_v0), "resto_v": str(sp.factor(resto_m)), "Ev": str(sp.factor(Ev))},
          open(os.path.join(AQUI, "resultados", "nogo_beta.json"), "w"), indent=1)
