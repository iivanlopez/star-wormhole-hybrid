# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): el simbolo principal del monopolo (simbolo_general.py) escrito con los perfiles del fondo.
Las funciones de la teoria se sustituyen por la reconstruccion: b4 = beta/Nc, a4 = Q Nc/2, a2 = p/Nc con p de H = 0,
A3' = 2 v (Q r')'/r, f_N = f'/Nc', f_NN = (f_N)'/Nc', y beta'' de su EDO.  Quedan r, v, Q, Nc, beta y derivadas.
Se busca la estructura de P2 (coeficiente de c^2, menos el coeficiente cinetico), P1, P0 y D.
Uso: python simbolo_perfiles.py
"""
import os
import sys
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
x = sp.Symbol("rho", real=True)
F = lambda n: sp.Function(n)(x)                      # noqa: E731
r, v, Q, N, be = F("r"), F("v"), F("Q"), F("Nc"), F("beta")
D = lambda f, n=1: sp.diff(f, x, n)                  # noqa: E731
K = 1 + r * D(r, 2) - D(r) ** 2
S = D(Q * r ** 4 * D(v ** 2 / r ** 2)) / (4 * r ** 3)
b2 = D(r) / r * D(be) - K * be / r ** 2 - r * S      # beta'' de la EDO


def reduce_beta(ex):
    """elimina beta'' y superiores con la EDO."""
    for _ in range(3):
        ex = ex.subs(D(be, 3), D(b2)).subs(D(be, 2), b2)
    return ex


p = (2 * be * (D(r) ** 2 - 1) + 4 * D(be) * r * D(r) + Q * D(r) * D(v ** 2 * r)) / r ** 2
th = {"b4": be / N, "a4": Q * N / 2, "a2": p / N}
val = {"N_0": N, "N_1": D(N), "R_0": r, "R_1": D(r), "R_2": D(r, 2), "V_0": v, "V_1": D(v)}
for nom, f in th.items():
    fN = D(f) / D(N)
    fNN = D(fN) / D(N)
    val.update({f"{nom}_0": f, f"{nom}_1": D(f), f"{nom}N_0": fN, f"{nom}N_1": D(fN), f"{nom}NN_0": fNN})
a3p = 2 * v * D(Q * D(r)) / r
a3N = a3p / D(N)
val.update({"a3N_0": a3N, "a3N_1": D(a3N), "a3NN_0": D(a3N) / D(N)})
sg = pickle.load(open(os.path.join(AQUI, "resultados", "simbolo_general.pkl"), "rb"))
c = sp.Symbol("c")
Pp = sp.Poly(sg["P"], c)
res = {}
for j in (2, 1, 0):
    co = Pp.coeff_monomial(c ** j)
    e = co.xreplace({sp.Symbol(k): vv for k, vv in val.items()})
    e = sp.together(reduce_beta(e))
    num, den = sp.fraction(e)
    num = sp.factor(sp.expand(reduce_beta(num)))
    res[j] = (num, sp.factor(den))
    print(f"P{j}: numerador = {num}", flush=True)
    print(f"     denominador = {res[j][1]}", flush=True)
Dd = sp.factor(sg["D"].xreplace({sp.Symbol(k): vv for k, vv in val.items()}))
print("D =", Dd, flush=True)
pickle.dump(dict(P=res, D=Dd), open(os.path.join(AQUI, "resultados", "simbolo_perfiles.pkl"), "wb"))
