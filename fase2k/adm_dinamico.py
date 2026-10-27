# -*- coding: utf-8 -*-
"""
Fase 2k: accion cuadratica l = 0 dependiente del tiempo en las variables ADM del corte phi = cte, para un fondo sin
desplazamiento (v = 0, es decir phi = qt t con metrica estatica), en el gauge h_rhorho = 1:
    L = Nc r^2 [A2 + A3 K + A4 (4 K_r K_th + 2 K_th^2) + B4 R3],
    K_r = -v'/Nc,  K_th = (r_t - v r')/(r Nc),  R3 = (2 - 4 r r'' - 2 r'^2)/r^2.
Perturbaciones Nc = N + e nu, v = e w, r = R + e d (funciones de t y rho).  Se comprueba que, salvo derivadas totales,
    L2 = L_est,2(nu, d) + a3N [nu (2 R d_t - (R^2 w)') + 2 R N' d w] - R (Q R')' w^2 + 2 R (Q d_t)' w + Q d_t^2,     Q = 2 A4/N,
donde L_est,2 es la parte cuadratica de la lagrangiana estatica (adm_estatico.py) con w = 0.  Con a3N = 0 (obligatorio
en el fondo: E_v = R^2 a3N N' = 0 y N no constante) w es auxiliar y, eliminandolo,
    T[d_t] = int [ Q d_t^2 + R ((Q d_t)')^2/(Q R')' ],
sin terminos giroscopicos.  En la garganta (R' = 0, R'' > 0): Q y (Q R')' = Q R'' tienen el signo de A4(0).
Uso: python adm_dinamico.py
"""
import os
import json
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
t, rho, e = sp.symbols("t rho e", real=True)
G = lambda n: sp.Function(n)(rho)                       # noqa: E731
H = lambda n: sp.Function(n)(t, rho)                    # noqa: E731
N, R = G("N"), G("R")
a2, a2N, a2NN, a3, a3N, a3NN, a4, a4N, a4NN, b4, b4N, b4NN = (G(n) for n in
    ("a2", "a2N", "a2NN", "a3", "a3N", "a3NN", "a4", "a4N", "a4NN", "b4", "b4N", "b4NN"))
nu, w, d = H("nu"), H("w"), H("d")
Dr = lambda x, k=1: sp.diff(x, rho, k)                  # noqa: E731
Dt = lambda x: sp.diff(x, t)                            # noqa: E731


def taylor(f0, f1, f2, x):
    return f0 + f1 * x + f2 * x ** 2 / 2


Nc = N + e * nu
v = e * w
r = R + e * d
A2, A3, A4, B4 = (taylor(*fs, e * nu) for fs in ((a2, a2N, a2NN), (a3, a3N, a3NN), (a4, a4N, a4NN), (b4, b4N, b4NN)))
Kr = -Dr(v) / Nc
Kth = (Dt(r) - v * Dr(r)) / (r * Nc)
R3 = (2 - 4 * r * Dr(r, 2) - 2 * Dr(r) ** 2) / r ** 2
L = Nc * r ** 2 * (A2 + A3 * (Kr + 2 * Kth) + A4 * (4 * Kr * Kth + 2 * Kth ** 2) + B4 * R3)
L2 = sp.expand(sp.diff(L, e, 2).subs(e, 0) / 2)
# parte estatica: misma lagrangiana con w = 0 y sin derivadas temporales
Lest = Nc * r ** 2 * A2 + Nc * B4 * (2 - 4 * r * Dr(r, 2) - 2 * Dr(r) ** 2)
Lest2 = sp.expand(sp.diff(Lest, e, 2).subs(e, 0) / 2)
Q = 2 * a4 / N
pred = (Lest2 + a3N * nu * (2 * R * Dt(d) - Dr(R ** 2 * w)) + 2 * R * a3N * Dr(N) * d * w - R * Dr(Q * Dr(R)) * w ** 2
        + 2 * R * Dr(Q * Dt(d)) * w + Q * Dt(d) ** 2)
dif = sp.expand(L2 - pred)


def EL(expr, f):
    return sp.simplify(sp.euler_equations(expr, [f], [t, rho])[0].lhs)


# las funciones de N a lo largo del fondo: a_i' = a_iN N' (regla de la cadena)
cadena = {}
for f0, f1, f2 in ((a2, a2N, a2NN), (a3, a3N, a3NN), (a4, a4N, a4NN), (b4, b4N, b4NN)):
    cadena[sp.Derivative(f0, rho)] = f1 * Dr(N)
    cadena[sp.Derivative(f0, (rho, 2))] = Dr(f1 * Dr(N))
els = {}
for f in (nu, w, d):
    if not dif.has(f):
        els[str(f.func)] = 0
        continue
    x = sp.euler_equations(dif, [f], [t, rho])[0].lhs
    x = sp.expand(x.doit().subs(cadena).doit())
    els[str(f.func)] = sp.simplify(x)
print("Euler-Lagrange de L2 - prediccion:", els, flush=True)
assert all(x == 0 for x in els.values())
# fondo: E_v = R^2 a3' (con v = 0) -> a3N N' = 0
Ev = sp.euler_equations(L, [w], [t, rho])[0].lhs
Ev0 = sp.simplify(sp.diff(Ev, e).subs(e, 0).subs({nu: 0, w: 0, d: 0}).doit())
print("ligadura de momento del fondo (orden e):", Ev0, flush=True)
json.dump({"EL_diferencia": {k: str(x) for k, x in els.items()}, "Ev_fondo": str(Ev0)},
          open(os.path.join(AQUI, "resultados", "adm_dinamico.json"), "w"), indent=1)
