# -*- coding: utf-8 -*-
"""
Fase 2k: ?puede el termino cubico (G3 != 0, desplazamiento v != 0) dar beta = Nc G4 > 0 en una garganta?
Ecuacion del fondo (nogo_beta.py):
    (beta'/r)' + K beta/r^3 = -[Q r^4 (v^2/r^2)']'/(4 r^3),     K = 1 + r r'' - r'^2,  Q = 2 A4/Nc,
con A3 fijado despues por la ligadura de momento r A3' = 2 v (Q r')' y p por H = 0.  Es lineal e inhomogenea: dada
la geometria (r, v, Q) la solucion par con beta acotada lejos es unica (sin afinar).  Se toma
    r = sqrt(rho^2 + 1),  v = V0 (rho/w) e^{-rho^2/(2 w^2)},  Q = Q0 (constante; A4 < 0  <=>  Q0 < 0)
y se mira si beta > 0 en todo el eje con beta_inf > 0 (beta es proporcional a Q0 V0^2).
Uso: python beta_cubico.py
"""
import os
import json
import numpy as np
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
L = 400.0


def resolver(V0, w, Q0):
    def geo(x):
        r = np.sqrt(x ** 2 + 1)
        r1, r2 = x / r, 1 / r ** 3
        return r, r1, r2

    def fuente(x):
        # S = [Q r^4 (v^2/r^2)']'/(4 r^3) con Q constante, derivadas exactas de g = v^2/r^2
        e = np.exp(-x ** 2 / w ** 2)
        g = V0 ** 2 * (x / w) ** 2 * e / (x ** 2 + 1)
        # g = V0^2/w^2 * x^2 e /(x^2+1);  derivadas por diferencias finitas complejas no: usar formula numerica fina
        return g

    h = 1e-4

    def S(x):
        G = lambda y: (Q0 * (y ** 2 + 1) ** 2 * (fuente(y + h) - fuente(y - h)) / (2 * h))   # noqa: E731
        return (G(x + h) - G(x - h)) / (2 * h) / (4 * (x ** 2 + 1) ** 1.5)

    def f(x, y):
        b, c = y                     # c = beta'/r
        r, r1, r2 = geo(x)
        K = 1 + r * r2 - r1 ** 2
        return [c * r, -K * b / r ** 3 - S(x)]

    sh = solve_ivp(f if False else (lambda x, y: [y[1] * geo(x)[0], -(1 + geo(x)[0] * geo(x)[2] - geo(x)[1] ** 2) *
                                                  y[0] / geo(x)[0] ** 3]), (0, L), [1.0, 0.0], dense_output=True,
                   rtol=1e-10, atol=1e-12, max_step=0.5)
    sp_ = solve_ivp(f, (0, L), [0.0, 0.0], dense_output=True, rtol=1e-10, atol=1e-12, max_step=0.05)
    xa = np.linspace(L / 2, L, 300)
    Mx = np.vstack([xa ** 2, np.ones_like(xa), 1 / xa]).T
    ah = np.linalg.lstsq(Mx, sh.sol(xa)[0], rcond=None)[0]
    ap = np.linalg.lstsq(Mx, sp_.sol(xa)[0], rcond=None)[0]
    b0 = -ap[0] / ah[0]
    x = np.linspace(0, L, 40001)
    beta = b0 * sh.sol(x)[0] + sp_.sol(x)[0]
    binf = b0 * ah[1] + ap[1]
    return dict(V0=V0, w=w, Q0=Q0, beta0=float(b0), beta_inf=float(binf), beta_min=float(beta.min()),
                positiva=bool(beta.min() > 0 and binf > 0), negativa=bool(beta.max() < 0 and binf < 0))


if __name__ == "__main__":
    RES = []
    for Q0 in (-2.0, 2.0):
        for w in (0.5, 1.0, 2.0, 4.0):
            for V0 in (0.1, 0.3, 0.6, 1.0, 2.0):
                c = resolver(V0, w, Q0)
                RES.append(c)
                print({k: (float(f"{x:.4g}") if isinstance(x, float) else x) for k, x in c.items()}, flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", "beta_cubico.json"), "w"), indent=1)
