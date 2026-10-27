# -*- coding: utf-8 -*-
"""
Fase 2k: modo cero beta > 0 con garganta afilada (ver modo_cero_beta.py y el informe, §5.3).
    r' = tanh(rho/a) (1 + B e^{-(rho/b)^2}),   r(0) = r0 = 1      (r'' (0) = (1 + B)/a)
Integra r junto con  r^2 beta'' - r r' beta' + (1 + r r'' - r'^2) beta = 0,  beta(0) = 1, beta'(0) = 0, y busca
cambios de signo del coeficiente de rho^2 de beta lejos con beta sin ceros.
Uso: python modo_cero_beta2.py
"""
import os
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

AQUI = os.path.dirname(os.path.abspath(__file__))
L = 600.0


def sistema(a, B, b):
    def f(x, y):
        r, be, be1 = y
        t = np.tanh(x / a)
        e = np.exp(-(x / b) ** 2)
        g = t * (1 + B * e)
        g1 = (1 - t ** 2) / a * (1 + B * e) + t * B * e * (-2 * x / b ** 2)
        return [g, be1, (r * g * be1 - (1 + r * g1 - g ** 2) * be) / r ** 2]
    return f


def analiza(a, B, b):
    s = solve_ivp(sistema(a, B, b), (0, L), [1.0, 1.0, 0.0], dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.2)
    x = np.linspace(0, L, 30001)
    be = s.sol(x)[1]
    ceros = int(np.sum(np.diff(np.sign(be)) != 0))
    xa = np.linspace(L / 2, L, 300)
    coef = np.linalg.lstsq(np.vstack([xa ** 2, np.ones_like(xa), 1 / xa]).T, s.sol(xa)[1], rcond=None)[0][0]
    return ceros, coef, s


if __name__ == "__main__":
    RES = []
    for a in (0.1, 0.2, 0.4):
        for b in (0.5, 1.0, 2.0):
            Bs = np.linspace(0.0, 4.0, 41)
            filas = [analiza(a, B, b)[:2] for B in Bs]
            print(f"a = {a}, b = {b}: (ceros, signo coef) = {[(z, int(np.sign(c))) for z, c in filas[::4]]}",
                  flush=True)
            for i in range(len(Bs) - 1):
                (z0, c0), (z1, c1) = filas[i], filas[i + 1]
                if c0 * c1 < 0:
                    Br = brentq(lambda B: analiza(a, B, b)[1], Bs[i], Bs[i + 1], xtol=1e-12)
                    z, c, s = analiza(a, Br, b)
                    x = np.linspace(0, L, 30001)
                    r, be, _ = s.sol(x)
                    RES.append(dict(a=a, b=b, B=float(Br), ceros_vecinos=(z0, z1), beta_min=float(be.min()),
                                    beta_inf=float(be[-1]), rpp0=(1 + Br) / a, rp_max=float(np.max(np.gradient(r, x)))))
                    print("   modo cero:", RES[-1], flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", "modo_cero_beta2.json"), "w"), indent=1)
