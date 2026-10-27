# -*- coding: utf-8 -*-
"""
Fase 2k: ?existe un fondo phi = qt t con G4 > 0?  El fondo exige que beta = N G4 sea un modo cero par y acotado de
    r^3 (beta'/r)' + (1 + r r'' - r'^2) beta = 0                      (reconstruccion_qt.py)
y beta > 0 exige (integrando) que 1 + r r'' - r'^2 < 0 en alguna zona (r' > 1: masa de Misner-Sharp negativa).
Barrido de geometrias con esa zona:
    r = sqrt(rho^2 + r0^2) + A [1 + tanh((rho - rho1)/w)]/2 + A [1 + tanh((-rho - rho1)/w)]/2      (par)
Para cada forma: numero de ceros de la beta par (beta(0) = 1, beta'(0) = 0) y coeficiente de rho^2 lejos.  Un modo
cero sin nodos aparece donde el coeficiente cambia de signo con 0 ceros.
Uso: python modo_cero_beta.py
"""
import os
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

AQUI = os.path.dirname(os.path.abspath(__file__))
L = 400.0
r0 = 1.0


def geom(A, rho1, w):
    def rr(x):
        s = np.sqrt(x ** 2 + r0 ** 2)
        s1, s2 = x / s, r0 ** 2 / s ** 3
        u1, u2 = (x - rho1) / w, (-x - rho1) / w
        t1, t2 = np.tanh(u1), np.tanh(u2)
        c1, c2 = 1 - t1 ** 2, 1 - t2 ** 2
        R = s + A * (2 + t1 + t2) / 2
        R1 = s1 + A * (c1 - c2) / (2 * w)
        R2 = s2 + A * (-2 * t1 * c1 - 2 * t2 * c2) / (2 * w ** 2)
        return R, R1, R2
    return rr


def beta(rr):
    def f(x, y):
        R, R1, R2 = rr(x)
        return [y[1], (R * R1 * y[1] - (1 + R * R2 - R1 ** 2) * y[0]) / R ** 2]
    return solve_ivp(f, (0, L), [1.0, 0.0], dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.5)


def analiza(A, rho1, w):
    rr = geom(A, rho1, w)
    s = beta(rr)
    x = np.linspace(0, L, 20001)
    b = s.sol(x)[0]
    ceros = int(np.sum(np.diff(np.sign(b)) != 0))
    xa = np.linspace(L / 2, L, 300)
    coef = np.linalg.lstsq(np.vstack([xa ** 2, np.ones_like(xa), 1 / xa]).T, s.sol(xa)[0], rcond=None)[0][0]
    rmin = rr(x)[0].min()
    return ceros, coef, rmin


if __name__ == "__main__":
    RES = []
    for rho1 in (1.0, 2.0, 4.0):
        for w in (0.3, 0.7, 1.5):
            As = np.linspace(-0.5, 6.0, 66)
            filas = [analiza(A, rho1, w) for A in As]
            for i in range(len(As) - 1):
                (z0, c0, m0), (z1, c1, m1) = filas[i], filas[i + 1]
                if c0 * c1 < 0 and m0 > 0 and m1 > 0:
                    Ar = brentq(lambda A: analiza(A, rho1, w)[1], As[i], As[i + 1], xtol=1e-12)
                    z, c, m = analiza(Ar, rho1, w)
                    # ceros justo antes del cruce (el modo cero es el limite)
                    RES.append(dict(rho1=rho1, w=w, A=float(Ar), ceros_beta=min(z0, z1), r_min=float(m)))
                    print(RES[-1], flush=True)
            print(f"rho1 = {rho1}, w = {w}: ceros de beta en A = {As[0]:.2f}..{As[-1]:.2f}: "
                  f"{[f[0] for f in filas[::5]]}", flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", "modo_cero_beta.json"), "w"), indent=1)
