# -*- coding: utf-8 -*-
"""
Busqueda de un contraejemplo: agujero ultraestatico EsGB + NLED con coeficiente cinetico K > 0 en toda la linea.
Familia flexible de formas pares (r0 = 1), asintoticamente planas con cola de tipo Maxwell:
    B(l) = 1 - r'^2 = 1/(1 + l^2 a(l)),   a(l) = exp( sum_k theta_k exp(-(log(1 + l^2) - t_k)^2/(2 w^2)) ),
con 6 nodos t_k en [0, log(1 + 400^2)] y w = 0.5 (theta = 0: Ellis).  r' puede no ser monotono (r'' < 0 permitido).
Parametros libres: theta_k in [-4, 4], log10 D in [-3, 2], c = dF/dl(0) in [-2, 2].
Objetivo: maximizar min_l K/(A^2 r^2 + 16 F'^2 B |E|) en |l| <= 400 (evolucion diferencial, 8 procesos), varias semillas.
Salida: resultados/busqueda_cinetico.json.   Uso: python busqueda_cinetico.py
"""
import os
import sys
import json
import time
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import differential_evolution

LMAX = 400.0
NK = 6
TK = np.linspace(0.0, np.log1p(LMAX**2), NK)
W = 0.5
_x = np.linspace(0.0, 1.0, 20001)
LP = LMAX * np.sinh(8 * _x) / np.sinh(8.0)


def geometria(theta):
    l = LP
    t = np.log1p(l * l)
    bump = np.exp(-((t[:, None] - TK[None, :]) ** 2) / (2 * W * W))
    loga = bump @ theta
    a = np.exp(loga)
    dloga_dt = (bump * (-(t[:, None] - TK[None, :]) / W**2)) @ theta
    dt_dl = 2 * l / (1 + l * l)
    q = l * l * a                                  # B = 1/(1+q)
    dq = 2 * l * a + l * l * a * dloga_dt * dt_dl
    B = 1.0 / (1.0 + q)
    rp = np.sqrt(q / (1.0 + q))
    rpp = np.empty_like(l)
    rpp[1:] = dq[1:] / (2 * rp[1:] * (1 + q[1:]) ** 2)
    rpp[0] = np.sqrt(a[0])
    r = 1.0 + cumulative_trapezoid(rp, l, initial=0.0)
    return r, rp, rpp, B


def kmin(x):
    theta, D, c = np.array(x[:NK]), 10.0 ** x[NK], x[NK + 1]
    r, rp, rpp, B = geometria(theta)
    r2P = r * rpp + D**2 / (4 * r**2)
    I = cumulative_trapezoid(r2P, LP, initial=0.0)
    peor = np.inf
    for s in (1.0, -1.0):
        rps, Is = s * rp, s * I
        h = (c - Is / 2) / B
        hp = (-r2P / 2 + 2 * rps * rpp * h) / B
        A = r - 4 * h * rps
        C = 1 - 4 * hp
        F1 = h * r**2 / D
        E = 4 * A * rpp + B * C
        K = A**2 * r**2 - 16 * F1**2 * B * E
        nrm = A**2 * r**2 + 16 * F1**2 * B * np.abs(E)
        peor = min(peor, float(np.min(K / nrm)))
    return peor


def objetivo(x):
    v = kmin(x)
    return -v if np.isfinite(v) else 2.0


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    t0 = time.time()
    limites = [(-4.0, 4.0)] * NK + [(-3.0, 2.0), (-2.0, 2.0)]
    ellis = kmin([0.0] * NK + [np.log10(2.0), 0.0])
    print(f"control: Ellis con D = 2, c = 0 (Canate Qs = Q): min K/norma = {ellis:+.4f}")
    res = []
    for semilla in range(4):
        r = differential_evolution(objetivo, limites, seed=semilla, popsize=25, maxiter=150, tol=1e-8, polish=True,
                                   workers=8, updating="deferred")
        res.append(dict(semilla=semilla, max_min_K=float(-r.fun), x=[float(v) for v in r.x], nfev=int(r.nfev)))
        print(f"semilla {semilla}: max_(forma, D, c) min_l K/norma = {-r.fun:+.5f}  (D = {10**r.x[NK]:.3g}, "
              f"c = {r.x[NK+1]:+.3f}, theta = {np.array2string(r.x[:NK], precision=2)}; {r.nfev} evaluaciones)", flush=True)
    json.dump(dict(control_ellis=ellis, busquedas=res), open("resultados/busqueda_cinetico.json", "w", encoding="utf-8"),
              indent=1)
    print(f"tiempo: {time.time() - t0:.0f} s")
