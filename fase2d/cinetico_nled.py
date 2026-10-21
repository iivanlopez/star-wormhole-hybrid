# -*- coding: utf-8 -*-
"""
Coeficiente cinetico K de las perturbaciones radiales de agujeros de gusano ULTRAESTATICOS de EsGB + electrodinamica
(Maxwell o NLED) para formas de garganta generales (cinetico_simbolico.py para las formulas).

Dada una forma par r(l) (r0 = 1), la carga escalar D y la asimetria c = dF/dl(0), el fondo queda fijado:
  B = 1 - r'^2,   P = r''/r + D^2/(4 r^4),   eps = (B + D^2/(4 r^2))/r^2 > 0     (NLED reconstruida)
  h = dF/dl = (c - (1/2) int_0^l r^2 P) / B,   F' = h r^2/D,   phi' = D/r^2
  K = A^2 r^2 - 16 F'^2 B (4 A r'' + B C),  G = A^2 r^2 - 16 F'^2 B^2,  A = r - 4 h r',  C = 1 - 4 h'
Familias de formas:
  hiperbola     r = 1 + sqrt(l^2 + L^2) - L        (L = 1: Ellis; L > 1: garganta "gorda", r''(0) = 1/L)
  superelipse   r = (1 + |l|^{2n})^{1/(2n)}        (n = 1: Ellis; n >= 2: garganta plana, r''(0) = 0)
Con Maxwell la forma es necesariamente Ellis (r'^2 = 1 - (Qe^2 - D^2/4)/r^2); el resto exige NLED.
Salida: resultados/cinetico_nled.json.   Uso: python cinetico_nled.py
"""
import os
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.integrate import cumulative_trapezoid

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def forma(tipo, par, l):
    """r, r', r'', B = 1 - r'^2 para l >= 0 (B analitico, sin cancelaciones)."""
    if tipo == "hiperbola":
        L = par
        s = np.sqrt(l * l + L * L)
        return 1.0 + s - L, l / s, L * L / s**3, (L / s) ** 2
    if tipo == "superelipse":
        n = par
        lu = np.log1p(l ** (2 * n))                         # log u, u = 1 + l^{2n}
        r = np.exp(lu / (2 * n))
        rp = l ** (2 * n - 1) * np.exp((1.0 / (2 * n) - 1.0) * lu)
        rpp = (2 * n - 1) * l ** (2 * n - 2) * np.exp((1.0 / (2 * n) - 2.0) * lu)
        # r'^2 = (1 - 1/u)^{(2n-1)/n}  ->  B = -expm1(((2n-1)/n) log1p(-1/u))
        B = -np.expm1((2 * n - 1) / n * np.log1p(-np.exp(-lu)))
        B[0] = 1.0
        return r, rp, rpp, B
    raise ValueError(tipo)


def malla(lmax, n=40001):
    x = np.linspace(0.0, 1.0, n)
    return lmax * np.sinh(8 * x) / np.sinh(8.0)        # densa cerca de la garganta


def perfil(tipo, par, D, c, lmax):
    lp = malla(lmax)
    if callable(tipo):
        r, rp, rpp, B = tipo(par, lp)
    else:
        r, rp, rpp, B = forma(tipo, par, lp)
    r2P = r * rpp + D**2 / (4 * r**2)
    I = cumulative_trapezoid(r2P, lp, initial=0.0)
    res = {}
    for lado, s in (("der", 1.0), ("izq", -1.0)):
        # en l < 0 (s = -1): r par, r' impar, r'' par, I impar
        rps = s * rp
        Is = s * I
        h = (c - Is / 2) / B
        hp = (-r2P / 2 + 2 * rps * rpp * h) / B
        A = r - 4 * h * rps
        C = 1 - 4 * hp
        F1 = h * r**2 / D
        E = 4 * A * rpp + B * C
        K = A**2 * r**2 - 16 * F1**2 * B * E
        G = A**2 * r**2 - 16 * F1**2 * B**2
        nrm = A**2 * r**2 + 16 * F1**2 * B * np.abs(E)
        res[lado] = dict(l=s * lp, K=K, G=G, k=K / nrm, eps=(B + D**2 / (4 * r**2)) / r**2, P=r2P / r**2)
    return res


def resumen(tipo, par, D, c, lmax):
    pr = perfil(tipo, par, D, c, lmax)
    kmin, lkmin, cambios_K, lam_neg = np.inf, None, 0, False
    for lado in ("der", "izq"):
        d = pr[lado]
        i = int(np.argmin(d["k"]))
        if d["k"][i] < kmin:
            kmin, lkmin = float(d["k"][i]), float(d["l"][i])
        cambios_K += int(np.sum(np.diff(np.sign(d["K"])) != 0))
        lam = d["G"] / d["K"]
        lam_neg = lam_neg or bool(np.any(lam < 0))
    return dict(tipo=tipo, par=par, D=D, c=c, kmin=kmin, l_kmin=lkmin, cambios_K=cambios_K, lambda_negativo=lam_neg)


def tarea(args):
    return resumen(*args)


if __name__ == "__main__":
    tareas = []
    for c in (0.0, 0.1, 0.3, 1.0):
        for L in (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0):
            for D in np.logspace(-3, 2, 26):
                tareas.append(("hiperbola", L, float(D), c, 400.0 * max(1.0, L)))
        for n in (1, 2, 3, 5, 10):
            for D in np.logspace(-3, 2, 26):
                tareas.append(("superelipse", n, float(D), c, 400.0))
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(tarea, tareas, chunksize=8))
    for c in (0.0, 0.1, 0.3, 1.0):
        print(f"c = {c:g}:")
        for fam in ("hiperbola", "superelipse"):
            for par in sorted({o["par"] for o in res if o["tipo"] == fam}):
                sub = [o for o in res if o["tipo"] == fam and o["par"] == par and o["c"] == c]
                pos = [o for o in sub if o["kmin"] > 0]
                peor = max(sub, key=lambda o: o["kmin"])
                print(f"  {fam:11s} {par:7g}: K > 0 en toda la linea en {len(pos):2d}/{len(sub)} valores de D; "
                      f"max_D min_l K/norma = {peor['kmin']:+.3e} (D = {peor['D']:.3g}, l = {peor['l_kmin']:.3g}); "
                      f"lambda < 0 en {sum(o['lambda_negativo'] for o in sub)}/{len(sub)}")
    print(f"total: {len(res)} configuraciones; con K > 0 en toda la linea: {sum(o['kmin'] > 0 for o in res)}; "
          f"sin lambda < 0: {sum(not o['lambda_negativo'] for o in res)}")
    json.dump(res, open("resultados/cinetico_nled.json", "w", encoding="utf-8"), indent=1)
