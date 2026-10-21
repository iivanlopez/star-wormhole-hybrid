# -*- coding: utf-8 -*-
"""
Agujeros ultraestaticos de EsGB con electrodinamica no lineal concreta (los modelos de Canate, arXiv:2310.08758):
  Born-Infeld (campo magnetico, escala b):  eps(r) = 2 b^2 (sqrt(1 + q^2/(b^2 r^4)) - 1)      (b -> inf: Maxwell)
  Euler-Heisenberg (debil):                 eps(r) = q^2/r^4 + eta q^4/r^8
con 8 pi T = diag(-eps, -eps, P, P), P = -eps - r eps'/2 (conservacion).  En unidades r0 = 1, la garganta exige
eps(1) = 1 + D^2/4 (ligadura ll), que fija q^2; la forma sale de integrar r'' = r (P - D^2/(4 r^4)) desde r(0) = 1,
r'(0) = 0 (necesita r''(0) > 0).  Con la forma, h = dF/dl, K y G como en cinetico_nled.py (c = dF/dl(0)).
Salida: resultados/nled_concretas.json.   Uso: python nled_concretas.py
"""
import os
import sys
import json
import math
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import brentq

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def modelo(tipo, par, D):
    """devuelve eps(r), P(r) con q^2 fijado por la garganta (None si no hay solucion)."""
    E1 = 1 + D * D / 4
    if tipo == "BI":
        b = par
        q2 = E1 + E1 * E1 / (4 * b * b)

        def eps(r):                                # 2 b^2 (sqrt(1+x) - 1) sin cancelacion
            x = q2 / (b * b * r**4)
            return 2 * q2 / (r**4 * (math.sqrt(1 + x) + 1))

        def P(r):
            return -eps(r) + 2 * q2 / (r**4 * math.sqrt(1 + q2 / (b * b * r**4)))
        return eps, P
    if tipo == "EH":
        eta = par
        if eta == 0:
            q2 = E1
        else:
            disc = 1 + 4 * eta * E1
            if disc < 0:
                return None
            q2 = (-1 + math.sqrt(disc)) / (2 * eta)
            if q2 <= 0:
                return None

        def eps(r):
            return q2 / r**4 + eta * q2 * q2 / r**8

        def P(r):
            return q2 / r**4 + 3 * eta * q2 * q2 / r**8
        return eps, P
    raise ValueError(tipo)


def analiza(tipo, par, D, c=0.0, lmax=400.0):
    m = modelo(tipo, par, D)
    if m is None:
        return None
    eps, P = m
    if P(1.0) - D * D / 4 <= 0:
        return dict(tipo=tipo, par=par, D=D, c=c, garganta=False)
    sol = solve_ivp(lambda l, y: [y[1], y[0] * (P(y[0]) - D * D / (4 * y[0] ** 4))], (0, lmax), [1.0, 0.0],
                    method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True)
    x = np.linspace(0.0, 1.0, 20001)
    lp = lmax * np.sinh(8 * x) / np.sinh(8.0)
    r, rp = sol.sol(lp)
    rpp = r * (np.array([P(v) for v in r]) - D * D / (4 * r**4))
    B = r * r * np.array([eps(v) for v in r]) - D * D / (4 * r * r)        # ligadura ll (comprobacion: 1 - r'^2)
    err_B = float(np.max(np.abs(B - (1 - rp**2))))
    if np.min(B[1:]) <= 0:
        # r' > 1 en algun punto: B = 0 hace divergir F_l = u/B (fondo singular salvo ajuste fino): se descarta
        return dict(tipo=tipo, par=par, D=D, c=c, garganta=True, regular=False)
    r2P = r * rpp + D * D / (4 * r * r)
    I = cumulative_trapezoid(r2P, lp, initial=0.0)
    kmin, cambios = np.inf, 0
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
        kmin = min(kmin, float(np.min(K / nrm)))
        cambios += int(np.sum(np.diff(np.sign(K)) != 0))
    return dict(tipo=tipo, par=par, D=D, c=c, garganta=True, regular=True, rpp0=float(rpp[0]), rp_fin=float(rp[-1]),
                err_B=err_B, kmin=kmin, cambios_K=cambios)


if __name__ == "__main__":
    res = []
    for b in (0.5, 1.0, 2.0, 5.0, 20.0, 1e3):
        for D in np.logspace(-2, 1.5, 15):
            o = analiza("BI", b, float(D))
            if o:
                res.append(o)
    for eta in (-0.2, -0.1, -0.03, 0.0, 0.03, 0.1, 0.3, 1.0):
        for D in np.logspace(-2, 1.5, 15):
            o = analiza("EH", eta, float(D))
            if o:
                res.append(o)
    for tipo in ("BI", "EH"):
        for par in sorted({o["par"] for o in res if o["tipo"] == tipo}):
            sub = [o for o in res if o["tipo"] == tipo and o["par"] == par]
            gar = [o for o in sub if o["garganta"] and o.get("regular")]
            sing = sum(1 for o in sub if o["garganta"] and not o.get("regular"))
            if not gar:
                print(f"{tipo} {par:g}: ninguna garganta regular ({sum(1 for o in sub if not o['garganta'])} sin garganta, {sing} con B = 0)")
                continue
            print(f"{tipo} {par:6g}: {len(gar)}/{len(sub)} con garganta regular ({sing} con B = 0, descartados); K cambia de signo en "
                  f"{sum(o['cambios_K'] > 0 for o in gar)}/{len(gar)}; max_D min K/norma = {max(o['kmin'] for o in gar):+.4f}; "
                  f"r'(l_max) en [{min(o['rp_fin'] for o in gar):.4f}, {max(o['rp_fin'] for o in gar):.4f}]; "
                  f"|B - (1 - r'^2)| <= {max(o['err_B'] for o in gar):.1e}")
    gar = [o for o in res if o["garganta"] and o.get("regular")]
    print(f"total con garganta regular: {len(gar)}; con K > 0 en toda la linea: {sum(o['kmin'] > 0 for o in gar)}")
    json.dump(res, open("resultados/nled_concretas.json", "w", encoding="utf-8"), indent=1)
