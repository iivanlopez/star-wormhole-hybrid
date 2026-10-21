# -*- coding: utf-8 -*-
"""Convergencia numerica de los autovalores del test H1 (rho_bc = 1e14, xi0 = 0.1)."""
import sys
import numpy as np
from scipy.optimize import brentq
from hibrido import Politropo, Fondo, RHO_GEOM

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
eos = Politropo(k=0.1, gamma=2.0)
rb = 1e14 * RHO_GEOM


def raiz(f, modelo, a, b, **kw):
    """a, b: intervalo en omega^2 fisico (medido en el infinito)."""
    k = np.exp(2.0 * f.Phi_inf)
    return f.w_fisico(brentq(lambda w: f.det(modelo, w, **kw), a * k, b * k, xtol=1e-16, rtol=1e-12))


for th, lmax, rtol in ((1e-7, 2e4, 1e-10), (1e-9, 2e4, 1e-10), (1e-7, 5e4, 1e-10), (1e-7, 2e4, 1e-12)):
    f = Fondo(r0=0.1 * eos.inv_A(rb), eos=eos, rho_bc=rb, theta_stop=th).resolver(l_max=lmax)
    wl = raiz(f, "ligado", 1.0e-4, 3.5e-4, rtol=rtol)
    wg = raiz(f, "fantasma", -3.0, -0.3, rtol=rtol)
    print(f"theta_stop={th:g} l_max={lmax:g} rtol={rtol:g}:  omega0^2 L = {wl:.10e}   fantasma = {wg:.10f}")

# Modo de la garganta del Modelo L (xi0 = 0.6, rho = 3e13): con la cola 1/r restada, la raiz no depende de l_max
rb6 = 3e13 * RHO_GEOM
for lmax in (2e3, 5e3, 2e4, 5e4):
    f = Fondo(r0=0.6 * eos.inv_A(rb6), eos=eos, rho_bc=rb6).resolver(l_max=lmax)
    print(f"xi0=0.6 rho=3e13 l_max={lmax:g}:  modo de la garganta (cola restada) = {raiz(f, 'ligado', -0.55, -0.42):.8f}")
