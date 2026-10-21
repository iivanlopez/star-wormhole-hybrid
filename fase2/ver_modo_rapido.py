# -*- coding: utf-8 -*-
"""Autofuncion del modo rapido del Modelo L (xi0 = 0.6, rho = 3e13): comprueba el mecanismo
(andamio que desliza C, fluido en reposo xi ~ -C r', garganta dr(0))."""
import sys
import math
import numpy as np
from scipy.optimize import brentq
from hibrido import Politropo, Fondo, RHO_GEOM
from autofuncion import autofuncion_gs

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
eos = Politropo(0.1, 2.0)
rb = 3e13 * RHO_GEOM
f = Fondo(r0=0.6 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
k = math.exp(2 * f.Phi_inf)
w = brentq(lambda x: f.det("ligado", x), -0.55 * k, -0.42 * k, xtol=1e-14, rtol=1e-12)
print(f"omega^2 = {f.w_fisico(w):.8f} km^-2")
l, X, sv = autofuncion_gs(f, "ligado", w)
print("valores singulares de las condiciones:", sv)
xi, DP, dr = X[0], X[1], X[2]
s = np.array([f.y_bg(x)[1] for x in l])
Phi = np.array([f.y_bg(x)[2] for x in l])
C = dr[-1] * math.exp(-Phi[-1])                      # deslizamiento del andamio en el infinito
dr0 = dr[0]
print(f"dr(0)/C = {dr0 / C:.4f}")
for x0 in (0.0, 0.5 * f.r0, f.r0, 3 * f.r0, 0.5 * f.ls, 0.9 * f.ls, f.ls, 100.0, 1000.0):
    i = np.argmin(abs(l - x0))
    print(f"  l={l[i]:9.3f}  xi/C={xi[i] / C:+.4f}  -s={-s[i]:+.4f}  dr/(C e^Phi)={dr[i] / (C * math.exp(Phi[i])):+.4f}"
          f"  (dr + s xi)/C={(dr[i] + s[i] * xi[i]) / C:+.4f}")
np.savez("resultados/modo_rapido_autofuncion.npz", l=l, X=X, s=s, Phi=Phi, C=C, w=f.w_fisico(w), ls=f.ls, r0=f.r0)
