# -*- coding: utf-8 -*-
"""Kodama puro: autovalor del Modelo L (raiz de u_inf(w), lineal en w) vs fantasma dinamico,
en funcion de la tolerancia con la que se congela la cola del campo (conv_tol)."""
import sys
import math
import numpy as np
from scipy.optimize import brentq
from hibrido import Fondo, autovalores
from validacion_estabilidad import ajustar_cuartico

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
f0 = ajustar_cuartico()
for ct in (1e-6, 1e-7, 1e-8, 1e-9):
    f = Fondo(r0=f0.r0, V0=f0.V0, Lam=1.0).resolver(l_max=2.0e3, conv_tol=ct)
    u = lambda w: f.condiciones("ligado", w)[1][0]
    wl = brentq(u, -0.6, -0.05, xtol=1e-14, rtol=1e-12)
    rg, _ = autovalores(f, "fantasma", np.array([-0.6, -0.05]))
    print(f"conv_tol={ct:g}:  Modelo L omega^2 = {f.w_fisico(wl):.9f}   fantasma omega^2 = {f.w_fisico(rg[0]):.9f}"
          f"   l_conv = {f.l_conv:.2f} km", flush=True)
