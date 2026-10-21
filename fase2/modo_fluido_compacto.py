# -*- coding: utf-8 -*-
"""
Configuracion mas compacta (xi0 = 0.6, rho_bc = 1e15): el barrido de test_h1.py (w interno hasta 2e-3) no
encuentra el modo fundamental del fluido.  Aqui se buscan los modos pares positivos del Modelo L en
omega^2 (infinito) en [1e-4, 0.03] km^-2, junto con xi0 = 0.6 a 4e14 como referencia.
Uso: python modo_fluido_compacto.py
"""
import sys
import math
import numpy as np
from scipy.optimize import brentq
from hibrido import Politropo, Fondo, RHO_GEOM

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
eos = Politropo(k=0.1, gamma=2.0)
for xi0, rho in ((0.6, 1e15), (0.6, 4e14)):
    rb = rho * RHO_GEOM
    f = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    k = math.exp(2 * f.Phi_inf)
    ws = np.linspace(1e-4, 0.03, 40) * k
    D = [f.det("ligado", w) for w in ws]
    r = [f.w_fisico(brentq(lambda w: f.det("ligado", w), ws[i], ws[i + 1], rtol=1e-10))
         for i in range(len(ws) - 1) if D[i] * D[i + 1] < 0]
    print(f"xi0={xi0} rho={rho:.0e}: modos pares con omega^2 en [1e-4, 0.03] km^-2: {[f'{x:.4e}' for x in r]}")
