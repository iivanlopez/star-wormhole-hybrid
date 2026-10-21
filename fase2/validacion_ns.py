# -*- coding: utf-8 -*-
"""
(e) Teorema del punto de retorno: para estrellas ordinarias omega_0^2 cambia de signo
    exactamente en el maximo de M(rho_c).  Valida la parte de fluido + gravedad.
(f) Limite de garganta diminuta del Modelo L (ligado): omega_0^2 -> omega_0^2 de la estrella ordinaria.
"""
import sys
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, minimize_scalar
from hibrido import Politropo, Fondo, autovalores, RHO_GEOM, MSUN_KM
from estrella_ns import EstrellaNS

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

eos = Politropo(k=0.1, gamma=2.0)
rhos = np.array([3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0]) * 1e14
Ms, ws = [], []
for rho in rhos:
    s = EstrellaNS(eos, rho * RHO_GEOM)
    Ms.append(s.M / MSUN_KM)
    ws.append(s.modo_fundamental(-2e-3, 3e-3, 60))
    print(f"  rho_c={rho:.2e}  M={Ms[-1]:.5f} Msun  omega0^2={ws[-1]:.4e} km^-2")
x = np.log(rhos)
Msp, wsp = CubicSpline(x, Ms), CubicSpline(x, ws)
x_max = minimize_scalar(lambda t: -Msp(t), bounds=(x[0], x[-1]), method="bounded").x
x_0 = brentq(lambda t: wsp(t), x[0], x[-1])
print(f"(e) rho_c(M_max) = {np.exp(x_max):.4e}   rho_c(omega0^2=0) = {np.exp(x_0):.4e}"
      f"   diferencia relativa = {abs(np.exp(x_0) / np.exp(x_max) - 1):.1e}")

print("\n(f) Modelo L con garganta diminuta (xi0 = 1e-3) frente a la estrella ordinaria")
for rho in (1e14, 4e14):
    rb = rho * RHO_GEOM
    w_ns = EstrellaNS(eos, rb).modo_fundamental(-2e-3, 3e-3, 60)
    f = Fondo(r0=1e-3 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    r, _ = autovalores(f, "ligado", np.linspace(-2e-3, 3e-3, 21))
    print(f"  rho={rho:.0e}: NS omega0^2 = {w_ns:.6e}   Modelo L omega0^2 = {f.w_fisico(r[0]):.6e}")
