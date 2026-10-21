# -*- coding: utf-8 -*-
"""
Validacion del solver de estabilidad radial (modelo 'fantasma' = DFKK):
  (c) agujero de gusano puro con fantasma cuartico (Kodama; R3 = arXiv:1302.5217):
      Lambda = 1, phi~'(0) = 1, L = 1 km  ->  DFKK: omega_0^2 = -0.11610638 / L^2
  (d) Ellis-Bronnikov (fantasma sin masa): un unico modo inestable par.
  (e) estrella de neutrones sin agujero: omega_0^2 = 0 en el maximo de M(rho_c).
"""
import sys
import math
import numpy as np
from scipy.integrate import solve_ivp
from hibrido import Fondo, Politropo, autovalores, RHO_GEOM, MSUN_KM, RTOL, ATOL
from estrella_ns import EstrellaNS

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def clasificar_cuartico(V0, Lam=1.0, P8=0.0, eos=None, rho_bc=0.0):
    """'sobre' si phi~ cruza 1/Lam, 'bajo' si psi se anula antes (normalizacion psi0 = 1)."""
    r0 = (0.5 - V0 - P8) ** -0.5
    f = Fondo(r0=r0, eos=eos, rho_bc=rho_bc, V0=V0, Lam=Lam)
    op = dict(method="DOP853", rtol=RTOL, atol=ATOL)
    y, l = f.y0(), 0.0
    if f.fluido:
        sol = solve_ivp(lambda l, y: f.rhs_bg(y, True), (0, 1e4), y,
                        events=[f.ev_superficie()] + f.ev_campo(), **op)
        if sol.t_events[1].size:
            return "sobre", f
        if sol.t_events[2].size:
            return "bajo", f
        y, l = sol.y[:, -1], sol.t[-1]
    sol = solve_ivp(lambda l, y: f.rhs_bg(y, False), (l, l + 400.0), y, events=f.ev_campo(), **op)
    if sol.t_events[0].size:
        return "sobre", f
    if sol.t_events[1].size:
        return "bajo", f
    return "indefinido", f


def ajustar_cuartico(Lam=1.0, eos=None, rho_bc=0.0):
    P8 = 8 * math.pi * eos.sigma(rho_bc) * rho_bc if eos is not None else 0.0
    lo, hi = 1e-6, 0.5 - P8 - 1e-9
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        res, _ = clasificar_cuartico(mid, Lam, P8, eos, rho_bc)
        if res == "sobre":
            lo = mid
        else:
            hi = mid
    V0 = 0.5 * (lo + hi)
    return Fondo(r0=(0.5 - V0 - P8) ** -0.5, eos=eos, rho_bc=rho_bc, V0=V0, Lam=Lam).resolver(l_max=2.0e3)


if __name__ == "__main__":
    print("(c) Agujero de gusano puro con fantasma cuartico (Lambda = 1, L = 1 km)")
    f = ajustar_cuartico()
    print(f"    V~(0) = {f.V0:.10f}   Sigma_c = r0 = {f.r0:.6f} km   M = {f.M / MSUN_KM:.4f} Msun"
          f"   (DFKK tabla R3, rho_bc=1e13: R_th = 2.0075 km, M = 0.3823)")
    ws = np.linspace(-1.0, -0.005, 25)
    raices, D = autovalores(f, "fantasma", ws)
    print(f"    omega^2 (fantasma): {[round(f.w_fisico(w), 8) for w in raices]}   DFKK: -0.11610638")

    print("\n(d) Ellis-Bronnikov (fantasma sin masa), r0 = 1")
    e = Fondo(r0=1.0).resolver(l_max=2.0e3)
    ws = np.linspace(-3.0, -0.02, 30)
    raices, D = autovalores(e, "fantasma", ws)
    print(f"    omega^2 r0^2 (fantasma): {[round(e.w_fisico(w), 8) for w in raices]}")
