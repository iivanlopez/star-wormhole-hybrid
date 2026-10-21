# -*- coding: utf-8 -*-
"""Validacion del fondo estatico contra las tablas de DFKK (R2 = arXiv:1203.3615, Tabla 2, EOS2)."""
import sys
from hibrido import Politropo, Fondo, RHO_GEOM
from estrella_ns import EstrellaNS

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

eos = Politropo(k=0.1, gamma=2.0)          # EOS2 de DFKK
print(f"K = {eos.K:.4f} km^2,  1/A = {eos.inv_A(1e14 * RHO_GEOM):.5f} km  (n = 1: independiente de rho_bc)")

# --- (a) sin agujero de gusano: R2 Tabla 2  [rho_bc: (R, R_prop, M)]
ref_ns = {1e13: (35.1484, 35.3940, 0.2777), 1e14: (30.9169, 32.8113, 1.8760),
          6e14: (20.9771, 25.6722, 3.1564), 1e15: (18.0313, 23.2159, 2.9715),
          3e15: (14.1095, 19.7367, 2.2541)}
print("\n(a) Estrella de neutrones sin agujero (R, R_prop en km; M en Msun) -- nuestro / DFKK")
for rho, (R, Rp, M) in ref_ns.items():
    s = EstrellaNS(eos, rho * RHO_GEOM).resumen()
    print(f"  rho={rho:.1e}: R {s['R']:.4f}/{R}  R_prop {s['R_prop']:.4f}/{Rp}  M {s['M']:.4f}/{M}")

# --- (b) hibrido con fantasma sin masa, xi0 = 0.1  [rho_bc: (R_th, R, R_prop, M, M_fl, M_sfint, M_sfext)]
ref_wh = {1e13: (1.1374, 35.0540, 35.2774, 0.2717, 0.2708, -0.3727, -0.0121),
          1e14: (1.1374, 30.8813, 32.7260, 1.8360, 1.8276, -0.3669, -0.0103),
          6e14: (1.1374, 21.1631, 25.7582, 3.0959, 3.0518, -0.3363, -0.0052),
          1e15: (1.1374, 18.3622, 23.4391, 2.9230, 2.8551, -0.3141, -0.0036),
          3e15: (1.1374, 15.1024, 20.6263, 2.2768, 2.1276, -0.2353, -0.0011)}
print("\n(b) Hibrido + fantasma sin masa, xi0 = 0.1 -- nuestro / DFKK")
for rho, ref in ref_wh.items():
    rb = rho * RHO_GEOM
    f = Fondo(r0=0.1 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    s = f.resumen()
    print(f"  rho={rho:.1e}: R_th {s['R_th']:.4f}/{ref[0]}  R {s['R']:.4f}/{ref[1]}  R_prop {s['R_prop']:.4f}/{ref[2]}"
          f"  M {s['M']:.4f}/{ref[3]}  M_fl {s['M_fl']:.4f}/{ref[4]}  M_sfint {s['M_sfint']:.4f}/{ref[5]}"
          f"  M_sfext {s['M_sfext']:.4f}/{ref[6]}   [M_suma {s['M_suma']:.4f}]")
