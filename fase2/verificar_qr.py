# -*- coding: utf-8 -*-
"""Comprobacion independiente del signo del determinante: Plucker (det) vs ortonormalizacion QR (det_gs)."""
import sys
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from hibrido import Politropo, Fondo, RHO_GEOM

eos = Politropo(0.1, 2.0)


def tarea(args):
    xi0, rho, wfis = args
    rb = rho * RHO_GEOM
    f = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    k = np.exp(2 * f.Phi_inf)
    return xi0, rho, [(w, np.sign(f.det("ligado", w * k)), np.sign(f.det_gs("ligado", w * k))) for w in wfis]


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    casos = [(0.6, 3e13, [-0.52, -0.49, -0.47, -0.45]), (0.1, 1e14, [-20.0, -19.0, -18.5, -17.5]),
             (0.1, 1e14, [1.5e-4, 2.0e-4, 2.1e-4, 2.6e-4])]
    t0 = time.time()
    with ProcessPoolExecutor(3) as ex:
        for xi0, rho, filas in ex.map(tarea, casos):
            print(f"xi0={xi0} rho={rho:.0e}:  " + "  ".join(
                f"w={w:+.5g}: Plucker {int(a):+d} / QR {int(b):+d}" for w, a, b in filas), flush=True)
    print("t =", round(time.time() - t0), "s")
