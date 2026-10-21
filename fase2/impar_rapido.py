# -*- coding: utf-8 -*-
"""Sector impar del Modelo L: raices reales en [-300, -1e-3] (fisico) para algunas configuraciones."""
import sys
import math
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.optimize import brentq
from hibrido import Politropo, Fondo, RHO_GEOM

eos = Politropo(k=0.1, gamma=2.0)


def tarea(xi0, rho):
    rb = rho * RHO_GEOM
    f = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    f.paridad = "impar"
    k = math.exp(2 * f.Phi_inf)
    ws = np.sort(-np.logspace(np.log10(300.0), -3.0, 41) * k)
    D = [f.det("ligado", w) for w in ws]
    r = [f.w_fisico(brentq(lambda w: f.det("ligado", w), ws[i], ws[i + 1], rtol=1e-10))
         for i in range(len(ws) - 1) if D[i] * D[i + 1] < 0]
    return f"impar xi0={xi0} rho={rho:.0e}: omega^2 < 0 en [-300, -1e-3]: {r}"


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    with ProcessPoolExecutor(max_workers=3) as ex:
        fut = [ex.submit(tarea, *c) for c in ((0.1, 1e14), (0.6, 3e13), (0.6, 1e15))]
        for fu in as_completed(fut):
            print(fu.result(), flush=True)
