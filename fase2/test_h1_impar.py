# -*- coding: utf-8 -*-
"""
H1, sector impar (antisimetrico en l -> -l): el fluido oscila a traves de la garganta entre las dos
estrellas Z2.  Solo Modelo L (en el gauge unitario no hay modo de traslacion espurio).
Salida: resultados/h1_impar.json   Uso: python test_h1_impar.py
"""
import sys
import os
import json
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from hibrido import Politropo, Fondo, autovalores, RHO_GEOM

eos = Politropo(k=0.1, gamma=2.0)
GRID = np.concatenate([-np.logspace(0.0, -3.5, 15), np.linspace(-2e-4, 2e-3, 45)])


def analizar(args):
    xi0, rho = args
    rb = rho * RHO_GEOM
    f = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    f.paridad = "impar"
    r, _ = autovalores(f, "ligado", GRID)
    return dict(xi0=xi0, rho=rho, impar=sorted(f.w_fisico(w) for w in r))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    rhos = [1e13, 3e13, 1e14, 3e14, 4e14, 5e14, 6e14, 8e14, 1e15, 2e15, 3e15]
    tareas = [(0.1, r) for r in rhos] + [(x, r) for x in (0.3, 0.6) for r in (1e14, 1e15)]
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=min(8, os.cpu_count() or 2)) as ex:
        res = list(ex.map(analizar, tareas))
    os.makedirs("resultados", exist_ok=True)
    json.dump(res, open("resultados/h1_impar.json", "w", encoding="utf-8"), indent=1)
    for o in res:
        w = o["impar"] + [float("nan")] * 2
        print(f"xi0={o['xi0']:<4} rho={o['rho']:.1e}  omega^2 impar: {w[0]:.5e}  {w[1]:.5e}")
    print(f"tiempo: {time.time() - t0:.0f} s")
