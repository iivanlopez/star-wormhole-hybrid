# -*- coding: utf-8 -*-
"""
Modo rapido de la garganta en el Modelo L: raices reales negativas de D(w) en [-300, -1e-3] km^-2
(determinante corregido por la cola 1/r de u).  Tambien el fantasma dinamico en el mismo rango.
Salida: resultados/modo_rapido.json (omega^2 fisicos).  Uso: python modo_rapido.py
"""
import sys
import os
import json
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.optimize import brentq
from hibrido import Politropo, Fondo, RHO_GEOM, MSUN_KM

eos = Politropo(k=0.1, gamma=2.0)
WF = -np.logspace(np.log10(300.0), -3.0, 41)[::-1][::-1]


def tarea(xi0, rho, modelo):
    rb = rho * RHO_GEOM
    f = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    k = np.exp(2 * f.Phi_inf)
    ws = np.sort(WF * k)
    D = [f.det(modelo, w) for w in ws]
    r = [f.w_fisico(brentq(lambda w: f.det(modelo, w), ws[i], ws[i + 1], rtol=1e-10))
         for i in range(len(ws) - 1) if D[i] * D[i + 1] < 0]
    fl = f.bg(f.y0(), True)
    return dict(xi0=xi0, rho=rho, modelo=modelo, r0=f.r0, M=f.M / MSUN_KM, Mfl=f.Mfl / MSUN_KM,
                eps0=fl[12], P0=fl[13], raices=sorted(r))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    casos = [(x, r, "ligado") for x in (0.1, 0.3, 0.6) for r in (1e13, 3e13, 1e14, 4e14, 1e15)]
    casos += [(0.1, 1e14, "fantasma"), (0.6, 3e13, "fantasma")]
    res = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        fut = [ex.submit(tarea, *c) for c in casos]
        for fu in as_completed(fut):
            o = fu.result()
            res.append(o)
            print(f"{o['modelo']:9s} xi0={o['xi0']:<4} rho={o['rho']:.0e}  r0={o['r0']:.3f}  "
                  f"omega^2 < 0: {[round(x, 6) for x in o['raices']]}", flush=True)
    os.makedirs("resultados", exist_ok=True)
    json.dump(res, open("resultados/modo_rapido.json", "w", encoding="utf-8"), indent=1)
