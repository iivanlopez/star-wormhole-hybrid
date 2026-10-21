# -*- coding: utf-8 -*-
"""
Condicion en el infinito del Modelo L: andamio libre (u_inf = 0, deslizamiento C libre) frente a
andamio sujeto (C = 0, con fuerza u_inf != 0 en el infinito).
Criterio fisico: con una garganta diminuta el hibrido debe reducirse a la estrella ordinaria
(modo fundamental de la EN a 1e14: 2.080022e-4 km^-2).
Salida: resultados/condicion_infinito.json.   Uso: python condicion_infinito.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq
from hibrido import Politropo, Fondo, RHO_GEOM

eos = Politropo(k=0.1, gamma=2.0)
CASOS = [(x, 1e14, p, inf) for x in (0.001, 0.1) for p in ("par", "impar") for inf in ("libre", "sujeto")]
CASOS += [(0.6, 3e13, "par", "sujeto")]


def tarea(args):
    xi0, rho, paridad, infinito = args
    rb = rho * RHO_GEOM
    f = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    f.paridad = paridad
    k = math.exp(2 * f.Phi_inf)
    ws = np.sort(np.concatenate([-np.logspace(np.log10(300.0), -6.0, 50) * k, np.linspace(2e-6, 6e-4, 25) * k]))
    D = [f.det("ligado", w, infinito=infinito) for w in ws]
    r = [f.w_fisico(brentq(lambda w: f.det("ligado", w, infinito=infinito), ws[i], ws[i + 1], rtol=1e-10))
         for i in range(len(ws) - 1) if D[i] * D[i + 1] < 0]
    return dict(xi0=xi0, rho=rho, paridad=paridad, infinito=infinito, r0=f.r0, raices=r)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    with ProcessPoolExecutor(max_workers=3) as ex:
        res = list(ex.map(tarea, CASOS))
    for o in res:
        print(f"{o['paridad']:5s} {o['infinito']:6s} xi0={o['xi0']:<5} rho={o['rho']:.0e}: "
              f"omega^2 en [-300, 6e-4] km^-2 = {[f'{x:.5e}' for x in o['raices']]}")
    json.dump(res, open("resultados/condicion_infinito.json", "w", encoding="utf-8"), indent=1)
