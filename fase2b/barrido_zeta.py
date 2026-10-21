# -*- coding: utf-8 -*-
"""
Metodo antiguo (seguimiento de signo); la referencia es trayectorias_robustas.py.
Dependencia del modo inestable con el modelo de capa: zeta escala el termino de gravedad inducida de la capa respecto
de la eleccion de KKK (zeta = 1).  Fondo: alpha = 0.03, phi'(0+) = -0.7 (r0 = 1, phi0 = 0): alpha/r0^2 = 0.028, f0 = 1.03.
Raices reales w < 0 de D(w) en [-1e6, -1e-4] (reloj de la garganta).
Salida: resultados/barrido_zeta.json.   Uso: python barrido_zeta.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq
from kkk import FondoKKK
from estabilidad_kkk import CapaKKK, barrido, vector_garganta, D

ZETAS = [-1.0, -0.5, 0.0, 0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 3.0]


def tarea(zeta):
    f = FondoKKK(alpha=0.03, phi1=-0.7).resolver()
    c = CapaKKK(f, zeta=zeta)
    ws = -np.logspace(6, -4, 121)
    out = barrido(f, c, ws)
    vals = [v for v, _ in out]
    rs = []
    for i in range(len(ws) - 1):
        if vals[i] * vals[i + 1] < 0:
            ref = vector_garganta(f, c, ws[i])[0]
            rs.append(brentq(lambda w: D(f, c, w, u=vector_garganta(f, c, w, ref)[0]), ws[i], ws[i + 1],
                             xtol=1e-9, rtol=1e-10))
    return dict(zeta=zeta, w=rs, winf=[w * math.exp(-2 * f.Phi_inf) for w in rs], Phi_inf=f.Phi_inf,
                res_tt=max(abs(r) for _, r in out))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(tarea, ZETAS))
    for o in res:
        print(f"zeta={o['zeta']:+.2f}: w* = {[f'{x:.6g}' for x in o['w']]}   omega_inf^2 r0^2 = "
              f"{[f'{x:.4g}' for x in o['winf']]}")
    json.dump(res, open("resultados/barrido_zeta.json", "w", encoding="utf-8"), indent=1)
