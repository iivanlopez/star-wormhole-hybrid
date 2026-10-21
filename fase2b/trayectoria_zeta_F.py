# -*- coding: utf-8 -*-
"""
SUSTITUIDO por trayectorias_robustas.py (metodo de menores): este script pierde el modo con zeta <= 0.2.
Trayectoria de los modos inestables reales frente a la intensidad zeta de la gravedad inducida de la capa, para un fondo
hiperbolico de la familia cuadratica de Antoniou et al. (alpha = 0.09, phi0 = 0.5, phi'(0+) = 0.19566:
alpha/eta0^2 = 0.309, f0 = 1.014).  Con zeta = 1 (su capa) hay un modo rapido y uno lento; con zeta = 0, ninguno real.
Salida: resultados/trayectoria_zeta_F.json.   Uso: python trayectoria_zeta_F.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import CapaF, raices

ZETAS = [1.5, 1.25, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
WS = -np.logspace(6, -4, 121)


def tarea(z):
    f = FondoF("cuadratico", 0.09, 0.5, 0.19566).resolver()
    rs, res = raices(f, CapaF(f, zeta=z), WS)
    return dict(zeta=z, w=rs, winf=[w * math.exp(-2 * f.Phi_inf) for w in rs], res_tt=res, Phi_inf=f.Phi_inf)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    with ProcessPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(tarea, ZETAS))
    for o in res:
        print(f"zeta={o['zeta']:+.2f}: w* = {[f'{x:.6g}' for x in o['w']]}  omega_inf^2 r0^2 = {[f'{x:.4g}' for x in o['winf']]}")
    json.dump(res, open("resultados/trayectoria_zeta_F.json", "w", encoding="utf-8"), indent=1)
