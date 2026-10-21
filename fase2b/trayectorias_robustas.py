# -*- coding: utf-8 -*-
"""
Repite con el metodo robusto (estabilidad_F.raices_menores) las trayectorias en zeta de barrido_zeta.py (fondo de
referencia de KKK, F = alpha e^{-phi}: alpha = 0.03, phi'(0+) = -0.7) y de trayectoria_zeta_F.py (fondo cuadratico de
Antoniou: alpha = 0.09, phi0 = 0.5, phi'(0+) = 0.19566), y busca el modo del fondo de KKK mas extremo
(alpha/r0^2 = 0.0029, f0 = 1.0001) en |w| hasta 1e8, fuera del rango de los barridos.  Tambien la trayectoria de un
fondo cuadratico casi negro (alpha/eta0^2 = 0.3465) hasta cerca de zeta = 2.
Salida: resultados/trayectorias_robustas.json.   Uso: python trayectorias_robustas.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import CapaF, raices_menores

WS = -np.logspace(6, -4, 121)
WS_EXTREMO = -np.logspace(8, 5, 61)
REFERENCIA = ("exp", 0.03, 0.0, -0.7)
ZETAS_KKK = [-1.0, -0.5, 0.0, 0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0]
ZETAS_CUAD = [3.0, 2.5, 2.0, 1.75, 1.5, 1.25, 1.0, 0.8, 0.6, 0.4, 0.2, 0.0]
ZETAS_CASI = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 1.9]


def tarea(args):
    nombre, (tipo, al, p0, u), z, ws = args
    f = FondoF(tipo, al, p0, u).resolver()
    rs = raices_menores(f, CapaF(f, zeta=z), ws)
    return dict(nombre=nombre, tipo=tipo, alpha=al, phi0=p0, phi1=u, zeta=z, Phi_inf=f.Phi_inf, w=rs,
                winf=[x * math.exp(-2 * f.Phi_inf) for x in rs], w_min=float(ws[0]), w_max=float(ws[-1]))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    extremo = [o for o in json.load(open("resultados/estabilidad_dominio.json", encoding="utf-8"))
               if abs(o["alpha_kkk"] - 0.0029) < 5e-5 and abs(o["f0"] - 1.0001) < 5e-5][0]
    fam = {f"{o['alpha_eta']:.4f}": ("cuadratico", o["alpha"], o["phi0"], o["phi1"])
           for o in json.load(open("resultados/familia_cuadratica.json", encoding="utf-8"))}
    CUADRATICO = fam["0.3087"]                     # alpha = 0.09, phi0 = 0.5
    CASI_NEGRO = fam["0.3465"]                     # alpha = 0.09, phi0 = 0.3: e^{2 Phi_inf} = 1.65e4
    casos = [("KKK referencia", REFERENCIA, z, WS) for z in ZETAS_KKK]
    casos += [("cuadratico 0.3087", CUADRATICO, z, WS) for z in ZETAS_CUAD]
    casos += [("cuadratico 0.3465", CASI_NEGRO, z, WS) for z in ZETAS_CASI]
    casos += [("KKK 0.0029 1.0001 (|w| hasta 1e8)", ("exp", extremo["alpha"], 0.0, extremo["phi1"]), z, WS_EXTREMO)
              for z in (1.0, 0.0)]
    res = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for o in ex.map(tarea, casos):
            res.append(o)
            print(f"{o['nombre']:34s} zeta={o['zeta']:+.2f}: w* = {[f'{x:.6g}' for x in o['w']]}  "
                  f"omega_inf^2 r0^2 = {[f'{x:.4g}' for x in o['winf']]}", flush=True)
            json.dump(res, open("resultados/trayectorias_robustas.json", "w", encoding="utf-8"), indent=1)
