# -*- coding: utf-8 -*-
"""
SUSTITUIDO por estabilidad_robusta.py: el seguimiento de signo del vector nulo unitario puede perder raices
(o crear raices espurias); se conserva por trazabilidad.
Estabilidad radial (modos pares, garganta/ecuador libre, capa dinamica) de la familia de Antoniou et al. (2019) con
F = alpha phi^2 y phi_inf = 0 (resultados/familia_cuadratica.json).  Raices reales de D(w) en [-1e6, -1e-4]
(reloj de l = 0) con la capa de Antoniou/KKK (zeta = 1, lambda0 = lambda1) y sin gravedad inducida (zeta = 0).
Salida: resultados/estabilidad_cuadratica.json.   Uso: python estabilidad_cuadratica.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import CapaF, raices

WS = -np.logspace(6, -4, 121)


def tarea(o):
    f = FondoF("cuadratico", o["alpha"], o["phi0"], o["phi1"]).resolver()
    out = dict(o)
    for z in (1.0, 0.0):
        try:
            rs, res = raices(f, CapaF(f, zeta=z), WS)
        except Exception as e:
            rs, res = [], float("nan")
            out[f"error_zeta{z:g}"] = str(e)
        out[f"w_zeta{z:g}"] = rs
        out[f"winf_zeta{z:g}"] = [w * math.exp(-2 * f.Phi_inf) for w in rs]
        out[f"res_tt_zeta{z:g}"] = res
    return out


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    fam = json.load(open("resultados/familia_cuadratica.json", encoding="utf-8"))
    # los fondos no hiperbolicos (hiperbolicidad.py) tienen el problema de valores iniciales mal planteado: su espectro
    # de modos no tiene cota inferior y no se calcula
    hip = {o["nombre"]: o for o in json.load(open("resultados/hiperbolicidad.json", encoding="utf-8"))}
    try:
        previos = {(o["alpha"], o["phi0"], round(o["phi1"], 8)): o
                   for o in json.load(open("resultados/estabilidad_cuadratica.json", encoding="utf-8"))}
    except FileNotFoundError:
        previos = {}
    res, pend = [], []
    for o in fam:
        h = hip[f"cuadratico {o['alpha_eta']:.4f}"]
        clave = (o["alpha"], o["phi0"], round(o["phi1"], 8))
        if not h["hiperbolico"]:
            res.append(dict(o, hiperbolico=False, l_no_hip=h["l_no_hip"], lambda_min=h["lambda_min"]))
        elif clave in previos:
            res.append(dict(previos[clave], hiperbolico=True))
        else:
            pend.append(o)
    print(f"{len(res)} ya resueltos o no hiperbolicos; {len(pend)} por calcular", flush=True)
    with ProcessPoolExecutor(max_workers=8) as ex:
        for o in ex.map(tarea, pend):
            o["hiperbolico"] = True
            res.append(o)
            print(f"alpha/eta0^2={o['alpha_eta']:.4f} f0={o['f0']:+.4g} ecuador={o['ecuador']}: "
                  f"omega_inf^2 r0^2 (zeta=1) = {[f'{x:.4g}' for x in o['winf_zeta1']]};  (zeta=0) = "
                  f"{[f'{x:.4g}' for x in o['winf_zeta0']]}  res={max(o['res_tt_zeta1'], o['res_tt_zeta0']):.0e}",
                  flush=True)
            json.dump(res, open("resultados/estabilidad_cuadratica.json", "w", encoding="utf-8"), indent=1)
