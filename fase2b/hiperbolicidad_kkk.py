# -*- coding: utf-8 -*-
"""
Hiperbolicidad (estabilidad_F.velocidades, mismo criterio que hiperbolicidad.py) de los 19 fondos de KKK de §6
(resultados/estabilidad_dominio.json).   Salida: resultados/hiperbolicidad_kkk.json.   Uso: python hiperbolicidad_kkk.py
"""
import sys
import json
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import velocidades

LS = np.concatenate([[0.0], np.logspace(-4, 1.7, 160)])


def tarea(o):
    f = FondoF("exp", o["alpha"], 0.0, o["phi1"]).resolver()
    lam = [min(v) if v else float("nan") for v in (velocidades(f, l) for l in LS if l <= 0.9 * f.l_end)]
    return dict(alpha_kkk=o["alpha_kkk"], f0=o["f0"], lambda_min=float(np.nanmin(lam)), hiperbolico=bool(np.nanmin(lam) > 0))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    dom = [o for o in json.load(open("resultados/estabilidad_dominio.json", encoding="utf-8")) if "w_zeta1" in o]
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(tarea, dom))
    for o in res:
        print(f"KKK {o['alpha_kkk']:.4f} {o['f0']:.4f}: lambda_min = {o['lambda_min']:+.4g}  hiperbolico = {o['hiperbolico']}")
    print(f"{sum(o['hiperbolico'] for o in res)}/{len(res)} hiperbolicos")
    json.dump(res, open("resultados/hiperbolicidad_kkk.json", "w", encoding="utf-8"), indent=1)
