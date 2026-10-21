# -*- coding: utf-8 -*-
"""
Hiperbolicidad y modos radiales inestables (metodo robusto, zeta = 1 y 0) de la muestra gaussiana
F = alpha e^{-phi^2}, phi_inf = 0 (resultados/familia_gauss.json, familia_gauss.py).
Salida: resultados/estabilidad_gauss.json.   Uso: python estabilidad_gauss.py [n_procesos]
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import CapaF, raices_menores, velocidades

WS = -np.logspace(6, -4, 121)
LS = np.concatenate([[0.0], np.logspace(-4, 1.7, 160)])
KW = {"gamma": 1.0}


def tarea(o):
    f = FondoF("gauss", o["alpha"], o["phi0"], o["phi1"], kw=KW).resolver()
    lam = []
    for l in LS:
        if l > 0.9 * f.l_end:
            break
        v = velocidades(f, l)
        lam.append(min(v) if v else float("nan"))
    lam = np.array(lam)
    neg = LS[:len(lam)][lam < 0]
    out = dict(o, lambda_min=float(np.nanmin(lam)), hiperbolico=bool(len(neg) == 0),
               l_no_hip=[float(neg.min()), float(neg.max())] if len(neg) else None)
    if out["hiperbolico"]:
        for z in (1.0, 0.0):
            rs = raices_menores(f, CapaF(f, zeta=z), WS)
            out[f"w_zeta{z:g}"] = rs
            out[f"winf_zeta{z:g}"] = [x * math.exp(-2 * f.Phi_inf) for x in rs]
    return out


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    fam = json.load(open("resultados/familia_gauss.json", encoding="utf-8"))
    res = []
    with ProcessPoolExecutor(max_workers=int(sys.argv[1]) if len(sys.argv) > 1 else 8) as ex:
        for o in ex.map(tarea, fam):
            res.append(o)
            txt = (f"zeta=1: {[f'{x:.4g}' for x in o['winf_zeta1']]}  zeta=0: {[f'{x:.4g}' for x in o['winf_zeta0']]}"
                   if o["hiperbolico"] else f"no hiperbolico en l en [{o['l_no_hip'][0]:.3g}, {o['l_no_hip'][1]:.3g}]")
            print(f"gauss alpha/eta0^2={o['alpha_eta']:.4f} f0={o['f0']:+.4g} ecuador={int(o['ecuador'])} "
                  f"lambda_min={o['lambda_min']:+.3g}: {txt}", flush=True)
            json.dump(res, open("resultados/estabilidad_gauss.json", "w", encoding="utf-8"), indent=1)
