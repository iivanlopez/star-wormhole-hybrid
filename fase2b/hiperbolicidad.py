# -*- coding: utf-8 -*-
"""
Hiperbolicidad de las perturbaciones radiales: cuadrado minimo de la velocidad caracteristica local lambda = w/k^2
(estabilidad_F.velocidades) a lo largo de cada fondo.  lambda < 0 en algun punto: el problema de valores iniciales esta
mal planteado (equivale a q_sigma < 0 en KKK/CKZ) y aparecen modos inestables sin cota de velocidad de crecimiento.
Fondos: familia cuadratica (resultados/familia_cuadratica.json), dobles gargantas exponenciales (estabilidad_doble.py)
y el fondo de referencia de KKK.
Salida: resultados/hiperbolicidad.json.   Uso: python hiperbolicidad.py
"""
import sys
import json
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import velocidades
from estabilidad_doble import CASOS as DOBLES

LS = np.concatenate([[0.0], np.logspace(-4, 1.7, 160)])


def tarea(args):
    nombre, tipo, al, p0, u = args
    f = FondoF(tipo, al, p0, u).resolver()
    lam = []
    for l in LS:
        if l > 0.9 * f.l_end:
            break
        v = velocidades(f, l)
        lam.append(min(v) if v else float("nan"))
    lam = np.array(lam)
    neg = LS[:len(lam)][lam < 0]
    return dict(nombre=nombre, tipo=tipo, alpha=al, phi0=p0, phi1=u, ecuador=bool(f.ecuador), alpha_eta=f.alpha_eta,
                lambda_min=float(np.nanmin(lam)), hiperbolico=bool(len(neg) == 0),
                l_no_hip=[float(neg.min()), float(neg.max())] if len(neg) else None)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    fam = json.load(open("resultados/familia_cuadratica.json", encoding="utf-8"))
    casos = [("KKK referencia", "exp", 0.03, 0.0, -0.7)]
    casos += [(f"cuadratico {o['alpha_eta']:.4f}", "cuadratico", o["alpha"], o["phi0"], o["phi1"]) for o in fam]
    casos += [(f"doble exp a={a} u={u}", "exp", a, 0.0, u) for a, u in DOBLES]
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(tarea, casos))
    for o in res:
        print(f"{o['nombre']:28s} ecuador={int(o['ecuador'])} hiperbolico={o['hiperbolico']}  lambda_min={o['lambda_min']:+.3g}"
              + (f"  no hiperbolico en l en [{o['l_no_hip'][0]:.3g}, {o['l_no_hip'][1]:.3g}]" if o["l_no_hip"] else ""))
    json.dump(res, open("resultados/hiperbolicidad.json", "w", encoding="utf-8"), indent=1)
