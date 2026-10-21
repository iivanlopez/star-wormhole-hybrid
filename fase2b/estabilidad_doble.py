# -*- coding: utf-8 -*-
"""
SUSTITUIDO por estabilidad_robusta.py (metodo de menores); se conserva por trazabilidad y por CASOS.
Estabilidad radial de agujeros de DOBLE garganta de Antoniou et al. (2019) con acoplamiento exponencial F = alpha e^{-phi}:
la capa esta en el ecuador (l = 0, maximo de r) y las dos gargantas (l = +-l_t) son lisas, sin materia.
Mismo metodo que estabilidad_cuadratica.py (raices reales de D en [-1e6, -1e-4], zeta = 1 y 0).
Salida: resultados/estabilidad_doble.json.   Uso: python estabilidad_doble.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import CapaF, raices

WS = -np.logspace(6, -4, 121)
CASOS = [(0.02, -50.119), (0.02, -12.340), (0.04, -8.268), (0.04, -4.534), (0.06, -4.534), (0.08, -3.038)]
# solo dos de ellos son hiperbolicos (hiperbolicidad.py); en el resto el problema esta mal planteado
HIPERBOLICOS = [(0.02, -12.340), (0.04, -4.534)]


def tarea(args):
    al, u = args
    f = FondoF("exp", al, 0.0, u).resolver()
    o = dict(alpha=al, phi1=u, ecuador=bool(f.ecuador), regular=bool(f.regular), f0=f.f0, alpha_eta=f.alpha_eta,
             r_garganta=getattr(f, "r_garganta", None), l_garganta=f.l_garganta, Phi_inf=f.Phi_inf)
    for z in (1.0, 0.0):
        rs, res = raices(f, CapaF(f, zeta=z), WS)
        o[f"w_zeta{z:g}"] = rs
        o[f"winf_zeta{z:g}"] = [w * math.exp(-2 * f.Phi_inf) for w in rs]
        o[f"res_tt_zeta{z:g}"] = res
    return o


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    res = []
    with ProcessPoolExecutor(max_workers=2) as ex:
        for o in ex.map(tarea, HIPERBOLICOS):
            res.append(o)
            print(f"alpha={o['alpha']} phi1={o['phi1']}: ecuador={o['ecuador']} r_t={o['r_garganta']:.4f} "
                  f"alpha/eta0^2={o['alpha_eta']:.4f}: omega_inf^2 r0^2 (zeta=1) = {[f'{x:.4g}' for x in o['winf_zeta1']]}; "
                  f"(zeta=0) = {[f'{x:.4g}' for x in o['winf_zeta0']]}  res={max(o['res_tt_zeta1'], o['res_tt_zeta0']):.0e}",
                  flush=True)
    json.dump(res, open("resultados/estabilidad_doble.json", "w", encoding="utf-8"), indent=1)
