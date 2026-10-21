# -*- coding: utf-8 -*-
"""
Convergencia de los dos modos inestables (rapido y lento) de dos fondos cuadraticos de Antoniou et al. con la capa de
Antoniou (zeta = 1): alpha/eta0^2 = 0.3087 (e^{2 Phi_inf} = 825) y 0.3465 (1.65e4, el modo lento mas lento de la familia,
omega_inf^2 r0^2 = -2.6e-5, longitud de decaimiento ~200 r0).  Variantes: cola de integracion (l_extra), tolerancia de
la integracion de D, tolerancia y alcance del fondo.
Salida: resultados/convergencia_F.json.   Uso: python convergencia_F.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq
import numpy as np
from fondo_F import FondoF
from estabilidad_F import CapaF, D, vector_menores

# corchetes (reloj de l = 0) alrededor de las raices del barrido robusto
FONDOS = {"0.3087": {"rapido": (-560.0, -300.0), "lento": (-1.2, -0.3)},
          "0.3465": {"rapido": (-1.2e4, -5e3), "lento": (-0.8, -0.2)}}
VARIANTES = [("referencia", {}, {}), ("cola 60/kappa", {"l_extra": 60.0}, {}), ("rtol 1e-12", {"rtol": 1e-12}, {}),
             ("fondo rtol 1e-13", {}, {"rtol": 1e-13}), ("fondo l_max 4e4", {}, {"l_max": 4e4})]


def tarea(args):
    nombre, kw_d, kw_f, fondo, modo, par = args
    f = FondoF(*par, **kw_f).resolver()
    c = CapaF(f, zeta=1.0)

    def g(w):
        return D(f, c, w, u=np.real(vector_menores(f, c, float(w))), **kw_d)
    a, b = FONDOS[fondo][modo]
    w = brentq(g, a, b, xtol=1e-12 * abs(a), rtol=1e-12)
    return dict(variante=nombre, fondo=fondo, modo=modo, w=w, winf=w * math.exp(-2 * f.Phi_inf), Phi_inf=f.Phi_inf)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    fam = {f"{o['alpha_eta']:.4f}": ("cuadratico", o["alpha"], o["phi0"], o["phi1"])
           for o in json.load(open("resultados/familia_cuadratica.json", encoding="utf-8"))}
    casos = [(n, kd, kf, fo, m, fam[fo]) for fo in FONDOS for n, kd, kf in VARIANTES for m in FONDOS[fo]]
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(tarea, casos))
    for o in res:
        print(f"{o['fondo']} {o['variante']:18s} {o['modo']:7s} w* = {o['w']:.10g}   omega_inf^2 r0^2 = {o['winf']:.8g}")
    json.dump(res, open("resultados/convergencia_F.json", "w", encoding="utf-8"), indent=1)
