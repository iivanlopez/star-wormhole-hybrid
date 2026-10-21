# -*- coding: utf-8 -*-
"""
Metodo antiguo (seguimiento de signo): con zeta = 0 pierde raices y crea una espuria.  La referencia es
estabilidad_robusta.py; este script se conserva por el residuo del empalme redundante.
Modos pares inestables de los agujeros de KKK sobre su dominio (alpha/r0^2 <= 0.05, f0 > 1), con la capa de KKK
(zeta = 1) y sin gravedad inducida (zeta = 0).  Raices reales w < 0 de D(w) en [-1e6, -1e-4] (reloj de la garganta).
omega_inf^2 r0^2 = w e^{-2 Phi_inf} es invariante bajo la simetria de escala (r0 = 1, phi0 = 0 en la integracion).
Salida: resultados/estabilidad_dominio.json.   Uso: python estabilidad_dominio.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq
from kkk import FondoKKK
from estabilidad_kkk import CapaKKK, barrido, vector_garganta, D

WS = -np.logspace(6, -4, 121)


def raices(f, c):
    out = barrido(f, c, WS)
    vals = [v for v, _ in out]
    rs = []
    for i in range(len(WS) - 1):
        if vals[i] * vals[i + 1] < 0:
            ref = vector_garganta(f, c, WS[i])[0]
            try:
                rs.append(brentq(lambda w: D(f, c, w, u=vector_garganta(f, c, w, ref)[0]), WS[i], WS[i + 1],
                                 xtol=1e-9, rtol=1e-10))
            except ValueError:
                rs.append(float("nan"))
    return rs, max(abs(r) for _, r in out)


def tarea(args):
    al, u0 = args
    f = FondoKKK(alpha=al, phi1=u0).resolver()
    o = dict(alpha=al, phi1=u0, f0=f.f0, alpha_kkk=f.alpha_kkk, regular=bool(f.regular), Phi_inf=f.Phi_inf)
    if not f.regular or f.alpha_kkk > 0.05:
        return o
    for z in (1.0, 0.0):
        rs, res = raices(f, CapaKKK(f, zeta=z))
        o[f"w_zeta{z:g}"] = rs
        o[f"winf_zeta{z:g}"] = [w * math.exp(-2 * f.Phi_inf) for w in rs]
        o[f"res_tt_zeta{z:g}"] = res
    return o


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    casos = [(a, u) for a in (0.003, 0.01, 0.02, 0.03, 0.045) for u in (-0.1, -0.3, -1.0, -3.0, -10.0)]
    res = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for o in ex.map(tarea, casos):
            res.append(o)
            if "w_zeta1" in o:
                print(f"alpha_kkk={o['alpha_kkk']:.4f} f0={o['f0']:9.4f}: w*(zeta=1) = {[f'{x:.5g}' for x in o['w_zeta1']]}"
                      f"  omega_inf^2 r0^2 = {[f'{x:.4g}' for x in o['winf_zeta1']]};  zeta=0: {[f'{x:.4g}' for x in o['winf_zeta0']]}",
                      flush=True)
            else:
                print(f"alpha={o['alpha']} phi1={o['phi1']}: fuera de rango (alpha_kkk={o['alpha_kkk']:.4f}, regular={o['regular']})",
                      flush=True)
    json.dump(res, open("resultados/estabilidad_dominio.json", "w", encoding="utf-8"), indent=1)
