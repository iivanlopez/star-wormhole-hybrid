# -*- coding: utf-8 -*-
"""
Autofuncion del modo rapido (zeta = 0 y 1) en fondos de KKK y cuadraticos: escala de decaimiento de la perturbacion del
lapso a' frente a la longitud 1/a, con a = Phi'(0+) la aceleracion propia de la garganta (§12 del informe).
Lee las raices de resultados/estabilidad_robusta.json.   Salida: resultados/autofuncion_F.json.
Uso: python autofuncion_F.py
"""
import sys
import json
import math
import numpy as np
from scipy.integrate import solve_ivp
from fondo_F import FondoF
from estabilidad_F import CapaF, vector_menores, derivadas, ligaduras

NOMBRES = ["KKK 0.0174 1.0004", "KKK 0.0452 1.0980", "cuadratico 0.3465", "cuadratico 0.3087", "cuadratico 0.4717"]

if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    rob = {(o["nombre"], o["zeta"]): o for o in json.load(open("resultados/estabilidad_robusta.json", encoding="utf-8"))}
    res = []
    for nom in NOMBRES:
        for z in (0.0, 1.0):
            o = rob[(nom, z)]
            f = FondoF(o["tipo"], o["alpha"], o["phi0"], o["phi1"]).resolver()
            a = f.y0()[3]
            c = CapaF(f, zeta=z)
            w = min(o["w"])
            u = np.real(vector_menores(f, c, w))
            u = u / abs(u[0])
            ls = np.geomspace(1e-4 / a, 3.0, 60)
            s = solve_ivp(lambda l, y: derivadas(f, l, y, w), (0.0, ls[-1]), u, method="DOP853", rtol=1e-11,
                          atol=1e-14, t_eval=ls)
            ap = s.y[0] / s.y[0][0]
            # longitud a la que |a'| cae a 1/e de su valor en la garganta
            k = np.argmax(np.abs(ap) < math.exp(-1))
            l_e = float(ls[k]) if k > 0 else float("nan")
            lig = max(max(abs(x) for x in ligaduras(f, l, s.y[:, i], w)) for i, l in enumerate(ls[::6]))
            res.append(dict(nombre=nom, zeta=z, a=a, w=w, C=-w / a**2, l_e=l_e, l_e_por_a=l_e * a,
                            d0=float(u[1]), q0=float(u[3]), ligaduras_max=lig))
            print(f"{nom:20s} zeta={z:g}: a r0 = {a:7.3f}  C = -w/a^2 = {-w / a**2:.4f}  "
                  f"a' cae a 1/e en l = {l_e:.4g} r0 = {l_e * a:.3f}/a   d(0)/a'(0) = {u[1]:+.2e}  q(0)/a'(0) = {u[3]:+.2e}",
                  flush=True)
    json.dump(res, open("resultados/autofuncion_F.json", "w", encoding="utf-8"), indent=1)
