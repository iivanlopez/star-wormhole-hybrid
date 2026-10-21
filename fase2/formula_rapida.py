# -*- coding: utf-8 -*-
"""
Prediccion semi-analitica del modo rapido del Modelo L (fluido debil sobre fondo de Ellis):
    w* = -1 / (16 r0^3 Sigma),   Sigma = int_0^{l_s} (eps + P) dl     [w interno, G = c = 1]
Se compara con las raices numericas de resultados/modo_rapido.json.
"""
import sys
import json
import math
import numpy as np
from scipy.integrate import quad
from hibrido import Politropo, Fondo, RHO_GEOM

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
eos = Politropo(k=0.1, gamma=2.0)
res = [o for o in json.load(open("resultados/modo_rapido.json", encoding="utf-8")) if o["modelo"] == "ligado"]
print(f"{'xi0':>4} {'rho':>8} {'r0':>6} {'Sigma':>10} | {'w* formula':>11} {'w* numerico':>12} {'cociente':>9}")
salida = []
for o in sorted(res, key=lambda o: (o["xi0"], o["rho"])):
    rb = o["rho"] * RHO_GEOM
    f = Fondo(r0=o["xi0"] * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    Sigma = quad(lambda l: sum(f.bg(f.y_bg(l), True)[12:14]), 0.0, f.ls, limit=400)[0]
    w_form = -1.0 / (16.0 * f.r0**3 * Sigma)
    w_num = [w * math.exp(2 * f.Phi_inf) for w in o["raices"]]
    wn = min(w_num) if w_num else float("nan")
    print(f"{o['xi0']:>4} {o['rho']:.1e} {f.r0:6.3f} {Sigma:10.3e} | {w_form:11.4e} {wn:12.4e} {wn / w_form:9.3f}")
    salida.append(dict(xi0=o["xi0"], rho=o["rho"], r0=f.r0, Sigma=Sigma, Phi_inf=f.Phi_inf,
                       w_formula_fis=w_form * math.exp(-2 * f.Phi_inf), w_num_fis=min(o["raices"]) if o["raices"] else None))
json.dump(salida, open("resultados/formula_rapida.json", "w", encoding="utf-8"), indent=1)
