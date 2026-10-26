# -*- coding: utf-8 -*-
"""
Fase 2h: espectro del mejor fondo de la optimizacion con la criba corregida (resultados/optimiza_nu_<semilla>.json):
familia c = 0 (sin funciones libres) y sistema completo frente a la rigidez Lam de las direcciones c != 0, con la
superficie del fluido libre y fija (diseno theta = 0, eps = 0.05, cm2 = 1, b_in = 0.05; malla n_in 240, L 1500, n_out 250).
Uso: python espectro_optimo.py [semilla]   -> resultados/optimo_nu_<semilla>_c0.json
"""
import os
import sys
import json
import time
import warnings
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2f"), os.path.join(AQUI, "..", "fase2g"), AQUI):
    sys.path.insert(0, p_)
import fondo_nu as FN                 # noqa: E402
import teoria_nu as TN                # noqa: E402
import espectro2 as E2                # noqa: E402

warnings.filterwarnings("ignore")
SEM = int(sys.argv[1]) if len(sys.argv) > 1 else 3
T0 = time.time()
x = json.load(open(os.path.join(AQUI, "resultados", f"optimiza_nu_{SEM}.json")))["x"]
F = FN.FondoNu(**x)
mallas = (240, 1500.0, 250)
out = []
for Lam in (10.0, 100.0, 1000.0):
    T = TN.TeoriaNu(F, theta=0.0, eps=0.05, cm2=1.0, b_in=0.05, Lam=Lam)
    lc = T.lc()
    a, b = E2.mallas(lc, *mallas)
    cf = E2.Coefs(T, np.concatenate([a, b[1:]]))
    for c0 in ((False, True) if Lam == 10.0 else (False,)):
        for libre in (True, False):
            for par in ("par", "impar"):
                o = E2.espectro(T, par, lc, *mallas, k=3, coefs=cf, libre=libre, c_cero=c0)
                out.append(dict(Lam=Lam, c_cero=c0, libre=libre, par=par, w=o.get("w"), Tmin=o["Tmin"]))
                print(f"[{time.time() - T0:6.1f} s] Lam {Lam:g} {'c=0' if c0 else 'completo'} "
                      f"{'libre' if libre else 'fija'} {par}: w = {np.array2string(np.array(o.get('w')), precision=4)}",
                      flush=True)
json.dump(out, open(os.path.join(AQUI, "resultados", f"optimo_nu_{SEM}_c0.json"), "w"), indent=1, default=float)
