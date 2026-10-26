# -*- coding: utf-8 -*-
"""
Optimizacion (evolucion diferencial) del margen minimo de la criba de condiciones necesarias (criba.criba_normalizada)
para gargantas de GLPV cuartica llenas de fluido.  Variables: log M, log a, log b, n (entero 1..3), d0, log wd,
u (Ns = N0 + (1 - N0) u), m, log (rho~+p~)_0, log lam, tauT.  Se exige N0 >= 0.1 (lejos de un horizonte).
Margen > 0: la configuracion cumple las condiciones necesarias (familia c = 0 y presupuesto de fantasma).
Uso: python criba_optimiza.py [generaciones] [semilla]   -> resultados/criba_optimiza_<semilla>.json
"""
import os
import sys
import json
import math
import time
import numpy as np
from scipy.optimize import differential_evolution

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import criba as C          # noqa: E402

GEN = int(sys.argv[1]) if len(sys.argv) > 1 else 40
SEM = int(sys.argv[2]) if len(sys.argv) > 2 else 1
LIM = [(math.log(0.2), math.log(10)), (math.log(0.5), math.log(12)), (math.log(1), math.log(10)), (1, 3), (0.05, 2.0),
       (math.log(0.3), math.log(5)), (0.05, 0.95), (1.0, 2.0), (math.log(0.01), math.log(10)), (math.log(0.02), 0.0),
       (1.0, 3.0)]


def parametros(v):
    M, a, b = math.exp(v[0]), math.exp(v[1]), math.exp(v[2])
    n = int(round(v[3]))
    N0 = math.exp(-M / a)
    Ns = N0 + (1 - N0) * v[6]
    m = v[7]
    rp0 = math.exp(v[8])
    lam = math.exp(v[9])
    Y0, Ys = 1 / (2 * N0 ** 2), 1 / (2 * Ns ** 2)
    A = rp0 / (2 * m * Y0 * (Y0 - Ys) ** (m - 1))
    return dict(M=M, a=a, b=b, n=n, d0=v[4], wd=math.exp(v[5]), Ns=Ns, m=m, A=A * lam, lam=lam, tauT=v[10])


def objetivo(v):
    x = parametros(v)
    if math.exp(-x["M"] / x["a"]) < 0.1:
        return 2.0
    try:
        F = C.FondoFluido(**x)
        if not F.ok:
            return 1.5
        r = C.criba_normalizada(F, n=80)
        if not np.isfinite(r["margen"]):
            return 1.5
        return -r["margen"]
    except Exception:                    # noqa: BLE001
        return 1.5


if __name__ == "__main__":
    t0 = time.time()
    hist = []

    def cb(xk, convergence=None):
        hist.append(float(objetivo(xk)))
        print(f"  gen {len(hist)}: mejor margen {-hist[-1]:+.4f}  [{time.time() - t0:.0f} s]", flush=True)

    res = differential_evolution(objetivo, LIM, integrality=[False, False, False, True] + [False] * 7, maxiter=GEN,
                                 popsize=12, seed=SEM, workers=8, updating="deferred", polish=False, callback=cb,
                                 tol=1e-6)
    x = parametros(res.x)
    F = C.FondoFluido(**x)
    fin = C.criba_normalizada(F, n=160)
    print("mejor:", {k: round(v, 4) for k, v in fin.items()})
    print("parametros:", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in x.items()}, f"N0 = {F.N0:.3f}, ls = {F.ls:.3f}")
    json.dump(dict(x=x, margenes=fin, historia=hist, v=list(map(float, res.x))),
              open(os.path.join(AQUI, "resultados", f"criba_optimiza_{SEM}.json"), "w"), indent=1)
