# -*- coding: utf-8 -*-
"""
Fase 2i: optimizacion (evolucion diferencial) de las estrellas elasticas con nucleo hueco (fondo_solido.FondoSolido)
frente a la criba de la familia c = 0 con las dos superficies libres (criba_solido.criba), mas dos condiciones sobre el
material: la condicion de energia dominante tangencial (|p~_t| <= rho~) y la estabilidad del material (e_v'' >= 0).
Variables: log wt, log nu0, f_m, phi1, log b, n, d0, log wd, log lam, log lT, log A, log (Ns - 1), log la, log lchi,
c_L^2, c_T^2/c_L^2.   Margen = min(E1n, Epar, Eimpar, Tpar, Timpar, dec_t, ev_min).
Uso: python optimiza_solido.py [generaciones] [semilla]   -> resultados/optimiza_solido_<semilla>.json
"""
import os
import sys
import json
import math
import time
import warnings
import numpy as np
from scipy.optimize import differential_evolution

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import fondo_solido as FS         # noqa: E402
import criba_solido as CS         # noqa: E402

GEN = int(sys.argv[1]) if len(sys.argv) > 1 else 30
SEM = int(sys.argv[2]) if len(sys.argv) > 2 else 1
L = math.log
LIM = [(L(0.3), L(4.0)), (L(0.5), L(6.0)), (0.3, 2.0), (0.1, 1.5), (L(1.0), L(10.0)), (1, 3), (0.05, 2.0),
       (L(0.3), L(5.0)), (L(0.05), 0.0), (L(2.0), L(300.0)), (L(0.003), L(1.0)), (L(0.3), L(15.0)), (L(0.3), L(10.0)),
       (L(0.1), L(10.0)), (0.3, 1.0), (0.0, 0.7)]


def parametros(v):
    return dict(wt=math.exp(v[0]), nu0=math.exp(v[1]), f_m=v[2], phi1=v[3], b=math.exp(v[4]), n=int(round(v[5])),
                d0=v[6], wd=math.exp(v[7]), lam=math.exp(v[8]), lT=math.exp(v[9]), A=math.exp(v[10]),
                Ns=1 + math.exp(v[11]), la=math.exp(v[12]), lchi=math.exp(v[13]), cL2=v[14], cT2=v[14] * v[15])


def evaluar(x, n=100, tol=1e-15):
    F = FS.FondoSolido(**x, tol=tol)
    if not F.ok or F.reentra:
        return F, None
    r = CS.criba(F, n=n)
    fis = F.fisica()
    r.update(dec_t=fis["dec_t"], ev_min=fis["ev_min"], masa=fis["masa"])
    r["margen"] = min(r["margen"], fis["dec_t"], fis["ev_min"])
    return F, r


def objetivo(v):
    warnings.filterwarnings("ignore")
    try:
        F, r = evaluar(parametros(v))
        if r is None or not np.isfinite(r["margen"]):
            return 1.5
        return -max(r["margen"], -1.4)
    except Exception:                    # noqa: BLE001
        return 1.5


if __name__ == "__main__":
    t0 = time.time()
    hist = []

    def cb(xk, convergence=None):
        hist.append(float(objetivo(xk)))
        print(f"  gen {len(hist)}: mejor margen {-hist[-1]:+.4f}  [{time.time() - t0:.0f} s]", flush=True)

    # semilla: el candidato de la Fase 2h con un nucleo hueco (la = 3) y sin cizalla
    x0 = [L(2.0), L(2.5), 1.0, 0.5, L(3.0), 2, 0.5, 0.0, 0.0, L(30.0), L(0.03), L(6.0), L(3.0), L(2.0), 1.0, 0.0]
    res = differential_evolution(objetivo, LIM, integrality=[False] * 5 + [True] + [False] * 10, maxiter=GEN,
                                 popsize=6, seed=SEM, workers=8, updating="deferred", polish=False, callback=cb,
                                 tol=1e-6, x0=x0)
    x = parametros(res.x)
    F, fin = evaluar(x, n=200)
    print("mejor:", {k: round(float(v), 4) for k, v in fin.items()})
    print("parametros:", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in x.items()},
          f"kappa = {F.kappa:.6f}, N0 = {F.N0:.4f}, M_adm = {F.M_adm:.4f}, ls = {F.ls:.4f}")
    json.dump(dict(x=x, kappa=F.kappa, margenes={k: float(v) for k, v in fin.items()}, historia=hist,
                   v=list(map(float, res.x))),
              open(os.path.join(AQUI, "resultados", f"optimiza_solido_{SEM}.json"), "w"), indent=1)
