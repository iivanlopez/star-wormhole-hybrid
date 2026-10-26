# -*- coding: utf-8 -*-
"""
Fase 2h: optimizacion (evolucion diferencial) de los fondos fondo_nu.FondoNu (violacion de la condicion de energia nula
confinada a la garganta) frente a la criba de condiciones necesarias de la Fase 2g (criba.criba_normalizada) y la
condicion exterior r sig(ls) > 0 (superficie fuera de la esfera de fotones; ver diseno_ext.py: con sig > 0 una Z por
encima de la natural, decreciente, sube U sin tocar K mas que en su presupuesto natural).
Variables: log wt, log nu0, f_fl, phi1, log b, n, d0, log wd, log lam, log lT, m, log (rho~+p~)(0), log (Ns - 1).
Margen = min(margenes de la criba, r sig(ls)).  La prueba E usa la superficie libre (criba_normalizada(libre=True): h'(ls)
libre y energia de superficie del salto de densidad); con m = 1 fijo (la variable 10 no se usa).
Uso: python optimiza_nu.py [generaciones] [semilla]   -> resultados/optimiza_nu_<semilla>.json
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
sys.path.insert(0, os.path.join(AQUI, "..", "fase2g"))
import fondo_nu as FN      # noqa: E402
import criba as C          # noqa: E402

GEN = int(sys.argv[1]) if len(sys.argv) > 1 else 40
SEM = int(sys.argv[2]) if len(sys.argv) > 2 else 1
L = math.log
LIM = [(L(0.3), L(4.0)), (L(0.5), L(6.0)), (0.3, 1.5), (0.1, 1.5), (L(1.0), L(10.0)), (1, 3), (0.05, 2.0),
       (L(0.3), L(5.0)), (L(0.05), 0.0), (L(2.0), L(300.0)), (1.0, 2.0), (L(0.003), L(1.0)), (L(0.3), L(15.0))]


def parametros(v):
    m = 1.0
    rp0 = math.exp(v[11])
    Ns = 1 + math.exp(v[12])
    Y0, Ys = 0.5, 1 / (2 * Ns ** 2)
    A = rp0 / (2 * m * Y0 * (Y0 - Ys) ** (m - 1))
    return dict(wt=math.exp(v[0]), nu0=math.exp(v[1]), f_fl=v[2], phi1=v[3], b=math.exp(v[4]), n=int(round(v[5])),
                d0=v[6], wd=math.exp(v[7]), lam=math.exp(v[8]), lT=math.exp(v[9]), m=m, A=A, Ns=Ns)


def evaluar(x, n=80, tol=1e-15):
    F = FN.FondoNu(**x, tol=tol)
    if not F.ok:
        return F, None
    if F.ls <= 0 or F.reentra:
        F.ok, F.motivo = False, "sin estrella o con reentrada del fluido"
        return F, None
    r = C.criba_normalizada(F, n=n, libre=True)
    r["ext"] = F.margen_exterior()
    r["margen"] = min(r["margen"], r["ext"])
    return F, r


def objetivo(v):
    warnings.filterwarnings("ignore")
    try:
        F, r = evaluar(parametros(v))
        if r is None or not np.isfinite(r["margen"]):
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

    # semilla: el mejor punto de una rejilla de pozos profundos (wt = 2, nu0 = 2.5, Ns - 1 = 6, (rho~+p~)(0) = 0.03)
    x0 = [L(2.0), L(2.5), 1.0, 0.5, L(3.0), 2, 0.5, 0.0, 0.0, L(30.0), 1.0, L(0.03), L(6.0)]
    res = differential_evolution(objetivo, LIM, integrality=[False] * 5 + [True] + [False] * 7, maxiter=GEN,
                                 popsize=8, seed=SEM, workers=8, updating="deferred", polish=False, callback=cb,
                                 tol=1e-6, x0=x0)
    x = parametros(res.x)
    F, fin = evaluar(x, n=160)
    print("mejor:", {k: round(v, 4) for k, v in fin.items()})
    print("parametros:", {k: (round(v, 5) if isinstance(v, float) else v) for k, v in x.items()},
          f"kappa = {F.kappa:.6f}, N0 = {F.N0:.4f}, M_adm = {F.M_adm:.4f}, ls = {F.ls:.4f}")
    json.dump(dict(x=x, kappa=F.kappa, margenes=fin, historia=hist, v=list(map(float, res.x))),
              open(os.path.join(AQUI, "resultados", f"optimiza_nu_{SEM}.json"), "w"), indent=1)
