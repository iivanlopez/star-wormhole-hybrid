# -*- coding: utf-8 -*-
"""
Busqueda aleatoria con la criba de condiciones necesarias (criba.py) sobre gargantas de GLPV cuartica llenas de fluido.
Parametros: M log-uniforme en [0.3, 10], a log-uniforme en [0.5, 12], b log-uniforme en [1, 10], n en {1, 2, 3},
d0 en [0.05, 2], wd log-uniforme en [0.3, 5], Ns = N0 + (1 - N0) u con u en [0.05, 0.95], m en [1, 2.5],
(rho~ + p~) en la garganta log-uniforme en [0.01, 10] (fija A), multiplicada por lam: acoplo G4 = S G4b y Gcal reducidos por
S = lam en el nucleo (lam log-uniforme en [0.02, 1] en el 80 % de las muestras), con transicion en lT = tauT ls,
tauT en [1, 3].
Para cada muestra valida (garganta unica, r' > 0 fuera de ella): E1 = E[d = N, xi = 0], autovalores minimos de E (par e
impar) y de T (par e impar).  Pasa la criba si los cinco son >= 0.  Escribe resultados/criba_busqueda.json.
Uso: python criba_busqueda.py [muestras] [semilla]
"""
import os
import sys
import json
import math
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import criba as C          # noqa: E402

NS = int(sys.argv[1]) if len(sys.argv) > 1 else 400
SEM = int(sys.argv[2]) if len(sys.argv) > 2 else 1


def muestras(ns, sem):
    rng = np.random.default_rng(sem)
    out = []
    for _ in range(ns):
        M = float(np.exp(rng.uniform(np.log(0.3), np.log(10))))
        a = float(np.exp(rng.uniform(np.log(0.5), np.log(12))))
        N0 = math.exp(-M / a)
        u = float(rng.uniform(0.05, 0.95))
        m = float(rng.uniform(1.0, 2.5)) if rng.uniform() < 0.7 else 1.0
        rp0 = float(np.exp(rng.uniform(np.log(0.01), np.log(10))))
        Ns = N0 + (1 - N0) * u
        Y0, Ys = 1 / (2 * N0 ** 2), 1 / (2 * Ns ** 2)
        A = rp0 / (2 * m * Y0 * (Y0 - Ys) ** (m - 1)) if m != 1 else rp0 / (2 * Y0)
        lam = float(np.exp(rng.uniform(np.log(0.02), np.log(1.0)))) if rng.uniform() < 0.8 else 1.0
        out.append(dict(M=M, a=a, b=float(np.exp(rng.uniform(np.log(1), np.log(10)))), n=int(rng.choice([1, 2, 3])),
                        d0=float(rng.uniform(0.05, 2.0)), wd=float(np.exp(rng.uniform(np.log(0.3), np.log(5)))),
                        Ns=float(Ns), m=m, A=float(A * lam), lam=lam, tauT=float(rng.uniform(1.0, 3.0))))
    return out


def evalua(x):
    out = dict(x=x)
    try:
        F = C.FondoFluido(**x)
        if not F.ok:
            out["invalido"] = F.motivo
            return out
        f0 = F.campos(np.array([1e-6]))
        out.update(N0=F.N0, ls=F.ls, rp0=float(f0["rp"][0]), cs2_0=float(f0["cs2"][0]))
        out.update(C.criba(F, n=120))
        out.update(C.presupuesto_fantasma(F))
        out["pasa_c0"] = bool(min(out["E1"], out["Emin_par"], out["Emin_impar"], out["Tmin_par"], out["Tmin_impar"]) >= 0)
        out["pasa"] = bool(out["pasa_c0"] and out["margen_fantasma"] > 0)
    except Exception as e:                    # noqa: BLE001
        out["error"] = repr(e)[:200]
    return out


if __name__ == "__main__":
    t0 = time.time()
    xs = muestras(NS, SEM)
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(evalua, xs, chunksize=4))
    val = [r for r in res if "E1" in r]
    print(f"{NS} muestras: {len(val)} fondos validos; errores {sum('error' in r for r in res)}")
    claves = ("E1", "Emin_par", "Emin_impar", "Tmin_par", "Tmin_impar")
    for k in claves:
        print(f"  {k} >= 0: {sum(r[k] >= 0 for r in val)}")
    print(f"  E1 >= 0 y T >= 0: {sum(r['E1'] >= 0 and r['Tmin_par'] >= 0 and r['Tmin_impar'] >= 0 for r in val)}")
    print(f"  E1, Emin_par >= 0 y T >= 0: "
          f"{sum(r['E1'] >= 0 and r['Emin_par'] >= 0 and r['Tmin_par'] >= 0 and r['Tmin_impar'] >= 0 for r in val)}")
    print(f"  pasan las cinco de la familia c = 0: {sum(r['pasa_c0'] for r in val)}")
    print(f"  presupuesto de fantasma positivo: {sum(r['margen_fantasma'] > 0 for r in val)}")
    print(f"  E1 >= 0 y presupuesto positivo: {sum(r['E1'] >= 0 and r['margen_fantasma'] > 0 for r in val)}")
    print(f"  pasan todas: {sum(r['pasa'] for r in val)}")
    for r in sorted(val, key=lambda r: -min(min(r[k] for k in claves), r["margen_fantasma"]))[:10]:
        print("   " + ", ".join(f"{k}={r[k]:+.4f}" for k in claves) + f", fantasma={r['margen_fantasma']:+.3f} | N0={r['N0']:.2f} rp0={r['rp0']:.3g} "
              f"cs2_0={r['cs2_0']:.2f} ls={r['ls']:.2f} | " + ", ".join(f"{k}={v:.3g}" for k, v in r["x"].items()))
    print(f"{time.time() - t0:.0f} s")
    json.dump(dict(muestras=NS, semilla=SEM, resultados=res),
              open(os.path.join(AQUI, "resultados", f"criba_busqueda_{SEM}.json"), "w"), indent=1, default=float)
