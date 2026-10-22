# -*- coding: utf-8 -*-
"""
Busqueda aleatoria de gargantas masivas estables en la familia de reconstruccion.Garganta (busqueda.evalua).
Parametros (distribuciones): n en {1, 2, 2, 3}; b log-uniforme en [1, 10]; frac = Gcal(0)/G4(0) en [0.6, 1]
(sector impar subluminal); wd/b en [0.3, 1]/sqrt(2n); M log-uniforme en [0.3, 10]; a log-uniforme en [0.5, 20];
ventana de K wK en [4, 25] y de P wP/wK en [1, 2]; cm2 (velocidad del monopolo en la ventana) en [0.5, 1];
F4_XX: cf2 = 0 o normal(0, 1) con anchura y centro aleatorios.
Para cada muestra: K > 0 y P > 0 (bien planteado y sin fantasmas; en la malla l <= 400 y en la zona lejana hasta
l = 4000), maxima velocidad radial fisica del monopolo y autovalores negativos (taquiones) por paridad.  Escribe resultados/busqueda_aleatoria.json.
Uso: python busqueda_aleatoria.py [n_muestras] [semilla]
"""
import os
import sys
import json
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import busqueda as B           # noqa: E402

NS = int(sys.argv[1]) if len(sys.argv) > 1 else 160
SEM = int(sys.argv[2]) if len(sys.argv) > 2 else 1


def muestras(ns, sem):
    rng = np.random.default_rng(sem)
    out = []
    for _ in range(ns):
        n = int(rng.choice([1, 2, 2, 3]))
        cf2 = float(rng.choice([0.0, rng.normal(0, 1)]))
        out.append(dict(n=n, b=float(np.exp(rng.uniform(np.log(1.0), np.log(10.0)))), frac=float(rng.uniform(0.6, 1.0)),
                        wdb=float(rng.uniform(0.3, 1.0) / np.sqrt(2 * n)),
                        M=float(np.exp(rng.uniform(np.log(0.3), np.log(10)))),
                        a=float(np.exp(rng.uniform(np.log(0.5), np.log(20)))), wK=float(rng.uniform(4, 25)),
                        wPK=float(rng.uniform(1.0, 2.0)), cm2=float(rng.uniform(0.5, 1.0)), cf2=cf2,
                        wf2=float(rng.uniform(0.5, 5)), lf2=float(rng.uniform(0, 5))))
    return out


def evalua(x):
    r = B.evalua(x)
    r["objetivo"] = B.objetivo(r)
    return r


if __name__ == "__main__":
    t0 = time.time()
    xs = muestras(NS, SEM)
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(evalua, xs))
    validos = [r for r in res if "n_par" in r]
    sublum = [r for r in validos if r["maxc2"] <= 1.0 + 1e-3]
    ntaq = [r["n_par"] + r["n_impar"] for r in validos]
    print(f"{NS} muestras: {len(validos)} con K > 0 y P > 0 ({len(sublum)} ademas subluminales); "
          f"taquiones del monopolo: minimo {min(ntaq) if ntaq else None}, mediana {int(np.median(ntaq)) if ntaq else None}")
    for r in sorted(validos, key=lambda r: r["objetivo"])[:8]:
        w = r["w_par"] + r["w_impar"]
        print(f"  objetivo {r['objetivo']:.3f}: ({r['n_par']}, {r['n_impar']}) taquiones, w0 = {min(w):+.4f}, "
              f"|w| minimo = {max(w):+.2e}, N0 = {r['N0']:.2f}, max c^2 = {r['maxc2']:.2f};  "
              + ", ".join(f"{k}={v:.3g}" for k, v in r["x"].items()))
    print(f"{time.time() - t0:.0f} s")
    json.dump(dict(muestras=NS, semilla=SEM, resultados=res), open(os.path.join(AQUI, "resultados", "busqueda_aleatoria.json"), "w"),
              indent=1, default=float)
