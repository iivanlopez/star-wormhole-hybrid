# -*- coding: utf-8 -*-
"""
Robustez de la clasificacion de las mitades repulsivas (a < 0) de los mapas.
Muestra estratificada (por acoplamiento y tipo de final; todas las de finales poco frecuentes) reintegrada con
rtol = 1e-12 (frente a 1e-9), l_max = 2000 (frente a 400) y umbral |Phi'|, |phi'| < 1e9 (frente a 1e7).
Se comparan el tipo de final y la distancia l_s.
Salida: resultados/robustez.json.   Uso: python robustez.py
"""
import os
import sys
import json
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fondo_doble      # noqa: E402
from fondo_doble import Lado      # noqa: E402


def rehacer(o):
    fondo_doble.PHIP_MAX = 1e9
    L = Lado(o["tipo"], o["alpha"], o["phi0"], o["phi1"], l_max=2000.0, rtol=1e-12).resolver()
    return dict(tipo=o["tipo"], alpha=o["alpha"], phi0=o["phi0"], phi1=o["phi1"], fin0=o["fin"], l0=o["l_end"],
                fin1=L.fin, l1=float(L.l_end))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    t0 = time.time()
    res = []
    for f in ("resultados/mapa_lados.json", "resultados/mapa_gauss.json", "resultados/mapa_lineal.json"):
        if os.path.exists(f):
            res += json.load(open(f, encoding="utf-8"))
    rep = [o for o in res if o["a"] < 0]
    rng = np.random.default_rng(2026)
    muestra = []
    grupos = {}
    for o in rep:
        grupos.setdefault((o["tipo"], o["fin"]), []).append(o)
    for k, lst in sorted(grupos.items()):
        n = len(lst) if len(lst) <= 150 else 150
        muestra += [lst[i] for i in rng.choice(len(lst), n, replace=False)]
        print(f"   {k}: {len(lst)} mitades, {n} en la muestra")
    with ProcessPoolExecutor(max_workers=8) as ex:
        out = list(ex.map(rehacer, muestra, chunksize=4))
    igual = sum(o["fin0"] == o["fin1"] for o in out)
    af = sum(o["fin1"] == "AF" for o in out)
    print(f"muestra: {len(out)}; mismo final: {igual}; AF con la integracion fina: {af}")
    cambios = {}
    for o in out:
        if o["fin0"] != o["fin1"]:
            cambios[(o["fin0"], o["fin1"])] = cambios.get((o["fin0"], o["fin1"]), 0) + 1
    print("cambios de clase:", cambios)
    hz = [o for o in out if o["fin0"] == "horizonte" and o["fin1"] == "horizonte"]
    if hz:
        d = np.array([abs(o["l1"] - o["l0"]) / o["l0"] for o in hz])
        print(f"horizonte -> horizonte: {len(hz)}; |dl_s|/l_s mediana {np.median(d):.2e}, maximo {d.max():.2e}")
    json.dump(out, open("resultados/robustez.json", "w", encoding="utf-8"), indent=1)
    print(f"tiempo: {time.time() - t0:.0f} s")
