# -*- coding: utf-8 -*-
"""
Convergencia de las raices reales de Delta_par y Delta_impar (modos_canate.py) frente a la longitud de la cola
(L = l2 + 15 Q -> l2 + 40 Q) y la tolerancia (rtol 1e-10 -> 1e-12), para Qs = Q: primeras raices y la mayor de cada
sector.  Entrada: resultados/modos_canate.json.  Salida: resultados/convergencia_canate.json.
Uso: python convergencia_canate.py
"""
import os
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import canate as C          # noqa: E402
import modos_canate as Mo   # noqa: E402


def refina(args):
    Qs, sector, w0, extra, rtol = args
    f = C.FondoCanate(Qs)
    l1, l2, I, _ = Mo.capa(f)
    k = 0 if sector == "par" else 1

    def g(w):
        return Mo.deltas(f, w, max(l2 + extra, 10.0 + 30.0 / math.sqrt(-w)), rtol)[k]
    a, b = w0 * (1 + 2e-3), w0 * (1 - 2e-3)
    return dict(Qs=Qs, sector=sector, w0=w0, extra=extra, rtol=rtol, w=brentq(g, a, b, xtol=1e-13, rtol=1e-13))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    mod = json.load(open("resultados/modos_canate.json", encoding="utf-8"))
    o = [m for m in mod if m["Qs"] == 1.0][0]
    tareas = []
    for sector in ("par", "impar"):
        rs = [r["w"] for r in o["raices"][sector]]
        for w0 in rs[:3] + rs[-1:]:
            for extra, rtol in ((15.0, 1e-10), (40.0, 1e-10), (15.0, 1e-12)):
                tareas.append((1.0, sector, w0, extra, rtol))
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(refina, tareas))
    for i in range(0, len(res), 3):
        a, b, c = res[i:i + 3]
        print(f"{a['sector']:5s} w0 = {a['w0']:+.10e}: L+15 -> {a['w']:+.12e}; L+40 -> {b['w']:+.12e} "
              f"(rel {abs(b['w'] - a['w']) / abs(a['w']):.1e}); rtol 1e-12 -> {c['w']:+.12e} (rel {abs(c['w'] - a['w']) / abs(a['w']):.1e})")
    json.dump(res, open("resultados/convergencia_canate.json", "w", encoding="utf-8"), indent=1)
