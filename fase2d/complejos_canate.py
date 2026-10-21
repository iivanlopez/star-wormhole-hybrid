# -*- coding: utf-8 -*-
"""
Principio del argumento para Delta_par(w) y Delta_impar(w) (modos_canate.py) del agujero de Canate et al. con Qs = Q:
numero de ceros en el rectangulo Re w in [W1, W2], Im w in [-H, H] (lejos del corte w > 0), con cola L fija (para que
Delta sea una unica funcion analitica en todo el contorno) y muestreo adaptativo (salto de fase < 0.35 rad).
Se compara con las raices reales del barrido (modos_canate.py) en [W1, W2]: si coinciden, no hay modos complejos.
Salida: resultados/complejos_canate.json.   Uso: python complejos_canate.py
"""
import os
import sys
import json
import math
import time
import cmath
from concurrent.futures import ProcessPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import canate as C          # noqa: E402
import modos_canate as Mo   # noqa: E402

QS = 1.0
W1, W2, H = -7.0, -0.02, 3.0
_F = None
_L = None


def init():
    global _F, _L
    _F = C.FondoCanate(QS)
    l1, l2, I, _ = Mo.capa(_F)
    _L = l2 + 40.0


def evalua(w):
    Dp, Di, lg, _ = Mo.deltas(_F, complex(w), _L)
    return complex(Dp), complex(Di)


def contorno():
    lados = [(complex(W1, -H), complex(W2, -H)), (complex(W2, -H), complex(W2, H)),
             (complex(W2, H), complex(W1, H)), (complex(W1, H), complex(W1, -H))]
    pts = []
    for a, b in lados:
        n = 60
        pts += [a + (b - a) * t for t in np.linspace(0, 1, n, endpoint=False)]
    return pts


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    t0 = time.time()
    pts = contorno()
    with ProcessPoolExecutor(max_workers=4, initializer=init) as ex:
        vals = dict(zip(pts, ex.map(evalua, pts, chunksize=4)))
        for it in range(12):
            orden = pts + [pts[0]]
            nuevos = []
            for k in range(2):
                for a, b in zip(orden[:-1], orden[1:]):
                    da = cmath.phase(vals[b][k] / vals[a][k])
                    if abs(da) > 0.35 and abs(b - a) > 1e-6:
                        nuevos.append((a + b) / 2)
            nuevos = sorted(set(nuevos), key=lambda z: pts.index(min(pts, key=lambda p: abs(p - z))))
            if not nuevos:
                break
            for z, v in zip(nuevos, ex.map(evalua, nuevos, chunksize=4)):
                vals[z] = v
            # insertar en orden a lo largo del contorno
            todos = list(vals.keys())
            def s_param(z):
                # parametro de longitud de arco a lo largo del rectangulo (sentido antihorario)
                x, y = z.real, z.imag
                if abs(y + H) < 1e-12:
                    return (x - W1)
                if abs(x - W2) < 1e-12:
                    return (W2 - W1) + (y + H)
                if abs(y - H) < 1e-12:
                    return (W2 - W1) + 2 * H + (W2 - x)
                return 2 * (W2 - W1) + 2 * H + (H - y)
            pts = sorted(todos, key=s_param)
            print(f"iteracion {it}: {len(pts)} puntos", flush=True)
    orden = pts + [pts[0]]
    res = {}
    for k, nombre in ((0, "par"), (1, "impar")):
        tot = sum(cmath.phase(vals[b][k] / vals[a][k]) for a, b in zip(orden[:-1], orden[1:]))
        salto = max(abs(cmath.phase(vals[b][k] / vals[a][k])) for a, b in zip(orden[:-1], orden[1:]))
        res[nombre] = dict(N=tot / (2 * math.pi), salto_max=salto)
        print(f"{nombre}: N = {tot / (2 * math.pi):.4f} ceros en Re w in [{W1}, {W2}], |Im w| <= {H}; salto max {salto:.3f}")
    json.dump(dict(Qs=QS, W1=W1, W2=W2, H=H, puntos=len(pts), **res), open("resultados/complejos_canate.json", "w",
                                                                             encoding="utf-8"), indent=1)
    print(f"tiempo: {time.time() - t0:.0f} s")
