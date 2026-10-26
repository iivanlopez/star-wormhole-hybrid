# -*- coding: utf-8 -*-
"""
Fase 2g, parte I: comprobacion numerica del teorema del taquion del monopolo (prueba d = N) sobre todas las gargantas
reconstruidas de la Fase 2f con K, P > 0: la de referencia, la sucesion de campo debil, las vias de escape (campo
fuerte y superluminalidad) y las 19 configuraciones validas de la busqueda aleatoria.

Para cada una, con los coeficientes K, P, U del monopolo en l in (0, L], L = 3000:
  E[N]    = int (P N'^2 + U N^2)  frente a  -2 int N^3 Gcal (r'/r - Phi')^2 (la diferencia es el borde -2 N^3 r' Gcal/r en L);
  qN      = E[N]/int K N^2   (cota de Rayleigh con d = N: el autovalor mas bajo w0 <= qN);
  q1      = int U/int K      (cota de la Fase 2f con d = cte; puede ser positiva en campo fuerte);
  w0      = autovalor mas bajo calculado en la Fase 2f (resultados/gargantas_bH.json y busqueda_aleatoria.json).
Escribe resultados/verificar_N.json.
Uso: python verificar_N.py
"""
import os
import sys
import json
import math
import time
import numpy as np
from scipy.integrate import simpson
from concurrent.futures import ProcessPoolExecutor

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
import busqueda as B          # noqa: E402

BASE = dict(n=2, b=4.0, frac=0.9, wdb=0.5, M=2.0, a=2.0, wK=10.0, wPK=1.5, cm2=0.8)


def malla(L=3000.0):
    return np.concatenate([np.linspace(1e-4, 30, 6001), np.geomspace(30, L, 2001)[1:]])


def evalua(item):
    nombre, x, w0 = item
    g = B.construir(x)
    ls = malla()
    acc = {k: [] for k in ("K", "P", "U", "N", "Phi1", "r", "r1", "Gc")}
    for i in range(0, len(ls), 2000):
        o = g.monopolo(ls[i:i + 2000])
        V = o["V"]
        for k in ("K", "P", "U"):
            acc[k].append(o[k].deriv(0))
        for k, n in (("N", "N"), ("Phi1", "Phi1"), ("r", "r0"), ("r1", "r1"), ("Gc", "Gc")):
            acc[k].append(V[n].deriv(0))
    f = {k: np.concatenate(v) for k, v in acc.items()}
    N, Np = f["N"], f["N"] * f["Phi1"]
    EN = simpson(f["P"] * Np ** 2 + f["U"] * N ** 2, x=ls)
    teo = -2 * simpson(N ** 3 * f["Gc"] * (f["r1"] / f["r"] - f["Phi1"]) ** 2, x=ls)
    borde = -2 * (N ** 3 * f["r1"] * f["Gc"] / f["r"])[-1]
    KN = simpson(f["K"] * N ** 2, x=ls)
    return dict(nombre=nombre, x=x, N0=float(g.N0), EN=float(EN), teo=float(teo), borde_L=float(borde),
                qN=float(EN / KN), q1=float(simpson(f["U"], x=ls) / simpson(f["K"], x=ls)), w0=w0,
                minK=float(f["K"].min()), minP=float(f["P"].min()))


def casos():
    out = []
    gb = json.load(open(os.path.join(F2F, "resultados", "gargantas_bH.json")))
    m = gb["masiva"]
    out.append(("referencia", BASE, min(m["w_par"] + m["w_impar"])))
    for M, r in zip([2.0, 1.0, 0.5, 0.25, 0.1], gb["campo_debil"]):
        if "n_par" in r:
            out.append((f"debil M={M}", dict(BASE, M=M, phi1=math.sqrt(M / 2)), min(r["w_par"] + r["w_impar"])))
    for e in gb["escape"]:
        if "n_par" in e:
            out.append((f"escape M={e['M']} cm2={e['cm2']}", dict(BASE, M=e["M"], a=e["a"], wK=e["wK"], cm2=e["cm2"]),
                        min(e["w_par"] + e["w_impar"])))
    ba = json.load(open(os.path.join(F2F, "resultados", "busqueda_aleatoria.json")))
    k = 0
    for r in ba["resultados"]:
        if "n_par" in r:
            k += 1
            out.append((f"aleatoria {k}", r["x"], min(r["w_par"] + r["w_impar"])))
    return out


if __name__ == "__main__":
    t0 = time.time()
    cs = casos()
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(evalua, cs))
    print(f"{len(res)} gargantas (K, P > 0) de la Fase 2f")
    print(f"{'caso':28s} {'N0':>6s} {'E[N]':>10s} {'-2intN3G s2':>12s} {'borde':>9s} {'qN':>10s} {'q1':>10s} {'w0':>10s}")
    for r in res:
        print(f"{r['nombre']:28s} {r['N0']:6.3f} {r['EN']:+10.5f} {r['teo']:+12.5f} {r['borde_L']:+9.5f} "
              f"{r['qN']:+10.5f} {r['q1']:+10.5f} {r['w0']:+10.5f}")
    print(f"max |E[N] - teo - borde| = {max(abs(r['EN'] - r['teo'] - r['borde_L']) for r in res):.2e}")
    print(f"todas con qN < 0: {all(r['qN'] < 0 for r in res)};  w0 <= qN en todas: {all(r['w0'] <= r['qN'] + 1e-9 for r in res)};"
          f"  q1 > 0 en {sum(r['q1'] > 0 for r in res)}")
    print(f"{time.time() - t0:.0f} s")
    json.dump(res, open(os.path.join(AQUI, "resultados", "verificar_N.json"), "w"), indent=1, default=float)
