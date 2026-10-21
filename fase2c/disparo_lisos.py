# -*- coding: utf-8 -*-
"""
Disparo dirigido de agujeros lisos.  Un agujero liso (alpha, phi0, phi1) necesita que la mitad "repulsiva" (a < 0: lapso
decreciente hacia fuera) sea AF.  Si existiese una region AF de mitades repulsivas, al acercarse a su frontera desde la
region singular la distancia al final singular l_s divergeria (o cambiaria el tipo de final).  Se usan dos diagnosticos:
  1. a lo largo de rectas de la malla, la distancia l_s (normalizada por r0 y por 1/|a|) y el minimo de N = e^Phi;
  2. biseccion en phi1 entre puntos vecinos con distinto tipo de final de la mitad repulsiva, hasta |dphi1| < 1e-10,
     registrando l_s y el tipo de final a cada lado de la frontera: una frontera con l_s acotado separa dos finales
     singulares y no esconde soluciones regulares.
Entrada: resultados/mapa_lados.json, mapa_gauss.json, mapa_lineal.json.  Salida: resultados/disparo_lisos.json.   Uso: python disparo_lisos.py
"""
import os
import sys
import json
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fondo_doble import Lado      # noqa: E402

L_MAX = 400.0


def mitad(tipo, al, p0, u):
    L = Lado(tipo, al, p0, u, l_max=L_MAX, rtol=1e-10).resolver()
    N = np.exp(L.sol.y[2])
    return dict(fin=L.fin, l_end=float(L.l_end), a=float(L.a), Nmin=float(N.min()), cruces=len(L.cruces))


def biseccion(args):
    tipo, al, p0, u_a, u_b = args
    A, B = mitad(tipo, al, p0, u_a), mitad(tipo, al, p0, u_b)
    traza = []
    for _ in range(60):
        if abs(u_b - u_a) < 1e-10 * max(1.0, abs(u_a)):
            break
        um = 0.5 * (u_a + u_b)
        Mm = mitad(tipo, al, p0, um)
        traza.append((um, Mm["fin"], Mm["l_end"]))
        if Mm["fin"] == A["fin"]:
            u_a, A = um, Mm
        else:
            u_b, B = um, Mm
        if Mm["fin"] == "AF":
            break
    return dict(tipo=tipo, alpha=al, phi0=p0, u_a=u_a, u_b=u_b, fin_a=A["fin"], fin_b=B["fin"], l_a=A["l_end"],
                l_b=B["l_end"], a_a=A["a"], a_b=B["a"], Nmin_a=A["Nmin"], Nmin_b=B["Nmin"],
                l_max_traza=max([t[2] for t in traza], default=0.0), af=any(t[1] == "AF" for t in traza))


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
    # fronteras entre tipos de final de la mitad repulsiva, en phi1 (mismo signo de phi1)
    por_linea = {}
    for o in rep:
        por_linea.setdefault((o["tipo"], o["alpha"], o["phi0"], np.sign(o["phi1"])), []).append(o)
    tareas = []
    for k, lst in por_linea.items():
        lst.sort(key=lambda o: o["phi1"])
        for o1, o2 in zip(lst[:-1], lst[1:]):
            if o1["fin"] != o2["fin"]:
                tareas.append((k[0], k[1], k[2], o1["phi1"], o2["phi1"]))
    print(f"mitades repulsivas: {len(rep)}; fronteras de tipo de final en phi1: {len(tareas)}", flush=True)
    with ProcessPoolExecutor(max_workers=8) as ex:
        bis = list(ex.map(biseccion, tareas, chunksize=2))
    n_af = sum(b["af"] for b in bis)
    print(f"biseccion: {len(bis)} fronteras; con alguna mitad AF en la traza: {n_af}")
    tipos = {}
    for b in bis:
        key = tuple(sorted((b["fin_a"], b["fin_b"])))
        tipos.setdefault(key, []).append(max(b["l_a"], b["l_b"]))
    for k, v in sorted(tipos.items()):
        print(f"   {k[0]:>9s} | {k[1]:<9s}: {len(v):4d} fronteras, l_s maximo en la frontera = {max(v):.4g} r0")
    # l_s de todas las mitades repulsivas singulares
    ls = np.array([o["l_end"] for o in rep if o["fin"] != "AF"])
    las = np.array([o["l_end"] * abs(o["a"]) for o in rep if o["fin"] != "AF"])
    print(f"mitades repulsivas singulares: {len(ls)}; l_s/r0: max {ls.max():.4g}, percentil 99 {np.percentile(ls, 99):.4g}; "
          f"l_s|a|: max {las.max():.4g}, mediana {np.median(las):.4g}")
    json.dump(dict(biseccion=bis, ls_max=float(ls.max()), ls_a_max=float(las.max())),
              open("resultados/disparo_lisos.json", "w", encoding="utf-8"), indent=1)
    print(f"tiempo: {time.time() - t0:.0f} s")
