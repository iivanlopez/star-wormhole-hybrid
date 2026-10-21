# -*- coding: utf-8 -*-
"""
Mapa de regularidad de las MITADES lisas (fondo_doble.Lado): para cada (F, alpha, phi0, phi1) se integra desde la
garganta (r0 = 1, datos lisos) y se clasifica el final: AF, noAF, horizonte (e^Phi -> 0, Phi' -> -inf), r0 (r -> 0),
det (parte principal degenerada), phi (phi' diverge), Phi+ (Phi' -> +inf), fallo.
Un agujero liso (alpha, phi0, phi1) existe si y solo si las mitades (phi0, phi1) y (phi0, -phi1) son ambas AF.
Mallas (phi1 en +-[0.01, 20], 34 valores logaritmicos por signo):
  cuadratico  F = alpha phi^2:     alpha en +-{0.02,...,2}, phi0 en [0.01, 5] (phi0 > 0 basta: phi -> -phi con l -> -l
              identifica (phi0, phi1) con (-phi0, phi1))
  exp         F = alpha e^{-phi}:  phi0 = 0 (simetria de traslacion de phi), alpha en +-[0.01, 30]
  gauss       F = alpha e^{-phi^2}: alpha en +-{0.05,...,5}, phi0 en [0.01, 3]
  lineal      F = alpha phi:       phi0 = 0 (simetria de traslacion), alpha en +-[0.001, 100]
Salida: resultados/mapa_lados.json (cuadratico y exp) o resultados/mapa_<tipo>.json.
Uso: python mapa_lados.py            (cuadratico y exp)
     python mapa_lados.py gauss      (o lineal)
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
U1 = [float(x) for x in np.logspace(-2, np.log10(20.0), 34)]


def tarea(args):
    tipo, al, p0 = args
    out = []
    for u in U1:
        for s in (1.0, -1.0):
            L = Lado(tipo, al, p0, s * u, l_max=L_MAX, rtol=1e-9).resolver()
            out.append(dict(tipo=tipo, alpha=al, phi0=p0, phi1=s * u, a=L.a, rpp0=L.rpp0, fin=L.fin, l_end=L.l_end,
                            cruces=len(L.cruces), M=L.M, D=L.D, phi_inf=L.phi_inf, Phi_inf=L.Phi_inf,
                            det_fin=float(L.det_fin), rp_fin=float(L.sol.y[1, -1])))
    return out


def tareas_de(tipo):
    if tipo == "cuadratico":
        return [("cuadratico", sa * al, float(p0)) for al in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0) for sa in (1.0, -1.0)
                for p0 in np.logspace(-2, np.log10(5.0), 28)]
    if tipo == "exp":
        return [("exp", float(sa * al), 0.0) for al in np.logspace(-2, np.log10(30.0), 36) for sa in (1.0, -1.0)]
    if tipo == "gauss":
        return [("gauss", sa * al, float(p0)) for al in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0) for sa in (1.0, -1.0)
                for p0 in np.logspace(-2, np.log10(3.0), 24)]
    if tipo == "lineal":
        return [("lineal", float(sa * al), 0.0) for al in np.logspace(-3, 2, 36) for sa in (1.0, -1.0)]
    raise ValueError(tipo)


def resumen(res, tipo):
    sub = [o for o in res if o["tipo"] == tipo]
    print(f"{tipo}: {len(sub)} mitades")
    for signo, nombre in ((1, "a > 0 (lapso creciente hacia fuera)"), (-1, "a < 0 (lapso decreciente hacia fuera)")):
        ss = [o for o in sub if np.sign(o["a"]) == signo]
        cuenta = {}
        for o in ss:
            cuenta[o["fin"]] = cuenta.get(o["fin"], 0) + 1
        print(f"   {nombre}: {len(ss)} -> {dict(sorted(cuenta.items()))}")
    idx = {(o["alpha"], o["phi0"], o["phi1"]): o for o in sub}
    lisos = [k for k, o in idx.items() if k[2] > 0 and o["fin"] == "AF" and idx[(k[0], k[1], -k[2])]["fin"] == "AF"]
    print(f"   mitades AF con a < 0: {sum(1 for o in sub if o['a'] < 0 and o['fin'] == 'AF')}")
    print(f"   agujeros lisos (dos mitades AF): {len(lisos)}")
    for k in lisos[:40]:
        o, o2 = idx[k], idx[(k[0], k[1], -k[2])]
        print(f"      alpha={k[0]:+.4g} phi0={k[1]:.4g} phi1={k[2]:.4g}: a={o['a']:+.4g} M+={o['M']:.4g} M-={o2['M']:.4g} "
              f"phi_inf+={o['phi_inf']:.4g} phi_inf-={o2['phi_inf']:.4g}")


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    t0 = time.time()
    tipos = sys.argv[1:] or ["cuadratico", "exp"]
    tareas = [t for tipo in tipos for t in tareas_de(tipo)]
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = [o for lst in ex.map(tarea, tareas, chunksize=4) for o in lst]
    os.makedirs("resultados", exist_ok=True)
    nombre = "resultados/mapa_lados.json" if not sys.argv[1:] else f"resultados/mapa_{'_'.join(tipos)}.json"
    json.dump(res, open(nombre, "w", encoding="utf-8"), indent=0)
    for tipo in tipos:
        resumen(res, tipo)
    print(f"salida: {nombre};  tiempo: {time.time() - t0:.0f} s")
