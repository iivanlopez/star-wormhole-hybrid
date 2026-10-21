# -*- coding: utf-8 -*-
"""
Verificacion independiente de cada raiz real de los barridos robustos con la funcion analitica Dc (estabilidad_F.Dc:
una sola integracion, sin renormalizar, dato inicial por menores).  La D renormalizada de los barridos solo conserva el
signo del coeficiente creciente; Dc conserva tambien su modulo, asi que un cero simple de verdad se ve como un cambio de
signo de Dc con |Dc(w*)| << |Dc(w*(1 +- delta))|.
Lee resultados/estabilidad_robusta.json, estabilidad_gauss.json y trayectorias_robustas.json.
Salida: resultados/verificar_raices.json.   Uso: python verificar_raices.py
"""
import sys
import json
import os
from concurrent.futures import ProcessPoolExecutor
from fondo_F import FondoF
from estabilidad_F import CapaF, Dc

DELTA = 1e-4


def tarea(args):
    nombre, tipo, al, p0, u, kw, z, ws = args
    f = FondoF(tipo, al, p0, u, kw=kw).resolver()
    c = CapaF(f, zeta=z)
    out = []
    for w in ws:
        dm, d0, dp = (Dc(f, c, w * x).real for x in (1 - DELTA, 1.0, 1 + DELTA))
        out.append(dict(nombre=nombre, zeta=z, w=w, cambio_signo=bool(dm * dp < 0),
                        cociente=abs(d0) / max(abs(dm), abs(dp))))
    return out


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    casos = [(o["nombre"], o["tipo"], o["alpha"], o["phi0"], o["phi1"], {}, o["zeta"], o["w"])
             for o in json.load(open("resultados/estabilidad_robusta.json", encoding="utf-8")) if o["w"]]
    if os.path.exists("resultados/estabilidad_gauss.json"):
        for o in json.load(open("resultados/estabilidad_gauss.json", encoding="utf-8")):
            if o["hiperbolico"]:
                casos += [(f"gauss {o['alpha_eta']:.4f}", "gauss", o["alpha"], o["phi0"], o["phi1"], {"gamma": 1.0}, z,
                           o[f"w_zeta{z:g}"]) for z in (1.0, 0.0) if o[f"w_zeta{z:g}"]]
    if os.path.exists("resultados/trayectorias_robustas.json"):
        casos += [(o["nombre"], o["tipo"], o["alpha"], o["phi0"], o["phi1"], {}, o["zeta"], o["w"])
                  for o in json.load(open("resultados/trayectorias_robustas.json", encoding="utf-8")) if o["w"]]
    res = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for lst in ex.map(tarea, casos):
            res += lst
    malas = [o for o in res if not o["cambio_signo"] or o["cociente"] > 1e-2]
    print(f"{len(res)} raices verificadas con Dc: {len(res) - len(malas)} ceros simples "
          f"(cambio de signo y |Dc(w*)| < 1e-2 |Dc(w*(1 +- {DELTA:g}))|); cociente maximo "
          f"{max(o['cociente'] for o in res):.1e}")
    for o in malas:
        print("  NO VERIFICADA:", o)
    json.dump(res, open("resultados/verificar_raices.json", "w", encoding="utf-8"), indent=1)
