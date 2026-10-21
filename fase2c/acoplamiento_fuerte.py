# -*- coding: utf-8 -*-
"""
Limite de acoplamiento fuerte (|F'(phi0)| >> r0^2) de las mitades repulsivas, fuera de la malla de mapa_lados.py.
En ese limite a = Phi'(0) -> 0 y el colapso inmediato del lapso no se produce; se comprueba que la mitad tampoco es
asintoticamente plana: r' se satura en un valor < 1 (geometria conica, masa de Misner-Sharp m = r (1 - r'^2)/2
proporcional a r) y el lapso acaba anulandose a una distancia l_s proporcional a |alpha|.
Salida: resultados/acoplamiento_fuerte.json.   Uso: python acoplamiento_fuerte.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fondo_doble import Lado      # noqa: E402

if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    casos = []
    for al in (-10.0, -100.0, -1000.0, -1e4, 10.0, 100.0, 1000.0, 1e4):
        for u in (0.1, 1.0, 5.0, 20.0):
            casos.append(("exp", al, 0.0, u))
            casos.append(("exp", al, 0.0, -u))
    for al in (-10.0, -100.0, -1000.0, 10.0, 100.0, 1000.0):
        for p0 in (0.1, 1.0):
            for u in (0.1, 1.0, 5.0, 20.0):
                casos.append(("cuadratico", al, p0, u))
                casos.append(("cuadratico", al, p0, -u))
    for al in (10.0, 100.0, 1000.0, 1e4):
        for u in (0.1, 1.0, 5.0, 20.0):
            casos.append(("lineal", al, 0.0, u))
            casos.append(("lineal", al, 0.0, -u))
    res = []
    for c in casos:
        L = Lado(*c, l_max=2e6, rtol=1e-10).resolver()
        if L.a > 0:
            continue                      # solo mitades repulsivas
        ys = L.sol.y
        r, rp = ys[0], ys[1]
        m = r * (1 - rp**2) / 2
        k = len(r) // 2
        o = dict(tipo=c[0], alpha=c[1], phi0=c[2], phi1=c[3], a=float(L.a), fin=L.fin, l_end=float(L.l_end),
                 r_fin=float(r[-1]), rp_max=float(rp.max()), rp_mitad=float(rp[k]), m_sobre_r_fin=float(m[-1] / r[-1]),
                 ls_a=float(L.l_end * abs(L.a)), ls_sobre_alpha=float(L.l_end / abs(c[1])))
        res.append(o)
        print(f"{c[0]:10s} alpha={c[1]:+9.3g} phi0={c[2]:.2g} phi1={c[3]:+6.3g}: a={L.a:+.3e} fin={L.fin:9s} l_s={L.l_end:.4g} "
              f"(l_s|a| = {o['ls_a']:.3g}, l_s/|alpha| = {o['ls_sobre_alpha']:.3g}) r_fin={r[-1]:.4g} r'_max={rp.max():.4f} "
              f"m/r(fin) = {o['m_sobre_r_fin']:+.3f}", flush=True)
    os.makedirs("resultados", exist_ok=True)
    json.dump(res, open("resultados/acoplamiento_fuerte.json", "w", encoding="utf-8"), indent=1)
    print(f"mitades repulsivas: {len(res)}; AF: {sum(o['fin'] == 'AF' for o in res)}; r'_max > 0.9: "
          f"{sum(o['rp_max'] > 0.9 for o in res)}")
