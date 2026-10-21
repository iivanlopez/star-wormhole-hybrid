# -*- coding: utf-8 -*-
"""
Analisis de los mapas de mitades lisas (mapa_lados.py): resultados/mapa_lados.json (cuadratico, exp),
mapa_gauss.json y mapa_lineal.json.
  - recuento de finales por acoplamiento, por signo de a = Phi'(0) en la direccion de la mitad y por signo de F''(phi0);
  - mitades AF con a < 0 (mitad repulsiva) y agujeros lisos (las dos mitades AF);
  - ley del colapso: l_s |a| por intervalos de |a| r0 (mitades repulsivas con horizonte singular);
  - r' maximo al final de las mitades repulsivas (cercania a la planitud asintotica);
  - fraccion de finales 'phi'/'r0' con centro tipo ecuador (r''(0) < 0).
Salida: resultados/analisis_mapa.json.   Uso: python analisis_mapa.py
"""
import os
import sys
import json
import math
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def F2(o):
    al, p = o["alpha"], o["phi0"]
    return {"cuadratico": 2 * al, "exp": al * math.exp(-p), "gauss": al * (4 * p * p - 2) * math.exp(-p * p),
            "lineal": 0.0}[o["tipo"]]


res = []
for f in ("resultados/mapa_lados.json", "resultados/mapa_gauss.json", "resultados/mapa_lineal.json"):
    if os.path.exists(f):
        res += json.load(open(f, encoding="utf-8"))
salida = {}
for tipo in ("cuadratico", "exp", "gauss", "lineal"):
    sub = [o for o in res if o["tipo"] == tipo]
    if not sub:
        continue
    s = dict(n=len(sub))
    print(f"\n=== {tipo}: {len(sub)} mitades")
    for signo, nombre in ((1, "atr"), (-1, "rep")):
        for sF2, nF in ((1, "F2>0"), (-1, "F2<0"), (0, "F2=0")):
            ss = [o for o in sub if np.sign(o["a"]) == signo and np.sign(F2(o)) == sF2]
            if not ss:
                continue
            cuenta = {}
            for o in ss:
                cuenta[o["fin"]] = cuenta.get(o["fin"], 0) + 1
            s[f"{nombre}_{nF}"] = dict(n=len(ss), finales=dict(sorted(cuenta.items())))
            print(f"   {'a > 0' if signo > 0 else 'a < 0'}, {nF}: {len(ss)} -> {dict(sorted(cuenta.items()))}")
    rep = [o for o in sub if o["a"] < 0]
    s["rep_af"] = sum(o["fin"] == "AF" for o in rep)
    idx = {(o["alpha"], o["phi0"], o["phi1"]): o for o in sub}
    s["lisos"] = sum(1 for k, o in idx.items() if k[2] > 0 and o["fin"] == "AF"
                     and idx.get((k[0], k[1], -k[2]), {}).get("fin") == "AF")
    print(f"   mitades repulsivas AF: {s['rep_af']};  agujeros lisos: {s['lisos']}")
    # ley del colapso
    hz = [o for o in rep if o["fin"] == "horizonte"]
    x = np.array([abs(o["a"]) for o in hz])
    y = np.array([o["l_end"] * abs(o["a"]) for o in hz])
    ley = []
    for lo, hi in ((0, 1), (1, 10), (10, 30), (30, 100), (100, 1e9)):
        m = (x >= lo) & (x < hi)
        if m.any():
            fila = dict(intervalo=[lo, hi], n=int(m.sum()), mediana=float(np.median(y[m])),
                        p10=float(np.percentile(y[m], 10)), p90=float(np.percentile(y[m], 90)))
            ley.append(fila)
            print(f"   horizonte, |a| r0 en [{lo:g}, {hi:g}): {fila['n']:5d}; l_s|a| mediana {fila['mediana']:.4f}, "
                  f"p10-p90 {fila['p10']:.3f}-{fila['p90']:.3f}")
    s["ley"] = ley
    rpf = np.array([o.get("rp_fin", np.nan) for o in rep])
    s["rp_fin_max"] = float(np.nanmax(rpf))
    s["ls_max"] = float(max(o["l_end"] for o in rep))
    pr = [o for o in rep if o["fin"] in ("phi", "r0")]
    s["patada_n"] = len(pr)
    s["patada_ecuador"] = sum(o["rpp0"] < 0 for o in pr)
    print(f"   r'_fin maximo (repulsivas): {s['rp_fin_max']:.4f};  l_fin maximo: {s['ls_max']:.4g} r0;  "
          f"finales phi/r0: {len(pr)} (con r''(0) < 0: {s['patada_ecuador']})")
    att_af = [o for o in sub if o["a"] > 0 and o["fin"] == "AF"]
    if att_af:
        aa = np.array([o["a"] for o in att_af])
        s["atr_af_a"] = [float(aa.min()), float(np.median(aa)), float(aa.max())]
        print(f"   mitades atractivas AF: {len(att_af)}, a r0 en [{aa.min():.3g}, {aa.max():.3g}] (mediana {np.median(aa):.3g})")
    salida[tipo] = s
tot_rep = sum(1 for o in res if o["a"] < 0)
tot_af = sum(1 for o in res if o["a"] < 0 and o["fin"] == "AF")
print(f"\nTOTAL: {len(res)} mitades, {tot_rep} repulsivas, {tot_af} repulsivas AF, "
      f"{sum(v['lisos'] for v in salida.values())} agujeros lisos")
salida["total"] = dict(mitades=len(res), repulsivas=tot_rep, repulsivas_AF=tot_af)
json.dump(salida, open("resultados/analisis_mapa.json", "w", encoding="utf-8"), indent=1)
