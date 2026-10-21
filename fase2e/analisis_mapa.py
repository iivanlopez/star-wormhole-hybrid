# -*- coding: utf-8 -*-
"""
Analisis del mapa de existencia (resultados/mapa_potencial.json):
  - recuento de finales por familia y por signo de lambda(0) (garganta bien o mal planteada);
  - distancia del horizonte singular a la garganta (comparada con la de Nariai, pi r0/2);
  - separatrices: cuantas llegan al vacio, con que r', chi, det M y N; cuantas son sanas o AF;
  - ceros regulares de chi en todas las mitades: comprobacion de K = -16 F'^2 (1 - r'^2)^2 chi'/r'^2 y signo de
    lambda en los ascendentes.
Salida: resultados/analisis_mapa.json.   Uso: python analisis_mapa.py
"""
import os
import sys
import json
import math
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
res = json.load(open(os.path.join(AQUI, "resultados", "mapa_potencial.json"), encoding="utf-8"))
out = {"familias": {}}


def clave(o):
    return f"{o['fam']}-{o['F']}-{o['V']}"


familias = sorted({clave(o) for o in res})
tot = dict(mitades=0, sin_vacio=0, seps=0, llegan=0, sanas=0, AF=0)
for k in familias:
    sub = [o for o in res if clave(o) == k]
    pts = [p for o in sub for p in o["puntos"] if p["valido"]]
    fines = {}
    for p in pts:
        fines[p["fin"]] = fines.get(p["fin"], 0) + 1
    bien = [p for p in pts if p["lam0"] > 0]
    hz = np.array([p["l_end"] for p in pts if p["fin"] == "horizonte"])
    seps = [s for o in sub for s in o["separatrices"]]
    lleg = [s for s in seps if s["llegada_lo"]["dphi"] < 1e-3]
    rp = np.array([s["llegada_lo"]["rp"] for s in lleg])
    d = dict(barridos=len(sub), mitades=len(pts), sin_vacio=sum(not p["valido"] for o in sub for p in o["puntos"]),
             finales=dict(sorted(fines.items())), lambda0_pos=len(bien),
             horizonte_l=[float(np.percentile(hz, q)) for q in (0, 10, 50, 90, 100)] if len(hz) else None,
             separatrices=len(seps), llegan=len(lleg), sanas=sum(s["sana"] for s in seps), AF=sum(s["AF"] for s in seps),
             rp_llegada=[float(rp.min()), float(np.median(rp)), float(rp.max())] if len(rp) else None,
             chi_llegada_max=max((s["llegada_lo"]["chi"] for s in lleg), default=None))
    out["familias"][k] = d
    for c in ("mitades", "sin_vacio", "llegan", "sanas", "AF"):
        tot[c] += d[c]
    tot["seps"] += d["separatrices"]
    print(f"{k}: {d['barridos']} barridos, {d['mitades']} mitades (+{d['sin_vacio']} sin vacio), lambda(0) > 0 en "
          f"{d['lambda0_pos']}")
    print(f"   finales {d['finales']}")
    if d["horizonte_l"]:
        print(f"   horizonte a l/r0 (min, p10, mediana, p90, max) = {', '.join(f'{x:.2f}' for x in d['horizonte_l'])}"
              f"   (pi/2 = {math.pi / 2:.3f})")
    print(f"   separatrices {d['separatrices']}, llegan al vacio {d['llegan']}, sanas {d['sanas']}, AF {d['AF']}"
          + (f"; r' al llegar (min, mediana, max) = {', '.join(f'{x:+.3f}' for x in d['rp_llegada'])}"
             if d["rp_llegada"] else ""))
out["total"] = tot
print(f"TOTAL: {tot}")

# horizontes: todos los finales 'horizonte' del mapa
hz = np.array([p["l_end"] for o in res for p in o["puntos"] if p.get("valido") and p["fin"] == "horizonte"])
cerca = np.mean(np.abs(hz - math.pi / 2) < 0.15) if len(hz) else float("nan")
print(f"horizontes: {len(hz)}; l/r0 mediana {np.median(hz):.3f}; fraccion a menos de 0.15 r0 de pi/2: {cerca:.2f}")
out["horizontes"] = dict(n=int(len(hz)), mediana=float(np.median(hz)), frac_cerca_pi2=float(cerca),
                         percentiles=[float(np.percentile(hz, q)) for q in (0, 10, 50, 90, 100)])

# ceros regulares de chi
zs = [z for o in res for p in o["puntos"] if p.get("valido") for z in p["ceros_chi"]]
zs += [z for o in res for s in o["separatrices"] for z in s["ceros_chi"]]
sub = [z for z in zs if z["sube"]]
err = max((abs(z["K"] - z["K_teo"]) / (abs(z["K"]) + abs(z["K_teo"]) + 1e-300) for z in zs), default=0.0)
mnum = max((abs(z["num"]) for z in zs), default=0.0)
kneg = sum(z["K"] < 0 for z in sub)
gpos = sum(z["G"] > 0 for z in sub)
lneg = sum(z["lam"] < 0 for z in sub)
print(f"ceros regulares de chi: {len(zs)} ({len(sub)} ascendentes); identidad max rel {err:.1e}; max |num| {mnum:.1e}")
print(f"   ascendentes: K < 0 en {kneg}, G > 0 en {gpos}, lambda < 0 en {lneg}")
excep = [z for z in sub if not z["lam"] < 0]
for z in excep[:10]:
    print(f"   excepcion: l={z['l']:.4f} r'={z['rp']:+.4f} K={z['K']:+.3e} G={z['G']:+.3e} lambda={z['lam']:+.3e}")
out["ceros_chi"] = dict(n=len(zs), ascendentes=len(sub), err_identidad=err, max_num=mnum, K_neg=kneg, G_pos=gpos,
                        lam_neg=lneg)
# mitades con lambda(0) < 0 (mal planteadas en la propia garganta)
pts = [p for o in res for p in o["puntos"] if p.get("valido")]
print(f"mitades con lambda(0) < 0 (garganta mal planteada): {sum(p['lam0'] < 0 for p in pts)} de {len(pts)}")
out["lambda0_neg"] = sum(p["lam0"] < 0 for p in pts)
json.dump(out, open(os.path.join(AQUI, "resultados", "analisis_mapa.json"), "w", encoding="utf-8"), indent=1)
