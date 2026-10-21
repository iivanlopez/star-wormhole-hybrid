# -*- coding: utf-8 -*-
"""
Ley empirica del modo rapido: con el reloj de l = 0, su tasa de crecimiento es la aceleracion propia a = Phi'(0+) de los
observadores estaticos en la garganta (o el ecuador),
    w = -C a^2,   C ~ 1 (zeta = 0),   C -> 2 (zeta = 1) en el limite de agujero negro;
equivalentemente omega_inf = -i sqrt(C) kappa_t, con kappa_t = e^{-Phi_inf} Phi'(0+) ("gravedad superficial" de la
garganta), que tiende a 1/(2 r0) en el limite de agujero negro.
Lee resultados/estabilidad_robusta.json, estabilidad_gauss.json y trayectorias_robustas.json.
Salida: resultados/ley_aceleracion.json.   Uso: python ley_aceleracion.py
"""
import sys
import json
import math
import os
from fondo_F import FondoF

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def aceleracion(tipo, al, p0, u, kw=None):
    return FondoF(tipo, al, p0, u, kw=kw or {}).y0()[3]


filas = []
for o in json.load(open("resultados/estabilidad_robusta.json", encoding="utf-8")):
    fam = "KKK" if o["nombre"].startswith("KKK") else ("doble" if o["nombre"].startswith("doble") else "cuadratico")
    if fam == "cuadratico" and o["ecuador"]:
        fam = "doble"
    a = aceleracion(o["tipo"], o["alpha"], o["phi0"], o["phi1"])
    filas.append(dict(familia=fam, nombre=o["nombre"], zeta=o["zeta"], a=a, e2Phi=o["e2Phi"],
                      kappa_t=a / math.sqrt(o["e2Phi"]), C=[-w / a**2 for w in sorted(o["w"])]))
if os.path.exists("resultados/estabilidad_gauss.json"):
    for o in json.load(open("resultados/estabilidad_gauss.json", encoding="utf-8")):
        if not o["hiperbolico"]:
            continue
        a = aceleracion("gauss", o["alpha"], o["phi0"], o["phi1"], {"gamma": 1.0})
        for z in (1.0, 0.0):
            filas.append(dict(familia="gauss" + (" doble" if o["ecuador"] else ""), nombre=f"gauss {o['alpha_eta']:.4f}",
                              zeta=z, a=a, e2Phi=math.exp(2 * o["Phi_inf"]), kappa_t=a * math.exp(-o["Phi_inf"]),
                              C=[-w / a**2 for w in sorted(o[f"w_zeta{z:g}"])]))

print("zeta = 1: C = -w/a^2 del modo mas rapido.  zeta = 0: C del modo mas proximo a C = 1 (y C de los demas modos):")
resumen = {}
for fam in ("KKK", "cuadratico", "doble", "gauss", "gauss doble"):
    for z in (1.0, 0.0):
        sub = [f for f in filas if f["familia"] == fam and f["zeta"] == z]
        if not sub:
            continue
        if z == 1.0:
            cs = [max(f["C"]) for f in sub if f["C"]]
            otros = [c for f in sub for c in f["C"] if c != max(f["C"])]
        else:
            cs = [min(f["C"], key=lambda c: abs(c - 1)) for f in sub if f["C"]]
            otros = [c for f in sub for c in f["C"] if c != min(f["C"], key=lambda c: abs(c - 1))]
        resumen[f"{fam} {z:g}"] = dict(n=len(sub), con_modo=len(cs), C=[min(cs), max(cs)] if cs else None,
                                      otros=[min(otros), max(otros)] if otros else None, n_otros=len(otros),
                                      kappa_t=[min(f["kappa_t"] for f in sub), max(f["kappa_t"] for f in sub)],
                                      a=[min(f["a"] for f in sub), max(f["a"] for f in sub)])
        print(f"  {fam:12s} zeta={z:g}: {len(cs)}/{len(sub)} con modo; C en "
              + (f"[{min(cs):.4f}, {max(cs):.4f}]" if cs else "-")
              + (f"; otros {len(otros)} modos con C en [{min(otros):.3g}, {max(otros):.3g}]" if otros else "")
              + f";  kappa_t r0 en [{min(f['kappa_t'] for f in sub):.4f}, {max(f['kappa_t'] for f in sub):.4f}]"
              + f";  a r0 en [{min(f['a'] for f in sub):.3g}, {max(f['a'] for f in sub):.3g}]")
out = dict(filas=filas, resumen=resumen)
if os.path.exists("resultados/trayectorias_robustas.json"):
    tr = json.load(open("resultados/trayectorias_robustas.json", encoding="utf-8"))
    print("\nTrayectorias: C(zeta) del modo rapido frente a 1/(1 - zeta/2):")
    out["trayectorias"] = []
    for nom in sorted({o["nombre"] for o in tr if not o["nombre"].startswith("KKK 0.0029")}):
        sub = sorted((o for o in tr if o["nombre"] == nom), key=lambda o: o["zeta"])
        a = aceleracion(sub[0]["tipo"], sub[0]["alpha"], sub[0]["phi0"], sub[0]["phi1"])
        e2 = math.exp(2 * sub[0]["Phi_inf"])
        print(f"  {nom} (a r0 = {a:.4g}, e^(2 Phi_inf) = {e2:.3g}):")
        for o in sub:
            cs = [-w / a**2 for w in sorted(o["w"])]
            ref = 1 / (1 - o["zeta"] / 2) if o["zeta"] < 2 else float("nan")
            print(f"    zeta={o['zeta']:+.2f}: C = {[round(c, 4) for c in cs]}   1/(1-zeta/2) = {ref:.4f}")
            out["trayectorias"].append(dict(nombre=nom, zeta=o["zeta"], a=a, C=cs, ref=ref))
json.dump(out, open("resultados/ley_aceleracion.json", "w", encoding="utf-8"), indent=1)
