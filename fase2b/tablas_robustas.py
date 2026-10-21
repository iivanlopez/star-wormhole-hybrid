# -*- coding: utf-8 -*-
"""
Tablas Markdown del informe (Parte I §6 y Parte II §11) a partir de resultados/estabilidad_robusta.json,
hiperbolicidad.json, familia_cuadratica.json y, si existen, trayectorias_robustas.json y estabilidad_gauss.json.
Salida: resultados/tablas_robustas.md.   Uso: python tablas_robustas.py
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
rob = json.load(open("resultados/estabilidad_robusta.json", encoding="utf-8"))
hip = {o["nombre"]: o for o in json.load(open("resultados/hiperbolicidad.json", encoding="utf-8"))}
fam = json.load(open("resultados/familia_cuadratica.json", encoding="utf-8"))
por = {}
for o in rob:
    por.setdefault(o["nombre"], {})[o["zeta"]] = o


def num(x):
    if abs(x) >= 0.01:
        return f"{x:.3f}".replace("-", "−")
    e = math.floor(math.log10(abs(x)))
    m = x / 10**e
    return f"${m:.1f}\\times10^{{{e}}}$"


def sg(x, fmt):
    return format(x, fmt).replace("-", "−")


def lista(o):
    if o is None:
        return "?"
    return ", ".join(num(x) for x in sorted(o["winf"])) or "sin raíz real"


def pot(x):
    if x < 100:
        return f"{x:.3g}"
    if x < 1e4:
        return f"{float(f'{x:.3g}'):.0f}"
    e = math.floor(math.log10(x))
    return f"${x / 10**e:.1f}\\times10^{{{e}}}$"


L = ["<!-- TABLA_KKK_ROBUSTA -->",
     "| $\\alpha/r_0^2$ | $f_0$ | $e^{2\\Phi_\\infty}$ | $\\omega^2_\\infty r_0^2$, $\\zeta=1$ (KKK) | $\\omega^2_\\infty r_0^2$, $\\zeta=0$ |",
     "|---|---|---|---|---|"]
kkk = sorted((v for k, v in por.items() if k.startswith("KKK")), key=lambda d: (d[1.0]["alpha_kkk"], d[1.0]["f0"]))
for d in kkk:
    o = d[1.0]
    L.append(f"| {o['alpha_kkk']:.4f} | {o['f0']:.4f} | {pot(o['e2Phi'])} | {lista(d.get(1.0))} | {lista(d.get(0.0))} |")

L += ["", "<!-- TABLA_UNA -->",
      "| $\\alpha/\\eta_0^2$ | $\\phi_0$ | $f_0$ | $e^{2\\Phi_\\infty}$ | $\\omega^2_\\infty r_0^2$, $\\zeta=1$ (Antoniou) "
      "| $\\omega^2_\\infty r_0^2$, $\\zeta=0$ |",
      "|---|---|---|---|---|---|"]
for o in sorted((o for o in fam if not o["ecuador"]), key=lambda o: o["alpha_eta"]):
    nom = f"cuadratico {o['alpha_eta']:.4f}"
    h = hip[nom]
    e2 = math.exp(2 * o["Phi_inf"])
    if h["hiperbolico"]:
        d = por.get(nom, {})
        L.append(f"| {o['alpha_eta']:.3f} | {o['phi0']:g} | {sg(o['f0'], '.3f')} | {pot(e2)} | {lista(d.get(1.0))} | {lista(d.get(0.0))} |")
    else:
        L.append(f"| {o['alpha_eta']:.3f} | {o['phi0']:g} | {sg(o['f0'], '.3f')} | {pot(e2)} | mal planteado | mal planteado |")

L += ["", "<!-- TABLA_DOBLE -->",
      "| $F$ | $\\alpha/\\eta_0^2$ | $f_0$ | $e^{2\\Phi_\\infty}$ | $\\omega^2_\\infty r_0^2$, $\\zeta=1$ (Antoniou) "
      "| $\\omega^2_\\infty r_0^2$, $\\zeta=0$ |",
      "|---|---|---|---|---|---|"]
for o in sorted((o for o in fam if o["ecuador"]), key=lambda o: o["alpha_eta"]):
    nom = f"cuadratico {o['alpha_eta']:.4f}"
    e2 = math.exp(2 * o["Phi_inf"])
    if hip[nom]["hiperbolico"]:
        d = por.get(nom, {})
        L.append(f"| $\\alpha\\phi^2$ | {o['alpha_eta']:.3f} | {sg(o['f0'], '.3f')} | {pot(e2)} | {lista(d.get(1.0))} | {lista(d.get(0.0))} |")
    else:
        L.append(f"| $\\alpha\\phi^2$ | {o['alpha_eta']:.3f} | {sg(o['f0'], '.3f')} | {pot(e2)} | mal planteado | mal planteado |")
for k, h in sorted(((k, h) for k, h in hip.items() if k.startswith("doble")), key=lambda kh: kh[1]["alpha_eta"]):
    f = FondoF("exp", h["alpha"], 0.0, h["phi1"]).resolver()
    e2 = math.exp(2 * f.Phi_inf)
    if h["hiperbolico"]:
        d = por.get(k, {})
        L.append(f"| $\\alpha e^{{-\\phi}}$ | {h['alpha_eta']:.3f} | {sg(f.f0, '.3f')} | {pot(e2)} | {lista(d.get(1.0))} | {lista(d.get(0.0))} |")
    else:
        L.append(f"| $\\alpha e^{{-\\phi}}$ | {h['alpha_eta']:.3f} | {sg(f.f0, '.3f')} | {pot(e2)} | mal planteado | mal planteado |")

if os.path.exists("resultados/trayectorias_robustas.json"):
    import numpy as np
    tr = json.load(open("resultados/trayectorias_robustas.json", encoding="utf-8"))
    L += ["", "<!-- TABLA_KKK_ZETA (Parte I, §5) -->",
          r"| $\zeta$ | $w^*$ (reloj de la garganta) | $\omega^2_\infty r_0^2$ |", "|---|---|---|"]
    for o in sorted((o for o in tr if o["nombre"] == "KKK referencia"), key=lambda o: o["zeta"]):
        ws = sorted(o["w"])
        L.append(f"| {sg(o['zeta'], 'g')} | " + (", ".join(sg(w, '.6g') for w in ws) or "sin raíz real") + " | "
                 + (", ".join(num(w) for w in sorted(o["winf"])) or "—") + " |")
    L += ["", "<!-- TABLA_TRAYECTORIAS -->",
          r"| Fondo | $a\,r_0$ | $e^{2\Phi_\infty}$ | $C(0)$ | $C(1)$ | $\zeta_c$ (ajuste lineal de $1/C$) | último $\zeta$ con modo rápido real |",
          "|---|---|---|---|---|---|---|"]
    for nom, etiqueta in (("KKK referencia", r"KKK, $\alpha/r_0^2=0.028$"), ("cuadratico 0.3087", r"cuadrático, $\alpha/\eta_0^2=0.309$"),
                          ("cuadratico 0.3465", r"cuadrático, $\alpha/\eta_0^2=0.346$")):
        sub = sorted((o for o in tr if o["nombre"] == nom), key=lambda o: o["zeta"])
        a = FondoF(sub[0]["tipo"], sub[0]["alpha"], sub[0]["phi0"], sub[0]["phi1"]).y0()[3]
        e2 = math.exp(2 * sub[0]["Phi_inf"])
        rap = [(o["zeta"], -min(o["w"]) / a**2) for o in sub if o["w"] and -min(o["w"]) / a**2 > 0.3 and o["zeta"] >= 0]
        z_fit = [(z, 1 / c) for z, c in rap if z >= 0.5]
        p = np.polyfit([z for z, _ in z_fit], [y for _, y in z_fit], 1)
        C = dict(rap)
        L.append(f"| {etiqueta} | {a:.1f} | {pot(e2)} | {C.get(0.0, float('nan')):.3f} | {C.get(1.0, float('nan')):.3f} | "
                 f"{-p[1] / p[0]:.2f} | {max(z for z, _ in rap):g} |")
    L += ["", "<!-- C_ZETA -->"]
    for nom in ("KKK referencia", "cuadratico 0.3087", "cuadratico 0.3465"):
        sub = sorted((o for o in tr if o["nombre"] == nom), key=lambda o: o["zeta"])
        a = FondoF(sub[0]["tipo"], sub[0]["alpha"], sub[0]["phi0"], sub[0]["phi1"]).y0()[3]
        L.append(nom + ": " + "; ".join(f"zeta={o['zeta']:g}: C={[round(float(-w / a**2), 4) for w in sorted(o['w'])]}" for o in sub))
    L += ["", "<!-- EXTREMO -->"]
    for o in tr:
        if o["nombre"].startswith("KKK 0.0029"):
            L.append(f"zeta={o['zeta']:g}: w* = {o['w']}  omega_inf^2 r0^2 = {lista(o)}")

if os.path.exists("resultados/estabilidad_gauss.json"):
    ga = json.load(open("resultados/estabilidad_gauss.json", encoding="utf-8"))
    L += ["", "<!-- TABLA_GAUSS -->",
          "| $\\alpha/\\eta_0^2$ | $\\phi_0$ | centro | $f_0$ | $e^{2\\Phi_\\infty}$ | $\\omega^2_\\infty r_0^2$, $\\zeta=1$ "
          "| $\\omega^2_\\infty r_0^2$, $\\zeta=0$ |",
          "|---|---|---|---|---|---|---|"]
    for o in sorted(ga, key=lambda o: (o["ecuador"], o["alpha_eta"])):
        e2 = math.exp(2 * o["Phi_inf"])
        centro = "ecuador" if o["ecuador"] else "garganta"
        if o["hiperbolico"]:
            z1 = ", ".join(num(x) for x in sorted(o["winf_zeta1"])) or "sin raíz real"
            z0 = ", ".join(num(x) for x in sorted(o["winf_zeta0"])) or "sin raíz real"
            L.append(f"| {o['alpha_eta']:.3f} | {o['phi0']:g} | {centro} | {sg(o['f0'], '.3f')} | {pot(e2)} | {z1} | {z0} |")
        else:
            L.append(f"| {o['alpha_eta']:.3f} | {o['phi0']:g} | {centro} | {sg(o['f0'], '.3f')} | {pot(e2)} | mal planteado "
                     f"| mal planteado |")

open("resultados/tablas_robustas.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
cu = {k: d for k, d in por.items() if k.startswith("cuadratico")}
for z in (1.0, 0.0):
    con = [d[z]["winf"] for d in cu.values() if z in d and d[z]["winf"]]
    rap = [min(w) for w in con]
    len_ = [max(w) for w in con if len(w) > 1]
    print(f"\ncuadraticos hiperbolicos: {len(cu)}; zeta={z:g}: con raiz real {len(con)}; rapido en "
          f"[{min(rap):.4f}, {max(rap):.4f}]" + (f"; lento en [{min(len_):.3g}, {max(len_):.3g}] ({len(len_)})" if len_ else ""))
