# -*- coding: utf-8 -*-
"""Genera las tablas Markdown del informe a partir de resultados/*.json (salida: resultados/tablas.md)."""
import sys
import json
import math

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
res = json.load(open("resultados/h1.json", encoding="utf-8"))
rap = json.load(open("resultados/modo_rapido.json", encoding="utf-8"))
form = {(o["xi0"], o["rho"]): o for o in json.load(open("resultados/formula_rapida.json", encoding="utf-8"))}
C_KMS = 3.0e5
tau = lambda w: 1.0 / (C_KMS * math.sqrt(abs(w))) * 1e6   # microsegundos


def pot(x):
    """3e13 -> 3×10^13 (LaTeX)."""
    e = math.floor(math.log10(x) + 1e-9)
    m = x / 10**e
    return f"$10^{{{e}}}$" if abs(m - 1.0) < 1e-9 else f"${m:g}\\times10^{{{e}}}$"


def num(x, sig=4):
    """sig cifras significativas; notacion cientifica (LaTeX) por debajo de 10^-3."""
    a = abs(x)
    if a >= 1e-3:
        return f"{x:#.{sig}g}".rstrip(".").replace("-", "−")
    e = math.floor(math.log10(a))
    return f"${x / 10**e:.{sig - 1}f}\\times10^{{{e}}}$"


L = ["<!-- TABLA_FANTASMA -->",
     "| $\\xi_0$ | $r_0$ [km] | $\\rho_{bc}$ [g/cm³] | $M$ [$M_\\odot$] | $R/r_0$ | $\\omega_0^2$ [km⁻²] "
     "| $\\omega_0^2r_0^2e^{2\\Phi_\\infty}$ (reloj de la garganta) | $\\tau$ [μs] |",
     "|---|---|---|---|---|---|---|---|"]
for o in sorted((o for o in res if o.get("fantasma")), key=lambda o: (o["xi0"], o["rho"])):
    w = o["fantasma"][0]
    L.append(f"| {o['xi0']} | {o['r0']:.3f} | {pot(o['rho'])} | {o['M']:.4f} | {o['R'] / o['r0']:.1f} | {num(w)} "
             f"| {num(w * o['r0']**2 * math.exp(2 * o['Phi_inf']))} | {tau(w):.2f} |")

L += ["", "<!-- TABLA_LENTOS -->",
      "| $\\rho_{bc}$ [g/cm³] | $M$ híbrido [$M_\\odot$] | $M$ EN [$M_\\odot$] | $\\omega^2_{\\rm fluido}$ Modelo L [km⁻²] "
      "| $\\omega_0^2$ EN [km⁻²] | $f$ Modelo L [Hz] |",
      "|---|---|---|---|---|---|"]
for o in sorted((o for o in res if o.get("xi0") == 0.1), key=lambda o: o["rho"]):
    lento = min(o["ligado"], key=abs)
    fz = math.sqrt(lento) * C_KMS / (2 * math.pi) if lento > 0 else float("nan")
    L.append(f"| {pot(o['rho'])} | {o['M']:.4f} | {o['M_ns']:.4f} | {num(lento)} | {num(o['ns'])} | "
             + (f"{fz:.0f}" if lento > 0 else "inestable") + " |")

L += ["", "<!-- TABLA_RAPIDO -->",
      "| $\\xi_0$ | $\\rho_{bc}$ [g/cm³] | $\\Sigma$ [km⁻¹] | $e^{2\\Phi_\\infty}$ | $\\omega^2$ numérico [km⁻²] "
      "| fórmula (5) [km⁻²] | cociente | $\\tau$ [μs] | $c\\tau/r_0$ |",
      "|---|---|---|---|---|---|---|---|---|"]
for o in sorted((o for o in rap if o["modelo"] == "ligado"), key=lambda o: (o["xi0"], o["rho"])):
    fo = form[(o["xi0"], o["rho"])]
    w = min(o["raices"])
    L.append(f"| {o['xi0']} | {pot(o['rho'])} | {num(fo['Sigma'], 3)} | {math.exp(2 * fo['Phi_inf']):.3f} | {num(w)} "
             f"| {num(fo['w_formula_fis'])} | {w / fo['w_formula_fis']:.3f} | {tau(w):.3f} "
             f"| {1.0 / math.sqrt(abs(w)) / o['r0']:.3f} |")
try:
    imp = json.load(open("resultados/impar_garganta.json", encoding="utf-8"))["secuencia"]
except FileNotFoundError:
    imp = []
par = {(o["xi0"], o["rho"]): o for o in res}
L += ["", "<!-- TABLA_IMPAR -->",
      "| $\\xi_0$ | $\\rho_{bc}$ [g/cm³] | $\\omega_1^2$ impar [km⁻²] | $\\tau_1$ [μs] "
      "| $\\lvert w_1\\rvert/4\\pi(\\varepsilon_0+3P_0)$ | $\\omega_2^2$ impar [km⁻²] | $\\omega^2$ par del fluido [km⁻²] |",
      "|---|---|---|---|---|---|---|"]
for o in sorted(imp, key=lambda o: (o["xi0"], o["rho"])):
    w1 = o["impar"][0]
    w2 = num(o["impar"][1]) if len(o["impar"]) > 1 else "—"
    lentos = [w for w in par[(o["xi0"], o["rho"])]["ligado"] if abs(w) < 3e-3]
    wp = num(min(lentos, key=abs)) if lentos else "—"
    L.append(f"| {o['xi0']} | {pot(o['rho'])} | {num(w1)} | {tau(w1):.0f} | {o['cociente']:.2f} | {w2} | {wp} |")
try:
    con = json.load(open("resultados/contorno.json", encoding="utf-8"))
except FileNotFoundError:
    con = []
L += ["", "<!-- TABLA_CONTORNO -->",
      "| $\\xi_0$ | $\\rho_{bc}$ [g/cm³] | ${\\rm Re}\\,w$ (reloj de la garganta) | $\\lvert{\\rm Im}\\,w\\rvert<$ | $d$ "
      "| $N=N_{\\rm garganta}+N_{\\rm lenta}+2N_{\\rm franja}$ | raíces reales | puntos |",
      "|---|---|---|---|---|---|---|---|"]
for c in con:
    n = [round(c[s]["N"]) + 0 for s in ("garganta", "lenta", "franja")]
    fiable = min(n) >= 0 and all(abs(c[s]["N"] - round(c[s]["N"])) < 0.01 for s in ("garganta", "lenta", "franja"))
    tot = f"**{c['N']:.3f}**" if fiable else f"{c['N']:.3f} (no concluyente: aliasing)"
    L.append(f"| {c['xi0']} | {pot(c['rho'])} | $[{-c['A']:.3g},\\ {c['B']:.2g}]$ | {c['C']:.3g} | {c['d']:.2g} "
             f"| {n[0]} + {n[1]} + 2·{n[2]} = {tot} | {c['reales']} | {c['puntos']} |")
open("resultados/tablas.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
print("\n".join(L))
