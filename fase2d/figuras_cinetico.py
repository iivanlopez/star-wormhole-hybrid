# -*- coding: utf-8 -*-
"""
Figura 4 de la Fase 2d (Parte II): coeficiente cinetico K de agujeros ultraestaticos EsGB + electrodinamica.
  izquierda: forma cerrada de Maxwell (Ellis), K/(A^2 r^2 + 16 F'^2 B |E|) frente a l para varias cargas D y asimetrias c
  derecha:   max_D min_l K/norma frente a la agudeza de la garganta r''(0) (familia hiperbola, c = 0 y 0.3), con los
             modelos Born-Infeld / Euler-Heisenberg y el supremo de la optimizacion global (busqueda_cinetico.py)
Uso: python figuras_cinetico.py
"""
import os
import sys
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cinetico_nled as CN      # noqa: E402

os.makedirs("figuras", exist_ok=True)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))

# ------------------------------------------------------------------ izquierda: Maxwell (forma cerrada)
l = np.linspace(-6, 6, 4001)
for D, c, col, ls in ((0.2, 0.0, "#8d99ae", "-"), (1.0, 0.0, "#e9c46a", "-"), (2.0, 0.0, "#e76f51", "-"),
                      (10.0, 0.0, "#264653", "-"), (2.0, 0.5, "#e76f51", "--")):
    Qe2 = 1 + D * D / 4
    g = c - (Qe2 / 2) * np.arctan(l)
    a2 = (1 - 4 * l * g) ** 2
    neg = (16 * g * g / D**2) * (5 + 2 * Qe2 - 24 * l * g)
    axs[0].plot(l, (a2 - neg) / (a2 + np.abs(neg)), color=col, ls=ls, lw=1.2,
                label=fr"$D={D:g}\,Q$, $c={c:g}$ ($Q_s/Q={D / 2:g}$)")
axs[0].axhline(0, color="k", lw=0.6)
axs[0].set_xlabel(r"$\ell/Q$")
axs[0].set_ylabel(r"$K/(A^2r^2+16F'^2B|E|)$")
axs[0].set_title(r"Maxwell (Ellis), forma cerrada: $K<0$ cerca de la garganta para todo $D$ y $c$", fontsize=9)
axs[0].legend(fontsize=7)

# ------------------------------------------------------------------ derecha: max_D min_l K frente a r''(0)
Ls = np.logspace(-3, 3, 19)
Ds = np.logspace(-2.5, 2.5, 41)
for c, col in ((0.0, "#2a9d8f"), (0.3, "#e76f51")):
    mx = []
    for L in Ls:
        mx.append(max(CN.resumen("hiperbola", float(L), float(D), c, 400.0 * max(1.0, L))["kmin"] for D in Ds))
    axs[1].plot(1 / Ls, mx, "o-", color=col, ms=3, lw=1.0, label=fr"hipérbola, $c={c:g}$")
nl = json.load(open("resultados/nled_concretas.json", encoding="utf-8"))
for tipo, mk in (("BI", "s"), ("EH", "^")):
    for par in sorted({o["par"] for o in nl if o["tipo"] == tipo}):
        sub = [o for o in nl if o["tipo"] == tipo and o["par"] == par and o.get("regular")]
        if sub:
            mejor = max(sub, key=lambda o: o["kmin"])
            axs[1].plot(mejor["rpp0"], mejor["kmin"], mk, color="#264653", ms=5, mfc="none")
axs[1].plot([], [], "s", color="#264653", mfc="none", label="Born–Infeld")
axs[1].plot([], [], "^", color="#264653", mfc="none", label="Euler–Heisenberg")
bus = json.load(open("resultados/busqueda_cinetico.json", encoding="utf-8"))
sup = max(b["max_min_K"] for b in bus["busquedas"])
axs[1].axhline(sup, color="k", lw=0.8, ls="--", label=f"supremo de la optimización global ({sup:.3f})")
axs[1].axhline(0, color="k", lw=0.6)
axs[1].set_xscale("log")
axs[1].set_ylim(-1.05, 0.15)
axs[1].set_xlabel(r"$r''(0)\,r_0$ (agudeza de la garganta)")
axs[1].set_ylabel(r"$\max_D\ \min_\ell\ K/\mathrm{norma}$")
axs[1].set_title("Ningún agujero ultraestático encontrado tiene $K>0$ en toda la línea", fontsize=9)
axs[1].legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig("figuras/fig4_cinetico.png", dpi=130)
print("figura: figuras/fig4_cinetico.png")
