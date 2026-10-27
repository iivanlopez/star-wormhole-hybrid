# -*- coding: utf-8 -*-
"""Figuras de la Fase 2l a partir de resultados/estimacion_rotacion.json.
Uso: python figuras_fase2l.py   -> figuras/fig1_rotacion.png"""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(AQUI, "resultados", "estimacion_rotacion.json")))
g = d["gargantas"]
os.makedirs(os.path.join(AQUI, "figuras"), exist_ok=True)
N0 = np.array([x["N0"] for x in g])
vc = np.array([x["v_c"] for x in g])
ka = np.array([x["kappa"] for x in g])

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].scatter(N0, vc, c="#2a6f97", s=28, label="gargantas de la Fase 2f")
ax[0].scatter([1.0], [d["ellis"]["v_c"]], marker="*", c="#c1121f", s=120, label="Ellis–Bronnikov")
ax[0].axhline(1.0, color="k", lw=0.8, ls="--")
ax[0].axhspan(0, 0.3, color="#90be6d", alpha=0.25, label="rotación lenta controlada ($v_0\\lesssim0.3$)")
ax[0].set_xlabel("lapso en la garganta $N_0$")
ax[0].set_ylabel("$v_c=\\kappa^{-1/2}$")
ax[0].set_title("(a) velocidad crítica de la prueba $d=N$")
ax[0].set_ylim(0, 1.9)
ax[0].legend(fontsize=8, loc="upper right")
v = np.linspace(0, 1, 200)
for k, c in ((ka.min(), "#577590"), (np.median(ka), "#2a6f97"), (ka.max(), "#003049"), (d["ellis"]["kappa"], "#c1121f")):
    ax[1].plot(v, -1 + k * v ** 2, color=c, label=f"$\\kappa={k:.2f}$")
ax[1].axhline(0, color="k", lw=0.8)
ax[1].axvspan(0, 0.3, color="#90be6d", alpha=0.25)
ax[1].set_xlabel("$v_0=r_0|\\omega(0)-\\omega(\\infty)|/N_0$")
ax[1].set_ylabel("$E_J[N]/|E_0|$  (orden $\\Omega^2$)")
ax[1].set_title("(b) energía de la prueba con rotación")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figuras", "fig1_rotacion.png"), dpi=150)
print("figuras/fig1_rotacion.png")
