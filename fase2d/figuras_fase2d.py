# -*- coding: utf-8 -*-
"""
Figuras de la Fase 2d (figuras/):
  fig1_hiperbolicidad.png   lambda(l) (cuadrado de la velocidad caracteristica) para varios Qs/Q: capa mal planteada
  fig2_torre.png            raices reales w_n frente a la prediccion WKB sqrt|w_n| I = (n + 1/2) pi
  fig3_autofunciones.png    autofunciones (q, d) de la torre: localizadas en la capa mal planteada
Entradas: resultados/hiperbolicidad_canate.json, resultados/modos_canate.json.   Uso: python figuras_fase2d.py
"""
import os
import sys
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import canate as C                       # noqa: E402
import estabilidad_lisa as E             # noqa: E402
from hiperbolicidad_canate import coeficientes   # noqa: E402

os.makedirs("figuras", exist_ok=True)
COL = {0.1: "#8d99ae", 0.5: "#e9c46a", 1.0: "#e76f51", 2.0: "#2a9d8f", 5.0: "#264653"}

# ------------------------------------------------------------------ figura 1
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
ls = np.unique(np.r_[np.logspace(-3, 0, 500), np.logspace(0, 4, 800)])
for Qs, c in COL.items():
    f = C.FondoCanate(Qs)
    cs = np.array([coeficientes(f, l) for l in ls])
    lam = -cs[:, 0] / cs[:, 1]
    axs[0].plot(ls, lam, color=c, lw=1.2, label=fr"$Q_s/Q={Qs:g}$")
    axs[1].plot(ls, -cs[:, 1] / cs[:, 0], color=c, lw=1.2)
axs[0].set_xscale("log")
axs[0].set_yscale("symlog", linthresh=1.0)
axs[0].axhline(0, color="k", lw=0.6)
axs[0].axhline(1, color="k", lw=0.6, ls=":")
axs[0].set_xlabel(r"$\ell/Q$")
axs[0].set_ylabel(r"$\lambda=w/k^2$ (luz: 1)")
axs[0].set_title(r"Cuadrado de la velocidad característica: $\lambda<0$ en una capa (mal planteado)", fontsize=9)
axs[0].legend(fontsize=7)
axs[1].set_xscale("log")
axs[1].axhline(0, color="k", lw=0.6)
axs[1].axhline(1, color="k", lw=0.6, ls=":")
axs[1].set_ylim(-3, 3)
axs[1].set_xlabel(r"$\ell/Q$")
axs[1].set_ylabel(r"$1/\lambda=-c_1/c_0$")
axs[1].set_title(r"Coeficiente cinético relativo $1/\lambda$: se anula en los polos y es negativo en la capa", fontsize=9)
fig.tight_layout()
fig.savefig("figuras/fig1_hiperbolicidad.png", dpi=130)
plt.close(fig)

# ------------------------------------------------------------------ figura 2
mod = json.load(open("resultados/modos_canate.json", encoding="utf-8"))
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
for o in mod:
    c = COL.get(o["Qs"], "k")
    for nombre, mk in (("par", "o"), ("impar", "s")):
        ws = np.array([r["w"] for r in o["raices"][nombre]])
        if len(ws) == 0:
            continue
        x = np.sqrt(-ws) * o["I"] / math.pi
        axs[0].plot(np.arange(len(x)), x, mk, color=c, ms=4, mfc="none" if nombre == "impar" else c,
                    label=fr"$Q_s/Q={o['Qs']:g}$, {nombre}")
        axs[1].plot(np.arange(1, len(x)), np.diff(x), mk, color=c, ms=4, mfc="none" if nombre == "impar" else c)
n = np.arange(0, 64)
axs[0].plot(n, n + 1.5, "k--", lw=0.8, label=r"pendiente 1: $n+3/2$")
axs[0].set_xlabel("n (orden de la raíz)")
axs[0].set_ylabel(r"$\sqrt{|w_n|}\,I/\pi$")
axs[0].set_title(r"Torre de modos inestables: $\omega_n^2=w_n\to-\infty$ (sin cota)", fontsize=9)
axs[0].legend(fontsize=6, ncol=2)
axs[1].axhline(1, color="k", lw=0.8, ls="--")
axs[1].set_ylim(0.6, 1.4)
axs[1].set_xlabel("n")
axs[1].set_ylabel(r"$(\sqrt{|w_n|}-\sqrt{|w_{n-1}|})\,I/\pi$")
axs[1].set_title(r"Espaciado entre raíces consecutivas: $\to1$ (WKB en la capa con $\lambda<0$; llenos: pares)",
                 fontsize=9)
fig.tight_layout()
fig.savefig("figuras/fig2_torre.png", dpi=130)
plt.close(fig)

# ------------------------------------------------------------------ figura 3
o = [m for m in mod if m["Qs"] == 1.0][0]
f = C.FondoCanate(1.0)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
for ax, nombre in zip(axs, ("par", "impar")):
    for k, c in zip((0, 1, 5, 12), ("#264653", "#2a9d8f", "#e9c46a", "#e76f51")):
        rs = o["raices"][nombre]
        if k >= len(rs):
            continue
        w = rs[k]["w"]
        L = max(o["l2"] + 15.0, 10.0 + 30.0 / math.sqrt(-w))
        y0 = E.dato_asintotico(f, L, w)
        s = solve_ivp(lambda x, yy: C.derivadas_Q(f, x, yy, w), (L, 0.0), y0, method="DOP853", rtol=1e-10, atol=1e-30,
                      dense_output=True)
        lz = np.linspace(0, min(L, o["l2"] * 2.0), 1500)
        # gauge d'(0) = 0 (la garganta perturbada sigue en l = 0): U_+ + c g con c = -d'_+(0)/r''(0)
        cg = -s.sol(0.0)[2] / f.rpp0
        q = np.array([s.sol(x)[3] + cg * E.modo_gauge(f, w, x)[3] for x in lz])
        ax.plot(lz, q / np.max(np.abs(q)), color=c, lw=1.0, label=fr"$n={k}$: $w={w:.4g}$")
    ax.axvspan(o["l1"], o["l2"], color="0.9", zorder=-1)
    ax.set_xlabel(r"$\ell/Q$")
    ax.set_ylabel(r"$q/\max|q|$")
    ax.set_title(f"Autofunciones del sector {nombre} ($Q_s=Q$); sombreado: capa mal planteada", fontsize=9)
    ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig("figuras/fig3_autofunciones.png", dpi=130)
plt.close(fig)
print("figuras: fig1_hiperbolicidad.png, fig2_torre.png, fig3_autofunciones.png")
