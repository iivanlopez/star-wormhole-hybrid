# -*- coding: utf-8 -*-
"""
Figuras de la Fase 2b (PNG, modo claro; paleta de la Fase 2).
Requiere resultados/dominio_kkk.json, estabilidad_dominio.json, trayectorias_robustas.json y modos_complejos.json.
Uso: python figuras.py
"""
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from kkk import FondoKKK
from estabilidad_kkk import CapaKKK, D, vector_garganta

AZUL, NARANJA, AGUA = "#2a78d6", "#eb6834", "#1baf7a"
TINTA, TINTA2, MUDO, REJILLA, EJE, FONDO = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": FONDO, "axes.facecolor": FONDO, "savefig.facecolor": FONDO,
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"], "font.size": 10,
    "axes.edgecolor": EJE, "axes.labelcolor": TINTA2, "xtick.color": MUDO, "ytick.color": MUDO,
    "axes.grid": True, "grid.color": REJILLA, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "lines.linewidth": 2.0, "legend.frameon": False, "legend.labelcolor": TINTA2,
    "axes.titlecolor": TINTA, "axes.titlesize": 11,
})


def fig_dominio():
    dom = json.load(open("resultados/dominio_kkk.json", encoding="utf-8"))
    est = json.load(open("resultados/estabilidad_dominio.json", encoding="utf-8"))
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    reg = [o for o in dom if o.get("regular")]
    sing = [o for o in dom if "f0" in o and not o.get("regular") and o["f0"] < 1e3]
    ax.scatter([o["alpha_kkk"] for o in sing], [o["f0"] for o in sing], s=6, color=EJE, label="sin solucion regular")
    ax.scatter([o["alpha_kkk"] for o in reg], [o["f0"] for o in reg], s=6, color=MUDO, label="agujero regular (KKK)")
    ins = [o for o in est if o.get("winf_zeta1")]
    if ins:
        tasa = [math.sqrt(-min(o["winf_zeta1"])) for o in ins]
        sc = ax.scatter([o["alpha_kkk"] for o in ins], [o["f0"] for o in ins], c=tasa, s=70, cmap="viridis",
                        edgecolor=TINTA, linewidth=0.6, label="modo inestable (capa de KKK)", zorder=3)
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label("tasa de crecimiento |omega_inf| r0 / c")
    sin = [o for o in est if "w_zeta1" in o and not o["w_zeta1"]]
    if sin:
        ax.scatter([o["alpha_kkk"] for o in sin], [o["f0"] for o in sin], s=70, marker="x", color=NARANJA,
                   label="sin raiz real en |w| < 1e6 (modo complejo o fuera de rango)", zorder=3)
    ax.axvline(0.0507, color=TINTA2, lw=0.8, ls="--")
    ax.text(0.052, 60, "0.0507", color=TINTA2, fontsize=8)
    ax.set_yscale("log")
    ax.set_xlim(0, 0.14)
    ax.set(xlabel="alpha / r0^2 (phi_inf = 0)", ylabel="f0 = 1/(r0 r''(0))",
           title="Dominio de KKK y su modo radial inestable (garganta libre, capa dinamica)")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig("figuras/fig1_dominio.png", dpi=150)


def fig_zeta():
    bz = [o for o in json.load(open("resultados/trayectorias_robustas.json", encoding="utf-8"))
          if o["nombre"] == "KKK referencia"]
    mc = json.load(open("resultados/modos_complejos.json", encoding="utf-8"))
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.9))
    # rama que diverge en zeta_c (el modo lento de zeta = 3 es otra rama: se muestra en el panel derecho)
    z = [o["zeta"] for o in bz if o["w"] and 0 <= o["zeta"] <= 2]
    inv = [1.0 / abs(min(o["w"])) for o in bz if o["w"] and 0 <= o["zeta"] <= 2]
    sel = [(a, b) for a, b in zip(z, inv) if 0.5 <= a <= 1.6]
    p = np.polyfit([a for a, _ in sel], [b for _, b in sel], 1)
    zc = -p[1] / p[0]
    ax[0].plot(z, inv, "o", color=AZUL, ms=6, label="modo real inestable")
    zz = np.linspace(0, zc, 50)
    ax[0].plot(zz, np.polyval(p, zz), color=AZUL, lw=1, ls="--", label=f"ajuste lineal: zeta_c = {zc:.2f}")
    ax[0].axvline(1.0, color=MUDO, lw=0.8)
    ax[0].text(1.02, max(inv) * 0.9, "capa de KKK", color=TINTA2, fontsize=8)
    ax[0].set_xlim(-0.05, 2.0)
    ax[0].set_ylim(0, max(inv) * 1.05)
    ax[0].set(xlabel="zeta (gravedad inducida de la capa / KKK)", ylabel="1 / |w*|  (reloj de la garganta)",
              title="La tasa diverge en zeta_c: la inercia efectiva se anula")
    ax[0].legend(fontsize=8)
    for o in mc:
        if o["convergido"] and o["w"][1] > 0:
            ax[1].plot(o["w"][0], o["w"][1], "o", color=NARANJA, ms=8)
            ax[1].plot(o["w"][0], -o["w"][1], "o", color=NARANJA, ms=8, mfc="none")
            ax[1].annotate(f"zeta = {o['zeta']}: w = {o['w'][0]:.0f} +/- {o['w'][1]:.0f} i", (o["w"][0], o["w"][1]),
                           (o["w"][0] * 0.35, o["w"][1] * 0.6), color=TINTA2, fontsize=8,
                           arrowprops=dict(arrowstyle="-", color=MUDO, lw=0.8))
    for o in bz:
        if o["zeta"] in (1.0, 3.0):
            for w in o["w"]:
                ax[1].plot(w, 0, "s", color=AZUL, ms=6)
                ax[1].annotate(f"zeta = {o['zeta']:g}", (w, 0), (w, 60), color=TINTA2, fontsize=8)
    ax[1].axhline(0, color=EJE, lw=0.8)
    ax[1].plot([0, 2500], [0, 0], color=TINTA2, lw=3, alpha=0.3)
    ax[1].text(1500, -40, "corte (continuo)", color=TINTA2, fontsize=8)
    ax[1].set(xlabel="Re w", ylabel="Im w", title="Modos inestables en el plano w (reloj de la garganta)")
    fig.tight_layout()
    fig.savefig("figuras/fig2_zeta.png", dpi=150)
    return zc


def fig_autofuncion():
    f = FondoKKK(alpha=0.03, phi1=-0.7).resolver()
    c = CapaKKK(f)
    w0 = -391.88143094
    v = vector_garganta(f, c, w0)[0]
    _, sol, _ = D(f, c, w0, u=v, devolver=True)
    l = np.logspace(-4, math.log10(6.0), 400)
    Y = sol.sol(l)
    s = np.sign(Y[0][0])
    fig, ax = plt.subplots(figsize=(7.0, 3.9))
    ax.plot(l, s * Y[0] / abs(Y[0][0]), color=AZUL, label="a' = d(delta Phi)/dl  (lapso)")
    ax.plot(l, s * Y[3] / abs(Y[0][0]) * 100, color=NARANJA, label="100 x q  (dilaton)")
    ax.plot(l, s * Y[1] / abs(Y[0][0]) * 100, color=AGUA, label="100 x d  (radio areal)")
    ax.axhline(0, color=EJE, lw=0.8)
    ax.set_xscale("log")
    ax.set(xlabel="longitud propia l / r0", ylabel="amplitud (a'(0) = 1)",
           title="Modo inestable de KKK (alpha/r0^2 = 0.028, f0 = 1.03, capa de KKK): localizado en la garganta")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figuras/fig3_autofuncion.png", dpi=150)


if __name__ == "__main__":
    import os
    os.makedirs("figuras", exist_ok=True)
    fig_dominio()
    print("zeta_c (ajuste) =", fig_zeta())
    fig_autofuncion()
