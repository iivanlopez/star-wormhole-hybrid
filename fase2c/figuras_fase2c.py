# -*- coding: utf-8 -*-
"""
Figuras de la Fase 2c (figuras/):
  fig1_mapa_mitades.png   final de las mitades lisas en el espacio de datos de la garganta (cuatro acoplamientos)
  fig2_colapso.png        ley del colapso del lapso: l_s |a| frente a |a| r0 en las mitades repulsivas
  fig3_perfiles.png       continuacion lisa de una garganta de Antoniou (colapso a l = -1/a) y atractor conico
Entradas: resultados/mapa_lados.json, mapa_gauss.json, mapa_lineal.json.   Uso: python figuras_fase2c.py
"""
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fondo_doble import Lado, segundas      # noqa: E402

os.makedirs("figuras", exist_ok=True)
res = []
for f in ("resultados/mapa_lados.json", "resultados/mapa_gauss.json", "resultados/mapa_lineal.json"):
    if os.path.exists(f):
        res += json.load(open(f, encoding="utf-8"))

CLASES = ["AF", "noAF", "horizonte", "phi", "r0", "det", "Phi+", "fallo"]
NOMBRES = ["AF", "no AF", "horizonte singular", "divergencia de φ'", "r → 0", "det M → 0", "Φ' → +∞", "fallo"]
COL = ["#2a9d8f", "#c9ada7", "#e76f51", "#8d99ae", "#264653", "#e9c46a", "#b5838d", "#000000"]
cmap = ListedColormap(COL)


# ------------------------------------------------------------------ figura 1
def panel(ax, pts, xkey, titulo, xlab):
    """Filas: phi1 (arriba la mitad con a > 0, abajo la mitad repulsiva a < 0, reflejada: misma |phi1| a la misma
    distancia de la linea central, de modo que las dos mitades de un mismo agujero liso quedan enfrentadas)."""
    xs = sorted({o[xkey] for o in pts}, key=abs)
    us = sorted({abs(o["phi1"]) for o in pts})
    Z = np.full((2 * len(us), len(xs)), np.nan)
    ix = {v: i for i, v in enumerate(xs)}
    iu = {v: i for i, v in enumerate(us)}
    for o in pts:
        j = iu[abs(o["phi1"])]
        fila = len(us) + j if o["a"] > 0 else len(us) - 1 - j
        Z[fila, ix[o[xkey]]] = CLASES.index(o["fin"]) if o["fin"] in CLASES else 7
    ax.pcolormesh(np.arange(len(xs) + 1), np.arange(2 * len(us) + 1), Z, cmap=cmap, vmin=-0.5,
                  vmax=len(CLASES) - 0.5)
    tx = np.linspace(0, len(xs) - 1, 5).astype(int)
    ax.set_xticks(tx + 0.5)
    ax.set_xticklabels([f"{abs(xs[i]):.2g}" for i in tx], fontsize=7)
    ty = [len(us) // 2, len(us) - 1]
    ax.set_yticks([len(us)] + [len(us) + t + 0.5 for t in ty] + [len(us) - 1 - t + 0.5 for t in ty])
    ax.set_yticklabels([f"{us[0]:.2g}"] + [f"{us[t]:.2g}" for t in ty] * 2, fontsize=7)
    ax.axhline(len(us), color="w", lw=1.2)
    ax.text(0.02, 0.97, "a > 0", transform=ax.transAxes, fontsize=7, va="top", color="w")
    ax.text(0.02, 0.03, "a < 0 (repulsiva)", transform=ax.transAxes, fontsize=7, va="bottom", color="w")
    ax.set_title(titulo, fontsize=8)
    ax.set_xlabel(xlab, fontsize=8)
    ax.set_ylabel(r"$|\phi_1|\,r_0$", fontsize=8)


paneles = []
for a in (-2.0, -0.1, 0.1, 2.0):
    paneles.append(([o for o in res if o["tipo"] == "cuadratico" and o["alpha"] == a], "phi0",
                    fr"$F=\alpha\phi^2$, $\alpha/r_0^2={a:+g}$", r"$\phi_0$"))
for sa in (1, -1):
    paneles.append(([o for o in res if o["tipo"] == "exp" and np.sign(o["alpha"]) == sa], "alpha",
                    fr"$F=\alpha e^{{-\phi}}$, $\phi_0=0$, $\alpha{'>' if sa > 0 else '<'}0$", r"$|\alpha|/r_0^2$"))
for a in (0.2, -0.2):
    paneles.append(([o for o in res if o["tipo"] == "gauss" and o["alpha"] == a], "phi0",
                    fr"$F=\alpha e^{{-\phi^2}}$, $\alpha/r_0^2={a:+g}$", r"$\phi_0$"))
for sa in (1, -1):
    paneles.append(([o for o in res if o["tipo"] == "lineal" and np.sign(o["alpha"]) == sa], "alpha",
                    fr"$F=\alpha\phi$, $\alpha{'>' if sa > 0 else '<'}0$", r"$|\alpha|/r_0^2$"))
paneles = [p for p in paneles if p[0]]
nc = 5
nf = int(np.ceil(len(paneles) / nc))
fig, axs = plt.subplots(nf, nc, figsize=(3.3 * nc, 3.2 * nf + 0.6), squeeze=False)
for ax, p in zip(axs.flat, paneles):
    panel(ax, *p)
for ax in list(axs.flat)[len(paneles):]:
    ax.axis("off")
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in COL]
fig.legend(handles, NOMBRES, loc="lower center", ncol=8, fontsize=8)
fig.suptitle("Final de las mitades lisas (r0 = 1). Un agujero liso necesita las dos mitades enfrentadas "
             "(misma columna y misma |φ1|) en verde: no ocurre en ningún punto", fontsize=9)
fig.tight_layout(rect=(0, 0.05, 1, 0.96))
fig.savefig("figuras/fig1_mapa_mitades.png", dpi=130)
plt.close(fig)

# ------------------------------------------------------------------ figura 2
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
mk = {"cuadratico": ("o", "#e76f51", r"$\alpha\phi^2$"), "exp": ("s", "#264653", r"$\alpha e^{-\phi}$"),
      "gauss": ("^", "#2a9d8f", r"$\alpha e^{-\phi^2}$"), "lineal": ("D", "#e9c46a", r"$\alpha\phi$")}
for tipo, (m, c, lab) in mk.items():
    hz = [o for o in res if o["tipo"] == tipo and o["a"] < 0 and o["fin"] == "horizonte"]
    if not hz:
        continue
    x = np.array([abs(o["a"]) for o in hz])
    y = np.array([o["l_end"] * abs(o["a"]) for o in hz])
    axs[0].scatter(x, y, s=3, marker=m, color=c, alpha=0.25, label=f"{lab} ({len(hz)})", rasterized=True)
    bordes = np.logspace(-3, 6, 28)
    xm, ym = [], []
    for lo, hi in zip(bordes[:-1], bordes[1:]):
        sel = (x >= lo) & (x < hi)
        if sel.sum() >= 10:
            xm.append(np.sqrt(lo * hi))
            ym.append(np.median(y[sel]))
    axs[0].plot(xm, ym, color=c, lw=2.2, marker=m, ms=4, mec="k", mew=0.4)
axs[0].axhline(1.0, color="k", lw=0.8, ls="--")
axs[0].set_xscale("log")
axs[0].set_yscale("log")
axs[0].set_xlabel(r"$|a|\,r_0=|\Phi'(0)|\,r_0$")
axs[0].set_ylabel(r"$\ell_s\,|a|$")
axs[0].set_title("Mitades repulsivas con horizonte singular (líneas: mediana por intervalo)", fontsize=9)
axs[0].legend(fontsize=7, markerscale=3)
otros = [o for o in res if o["a"] < 0 and o["fin"] != "horizonte"]
for fin, c in (("phi", "#8d99ae"), ("r0", "#264653"), ("det", "#e9c46a"), ("fallo", "k"), ("noAF", "#c9ada7")):
    ss = [o for o in otros if o["fin"] == fin]
    if ss:
        axs[1].scatter([abs(o["a"]) for o in ss], [o["l_end"] for o in ss], s=4, color=c, label=f"{fin} ({len(ss)})",
                       rasterized=True)
axs[1].set_xscale("log")
axs[1].set_yscale("log")
axs[1].set_xlabel(r"$|a|\,r_0$")
axs[1].set_ylabel(r"$\ell_{\rm fin}/r_0$")
axs[1].set_title("Mitades repulsivas con otros finales singulares", fontsize=9)
axs[1].legend(fontsize=7, markerscale=3)
fig.tight_layout()
fig.savefig("figuras/fig2_colapso.png", dpi=130)
plt.close(fig)

# ------------------------------------------------------------------ figura 3
fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
D = Lado("cuadratico", 0.1, 0.1, 0.32989, l_max=60, rtol=1e-12).resolver()      # mitad AF (Antoniou, phi_inf ~ 0)
I = Lado("cuadratico", 0.1, 0.1, -0.32989, l_max=60, rtol=1e-12).resolver()     # mitad repulsiva
a = D.a
ld = np.linspace(0, 6 / a, 400)
li = np.linspace(0, I.l_end * (1 - 1e-9), 400)
ax = axs[0]
ax.plot(ld * a, np.exp(D.sol.sol(ld)[2]), color="#2a9d8f", label=r"$N=e^{\Phi}$, liso, $\ell>0$")
ax.plot(-li * a, np.exp(I.sol.sol(li)[2]), color="#e76f51", label=r"$N$, liso, $\ell<0$ (colapso)")
ax.plot(-ld * a, np.exp(D.sol.sol(ld)[2]), color="#2a9d8f", ls=":", label=r"$N$, $\mathbb{Z}_2$ con cáscara")
ax.plot(np.r_[-1.5, 6], 1 + a * np.r_[-1.5, 6] / a, color="k", lw=0.6, ls="--", label=r"$1+a\ell$ (Rindler)")
ax.axvline(-1, color="0.6", lw=0.6)
ax.set_xlim(-1.5, 6)
ax.set_ylim(0, 4)
ax.set_xlabel(r"$a\,\ell$")
ax.set_title(fr"Garganta de Antoniou ($\alpha=0.1$, $\phi_0=0.1$): $a\,r_0={a:.1f}$", fontsize=9)
ax.legend(fontsize=7)
ax = axs[1]
lz = I.l_end * (1 - np.logspace(-7, 0, 600))          # muestreo logaritmico en la distancia al final
K = []
for x in lz:
    y = I.sol.sol(x)
    s2, _ = segundas(y, I.Fd)
    r, rp, P = y[0], y[1], y[3]
    K.append(4 * (s2[1] + P * P) ** 2 + 8 * (rp * P / r) ** 2 + 8 * (s2[0] / r) ** 2 + 4 * ((1 - rp * rp) / r**2) ** 2)
l_s = I.l_end - 1.0 / I.sol.y[3, -1]          # Phi' ~ -1/(l_s - l): posicion del punto singular
xs = (l_s - lz) * a
ax.loglog(xs, np.exp(I.sol.sol(lz)[2]), color="#e76f51", label=r"$N$")
ax.loglog(xs, np.array(K) / a**4, color="#264653", label=r"$K/a^4$ (Kretschmann)")
ax.loglog(xs, xs, color="k", lw=0.6, ls="--", label=r"$\propto(\ell_s-\ell)$")
ax.loglog(xs, 1 / xs**2, color="k", lw=0.6, ls=":", label=r"$\propto(\ell_s-\ell)^{-2}$")
ax.set_xlabel(r"$a\,(\ell_s-|\ell|)$")
ax.set_title(fr"Mitad repulsiva: $N\to0$ lineal y $K\to\infty$ en $\ell_s={l_s * a:.3f}/a$", fontsize=9)
ax.legend(fontsize=7)
ax = axs[2]
for al, c in ((-10.0, "#8d99ae"), (-100.0, "#2a9d8f"), (-1000.0, "#264653")):
    L = Lado("cuadratico", al, 1.0, 20.0 if al < 0 else -20.0, l_max=1e5, rtol=1e-10).resolver()
    lz = np.linspace(0, L.l_end * (1 - 1e-9), 3000)
    y = L.sol.sol(lz)
    ax.plot(lz / L.l_end, y[1], color=c, label=fr"$r'$, $\alpha={al:g}$ ($\ell_s={L.l_end:.0f}\,r_0$)")
    ax.plot(lz / L.l_end, np.exp(y[2]), color=c, ls="--")
ax.axhline(1 / np.sqrt(3), color="k", lw=0.6, ls=":", label=r"$1/\sqrt{3}$")
ax.set_xlabel(r"$\ell/\ell_s$")
ax.set_title(r"Acoplamiento fuerte ($\alpha\phi^2$, $\phi_0=1$): $r'\to1/\sqrt{3}$; trazos: $N$", fontsize=9)
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig("figuras/fig3_perfiles.png", dpi=130)
plt.close(fig)
print("figuras: fig1_mapa_mitades.png, fig2_colapso.png, fig3_perfiles.png")
