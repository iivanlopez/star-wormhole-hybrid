# -*- coding: utf-8 -*-
"""
Figuras de la Fase 2e (r0 = 1):
  fig1_perfiles.png  izquierda: kink grueso (F = -0.07 phi^2, doble pozo v = 1, u = 1): el lapso sigue al de Nariai,
                     cos(l), y se anula en l ~ pi/2 donde chi -> 0 con num != 0 (horizonte);
                     centro: separatrices del disparo (el escalar llega al vacio): r'(l) para paredes finas (v = 0.3)
                     y gruesas (v = 1), siempre con r' < 0 al llegar; derecha: cruce regular de chi (F = -0.5 phi^2,
                     v = 1, u = 1): chi, K/norma y G/norma; lambda < 0 alrededor del cero ascendente.
  fig2_teorema.png   K frente a -16 F'^2 (1 - r'^2)^2 chi'/r'^2 en todos los ceros regulares de chi del mapa, y
                     lambda en los ascendentes (resultados/mapa_potencial.json).
  fig3_mapa.png      r' con que llegan al vacio las separatrices del mapa, por familia y pendiente del potencial.
Uso: python figuras_fase2e.py
"""
import os
import sys
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import potencial as P               # noqa: E402
import mapa_potencial as MP         # noqa: E402

AZUL, NARANJA, AGUA = "#2a78d6", "#eb6834", "#1baf7a"       # paleta categorica validada (3 series)
TINTA, TINTA2, REJILLA = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2,
                     "ytick.color": TINTA2, "axes.grid": True, "grid.color": REJILLA, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False})
os.makedirs(os.path.join(AQUI, "figuras"), exist_ok=True)


def kink(alpha, v, u, rtol=1e-12):
    Fd = P.acoplamiento("cuadratico", alpha)
    amp = (2 + u * u / 2) / v**4
    return P.LadoV(Fd, P.potencial("doble_pozo", amp, v=v), 0.0, u, 0.0, phi_v=v, l_max=300.0, rtol=rtol).resolver()


def serie(L, n=1500, frac=1.0):
    ls = np.linspace(1e-4, L.l_end * frac, n)
    return ls, [MP.magnitudes(L, x) for x in ls]


# ------------------------------------------------------------------ figura 1
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.6))

L = kink(-0.07, 1.0, 1.0)
ls, ms = serie(L, frac=1 - 1e-7)
N = np.array([m["N"] for m in ms])
chi = np.array([m["chi"] for m in ms])
num = np.array([m["num"] for m in ms])
a1.plot(ls, N, color=AZUL, lw=1.8, label=r"lapso $N=e^\Phi$")
lc = np.linspace(0, math.pi / 2, 200)
a1.plot(lc, np.cos(lc), color=TINTA2, lw=1.1, ls=":", label=r"Nariai: $\cos(\ell/r_0)$")
a1.plot(ls, chi, color=NARANJA, lw=1.6, ls="--", label=r"$\chi=rr'+2F_\ell(1-3r'^2)$")
a1.plot(ls, num, color=AGUA, lw=1.6, ls="-.", label=r"num $=(1-r'^2)/r^2+\phi'^2/4-V/2$")
a1.axvline(L.l_end, color=TINTA2, lw=0.8)
a1.text(L.l_end - 0.03, 1.25, f"horizonte\n$\\ell_s={L.l_end:.3f}$", ha="right", fontsize=8, color=TINTA2)
a1.set_xlim(0, 1.62)
a1.set_ylim(-0.4, 1.6)
a1.set_xlabel(r"$\ell/r_0$")
a1.set_title("Kink grueso: el lapso de Nariai colapsa donde\n"
             r"$\chi\to0$ con num$\neq0$ ($2\chi\Phi'=r^2$num)", fontsize=9, color=TINTA)
a1.legend(fontsize=7.5, frameon=False, loc="center left", bbox_to_anchor=(0.0, 0.42))

res = json.load(open(os.path.join(AQUI, "resultados", "mapa_potencial.json"), encoding="utf-8"))
casos = []
for o in res:
    if o["fam"] == "B" and o["F"] == "cuadratico" and o["alpha"] == -0.07 and o["v"] in (0.3, 0.5, 1.0):
        for s in o["separatrices"]:
            casos.append((o["v"], s["x"]))
estilos = {0.3: (NARANJA, "-", "pared fina"), 0.5: (AGUA, "-.", "intermedia"), 1.0: (AZUL, "--", "pared gruesa")}
for v, x in sorted(casos):
    Ls = kink(-0.07, v, x)
    ll = MP.llegada(Ls)
    l2 = np.linspace(1e-4, ll["l"], 600)                    # hasta la maxima aproximacion al vacio (despues: fuga)
    m2 = [MP.magnitudes(Ls, y) for y in l2]
    rp = np.array([m["rp"] for m in m2])
    Nl = np.array([m["N"] for m in m2])
    col, est, nombre = estilos[v]
    a2.plot(l2, rp, color=col, ls=est, lw=1.5, label=f"$v={v}$ ({nombre}): $r'$")
    a2.plot(l2, Nl, color=col, ls=":", lw=1.1)
    a2.plot(ll["l"], ll["rp"], "o", color=col, ms=7, mec="white", mew=1.3)
    a2.plot(ll["l"], ll["N"], "o", color=col, ms=5, mfc="white", mew=1.1)
a2.plot([], [], color=TINTA2, ls=":", lw=1.1, label="lapso $N$ (punteado)")
a2.axhline(0, color=TINTA2, lw=0.8)
a2.set_xlim(0, 0.8)
a2.set_ylim(-0.6, 1.05)
a2.set_xlabel(r"$\ell/r_0$")
a2.set_title("Separatrices del disparo (kink, $F=-0.07\\,\\phi^2$): $r'$ al acercarse\n"
             r"al vacío ($\bullet$) es negativo; con pared gruesa $N$ cae hacia Nariai", fontsize=9, color=TINTA)
a2.legend(fontsize=7.5, frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.36))

L = kink(-0.5, 1.0, 1.0)
ls, ms = serie(L, n=2500, frac=0.95)
chi = np.array([m["chi"] for m in ms])
K = np.array([m["K"] for m in ms])
G = np.array([m["G"] for m in ms])
lam = G * np.array([m["N"] for m in ms]) ** 2 / K
norma = np.abs(K) + np.abs(G)
a3.plot(ls, 10 * chi, color=NARANJA, lw=1.6, ls="--", label=r"$10\,\chi$")
a3.plot(ls, K / norma, color=AZUL, lw=1.8, label=r"$K/(|K|+|G|)$ (cinético)")
a3.plot(ls, G / norma, color=AGUA, lw=1.6, ls="-.", label=r"$G/(|K|+|G|)$ (gradiente)")
mal = lam < 0
a3.fill_between(ls, -1.05, 1.05, where=mal, color=REJILLA, alpha=0.9, lw=0, label=r"$\lambda=N^2G/K<0$")
for z in MP.ceros_chi(L):
    a3.axvline(z["l"], color=TINTA2, lw=0.8, ls=":")
a3.axhline(0, color=TINTA2, lw=0.8)
a3.set_ylim(-1.05, 1.05)
a3.set_xlabel(r"$\ell/r_0$")
a3.set_title("Cruce regular de $\\chi$ (kink, $\\alpha=-0.5$): en el cero ascendente\n"
             r"$K=-16F'^2(1-r'^2)^2\chi'/r'^2<0$ y $G>0$", fontsize=9, color=TINTA)
a3.legend(fontsize=7.5, frameon=False, loc="lower right")
fig.savefig(os.path.join(AQUI, "figuras", "fig1_perfiles.png"), dpi=130, bbox_inches="tight")
print("figura: figuras/fig1_perfiles.png")

# ------------------------------------------------------------------ figura 2
zs = [z for o in res for p in o["puntos"] if p.get("valido") for z in p["ceros_chi"]]
zs += [z for o in res for s in o["separatrices"] for z in s["ceros_chi"]]
fig, (b1, b2) = plt.subplots(1, 2, figsize=(11.5, 4.6))
for sube, col, mk, nombre in ((True, AZUL, "o", "cero ascendente ($\\chi'>0$)"),
                              (False, NARANJA, "^", "cero descendente ($\\chi'<0$)")):
    sel = [z for z in zs if z["sube"] == sube]
    x = np.array([z["K_teo"] for z in sel])
    y = np.array([z["K"] for z in sel])
    b1.scatter(np.sign(x) * np.log10(1 + np.abs(x)), np.sign(y) * np.log10(1 + np.abs(y)), s=14, color=col,
               marker=mk, edgecolors="white", linewidths=0.4, label=f"{nombre}: {len(sel)}")
lim = 6
b1.plot([-lim, lim], [-lim, lim], color=TINTA2, lw=0.8, ls=":")
b1.axhline(0, color=TINTA2, lw=0.8)
b1.axvline(0, color=TINTA2, lw=0.8)
b1.set_xlabel(r"$-16F'^2(1-r'^2)^2\chi'/r'^2$   [escala $\mathrm{sgn}\cdot\log_{10}(1+|x|)$]")
b1.set_ylabel(r"$K$ calculado   [misma escala]")
b1.set_title("Identidad en los ceros regulares de $\\chi$ (todas las mitades del mapa)", fontsize=9, color=TINTA)
b1.legend(fontsize=8, frameon=False, loc="upper left")
sel = [z for z in zs if z["sube"]]
lam = np.array([z["lam"] for z in sel])
b2.hist(np.sign(lam) * np.log10(1 + np.abs(lam)), bins=40, color=AZUL, edgecolor="white", linewidth=0.6)
b2.axvline(0, color=TINTA2, lw=0.8)
b2.set_xlabel(r"$\lambda=N^2G/K$ en los ceros ascendentes   [escala $\mathrm{sgn}\cdot\log_{10}(1+|\lambda|)$]")
b2.set_ylabel("ceros")
b2.set_title(f"$\\lambda<0$ en {int((lam < 0).sum())} de {len(lam)} ceros ascendentes: problema mal planteado",
             fontsize=9, color=TINTA)
fig.savefig(os.path.join(AQUI, "figuras", "fig2_teorema.png"), dpi=130, bbox_inches="tight")
print("figura: figuras/fig2_teorema.png")

# ------------------------------------------------------------------ figura 3
# estado de cada separatriz en su punto de maxima aproximacion al vacio: (N, r'); r' > 0 seria una mitad que se abre
fig, ax = plt.subplots(figsize=(8.5, 5.0))
colores = {0.3: NARANJA, 0.5: AGUA}
seps = [(o, s) for o in res for s in o["separatrices"]]
n_pos = sum(s["llegada_lo"]["rp"] > 0 for o, s in seps)
for (fam, mk) in (("A", "o"), ("B", "^")):
    for vv, col, nombre in ((0.3, NARANJA, "$v=0.3$ (empinado)"), (0.5, AGUA, "$v=0.5$"),
                            (1.0, AZUL, r"$v\geq1$ (suave)")):
        sel = [s["llegada_lo"] for o, s in seps if o["fam"] == fam and (o["v"] == vv or (vv == 1.0 and o["v"] >= 1.0))]
        if not sel:
            continue
        N_ = np.array([d["N"] for d in sel])
        rp_ = np.array([d["rp"] for d in sel])
        llega = np.array([d["dphi"] < 1e-3 for d in sel])
        ax.scatter(N_[llega], rp_[llega], s=40, color=col, marker=mk, edgecolors="white", linewidths=0.6,
                   label=f"clase {fam}, {nombre}")
        ax.scatter(N_[~llega], rp_[~llega], s=40, facecolors="none", edgecolors=col, marker=mk, linewidths=1.1)
ax.scatter([], [], s=40, facecolors="none", edgecolors=TINTA2, marker="o", label="no llega al vacío (corte singular)")
xs = np.linspace(0, 1, 50)
ax.plot(xs, -0.5 * xs, color=TINTA2, lw=0.9, ls=":")
ax.text(0.47, -0.31, r"$r'=-N/2$", fontsize=8, color=TINTA2)
ax.axhline(0, color=TINTA, lw=0.9)
ax.set_xlim(0, 1.02)
ax.set_ylim(-0.62, 0.12)
ax.set_xlabel(r"lapso $N$ en el punto de máxima aproximación al vacío")
ax.set_ylabel(r"$r'$ en ese punto")
ax.set_title(f"Las {len(seps)} separatrices del mapa: {n_pos} con $r'>0$. Empinado: pared fina que cierra la\n"
             r"garganta ($r'\approx-0.5$); suave: $r'\to0^-$ a la vez que $N\to0$ (horizonte de Nariai)", fontsize=9,
             color=TINTA)
ax.legend(fontsize=7.5, frameon=False, loc="upper right", ncol=2)
fig.savefig(os.path.join(AQUI, "figuras", "fig3_mapa.png"), dpi=130, bbox_inches="tight")
print("figura: figuras/fig3_mapa.png")
