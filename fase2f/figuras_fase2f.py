# -*- coding: utf-8 -*-
"""
Figuras de la Fase 2f (r0 = 1), a partir de resultados/perfiles_bH.npz y resultados/gargantas_bH.json
(gargantas_bH.py):
  fig1_fondo.png      garganta masiva de referencia (M = 2, a = 2, b = 4, n = 2): (a) r', lapso N y phi'/phi1;
                      (b) G4, Gcal = G4 - 2X G4_X - 4X^2 F4 y H = G4 - 2X G4_X sobre el fondo (cero doble de H en la
                      garganta; Gcal <= G4: sector impar subluminal); (c) W(l), cuya derivada es N^2 K/2: tiende a
                      3M en la masiva y a 0 en una ultraestatica (b = 4, n = 2).
  fig2_monopolo.png   (a) K y P de la masiva (positivos: sin fantasmas ni inestabilidad de gradiente); (b) K de la
                      ultraestatica: cambia de signo con integral nula (fantasma obligado); (c) U de la masiva y su
                      parte -2 N Gcal r'^2/r^2, que fija el signo de int U (taquion por la funcion de prueba d = cte).
  fig3_espectro.png   (a) autovalores negativos del monopolo de la masiva por paridad (taquiones rapido y lentos);
                      (b) campo debil: |w0| y la cota variacional |int U/int K| frente a M; (c) numero de taquiones
                      frente a M.
Uso: python figuras_fase2f.py
"""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
AZUL, NARANJA, AGUA = "#2a78d6", "#eb6834", "#1baf7a"       # paleta categorica validada (3 series)
TINTA, TINTA2, REJILLA = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2,
                     "ytick.color": TINTA2, "axes.grid": True, "grid.color": REJILLA, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "font.size": 10.5,
                     "legend.frameon": False})
os.makedirs(os.path.join(AQUI, "figuras"), exist_ok=True)
D = np.load(os.path.join(AQUI, "resultados", "perfiles_bH.npz"))
R = json.load(open(os.path.join(AQUI, "resultados", "gargantas_bH.json")))
l = D["l"]
M = R["masiva"]["x"]["M"]

# ------------------------------------------------------------------ figura 1
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.4))
m = l <= 20
a1.plot(l[m], D["m_r1"][m], color=AZUL, lw=2, label=r"$r'$")
a1.plot(l[m], D["m_N"][m], color=NARANJA, lw=2, label=r"lapso $N=e^\Phi$")
a1.plot(l[m], D["m_P1"][m] / D["m_P1"][0], color=AGUA, lw=2, label=r"$\phi'/\phi'(0)$")
a1.set_xlabel(r"$\ell/r_0$")
a1.set_title("(a) Garganta masiva de referencia", fontsize=11, color=TINTA)
a1.legend(loc="center right")
a1.set_ylim(0, 1.08 * max(1.0, float(np.max(D["m_r1"][l <= 20]))))
m = l <= 12
a2.plot(l[m], D["m_G4"][m], color=AZUL, lw=2, label=r"$G_4$")
a2.plot(l[m], D["m_Gc"][m], color=NARANJA, lw=2, label=r"$\mathcal{G}=G_4-2XG_{4X}-4X^2F_4$")
a2.plot(l[m], D["m_H"][m], color=AGUA, lw=2, label=r"$\mathcal{H}=G_4-2XG_{4X}$")
a2.set_xlabel(r"$\ell/r_0$")
a2.set_title("(b) Funciones de la teoría sobre el fondo", fontsize=11, color=TINTA)
a2.legend(loc="center right")
a2.annotate(r"$\mathcal{H}$: cero doble en la garganta", xy=(0.08, 0.01), xytext=(3.4, 0.1), color=TINTA2,
            arrowprops=dict(arrowstyle="->", color=TINTA2, lw=0.8))
m = l <= 40
a3.plot(l[m], D["m_W"][m], color=AZUL, lw=2, label=f"masiva ($M={M:g}\\,r_0$)")
a3.plot(l[m], D["u_W"][m], color=NARANJA, lw=2, label="ultraestática ($M=0$)")
a3.axhline(3 * M, color=TINTA2, lw=0.9, ls=":")
a3.text(39, 3 * M + 0.12, r"$3M$", ha="right", va="bottom", color=TINTA2)
a3.axhline(0, color=TINTA2, lw=0.8)
a3.set_xlabel(r"$\ell/r_0$")
a3.set_ylabel(r"$W$")
a3.set_title(r"(c) $N^2K/2=W'$: $\int N^2K=6(M_++M_-)$", fontsize=11, color=TINTA)
a3.legend(loc="center right")
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figuras", "fig1_fondo.png"), dpi=150)
plt.close(fig)

# ------------------------------------------------------------------ figura 2
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.4))
m = l <= 40
a1.semilogy(l[m], D["m_K"][m], color=AZUL, lw=2, label=r"$K$ (cinético)")
a1.semilogy(l[m], D["m_P"][m], color=NARANJA, lw=2, ls="--", label=r"$P$ (gradiente)")
a1.set_xlabel(r"$\ell/r_0$")
a1.set_title("(a) Monopolo de la masiva: $K>0$, $P>0$", fontsize=11, color=TINTA)
a1.legend(loc="upper right")
m = l <= 14
K = D["u_K"][m]
a2.plot(l[m], K, color=AZUL, lw=2, label=r"$K$ (ultraestática)")
a2.fill_between(l[m], 0, K, where=K > 0, color=AZUL, alpha=0.18, lw=0)
a2.fill_between(l[m], 0, K, where=K < 0, color=NARANJA, alpha=0.25, lw=0)
a2.axhline(0, color=TINTA2, lw=0.8)
a2.set_xlabel(r"$\ell/r_0$")
a2.set_title(r"(b) Ultraestática: $\int_0^\infty K\,d\ell=0$ (fantasma)", fontsize=11, color=TINTA)
a2.text(0.97, 0.06, "zona sombreada naranja: $K<0$", transform=a2.transAxes, ha="right", color=TINTA2)
a2.set_ylim(0.025, 400)
a2.legend(loc="upper right")
m = l <= 15
a3.plot(l[m], D["m_U"][m], color=AZUL, lw=2, label=r"$U$")
lead = -2 * D["m_N"] * D["m_Gc"] * D["m_r1"] ** 2 / D["m_r"] ** 2
a3.plot(l[m], lead[m], color=NARANJA, lw=2, ls="--", label=r"$-2N\mathcal{G}r'^2/r^2$")
a3.axhline(0, color=TINTA2, lw=0.8)
a3.set_xlabel(r"$\ell/r_0$")
r_m = R["masiva"]
a3.set_title(r"(c) $\int U\,d\ell<0$ $\Rightarrow$ taquión (prueba $d=$ cte)", fontsize=11, color=TINTA)
a3.text(0.97, 0.06, f"$\\int_0^\\infty U={r_m['int_U']:.2f}$,  $-2\\int N\\mathcal{{G}}r'^2/r^2={r_m['I1']:.2f}$",
        transform=a3.transAxes, ha="right", color=TINTA2)
a3.legend(loc="center right")
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figuras", "fig2_monopolo.png"), dpi=150)
plt.close(fig)

# ------------------------------------------------------------------ figura 3
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 4.4))
wp = [-w for w in r_m["w_par"]]
wi = [-w for w in r_m["w_impar"]]
a1.scatter(np.zeros(len(wp)), wp, s=60, color=AZUL, label="paridad par", zorder=3, edgecolor="white", linewidth=1.5)
a1.scatter(np.ones(len(wi)), wi, s=60, color=NARANJA, label="paridad impar", zorder=3, edgecolor="white", linewidth=1.5)
for x0, ws in ((0, wp), (1, wi)):
    for w in ws:
        a1.text(x0 + 0.08, w, f"$\\omega^2r_0^2={-w:.2g}$", va="center", color=TINTA2)
a1.set_yscale("log")
a1.set_xlim(-0.4, 1.9)
a1.set_xticks([0, 1])
a1.set_xticklabels(["par", "impar"])
a1.set_ylabel(r"$|\omega^2|\,r_0^2$ (taquiones)")
a1.set_title("(a) Monopolo de la masiva: 1 rápido, 2 lentos", fontsize=11, color=TINTA)
o2 = r_m["impar_l2"]["nodos"]
o3 = r_m["impar_l3"]["nodos"]
a1.set_ylim(2e-3, 0.4)
a1.text(0.97, 0.96, f"sector impar $l=2,3$: {sum(o2)} y {sum(o3)} taquiones", transform=a1.transAxes, ha="right",
        va="top", color=TINTA2)
CD = [r for r in R["campo_debil"]]
Ms = np.array([r["x"]["M"] for r in CD])
ok = np.array(["n_par" in r for r in CD])
w0 = np.array([-min(r["w_par"] + r["w_impar"]) if "n_par" in r else np.nan for r in CD])
cota = np.array([-r["int_U_mitad"] / r["int_K_mitad"] for r in CD])
a2.loglog(Ms[ok], w0[ok], "o-", color=AZUL, lw=2, ms=8, mec="white", mew=1.5, label=r"$|\omega_0^2|$ (el más rápido)")
a2.loglog(Ms, cota, "s--", color=NARANJA, lw=2, ms=8, mec="white", mew=1.5, label=r"cota $|\int U/\int K|$")
mm = np.geomspace(Ms.min(), Ms.max(), 10)
a2.loglog(mm, cota[ok][0] * Ms[ok][0] / mm, ":", color=TINTA2, lw=1, label=r"$\propto 1/M$")
a2.set_xlabel(r"$M/r_0$")
a2.set_ylabel(r"$|\omega^2|\,r_0^2$")
a2.set_title("(b) Campo débil: el taquión se acelera", fontsize=11, color=TINTA)
a2.set_ylim(0.025, 400)
a2.legend(loc="upper right")
npar = np.array([r.get("n_par", np.nan) for r in CD], dtype=float)
nimp = np.array([r.get("n_impar", np.nan) for r in CD], dtype=float)
x = np.arange(len(Ms))
a3.bar(x[ok], npar[ok], color=AZUL, width=0.6, label="par", edgecolor="white", linewidth=2)
a3.bar(x[ok], nimp[ok], bottom=npar[ok], color=NARANJA, width=0.6, label="impar", edgecolor="white", linewidth=2)
for i in x[~ok]:
    a3.text(i, 0.2, "fantasma\n(K<0)", ha="center", va="bottom", color=TINTA2, fontsize=9)
a3.set_xticks(x)
a3.set_xticklabels([f"{v:g}" for v in Ms])
a3.set_xlabel(r"$M/r_0$")
a3.set_ylabel("número de taquiones")
a3.set_title("(c) Taquiones del monopolo frente a $M$", fontsize=11, color=TINTA)
a3.legend(loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(AQUI, "figuras", "fig3_espectro.png"), dpi=150)
plt.close(fig)
print("figuras: fig1_fondo.png, fig2_monopolo.png, fig3_espectro.png")
