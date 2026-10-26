# -*- coding: utf-8 -*-
"""
Figuras de la Fase 2h (figuras/), a partir de resultados/analisis_final.json y resultados/perfiles_fase2h.npz:
  fig1_exterior.png    (a) asintotica de la teoria natural fuera de la estrella: U l^4 -> -2Q^2 y P l^2, K l^2 -> Q^2
                       (escalar sin masa, estable por Hardy aunque U < 0); (b) autovalor de Dirichlet del exterior
                       frente a la fraccion theta del presupuesto de K movida fuera.
  fig2_candidato.png   fondo del candidato rejilla_pozo: perfiles y r/N (esfera de fotones dentro de la estrella).
  fig3_espectro.png    (a) autovalor mas bajo frente a la rigidez Lam (completo y familia c = 0, superficie fija y libre);
                       (b) el modo inestable y la prueba sin compresion y(l).
  fig4_criba.png       (a) criba con superficie libre frente a la escala del muelle de superficie; (b) partes de la
                       energia de la prueba sin compresion en el candidato y en las semillas de la Fase 2g.
  fig5_velocidades.png velocidades caracteristicas en la estrella (modo rapido superluminico, divergente en la garganta).
Uso: python figuras_fase2h.py
"""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt            # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RES = json.load(open(os.path.join(AQUI, "resultados", "analisis_final.json")))
P = np.load(os.path.join(AQUI, "resultados", "perfiles_fase2h.npz"))
FIG = os.path.join(AQUI, "figuras")
os.makedirs(FIG, exist_ok=True)

# ------------------------------------------------------------------ 1. exterior
Q2 = RES["A"]["Q2"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
l = P["A_l"]
ax[0].plot(l, P["A_Ul4"] / Q2, label=r"$U\ell^4/Q^2$")
ax[0].plot(l, P["A_U0l4"] / Q2, "--", label=r"$U_0\ell^4/Q^2$ (sin el término de $P$)")
ax[0].plot(l, P["A_Pl2"] / Q2, label=r"$P\ell^2/Q^2$")
ax[0].plot(l, P["A_Kl2"] / Q2, ":", label=r"$K\ell^2/Q^2$")
ax[0].axhline(-2, color="k", lw=0.6)
ax[0].axhline(1, color="k", lw=0.6)
ax[0].set_xscale("log")
ax[0].set_xlabel(r"$\ell/r_0$")
ax[0].set_title(r"(a) teoría natural fuera: $U<0$, pero $|U|\ell^2/P\to2<9/4$")
ax[0].legend(fontsize=8)
for cm2, mk in ((1.0, "o-"), (3.0, "s--")):
    d = [e for e in RES["A"]["dirichlet"] if e["cm2"] == cm2]
    ax[1].plot([e["theta"] for e in d], [e["w0"] for e in d], mk, label=rf"$P=c^2N^2K$, $c^2={cm2:g}$")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_yscale("symlog", linthresh=1e-6)
ax[1].set_xlabel(r"$\theta$ (fracción del presupuesto natural de $K$ movida a $D>0$)")
ax[1].set_ylabel(r"$\omega_0^2$ del exterior (Dirichlet en $\ell_s$)")
ax[1].set_title("(b) mover $K$ fuera con $P$ lumínica desestabiliza")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig1_exterior.png"), dpi=140)

# ------------------------------------------------------------------ 2. candidato
B = RES["B"]
l, ls = P["B_l"], float(P["B_ls"][0])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(l, P["B_N"], label="$N$")
ax[0].plot(l, P["B_rp"] / P["B_rp"][0], label=r"$(\tilde\rho+\tilde p)/(\tilde\rho+\tilde p)_0$")
ax[0].plot(l, P["B_Gc"], label=r"$\mathcal{G}$")
ax[0].plot(l, P["B_G4"], label="$G_4$")
ax[0].plot(l, P["B_r1"], label="$r'$")
ax[0].axvline(ls, color="k", ls=":", lw=0.8)
ax[0].set_xlabel(r"$\ell/r_0$")
ax[0].set_title(rf"(a) candidato: $M={B['M_adm']:.1f}$, $\ell_s={ls:.1f}$, $r_s/M={B['rs_sobre_M']:.2f}$, $N_0={B['N0']:.3f}$")
ax[0].legend(fontsize=8)
ax[1].plot(l, P["B_r"] / P["B_N"], label=r"$r/N$")
ax[1].axvline(ls, color="k", ls=":", lw=0.8, label="superficie")
ax[1].set_xlabel(r"$\ell/r_0$")
ax[1].set_xscale("log")
ax[1].set_yscale("log")
ax[1].set_title(r"(b) $r/N$ crece desde la garganta (su única esfera de fotones)")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig2_candidato.png"), dpi=140)

# ------------------------------------------------------------------ 3. espectro
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
esp = B["espectro"]
for libre in (False, True):
    for c0 in (False, True):
        d = [e for e in esp if e["libre"] == libre and e["c_cero"] == c0 and e["par"] == "par"]
        Lam = [max(e["Lam"], 1e-3) for e in d]
        ax[0].plot(Lam, [e["w"][0] for e in d], ("o-" if libre else "s--") if not c0 else ("^-" if libre else "v--"),
                   label=f"{'completo' if not c0 else '$c=0$ dentro'}, superficie {'libre' if libre else 'fija'}")
ax[0].axhline(0, color="k", lw=0.6)
ax[0].set_xscale("log")
ax[0].set_yscale("symlog", linthresh=1e-6)
ax[0].set_xlabel(r"$\Lambda$ (rigidez de $c\neq0$ en la estrella; $10^{-3}$ = sin rigidez)")
ax[0].set_ylabel(r"$\omega_0^2 r_0^2$ (paridad par)")
ax[0].set_title("(a) el taquión residual es el de la familia $c=0$ con superficie libre")
ax[0].legend(fontsize=7)
ax[1].plot(P["B_lmodo"], P["B_hmodo"], label=r"modo inestable ($\Lambda=10$, superficie libre): $h/h(0)$")
ax[1].plot(P["B_ly"], P["B_yQ0"], "--", label=r"prueba sin compresión $y$ ($Q\equiv0$)")
ax[1].axvline(ls, color="k", ls=":", lw=0.8)
ax[1].set_xlim(0, 3 * ls)
ax[1].set_xlabel(r"$\ell/r_0$")
ax[1].set_title("(b) el fluido sigue al desplazamiento sin comprimirse")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig3_espectro.png"), dpi=140)

# ------------------------------------------------------------------ 4. criba con superficie libre
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
mu = B["muelle"]
ax[0].plot([max(m["factor"], 0.1) for m in mu], [m["Epar"] for m in mu], "o-")
ax[0].axhline(B["criba_fija"]["Epar"], color="C2", ls="--", lw=0.8, label="superficie fija (criba de la Fase 2g)")
ax[0].axhline(0, color="k", lw=0.6)
ax[0].axvline(1, color="C3", ls=":", lw=0.8, label="muelle físico")
ax[0].set_xscale("log")
ax[0].set_xlabel(r"escala del muelle de superficie ($0.1$ = sin muelle)")
ax[0].set_ylabel(r"$E_{\rm par}$ normalizada")
ax[0].set_title("(a) candidato: hace falta un muelle ~100 veces el físico")
ax[0].legend(fontsize=8)
casos = [("candidato", B["sin_compresion"])] + [(k.replace("semilla", "2g s.").replace("_m1", " m=1"), v["sin_compresion"])
                                                for k, v in RES["C"].items()]
x = np.arange(len(casos))
for j, (k, lab) in enumerate((("vacio", r"$-2\int N^3\mathcal{G}(y'-\sigma y)^2$"), ("exterior", r"$E_{\rm out}\,y(\ell_s)^2$"),
                               ("muelle", r"muelle $E_s y'(\ell_s)^2$"), ("E_Q0", r"$E_{Q0}$ (total)"))):
    ax[1].bar(x + (j - 1.5) * 0.2, [c[1][k] / abs(c[1]["E1"]) for c in casos], 0.2, label=lab)
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xticks(x)
ax[1].set_xticklabels([c[0] for c in casos], fontsize=8)
ax[1].set_ylabel(r"energía / $|E[1]|$")
ax[1].set_title(r"(b) prueba sin compresión: $E_{Q0}<0$ en todos los casos")
ax[1].legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig4_criba.png"), dpi=140)

# ------------------------------------------------------------------ 5. velocidades
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
lv = P["B_lv"]
for j, Lv in enumerate((0, 10)):
    ax[j].plot(lv, P[f"B_c2r_{Lv}"], label=r"$c^2_{\rm rápida}$")
    ax[j].plot(lv, P[f"B_c2l_{Lv}"], label=r"$c^2_{\rm lenta}$")
    ax[j].plot(lv, P[f"B_cphi2_{Lv}"], "--", label=r"$P/(N^2K)$ (escalar sin acoplo)")
    ax[j].plot(lv, P["B_cs2"], ":", label=r"$c_s^2$ (fluido sin acoplo)")
    ax[j].plot(lv, P[f"B_detK_{Lv}"], "-.", label=r"$\det\mathbb{K}/(\kappa_b K)$")
    ax[j].axhline(1, color="k", lw=0.6)
    ax[j].set_xscale("log")
    ax[j].set_yscale("log")
    ax[j].set_xlabel(r"$\ell/r_0$")
    ax[j].set_title(rf"(a) $\Lambda=0$: el acoplo fluido–escalar es superlumínico" if Lv == 0 else
                    rf"(b) $\Lambda=10$: la rigidez de $c\neq0$ dispara $P/(N^2K)$")
    ax[j].legend(fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig5_velocidades.png"), dpi=140)
print("figuras escritas en", FIG)
