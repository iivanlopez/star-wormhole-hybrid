# -*- coding: utf-8 -*-
"""
Figura 5 de la Fase 2d (Parte III): la demostracion general con chi = r r' - 2 F_l (3 r'^2 - 1).
  izquierda: chi/|chi_min| (arriba) y K/norma (abajo) frente a l/l_chi para tres agujeros (Maxwell/Ellis, hiperbola
             afilada, Born-Infeld); en el minimo de chi K = b chi < 0, y en el cero de chi K = -a chi' <= 0
  derecha:   K/norma en el cero de chi frente a la carga escalar D (hiperbolas): tiende a -1/2 cuando D -> 0
Uso: python figuras_chi.py
"""
import os
import sys
import io
import contextlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
with contextlib.redirect_stdout(io.StringIO()):
    import teorema_chi as T     # noqa: E402
import cinetico_nled as CN      # noqa: E402

AZUL, NARANJA, AGUA = "#2a78d6", "#eb6834", "#1baf7a"       # paleta categorica validada (3 series)
TINTA, TINTA2, REJILLA = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2,
                     "ytick.color": TINTA2, "axes.grid": True, "grid.color": REJILLA, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False})
os.makedirs("figuras", exist_ok=True)


def hiperbola(L, Dv):
    l = T.malla(400.0 * max(1.0, L))
    return (l, *CN.forma("hiperbola", L, l))


casos = [("Maxwell (Ellis), $D=2$", hiperbola(1.0, 2.0), 2.0, AZUL, "-"),
         ("hipérbola $L=0.03$, $D=2$", hiperbola(0.03, 2.0), 2.0, NARANJA, "--"),
         ("Born–Infeld $b=2$, $D=1$", T.forma_nled("BI", 2.0, 1.0), 1.0, AGUA, "-.")]

fig = plt.figure(figsize=(12, 5.2))
gs = fig.add_gridspec(2, 2, width_ratios=(1.15, 1), hspace=0.08, wspace=0.28)
ax_chi = fig.add_subplot(gs[0, 0])
ax_K = fig.add_subplot(gs[1, 0], sharex=ax_chi)
ax_D = fig.add_subplot(gs[:, 1])

for nombre, (l, rr, rp, rpp, Bv), Dv, col, ls in casos:
    f = T.campos(l, rr, rp, rpp, Bv, Dv, 0.0)
    lchi = T.cruces(l, f["chi"], sube=True)[0]
    lm = T.cruces(l, f["chip"], sube=True)[0]
    chim = float(np.interp(lm, l, f["chi"]))
    x = l / lchi
    k = (x > 1e-2) & (x < 1e2)
    ax_chi.plot(x[k], f["chi"][k] / abs(chim), color=col, ls=ls, lw=1.6, label=nombre)
    ax_chi.plot(lm / lchi, -1.0, "o", color=col, ms=6, mec="white", mew=1.2)
    kn = f["K"] / f["norma"]
    ax_K.plot(x[k], kn[k], color=col, ls=ls, lw=1.6)
    ax_K.plot(lm / lchi, np.interp(lm, l, kn), "o", color=col, ms=6, mec="white", mew=1.2)
    ax_K.plot(1.0, np.interp(lchi, l, kn), "s", color=col, ms=6, mec="white", mew=1.2)
for ax in (ax_chi, ax_K):
    ax.axhline(0, color=TINTA2, lw=0.8)
    ax.axvline(1.0, color=TINTA2, lw=0.8, ls=":")
    ax.set_xscale("log")
ax_chi.set_ylabel(r"$\chi/|\chi_{\min}|$")
ax_chi.set_ylim(-1.3, 3.0)
ax_chi.legend(fontsize=8, loc="upper left", frameon=False)
ax_chi.set_title(r"$\chi=rr'-2F_\ell(3r'^2-1)$: mínimo negativo ($\bullet$) y cero ascendente (línea de puntos)",
                 fontsize=9, color=TINTA)
plt.setp(ax_chi.get_xticklabels(), visible=False)
ax_K.set_ylabel(r"$K/(A^2r^2+16F'^2B|E|)$")
ax_K.set_xlabel(r"$\ell/\ell_\chi$  ($\ell_\chi$: primer cero ascendente de $\chi$)")
ax_K.set_ylim(-1.05, 1.05)
ax_K.text(0.012, -0.70, r"$\bullet$ mínimo de $\chi$: $K=b\,\chi_{\min}<0$", fontsize=8, color=TINTA2)
ax_K.text(0.012, -0.86, r"$\blacksquare$ cero de $\chi$: $K=-a\chi'\leq0$", fontsize=8, color=TINTA2)

# ------------------------------------------------------------------ derecha: K/norma en el cero de chi frente a D
Ds = np.logspace(-2, 2, 33)
for L, col, ls, mk in ((0.03, NARANJA, "--", "^"), (1.0, AZUL, "-", "o"), (30.0, AGUA, "-.", "s")):
    l = T.malla(400.0 * max(1.0, L))
    rr, rp, rpp, Bv = CN.forma("hiperbola", L, l)
    kz = [T.lado(l, rr, rp, rpp, Bv, float(Dv), 0.0)["K_norma_en_l_chi"] for Dv in Ds]
    etiqueta = f"hipérbola $L={L:g}$" + (" (Maxwell)" if L == 1.0 else "")
    ax_D.plot(Ds, kz, color=col, ls=ls, lw=1.6, marker=mk, ms=4, mfc="white", label=etiqueta)
ax_D.axhline(-0.5, color=TINTA, lw=0.9, ls=":")
ax_D.text(1.3e-2, -0.30, r"línea de puntos: $-1/2$, límite $D\to0$ ($q=4r^2\chi'/D^2\to2$)", fontsize=8, color=TINTA2)
ax_D.text(1.3e-2, -0.36, "(la optimización global de la Parte II se quedó en $-0.497$)", fontsize=8, color=TINTA2)
ax_D.set_xscale("log")
ax_D.set_ylim(-1.05, 0.05)
ax_D.set_xlabel(r"carga escalar $D\,r_0$")
ax_D.set_ylabel(r"$K/\mathrm{norma}$ en el cero de $\chi$")
ax_D.set_title(r"En el cero de $\chi$, $K/\mathrm{norma}=-q/(2+q)\leq0$; tiende a $-1/2$ si $D\to0$", fontsize=9,
               color=TINTA)
ax_D.legend(fontsize=8, loc="lower left", frameon=False)
fig.savefig("figuras/fig5_chi.png", dpi=130, bbox_inches="tight")
print("figura: figuras/fig5_chi.png")
