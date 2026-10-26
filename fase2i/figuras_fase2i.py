# -*- coding: utf-8 -*-
"""
Figuras de la Fase 2i (figuras/), a partir de resultados/analisis_fase2i.json, perfiles_fase2i.npz y
optimiza_solido_<semilla>.json:
  fig1_cizalla.png     (a) margenes de la criba c = 0 frente a c_T^2 en la estrella de la Fase 2h hecha solido (c_L^2 = 1):
                       la cizalla arregla el sector impar, no el par, y empeora E[1]; (b) el modo par (c_T^2 = 0 y 0.2).
  fig2_traslacion.png  (a) desviacion de la relacion de RG nu = r w/(4 Gcal) en el material; (b) energia acumulada de la
                       prueba h = r'/N (vacio, material y total): la traslacion de las capas exteriores no cuesta nada.
  fig3_cascaras.png    (a) cociente Q(x) de la cascara delgada (energia de la prueba uniforme frente a |E_out|);
                       (b) margenes de las cascaras de membrana frente a r_b/M.
  fig4_optimizacion.png historia de la evolucion diferencial sobre las estrellas elasticas con nucleo hueco.
Uso: python figuras_fase2i.py
"""
import os
import json
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt            # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RES = json.load(open(os.path.join(AQUI, "resultados", "analisis_fase2i.json")))
P = np.load(os.path.join(AQUI, "resultados", "perfiles_fase2i.npz"))
FIG = os.path.join(AQUI, "figuras")
os.makedirs(FIG, exist_ok=True)

# ------------------------------------------------------------------ 1. cizalla
B = RES["B"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
c = [e["cT2"] for e in B]
ax[0].plot(c, [e["E1n"] for e in B], "o-", label=r"$E[1]$ (prueba uniforme)")
ax[0].plot(c, [e["Epar"] for e in B], "s-", label=r"$E_{\rm par}$")
ax[0].plot(c, [e["Eimpar"] for e in B], "^-", label=r"$E_{\rm impar}$")
ax[0].axhline(0, color="k", lw=0.6)
ax[0].set_xlabel(r"$c_T^2$ (con $c_L^2=1$)")
ax[0].set_ylabel("margen normalizado")
ax[0].set_title("(a) la cizalla arregla el sector impar, no el par")
ax[0].legend(fontsize=8)
for k, est in (("0", "-"), ("20", "--")):
    ax[1].plot(P[f"B_l_{k}"], P[f"B_h_{k}"], est, label=rf"$c_T^2={int(k) / 100:g}$")
ax[1].set_xlabel(r"$\ell/r_0$")
ax[1].set_ylabel(r"$h$ del modo par (normalizado)")
ax[1].set_title("(b) el modo par se extiende por toda la estrella")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig1_cizalla.png"), dpi=140)

# ------------------------------------------------------------------ 2. traslacion
l = P["C_l"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].semilogy(l, np.abs(P["C_desv"]))
ax[0].axhline(0.1, color="k", lw=0.6, ls=":")
ax[0].set_xlabel(r"$\ell/r_0$")
ax[0].set_ylabel(r"$|\nu/(rw/4\mathcal{G})-1|$")
ax[0].set_title("(a) el material solo está anclado cerca de la garganta")
ev, em = P["C_ev"], P["C_em"]
cum = lambda y: np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(l))])   # noqa: E731
m = l >= 20
lm = l[m]
base = lambda y: cum(y)[m] - cum(y)[m][0]                                                 # noqa: E731
ax[1].plot(lm, base(ev), label="vacío")
ax[1].plot(lm, base(em), label="material")
ax[1].plot(lm, base(ev + em), "k", label="total")
ax[1].axhline(RES["C"]["E_out_h2"], color="r", ls="--", lw=0.8, label=r"exterior $E_{\rm out}h(\ell_b)^2$")
ax[1].set_xlabel(r"$\ell/r_0$")
ax[1].set_ylabel(r"energía acumulada desde $\ell=20$")
ax[1].set_title(r"(b) prueba $h=r'/N$: la traslación no cuesta nada")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig2_traslacion.png"), dpi=140)

# ------------------------------------------------------------------ 3. cascaras
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(P["D_x"], P["D_Q"])
ax[0].axhline(1, color="k", lw=0.6)
ax[0].set_xlabel(r"$x=M/R$")
ax[0].set_ylabel(r"$Q$ = energía de la cáscara / $|E_{\rm out}|$")
ax[0].set_title(r"(a) cáscara delgada ($c_L^2=1$, $c_T^2=0.1$): $Q\lesssim1.08$")
D = RES["D"]
ax[1].plot([e["rb_M"] for e in D], [e["E1n"] for e in D], "o", label=r"$E[1]$")
ax[1].plot([e["rb_M"] for e in D], [e["Epar"] for e in D], "s", label=r"$E_{\rm par}$")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel(r"$r_b/M$")
ax[1].set_ylabel("margen normalizado")
ax[1].set_title("(b) cáscaras de membrana de masa dominante")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig3_cascaras.png"), dpi=140)

# ------------------------------------------------------------------ 4. optimizacion
fig, ax = plt.subplots(figsize=(6, 4))
for f in sorted(glob.glob(os.path.join(AQUI, "resultados", "optimiza_solido_*.json"))):
    d = json.load(open(f))
    ax.plot(np.arange(1, len(d["historia"]) + 1), [-h for h in d["historia"]], "o-", ms=3,
            label=f"semilla {os.path.splitext(os.path.basename(f))[0].split('_')[-1]}")
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("generación")
ax.set_ylabel("mejor margen")
ax.set_title("Evolución diferencial: estrellas elásticas con núcleo hueco")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig4_optimizacion.png"), dpi=140)
print("figuras escritas en", FIG)
