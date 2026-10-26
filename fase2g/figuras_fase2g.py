# -*- coding: utf-8 -*-
"""
Figuras de la Fase 2g (figuras/):
  fig1_teorema.png   (a) cota de Rayleigh con d = N frente al autovalor mas bajo de la Fase 2f (todas las gargantas
                     con K, P > 0); (b) la cota de la Fase 2f con d = cte, que falla en campo fuerte.
  fig2_fluido.png    criba de los hibridos con fluido: (a) margenes de las muestras aleatorias; (b) perfiles de la
                     mejor configuracion (lapso, radio, Gcal, G4, densidad, velocidad del sonido).
  fig3_densidad.png  densidades de energia de la prueba d = N: deficit -2 N^3 Gcal sig^2 y rigidez del fluido
                     2 N^3 (rho~+p~) c_s^2 en la mejor configuracion.
Uso: python figuras_fase2g.py
"""
import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt            # noqa: E402
import sys                                 # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

AQUI = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(AQUI, "resultados")
FIG = os.path.join(AQUI, "figuras")
os.makedirs(FIG, exist_ok=True)

# ------------------------------------------------------------------ figura 1: teorema del vacio
v = [r for r in json.load(open(os.path.join(RES, "verificar_N.json"))) if r.get("valida", True)]
qN = np.array([r["qN"] for r in v])
q1 = np.array([r["q1"] for r in v])
w0 = np.array([r["w0"] for r in v])
N0 = np.array([r["N0"] for r in v])
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
sc = ax[0].scatter(-w0, -qN, c=N0, cmap="viridis", s=28)
x = np.logspace(np.log10(min(-w0.min(), -qN.min()) * 0.5), np.log10(max(-w0.max(), -qN.max()) * 2), 10)
ax[0].plot(x, x, "k--", lw=0.8, label=r"$|q_N| = |\omega_0^2|$")
ax[0].set_xscale("log")
ax[0].set_yscale("log")
ax[0].set_xlabel(r"$|\omega_0^2|\,r_0^2$ (autovalor más bajo, Fase 2f)")
ax[0].set_ylabel(r"$|q_N| = |E[N]|/\int K N^2$")
ax[0].set_title(r"(a) prueba $d=N$: $\omega_0^2 \leq q_N < 0$ siempre")
ax[0].legend(fontsize=8)
plt.colorbar(sc, ax=ax[0], label=r"$N(0)$")
ax[1].scatter(N0, q1, s=22, label=r"$q_1=\int U/\int K$ ($d$ = cte, Fase 2f)")
ax[1].scatter(N0, qN, s=22, marker="x", label=r"$q_N$ ($d=N$, este trabajo)")
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_xlabel(r"$N(0)$ (lapso en la garganta)")
ax[1].set_ylabel("cota variacional")
ax[1].set_yscale("symlog", linthresh=1e-3)
ax[1].set_title(r"(b) $d$ = cte falla en campo fuerte; $d=N$ no")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig1_teorema.png"), dpi=140)

# ------------------------------------------------------------------ figura 2: criba con fluido
claves = ("E1", "Emin_par", "Emin_impar", "Tmin_par", "Tmin_impar")
mues = []
for sem in (1, 2):
    p = os.path.join(RES, f"criba_busqueda_{sem}.json")
    if os.path.exists(p):
        mues += [r for r in json.load(open(p))["resultados"] if "E1" in r and r.get("N0", 0) > 0.05]
pf = np.load(os.path.join(RES, "perfiles_criba.npz"))
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
e1 = np.array([np.sign(r["E1"]) for r in mues])
fam = np.array([min(r["Emin_par"], r["Emin_impar"]) for r in mues])
fan = np.array([r.get("margen_fantasma", np.nan) for r in mues])
rp0 = np.array([r["rp0"] for r in mues])
ok = np.array([r["E1"] >= 0 and r["Emin_par"] >= 0 and r["Emin_impar"] >= 0 for r in mues])
ax[0].scatter(rp0[~ok], fan[~ok], s=8, c="0.6", label="falla la familia $c=0$")
ax[0].scatter(rp0[ok], fan[ok], s=30, c="C3", label="pasa la familia $c=0$")
mej = json.load(open(os.path.join(RES, "criba_mejores.json")))
for i, (k, m) in enumerate(mej.items()):
    F_ = m["fantasma"]
    from criba import FondoFluido                                          # noqa: E402
    rp0_ = FondoFluido(**m["x"]).campos(np.array([1e-6]))["rp"][0]
    ax[0].scatter([rp0_], [F_["presupuesto"] - F_["gasto"]], marker="*", s=220, c="C2", edgecolors="k",
                  label="optimización (pasa todas)" if i == 0 else None, zorder=5)
ax[0].axhline(0, color="k", lw=0.6)
ax[0].set_xscale("log")
ax[0].set_yscale("symlog", linthresh=0.1)
ax[0].set_xlabel(r"$(\tilde\rho+\tilde p)(0)\,r_0^2$")
ax[0].set_ylabel(r"margen del presupuesto de fantasma $3M-$gasto")
ax[0].set_title("(a) búsqueda aleatoria (1600 muestras)")
ax[0].legend(fontsize=8)
l = pf["s2_l"]
ls = float(pf["s2_ls"][0])
for k, lab in (("N", "$N$"), ("G4", "$G_4$"), ("Gc", r"$\mathcal{G}$"), ("cs2", "$c_s^2$")):
    ax[1].plot(l, pf[f"s2_{k}"], label=lab)
ax[1].plot(l, pf["s2_rp"] / pf["s2_rp"][0], label=r"$(\rho+p)/(\rho+p)_0$")
ax[1].axvline(ls, color="k", ls=":", lw=0.8)
ax[1].set_xlim(0, 1.6 * ls)
ax[1].set_xlabel(r"$\ell/r_0$")
ax[1].set_title("(b) mejor híbrido (semilla 2): perfiles")
ax[1].legend(fontsize=8)
ax[2].plot(l, pf["s2_r"], label="$r$")
ax[2].plot(l, l, "k--", lw=0.6)
ax[2].axvline(ls, color="k", ls=":", lw=0.8)
ax[2].set_xlim(0, 1.6 * ls)
ax[2].set_ylim(0, pf["s2_r"][np.searchsorted(l, 1.6 * ls)] * 1.1)
ax[2].set_xlabel(r"$\ell/r_0$")
ax[2].set_ylabel(r"$r/r_0$")
ax[2].set_title("(c) radio areal: garganta en el centro de la estrella")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig2_fluido.png"), dpi=140)

# ------------------------------------------------------------------ figura 3: densidades de la prueba d = N
fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
for i, sem in enumerate((1, 2)):
    l = pf[f"s{sem}_l"]
    N = pf[f"s{sem}_N"]
    defi = -2 * N ** 3 * pf[f"s{sem}_Gc"] * pf[f"s{sem}_sig"] ** 2
    rig = 2 * N ** 3 * pf[f"s{sem}_rp"] * pf[f"s{sem}_cs2"]
    ax[i].plot(l, defi, label=r"$-2N^3\mathcal{G}\sigma^2$ (déficit, vacío)")
    ax[i].plot(l, rig, label=r"$2N^3(\tilde\rho+\tilde p)c_s^2$ (rigidez del fluido)")
    ax[i].plot(l, defi + rig, "k", lw=1, label="suma")
    ax[i].axhline(0, color="k", lw=0.5)
    ax[i].axvline(float(pf[f"s{sem}_ls"][0]), color="k", ls=":", lw=0.8)
    ax[i].set_xlim(0, 1.6 * float(pf[f"s{sem}_ls"][0]))
    ax[i].set_xlabel(r"$\ell/r_0$")
    ax[i].set_title(f"semilla {sem}: densidad de $E[N,\\xi=0]$")
    ax[i].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig3_densidad.png"), dpi=140)
print("figuras escritas en", FIG)
