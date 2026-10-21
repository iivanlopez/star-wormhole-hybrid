# -*- coding: utf-8 -*-
"""
Figuras de la Fase 2 (PNG estaticos, modo claro).
Paleta de referencia del skill dataviz, ranuras 1-3 validadas (todas las parejas):
  Modelo L (ligado) = azul, fantasma dinamico (DFKK) = naranja, estrella ordinaria = aguamarina.
Uso: python figuras.py   (requiere resultados/h1.json de test_h1.py)
"""
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hibrido import Politropo, Fondo, RHO_GEOM, PI8, MSUN_KM

AZUL, NARANJA, AGUA = "#2a78d6", "#eb6834", "#1baf7a"
TINTA, TINTA2, MUDO, REJILLA, EJE, FONDO = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
C_KMS = 3.0e5
plt.rcParams.update({
    "figure.facecolor": FONDO, "axes.facecolor": FONDO, "savefig.facecolor": FONDO,
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"], "font.size": 10,
    "axes.edgecolor": EJE, "axes.labelcolor": TINTA2, "xtick.color": MUDO, "ytick.color": MUDO,
    "axes.grid": True, "grid.color": REJILLA, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2.0,
    "legend.frameon": False, "legend.labelcolor": TINTA2, "axes.titlecolor": TINTA, "axes.titlesize": 11,
})
eos = Politropo(k=0.1, gamma=2.0)


def fondo(xi0, rho):
    rb = rho * RHO_GEOM
    return Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()


def fig_perfil():
    f = fondo(0.1, 6e14)
    l = np.concatenate([np.linspace(0, f.ls, 600), np.linspace(f.ls, f.ls + 15, 150)])
    Y = np.array([f.y_bg(x) for x in l])
    eps = np.array([f.bg(y, x <= f.ls)[12] for x, y in zip(l, Y)])
    rho_x = -0.5 * Y[:, 5] ** 2                         # 8 pi rho del sector exotico (fantasma sin masa)
    L = np.concatenate([-l[::-1], l])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(L, np.concatenate([Y[::-1, 0], Y[:, 0]]), color=AZUL)
    ax[0].axvspan(-f.ls, f.ls, color=AZUL, alpha=0.06, lw=0)
    ax[0].annotate(f"garganta r0 = {f.r0:.3f} km", (0, f.r0), (6, f.r0 + 4), color=TINTA2,
                   arrowprops=dict(arrowstyle="-", color=MUDO, lw=0.8))
    ax[0].text(0, 31.5, "region con fluido (ambos lados)", color=TINTA2, ha="center", fontsize=9)
    ax[0].set(xlabel="longitud propia l [km]", ylabel="radio areal r [km]",
              title="Geometria Z2: dos estrellas unidas por la garganta")
    ax[1].plot(L, np.concatenate([PI8 * eps[::-1], PI8 * eps]), color=AZUL, label="fluido  8 pi eps")
    ax[1].plot(L, np.concatenate([rho_x[::-1], rho_x]), color=NARANJA, label="sector exotico  8 pi rho_x")
    ax[1].axhline(0, color=EJE, lw=0.8)
    ax[1].set_xlim(-12, 12)
    ax[1].set(xlabel="longitud propia l [km]", ylabel="densidad de energia x 8 pi [km^-2]",
              title="Nucleo: energia negativa concentrada en la garganta")
    ax[1].annotate(f"{rho_x[0]:.2f} km^-2 en la garganta", (0, rho_x[0]), (1.5, rho_x[0] + 0.05),
                   color=TINTA2, fontsize=9)
    ax[1].text(-11.5, -0.07, f"fluido en la garganta: {PI8 * eps[0]:.3f} km^-2", color=TINTA2, fontsize=9)
    ax[1].legend(loc="center right")
    fig.suptitle(f"Hibrido DFKK (EOS2, xi0 = 0.1, rho_bc = 6e14 g/cm3): M = {f.M / MSUN_KM:.3f} Msun,"
                 f" R = {f.Rs:.2f} km", color=TINTA, fontsize=11)
    fig.tight_layout()
    fig.savefig("figuras/fig1_perfil.png", dpi=150)


def fig_rapido(res, rap, form):
    """Modo inestable mas rapido frente a las leyes locales de la garganta, desplazadas al rojo:
    fantasma -> Ellis-Bronnikov, -1.397 e^{-2 Phi_inf} / r0^2;  Modelo L -> formula (5)."""
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.0), sharey=True)
    aro = dict(ls="none", ms=11, mfc="none", mec=TINTA, mew=1.0)
    for a, xi0 in zip(ax, (0.1, 0.3, 0.6)):
        g = sorted((o for o in res if o.get("xi0") == xi0 and o.get("fantasma")), key=lambda o: o["rho"])
        l = sorted((o for o in rap if o["xi0"] == xi0 and o["modelo"] == "ligado" and o["raices"]),
                   key=lambda o: o["rho"])
        fo = sorted((o for o in form if o["xi0"] == xi0), key=lambda o: o["rho"])
        a.plot([o["rho"] for o in g], [-o["fantasma"][0] for o in g], color=NARANJA, marker="s", ms=6,
               label="fantasma dinamico (DFKK)")
        a.plot([o["rho"] for o in l], [-min(o["raices"]) for o in l], color=AZUL, marker="o", ms=6,
               label="Modelo L: modo de la garganta")
        a.plot([o["rho"] for o in g], [1.39705 * math.exp(-2 * o["Phi_inf"]) / o["r0"] ** 2 for o in g],
               marker="s", label="prediccion: Ellis-Bronnikov x e^(-2 Phi_inf)", **aro)
        a.plot([o["rho"] for o in fo], [-o["w_formula_fis"] for o in fo], marker="o",
               label="prediccion: formula (5)", **aro)
        a.set_xscale("log")
        a.set_yscale("log")
        a.set_ylim(1e-3, 1e3)
        a.set_title(f"xi0 = {xi0}  (r0 = {fo[0]['r0']:.2f} km)" if fo else f"xi0 = {xi0}")
        a.set_xlabel("rho_bc [g/cm3]")
    ax[0].set_ylabel("|omega^2| del modo inestable [km^-2]")
    ax[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Ambas dinamicas son inestables; la ligada lo es mas cuanto menos fluido hay "
                 "(tiempo de e-pliegue = 1/(c sqrt|omega^2|))", color=TINTA, fontsize=11)
    fig.tight_layout()
    fig.savefig("figuras/fig2_modo_rapido.png", dpi=150)


def fig_h1(res):
    seq = sorted((o for o in res if o.get("xi0") == 0.1 and "error" not in o), key=lambda o: o["rho"])
    rho = np.array([o["rho"] for o in seq])
    wg = np.array([o["fantasma"][0] if o["fantasma"] else np.nan for o in seq])
    wl = np.array([o["ligado"][0] if o["ligado"] else np.nan for o in seq])
    wn = np.array([o["ns"] if o["ns"] is not None else np.nan for o in seq])
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))
    a = ax[0]
    a.set_yscale("symlog", linthresh=1e-5)
    a.plot(rho, wg, color=NARANJA, marker="s", ms=6, label="fantasma dinamico (DFKK)")
    a.plot(rho, wl, color=AZUL, marker="o", ms=6, label="Modelo L (sector ligado, NOF)")
    a.plot(rho, wn, color=AGUA, marker="^", ms=6, lw=1.5, label="estrella ordinaria (sin garganta)")
    a.axhline(0, color=TINTA2, lw=0.8)
    a.set_xscale("log")
    a.set(xlabel="densidad en la garganta rho_bc [g/cm3]", ylabel="omega^2 [km^-2]  (escala symlog)",
          title="xi0 = 0.1: modo lento del fluido (el de la garganta: Fig. 2)")
    a.text(1.3e13, -0.25, "fantasma: inestable (microsegundos)", color=TINTA2, fontsize=9)
    a.legend(loc="upper left", bbox_to_anchor=(0.0, 0.62), fontsize=8.5)
    b = ax[1]
    b.plot(rho, [o["M"] for o in seq], color=AZUL, marker="o", ms=6, label="hibrido (xi0 = 0.1)")
    b.plot(rho, [o["M_ns"] for o in seq], color=AGUA, marker="^", ms=6, lw=1.5, label="estrella ordinaria")
    b.set_xscale("log")
    b.set(xlabel="densidad en la garganta rho_bc [g/cm3]", ylabel="masa total M [Msun]",
          title="Secuencias de equilibrio (punto de retorno)")
    b.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig("figuras/fig3_modos_lentos.png", dpi=150)


def autofuncion(f, modelo, w, n=500):
    F1, F2, sols = f.condiciones(modelo, w, devolver=True)
    c = [1.0] if len(F2) == 1 else [F1[1], -F1[0]]
    ls, vals = [], []
    for sol in sols:
        x = np.linspace(sol.t[0], sol.t[-1], n)
        Y = sol.sol(x)
        p = sum(ck * Y[6 * k:6 * k + 6] for k, ck in enumerate(c))
        ls.append(x)
        vals.append(p)
    return np.concatenate(ls), np.concatenate(vals, axis=1)


def fig_modos(res):
    from hibrido import autovalores
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    e = Fondo(r0=1.0).resolver(l_max=200.0)
    r_, _ = autovalores(e, "fantasma", np.linspace(-1.6, -1.2, 5))
    l, p = autofuncion(e, "fantasma", r_[0])
    l, p = l[l < 12], p[:, l < 12]
    ax[0].plot(l, p[2] / p[2][0], color=AZUL, label="dr / dr(0)  (radio areal)")
    ax[0].plot(l, p[4] / np.max(np.abs(p[4])), color=NARANJA, lw=1.5, label="d phi~ (fantasma, normalizado)")
    ax[0].axhline(0, color=EJE, lw=0.8)
    ax[0].set(xlabel="l / r0", ylabel="amplitud relativa",
              title=f"Ellis-Bronnikov + fantasma: omega^2 r0^2 = {e.w_fisico(r_[0]):.4f}")
    ax[0].legend(fontsize=8.5)
    d = np.load("resultados/modo_rapido_autofuncion.npz")
    l, X, s, Phi, C = d["l"], d["X"], d["s"], d["Phi"], float(d["C"])
    m = l <= 0.97 * float(d["ls"])                     # se excluye la capa superficial evanescente
    l, xi, dr, s, Phi = l[m], X[0][m], X[2][m], s[m], Phi[m]
    ax[1].plot(l, dr / (C * np.exp(Phi)), color=AZUL, label="dr / (C e^Phi): deslizamiento del andamio")
    ax[1].plot(l, (dr + s * xi) / C, color=AGUA, label="Delta R / C: radio areal de cada capa de fluido")
    ax[1].plot(l, -xi / C, color=TINTA2, lw=1.3, label="-xi / C: fluido respecto del andamio")
    ax[1].axhline(0, color=EJE, lw=0.8)
    ax[1].annotate(f"garganta: dr(0) = {dr[0] / C:.2f} C", (0, dr[0] / C), (4, 1.25), color=TINTA2, fontsize=9,
                   arrowprops=dict(arrowstyle="-", color=MUDO, lw=0.8))
    ax[1].set_ylim(-0.1, 1.45)
    ax[1].set(xlabel="longitud propia l [km]", ylabel="amplitud / C",
              title=f"Modelo L + fluido, garganta: omega^2 = {float(d['w']):.3f} km^-2")
    ax[1].legend(fontsize=8, loc="center right")
    fig.tight_layout()
    fig.savefig("figuras/fig4_autofunciones.png", dpi=150)
    return dr[0] / C


if __name__ == "__main__":
    import os
    os.makedirs("figuras", exist_ok=True)
    res = json.load(open("resultados/h1.json", encoding="utf-8"))
    rap = json.load(open("resultados/modo_rapido.json", encoding="utf-8"))
    form = json.load(open("resultados/formula_rapida.json", encoding="utf-8"))
    fig_perfil()
    fig_rapido(res, rap, form)
    fig_h1(res)
    print("dr(0)/C del modo rapido:", fig_modos(res))
