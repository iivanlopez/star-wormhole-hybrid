# -*- coding: utf-8 -*-
"""
Figuras de la Parte II (familias de Antoniou et al. 2019), con la paleta de figuras.py.
Requiere resultados/familia_cuadratica.json, hiperbolicidad.json, estabilidad_robusta.json, trayectorias_robustas.json
y ley_aceleracion.json.
Uso: python figuras_antoniou.py
"""
import json
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from figuras import AZUL, NARANJA, AGUA, TINTA, TINTA2, MUDO, EJE  # noqa: F401  (tambien aplica el estilo)


def fig_familia():
    fam = json.load(open("resultados/familia_cuadratica.json", encoding="utf-8"))
    hip = {o["nombre"]: o for o in json.load(open("resultados/hiperbolicidad.json", encoding="utf-8"))}
    rob = {}
    for o in json.load(open("resultados/estabilidad_robusta.json", encoding="utf-8")):
        rob.setdefault(o["nombre"], {})[o["zeta"]] = o
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)
    for k, z in enumerate((1.0, 0.0)):
        a = ax[k]
        xs, ys, cs, ms = [], [], [], []
        for o in fam:
            nom = f"cuadratico {o['alpha_eta']:.4f}"
            x, y = o["alpha_eta"], o["D"] / o["M"]
            mk = "D" if o["ecuador"] else "o"
            if not hip[nom]["hiperbolico"]:
                a.scatter([x], [y], marker="x", s=45, color=MUDO, zorder=2)
                continue
            r = rob.get(nom, {}).get(z)
            if r is None or not r["winf"]:
                a.scatter([x], [y], marker=mk, s=60, facecolor="none", edgecolor=NARANJA, lw=1.4, zorder=3)
                continue
            xs.append(x)
            ys.append(y)
            cs.append(math.sqrt(-min(r["winf"])))
            ms.append(mk)
        sc = None
        for mk in ("o", "D"):
            sel = [i for i, m in enumerate(ms) if m == mk]
            if sel:
                sc = a.scatter([xs[i] for i in sel], [ys[i] for i in sel], c=[cs[i] for i in sel], marker=mk, s=60,
                               cmap="viridis", vmin=0.1, vmax=1.0, edgecolor=TINTA, linewidth=0.6, zorder=3)
        a.axvline(0.205, color=TINTA2, lw=0.8, ls="--")
        a.axvline(0.480, color=TINTA2, lw=0.8, ls="--")
        a.set(xlabel="alpha / eta0^2", title=f"zeta = {z:g}" + (" (capa de Antoniou et al.)" if z == 1 else
                                                                 " (sin gravedad inducida)"))
        if k == 0:
            a.set_ylabel("D / M (carga escalar / masa)")
        if sc is not None:
            cb = fig.colorbar(sc, ax=a)
            cb.set_label("tasa del modo rapido |omega_inf| r0 / c")
    ax[0].scatter([], [], marker="o", color=MUDO, label="una garganta")
    ax[0].scatter([], [], marker="D", color=MUDO, label="doble garganta (capa en el ecuador)")
    ax[0].scatter([], [], marker="x", color=MUDO, label="no hiperbolico (mal planteado)")
    h, lab = ax[0].get_legend_handles_labels()
    fig.legend(h, lab, fontsize=8, loc="lower center", ncol=4)
    for a in ax:
        a.annotate("", xy=(0.205, 0.97), xytext=(0.480, 0.97), xycoords=("data", "axes fraction"),
                   arrowprops=dict(arrowstyle="<->", color=TINTA2, lw=0.8))
        a.text(0.2575, 0.93, "intervalo publicado", color=TINTA2, fontsize=7.5, ha="center",
               transform=a.get_xaxis_transform())
    fig.suptitle("Familia cuadratica F = alpha phi^2, phi_inf = 0: todos los fondos bien planteados son inestables",
                 color=TINTA, fontsize=11)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig("figuras/fig4_antoniou.png", dpi=150)


def fig_trayectorias():
    """Izquierda: inercia efectiva 1/C = -a^2/w del modo rapido frente a zeta (lineal; se anula en zeta_c), con la recta
    1 - zeta/2 del limite de agujero negro.  Derecha: modo lento, -omega_inf^2 r0^2 frente a zeta."""
    from fondo_F import FondoF
    tr = json.load(open("resultados/trayectorias_robustas.json", encoding="utf-8"))
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))
    series = (("KKK referencia", AZUL, "s", "KKK (a r0 = 13.4)"),
              ("cuadratico 0.3087", NARANJA, "o", "cuadratico 0.309 (a r0 = 14.3)"),
              ("cuadratico 0.3465", AGUA, "^", "cuadratico 0.346 (a r0 = 64)"))
    for nom, color, mk, etiqueta in series:
        sub = sorted((o for o in tr if o["nombre"] == nom and o["zeta"] >= 0), key=lambda o: o["zeta"])
        a2 = FondoF(sub[0]["tipo"], sub[0]["alpha"], sub[0]["phi0"], sub[0]["phi1"]).y0()[3] ** 2
        rap = [(o["zeta"], a2 / -min(o["w"])) for o in sub if o["w"] and -min(o["w"]) / a2 > 0.3]
        ax[0].plot([z for z, _ in rap], [y for _, y in rap], mk + "-", color=color, ms=6, lw=1.2, label=etiqueta)
        lento = [(o["zeta"], -max(o["winf"])) for o in sub if o["w"] and -max(o["w"]) / a2 < 0.1]
        if lento:
            ax[1].plot([z for z, _ in lento], [y for _, y in lento], mk + "-", color=color, ms=6, lw=1.2, label=etiqueta)
    zz = [0, 2]
    ax[0].plot(zz, [1 - z / 2 for z in zz], color=TINTA, lw=0.8, ls="--", label="1 - zeta/2 (limite de agujero negro)")
    ax[0].axhline(0, color=EJE, lw=0.8)
    ax[0].axvline(1.0, color=MUDO, lw=0.8)
    ax[0].set(xlabel="zeta (gravedad inducida de la capa)", ylabel="1/C = -a^2/w  (inercia efectiva)", xlim=(-0.05, 2.05),
              title="Modo rapido: la inercia se anula en zeta_c -> 2")
    ax[0].legend(fontsize=7.5, loc="lower left")
    ax[1].set_yscale("log")
    ax[1].axvline(1.0, color=MUDO, lw=0.8)
    ax[1].set(xlabel="zeta (gravedad inducida de la capa)", ylabel="-omega_inf^2 r0^2 (modo lento)",
              title="Modo lento: aparece hacia zeta ~ 0.7-1 y sigue para zeta >= 2")
    ax[1].legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig("figuras/fig5_trayectorias.png", dpi=150)


def fig_limite():
    """Ley local (§12): C = -w/a^2 frente a a r0, con a = Phi'(0+) la aceleracion propia de la garganta.
    zeta = 0: modo con C mas proximo a 1 (y los demas modos, huecos); zeta = 1: modo mas rapido."""
    ley = json.load(open("resultados/ley_aceleracion.json", encoding="utf-8"))["filas"]
    estilos = {"KKK": (AZUL, "s", "KKK (exponencial)"), "cuadratico": (NARANJA, "o", "cuadratico, una garganta"),
               "doble": (TINTA2, "D", "doble garganta"), "gauss": (AGUA, "^", "gaussiano")}
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for k, z in enumerate((0.0, 1.0)):
        a = ax[k]
        for fam, (color, mk, etiqueta) in estilos.items():
            sub = [f for f in ley if f["familia"] == fam and f["zeta"] == z and f["C"]]
            if not sub:
                continue
            if z == 0.0:
                prin = [min(f["C"], key=lambda c: abs(c - 1)) for f in sub]
                otros = [(f["a"], c) for f in sub for c in f["C"] if c != min(f["C"], key=lambda c: abs(c - 1))]
            else:
                prin = [max(f["C"]) for f in sub]
                otros = [(f["a"], c) for f in sub for c in f["C"] if c != max(f["C"])]
            a.plot([f["a"] for f in sub], prin, mk, color=color, ms=6, label=etiqueta, alpha=0.9)
            if otros:
                a.plot([x for x, _ in otros], [c for _, c in otros], mk, color=color, ms=6, mfc="none", alpha=0.7)
        a.axhline(1.0 if z == 0 else 2.0, color=TINTA, lw=0.8, ls="--")
        a.set_xscale("log")
        if z == 0.0:
            a.set_ylim(0, 1.8)
            a.set_title("zeta = 0: C = 1 (tasa = gravedad superficial)")
        else:
            a.set_yscale("log")
            a.set_title("zeta = 1 (capa de Antoniou et al.): C -> 2")
        a.set(xlabel="a r0 = Phi'(0+) r0  (aceleracion propia de la garganta)", ylabel="C = -w / a^2")
        a.legend(fontsize=7.5, loc="lower left" if z == 0 else "lower right")
    fig.suptitle("Modo rapido: w = -C a^2 (reloj de la garganta), es decir omega_inf = -i sqrt(C) kappa_t;"
                 " huecos: otros modos", color=TINTA, fontsize=10.5)
    fig.tight_layout()
    fig.savefig("figuras/fig6_ley.png", dpi=150)


def fig_derivacion():
    """Deduccion de §13 (ley_analitica.py, comparar_ley_analitica.py).  Izquierda: inercia 1/C del modo rapido frente a
    zeta >= 0; puntos numericos, lineas de la relacion de zona cercana.  Centro: error relativo de la relacion frente a
    a r0.  Derecha: modos complejos (fase sobreestable) en el plano (|Re omega_inf|, |Im omega_inf|) r0."""
    import cmath
    import numpy as np
    import sympy as sp
    from fondo_F import FondoF
    ley = json.load(open("resultados/ley_analitica.json", encoding="utf-8"))
    cmp_ = json.load(open("resultados/comparar_ley_analitica.json", encoding="utf-8"))
    u_, F0_, F2_, z_ = sp.symbols("u F0 F2 zeta", real=True)
    s_ = sp.symbols("s", positive=True)
    coefs = sp.lambdify((u_, F0_, F2_, z_), sp.Poly(sp.sympify(ley["polinomio"], locals=dict(
        s=s_, u=u_, F0=F0_, F2=F2_, zeta=z_)), s_).all_coeffs())
    tr = json.load(open("resultados/trayectorias_robustas.json", encoding="utf-8"))
    par = {o["nombre"]: o for o in json.load(open("resultados/estabilidad_robusta.json", encoding="utf-8"))}
    par.update({o["nombre"]: o for o in tr})

    def fondo(nom):
        o = par[nom] if nom in par else next(v for k, v in par.items() if k.startswith(nom + " "))
        return FondoF(o["tipo"], o["alpha"], o["phi0"], o["phi1"]), o

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.5))
    # --- izquierda: 1/C(zeta), zeta >= 0
    series = (("KKK referencia", AZUL, "s", "KKK (a r0 = 13.4)"),
              ("cuadratico 0.3087", NARANJA, "o", "cuadratico 0.309 (a r0 = 14.3)"),
              ("cuadratico 0.3465", AGUA, "^", "cuadratico 0.346 (a r0 = 64)"))
    for nom, color, mk, etiqueta in series:
        f, o = fondo(nom)
        a = f.y0()[3]
        F0, _, F2, _ = f.Fd(o["phi0"])
        u = o["phi1"]
        zc = 2 * (1 - F0) / (1 + u * u / 4)
        zs = np.linspace(0.0, zc, 400)
        inv = []
        for z in zs:
            rts = np.roots(np.array(coefs(u, F0, F2, z), dtype=float))
            Cs = [float((x * x).real) for x in rts if x.real > 1e-9 and abs(x.imag) < 1e-9 and (x * x).real > 0.3]
            inv.append(1 / max(Cs) if Cs else np.nan)
        ax[0].plot(zs, inv, "-", color=color, lw=2, zorder=2)
        sub = sorted((x for x in tr if x["nombre"] == nom and x["zeta"] >= 0), key=lambda x: x["zeta"])
        pts = [(x["zeta"], a * a / -min(x["w"])) for x in sub if x["w"] and -min(x["w"]) / a**2 > 0.3]
        ax[0].plot([p[0] for p in pts], [p[1] for p in pts], mk, color=color, ms=7, mec=FONDO_FIG, mew=1.2, zorder=3,
                   label=f"{etiqueta}; zeta_c = {zc:.2f}")
    zz = np.array([0.0, 2.0])
    ax[0].plot(zz, 1 - zz / 2, color=TINTA, lw=0.8, ls="--", label="1 - zeta/2 (limite de agujero negro)", zorder=1)
    ax[0].axhline(0, color=EJE, lw=0.8)
    ax[0].set(xlabel="zeta (gravedad inducida de la capa)", ylabel="1/C = -a^2/w (inercia efectiva)", xlim=(-0.03, 2.08),
              ylim=(-0.03, 1.08), title="Modo rapido: puntos numericos, lineas de zona cercana")
    ax[0].legend(fontsize=7.5, loc="upper right")
    # --- centro: error frente a a r0
    for z, mk, etq in ((0.0, "o", "zeta = 0"), (1.0, "D", "zeta = 1")):
        for valido, estilo in ((True, dict(color=TINTA2, mfc=TINTA2)), (False, dict(color=MUDO, mfc="none"))):
            xs, ys = [], []
            for fila in cmp_["filas"]:
                if fila["zeta"] != z or (fila["a"] >= 10 and abs(fila["u"]) / fila["a"] <= 0.1) != valido:
                    continue
                for Cn, Cp, e in fila["pares"]:
                    if Cn > 0.3:
                        xs.append(fila["a"])
                        ys.append(max(abs(e), 1e-5))
            ax[1].plot(xs, ys, mk, ms=5.5 if mk == "o" else 5, mew=1.0, ls="none", alpha=0.85,
                       label=f"{etq}, " + ("a r0 >= 10 y |phi1|/a <= 0.1" if valido else "fuera de ese regimen"), **estilo)
    xa = np.array([4, 900])
    ax[1].plot(xa, 2 * xa**-2.0, color=TINTA, lw=0.8, ls="--", label="2 (a r0)^-2")
    ax[1].set(xscale="log", yscale="log", xlabel="a r0", ylabel="|C_zona cercana / C_numerico - 1|", ylim=(5e-6, 1.5),
              title="Error de la relacion: O((a r0)^-2)")
    ax[1].legend(fontsize=7, loc="upper right")
    # --- derecha: modos complejos (sin los dos gaussianos casi reales, cuyo W predicho es real)
    conocidos = {"KKK 0.0035 1.0571"}
    casos = [(c["nombre"], c["zeta"], c["W_num"], c["W_pred"], False) for c in cmp_["complejos"]]
    casos += [(c["nombre"], c["zeta"], c["W_num"], c["W_pred"], c["nombre"] not in conocidos) for c in cmp_["prediccion"]]
    cache = {}
    for nom, z, Wn, Wp, nuevo in casos:
        if abs(Wp[1]) < 1e-9:
            continue
        if nom not in cache:
            f, o = fondo(nom)
            f.resolver()
            cache[nom] = f.y0()[3] * math.exp(-f.Phi_inf)
        kt = cache[nom]
        on = kt * cmath.sqrt(-complex(*Wn))
        op = kt * cmath.sqrt(-complex(*Wp))
        ax[2].plot([abs(op.real), abs(on.real)], [abs(op.imag), abs(on.imag)], color=EJE, lw=0.8, zorder=1)
        ax[2].plot(abs(op.real), abs(op.imag), "o", ms=10, mfc="none", mec=TINTA2, mew=1.2, zorder=2)
        ax[2].plot(abs(on.real), abs(on.imag), "o", ms=5, color=NARANJA if nuevo else AZUL, zorder=3)
        corto = (nom.replace("cuadratico ", "cuad. ").replace("KKK referencia", "KKK ref.")
                 .replace("KKK 0.0035 1.0571", "KKK f0=1.06").replace("KKK 0.0141 1.2186", "KKK f0=1.22")
                 .replace("KKK 0.0333 1.5501", "KKK f0=1.55").replace("doble exp a=0.04", "doble exp."))
        if corto in ("KKK f0=1.06", "KKK f0=1.55"):
            continue
        if corto == "KKK f0=1.22":
            corto = "KKK f0 = 1.06, 1.22, 1.55"
        dy = -10 if corto.startswith("KKK f0") or (corto == "cuad. 0.3087" and z == 3.0) else 4
        ax[2].annotate(f"{corto}, zeta={z:g}", (abs(on.real), abs(on.imag)), xytext=(5, dy), textcoords="offset points",
                       fontsize=6.5, color=TINTA2)
    ax[2].plot([], [], "o", ms=10, mfc="none", mec=TINTA2, mew=1.2, label="prediccion de zona cercana")
    ax[2].plot([], [], "o", ms=5, color=AZUL, label="numerico, ya conocido")
    ax[2].plot([], [], "o", ms=5, color=NARANJA, label="numerico, buscado tras la prediccion")
    ax[2].set(xscale="log", yscale="log", xlabel="|Re omega_inf| r0 (oscilacion)",
              ylabel="|Im omega_inf| r0 (tasa de crecimiento)", xlim=(0.04, 6), ylim=(5e-4, 1.5),
              title="Fase compleja (zeta > zeta_c): resonancias crecientes")
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
    ax[2].xaxis.set_major_locator(FixedLocator([0.05, 0.1, 0.2, 0.5, 1, 2, 5]))
    ax[2].xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax[2].xaxis.set_minor_formatter(NullFormatter())
    ax[2].legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig("figuras/fig7_deduccion.png", dpi=150)


def fig_tension():
    """§13.7: cascara con tension dependiente del escalar (tension_cercana.py, ventanas_tension.py, tension_completa.py).
    Izquierda y centro: abscisa espectral max Re s de la relacion de zona cercana en el plano (t, s2) para el fondo de KKK
    con zeta = 0 y zeta = 1, con el contorno s = 1 y los casos del calculo completo.  Derecha: min sobre s2 de max Re s
    frente a t para varios zeta (KKK); nunca baja de 1."""
    import numpy as np
    import sympy as sp
    from matplotlib.colors import LinearSegmentedColormap
    d = json.load(open("resultados/tension_cercana.json", encoding="utf-8"))
    u_, F0_, F2_, z_, l1_, s2_ = sp.symbols("u F0 F2 zeta l1 s2", real=True)
    s_ = sp.symbols("s", positive=True)
    cf = sp.lambdify((u_, F0_, F2_, z_, l1_, s2_), sp.Poly(sp.sympify(d["polinomio"], locals=dict(
        s=s_, u=u_, F0=F0_, F2=F2_, zeta=z_, l1=l1_, s2=s2_)), s_).all_coeffs())
    u, F0, F2 = -0.7, 0.03, 0.03

    def alfa(z, t, s2):
        c = np.array(cf(u, F0, F2, z, t / u, s2), dtype=float)
        while len(c) > 1 and abs(c[0]) < 1e-13 * np.abs(c).max():
            c = c[1:]
        r = [x for x in np.roots(c) if abs(x) > 1e-9]
        return max(x.real for x in r)

    azul = LinearSegmentedColormap.from_list("azul", ["#e9f1fb", "#9cc3ef", "#2a78d6", "#123e73"])
    ts = np.linspace(-6, 10, 161)
    ys = np.linspace(-3, 4, 141)                      # s2 = 2 sinh(y ln 10): s2 ~ +-10^|y| para |y| >= 1
    s2s = 2 * np.sinh(ys * np.log(10.0))
    comp = None
    try:
        comp = json.load(open("resultados/tension_completa.json", encoding="utf-8"))["casos"]
    except OSError:
        pass
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))
    for k, z in enumerate((0.0, 1.0)):
        A = np.array([[alfa(z, t, s2) for t in ts] for s2 in s2s])
        im = ax[k].pcolormesh(ts, ys, np.log10(A), cmap=azul, shading="auto", vmin=0, vmax=2)
        ax[k].contour(ts, ys, A, levels=[1.0 + 1e-3], colors=[TINTA], linewidths=0.8)
        if comp:
            for o in comp:
                if o["nombre"] == "KKK referencia" and o["zeta"] == z:
                    y = np.arcsinh(o["s2"] / 2) / np.log(10.0)
                    ax[k].plot(o["t"], y, "o", ms=6, mfc=NARANJA, mec=FONDO_FIG, mew=1.2)
        ax[k].set(xlabel="t = Lambda_s'(phi0) phi1 / (a^2 r0)", ylabel="s2 (escala s2 = 2 sinh(y ln 10))",
                  title=f"KKK, zeta = {z:g}: log10 max Re s (zona cercana)")
        ax[k].set_yticks([-3, -2, -1, 0, 1, 2, 3, 4])
        ax[k].set_yticklabels(["-1e3", "-1e2", "-10", "0", "10", "1e2", "1e3", "1e4"])
        cb = fig.colorbar(im, ax=ax[k])
        cb.set_label("log10 max Re s   (0: tasa = gravedad superficial)")
    tt = np.linspace(-6, 10, 321)
    for s2, color in ((0.0, AZUL), (10.0, AGUA), (100.0, NARANJA), (1000.0, TINTA2)):
        for z, estilo in ((1.0, "-"), (0.0, ":")):
            ax[2].plot(tt, [alfa(z, t, s2) for t in tt], estilo, color=color, lw=2 if z == 1.0 else 1.4,
                       label=f"s2 = {s2:g}" if z == 1.0 else None)
    ax[2].plot([], [], "-", color=MUDO, lw=2, label="(trazo) zeta = 1")
    ax[2].plot([], [], ":", color=MUDO, lw=1.4, label="(trazo) zeta = 0")
    ax[2].axhline(1.0, color=TINTA, lw=0.8)
    ax[2].axvline(2 * (1 + u * u / 4), color=EJE, lw=0.8)
    ax[2].set(yscale="log", xlabel="t = Lambda_s'(phi0) phi1 / (a^2 r0)", ylabel="max Re s (zona cercana)",
              title="KKK: max Re s >= 1 para todo t (C >= 1)", ylim=(0.8, 60))
    ax[2].legend(fontsize=7.5, loc="upper left", ncol=2)
    fig.tight_layout()
    fig.savefig("figuras/fig8_tension.png", dpi=150)


FONDO_FIG = "#fcfcfb"

if __name__ == "__main__":
    import sys
    if "deduccion" in sys.argv[1:]:
        fig_derivacion()
        raise SystemExit
    if "tension" in sys.argv[1:]:
        fig_tension()
        raise SystemExit
    fig_familia()
    fig_trayectorias()
    fig_limite()
    fig_derivacion()
    fig_tension()
