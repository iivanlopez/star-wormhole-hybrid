# -*- coding: utf-8 -*-
"""
Fase 2l, parte III: cuanta rotacion haria falta para que la prueba d = N deje de detectar el taquion.

Con rotacion lenta (monopolo_rotacion.py), a O(Omega^2) y exactamente para la familia c = 0:
    E_J[N] = E_0 + pi^2 C_1,     E_0 = -2 int N^3 Gcal sig^2 < 0,     C_1 = 9 int N^3/(Gcal r^6) > 0,
con pi = J w' constante (pi = -4 J_Kerr en un exterior de Kerr lento).  La rotacion se mide con la velocidad del
observador de area cero (ZAMO) de la garganta respecto del sistema estatico de un infinito, prolongado rigidamente a la
garganta:
    v0 = r0 |w(0) - w(+inf)| / N0 = r0 |pi| I / N0,      I = int_0^inf dl / J = (3/2) int_0^inf N/(r^4 Gcal) dl
(en una garganta simetrica w(l) - w(0) es impar: los dos infinitos giran a +-pi I respecto de la garganta y la
velocidad es la misma vista desde cualquiera de ellos; v0 < 1 es la ausencia de ergorregion en la garganta).  Entonces
    E_J[N]/|E_0| = -1 + kappa v0^2,     kappa = C_1 N0^2/(r0^2 I^2 |E_0|),     v_c = kappa^{-1/2}.
Fondos: Ellis-Bronnikov (RG + escalar fantasma: N = 1, Gcal = 1, r^2 = l^2 + r0^2; analitico: kappa = 24/pi^2) y las 34
gargantas reconstruidas de la Fase 2f (las de fase2g/verificar_N.py).  Tambien chi_c = J_c/M^2 (M de Komar).
Uso: python estimacion_rotacion.py     -> resultados/estimacion_rotacion.json
"""
import os
import sys
import json
import math
import time
import numpy as np
from scipy.integrate import simpson, quad
from concurrent.futures import ProcessPoolExecutor

AQUI = os.path.dirname(os.path.abspath(__file__))
F2G = os.path.join(AQUI, "..", "fase2g")
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
sys.path.insert(0, F2G)


def malla(L=4000.0):
    return np.concatenate([np.linspace(0.0, 30, 30001), np.geomspace(30, L, 4001)[1:]])


def integrales(ls, N, r, r1, Phi1, Gc):
    sig = r1 / r - Phi1
    E0 = -2 * 2 * simpson(N ** 3 * Gc * sig ** 2, x=ls)            # linea completa = 2 x semirrecta
    cola_E0 = -2 * 2 * (N ** 3 * Gc * sig ** 2)[-1] * ls[-1]        # sig^2 ~ 1/l^2 -> cola ~ f(L) L
    C1 = 9 * 2 * simpson(N ** 3 / (Gc * r ** 6), x=ls)
    I = 1.5 * simpson(N / (r ** 4 * Gc), x=ls)
    N0, r0 = N[0], r[0]
    kappa = C1 * N0 ** 2 / (r0 ** 2 * I ** 2 * abs(E0 + cola_E0))
    return dict(E0=float(E0 + cola_E0), cola_E0=float(cola_E0), C1=float(C1), I=float(I), N0=float(N0), r0=float(r0),
                kappa=float(kappa), v_c=float(kappa ** -0.5))


def ellis():
    ls = malla()
    r = np.sqrt(ls ** 2 + 1)
    one = np.ones_like(ls)
    num = integrales(ls, one, r, ls / r, 0 * ls, one)
    num["kappa_analitico"] = 24 / math.pi ** 2
    num["v_c_analitico"] = math.pi / math.sqrt(24)
    # comprobacion independiente con quad
    E0 = -2 * quad(lambda x: x * x / (1 + x * x) ** 2, -np.inf, np.inf)[0]
    C1 = 9 * quad(lambda x: (1 + x * x) ** -3, -np.inf, np.inf)[0]
    I = 1.5 * quad(lambda x: (1 + x * x) ** -2, 0, np.inf)[0]
    num["kappa_quad"] = C1 / (I ** 2 * abs(E0))
    return num


def evalua(item):
    import busqueda as B
    nombre, x, w0 = item
    g = B.construir(x)
    ls = malla()
    acc = {k: [] for k in ("N", "r", "r1", "Phi1", "Gc")}
    for i in range(0, len(ls), 3000):
        V = g._fondo_basico(ls[i:i + 3000])
        acc["N"].append(V["N"].deriv(0))
        acc["r"].append(V["r0"].deriv(0))
        acc["r1"].append(V["r1"].deriv(0))
        acc["Phi1"].append(V["Phi"].deriv(1))
        acc["Gc"].append(V["p"]["Gc"].deriv(0))
    f = {k: np.concatenate(v) for k, v in acc.items()}
    out = integrales(ls, f["N"], f["r"], f["r1"], f["Phi1"], f["Gc"])
    Mk = float(g.masa_komar)
    pic = out["v_c"] * out["N0"] / (out["r0"] * out["I"])         # |pi| critico (v0 = v_c)
    out.update(nombre=nombre, M_komar=Mk, pi_c=float(pic), J_c=float(pic / 4), chi_c=float(pic / 4 / Mk ** 2),
               minGc=float(f["Gc"].min()), w0_fase2f=w0)
    return out


if __name__ == "__main__":
    t0 = time.time()
    import verificar_N as VN
    el = ellis()
    print(f"Ellis-Bronnikov: E0 = {el['E0']:+.6f} (analitico -pi = {-math.pi:+.6f}), kappa = {el['kappa']:.6f} "
          f"(analitico 24/pi^2 = {el['kappa_analitico']:.6f}, quad {el['kappa_quad']:.6f}), v_c = {el['v_c']:.4f}")
    cs = VN.casos()
    vistos, unicos = set(), []
    for c in cs:
        k = json.dumps(c[1], sort_keys=True)
        if k not in vistos:
            vistos.add(k)
            unicos.append(c)
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(evalua, unicos))
    print(f"{len(res)} gargantas distintas de la Fase 2f (de {len(cs)})")
    print(f"{'caso':28s} {'N0':>6s} {'r0':>6s} {'E0':>9s} {'C1':>9s} {'I':>8s} {'kappa':>8s} {'v_c':>7s} {'M':>6s} "
          f"{'chi_c':>8s}")
    for r in res:
        print(f"{r['nombre']:28s} {r['N0']:6.3f} {r['r0']:6.3f} {r['E0']:+9.5f} {r['C1']:9.5f} {r['I']:8.4f} "
              f"{r['kappa']:8.4f} {r['v_c']:7.4f} {r['M_komar']:6.3f} {r['chi_c']:8.4f}")
    vc = [r["v_c"] for r in res]
    print(f"v_c: min {min(vc):.4f}, mediana {float(np.median(vc)):.4f}, max {max(vc):.4f};  "
          f"con v_c < 1: {sum(v < 1 for v in vc)} de {len(vc)}")
    # fondos distintos (cm2 solo cambia G2_XX, y la referencia coincide con 'debil M=2')
    u = {(round(r["N0"], 6), round(r["E0"], 6), round(r["C1"], 6)): r for r in res}
    ku = np.array([r["kappa"] for r in u.values()])
    vu = np.array([r["v_c"] for r in u.values()])
    print(f"{len(u)} fondos distintos: kappa {ku.min():.3f}-{ku.max():.3f} (mediana {np.median(ku):.3f}), v_c "
          f"{vu.min():.3f}-{vu.max():.3f} (mediana {np.median(vu):.3f}), v_c < 1 en {int((vu < 1).sum())};  "
          f"reduccion maxima de |E0| a v0 = 0.1, 0.3: {ku.max() * 0.01:.3f}, {ku.max() * 0.09:.3f}")
    print(f"{time.time() - t0:.0f} s")
    json.dump(dict(ellis=el, gargantas=res), open(os.path.join(AQUI, "resultados", "estimacion_rotacion.json"), "w"),
              indent=1, default=float)
