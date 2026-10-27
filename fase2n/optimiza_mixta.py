# -*- coding: utf-8 -*-
"""
Fase 2n: evolucion diferencial sobre la clase mixta de criba_mixta.py (nucleo de Ellis, cola de energia positiva).
Objetivo (0 = monopolo localmente sano en todos los puntos muestreados):
    beta con ceros:  2 - min(beta)/max|beta|  (> 2)
    beta > 0:        fraccion de puntos con s2 <= 0  +  fraccion con s0 >= 0  +  fraccion no hiperbolica
Uso: python optimiza_mixta.py [generaciones] [semilla]
"""
import os
import sys
import json
import time
import numpy as np
from scipy.optimize import differential_evolution
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import criba_mixta as cx                                 # noqa: E402
cm = cx.cm
NOM = ("lM", "n", "fMN", "lw", "la", "lwq", "lw2", "lam", "N0")
CAJA = [(np.log(0.05), np.log(2.0)), (0.2, 0.95), (1.0, 3.0), (np.log(0.5), np.log(5)), (np.log(0.3), np.log(300)),
        (np.log(0.05), np.log(3)), (np.log(0.05), np.log(0.8)), (0.0, 0.6), (0.1, 0.95)]


def prm_de(z):
    d = dict(zip(NOM, z))
    M = float(np.exp(d["lM"]))
    MN = M * d["fMN"]
    return dict(M=M, n=float(d["n"]), MN=MN, w=float(np.exp(d["lw"])), a=float(np.exp(d["la"])),
                wq=float(np.exp(d["lwq"])), wn=MN / (1 - d["N0"]), w2=float(np.exp(d["lw2"])), lam=float(d["lam"]))


def medida(prm):
    """repite cm.evalua devolviendo fracciones continuas."""
    m, mr, mrr = cm.masa(prm["M"], prm["n"])
    if not cx.geometria_ok(prm["M"], prm["n"]):
        return 5.0, {}
    w1, w2, lam, a, wq, MN, wn = (prm[k] for k in ("w", "w2", "lam", "a", "wq", "MN", "wn"))
    out = {}
    # reutiliza evalua para beta y el simbolo: se intercepta el resultado con un envoltorio de solve_ivp
    c = cm.evalua(prm)
    if not c.get("B"):
        bm = beta_min(prm)
        return 2.0 + min(1.0, max(0.0, -bm)), dict(c, bmin=bm)
    pen = (1 - c["cin_pos"]) + (1 - c["grad_pos"]) + (0.0 if c["H"] else 0.5)
    return float(pen), c


def beta_min(prm):
    """min(beta)/max|beta| con la misma EDO y la misma eliminacion del modo creciente que cm.evalua."""
    m, mr, mrr = cm.masa(prm["M"], prm["n"])
    w1, w2, lam, a, wq = (prm[k] for k in ("w", "w2", "lam", "a", "wq"))
    s_ = np.sqrt
    wn, MN = prm["wn"], prm["MN"]

    def S(xx, r, r1, r2):
        z1, z2 = xx / w1, xx / w2
        e1, e2 = np.exp(-z1 ** 2 / 2), np.exp(-z2 ** 2 / 2)
        v = z1 * e1 + lam * z2 * e2
        v1 = (1 - z1 ** 2) * e1 / w1 + lam * (1 - z2 ** 2) * e2 / w2
        v2 = (z1 ** 3 - 3 * z1) * e1 / w1 ** 2 + lam * (z2 ** 3 - 3 * z2) * e2 / w2 ** 2
        g = np.exp(-xx ** 2 / wq ** 2)
        A = -2 * (1 + a * g)
        A1 = -2 * a * g * (-2 * xx / wq ** 2)
        q = s_(xx ** 2 + wn ** 2)
        N = 1 - MN / q
        N1 = MN * xx / q ** 3
        Q, Q1 = A / N, A1 / N - A * N1 / N ** 2
        u1 = 2 * v * v1 / r ** 2 - 2 * v ** 2 * r1 / r ** 3
        u2 = (2 * v1 ** 2 + 2 * v * v2) / r ** 2 - 8 * v * v1 * r1 / r ** 3 - 2 * v ** 2 * r2 / r ** 3             + 6 * v ** 2 * r1 ** 2 / r ** 4
        return (Q1 * r ** 4 * u1 + Q * 4 * r ** 3 * r1 * u1 + Q * r ** 4 * u2) / (4 * r ** 3)

    def rhs(xx, y):
        r, r1, bh, ch, bp, cp = y
        r2 = m(r) / r ** 2 - mr(r) / r
        K = 1 + r * r2 - r1 ** 2
        return [r1, r2, ch * r, -K * bh / r ** 3, cp * r, -K * bp / r ** 3 - S(xx, r, r1, r2)]

    sol = solve_ivp(rhs, (0, cm.L), [1.0, 0.0, 1.0, 0.0, 0.0, 0.0], dense_output=True, rtol=1e-10, atol=1e-13,
                    max_step=0.05)
    YL = sol.sol(cm.L)
    cc = -YL[5] / YL[3]
    Y = sol.sol(np.linspace(1e-3, cm.L, 30001))
    bet = cc * Y[2] + Y[4]
    return float(bet.min() / np.abs(bet).max())


def objetivo(z):
    try:
        f, _ = medida(prm_de(z))
    except Exception:                                    # noqa: BLE001
        f = 6.0
    return f


if __name__ == "__main__":
    gen = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    sem = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    t0 = time.time()
    hist = []

    def cb(xk, convergence=None):
        f = objetivo(xk)
        hist.append(dict(f=f, prm=prm_de(xk), t=time.time() - t0))
        print(f"[{time.time() - t0:6.0f} s] mejor {f:.4f}  {prm_de(xk)}", flush=True)

    r = differential_evolution(objetivo, CAJA, maxiter=gen, popsize=12, seed=sem, polish=False, callback=cb, workers=6, updating="deferred",
                               tol=1e-6)
    f, c = medida(prm_de(r.x))
    print("final:", f, c, flush=True)
    json.dump(dict(hist=hist, final=c, f=f), open(os.path.join(AQUI, "resultados", f"optimiza_mixta_{sem}.json"),
                                                  "w"), indent=1, default=float)
