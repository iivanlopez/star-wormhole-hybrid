# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): ?hay beta = Nc G4 > 0 con cola de energia positiva?  Evolucion diferencial sobre la familia de
criba_masa.py (geometria m(r) = M - (M - 1/2) r^-n, flujo de dos componentes, A4 = -[1 + a e^{-rho^2/wq^2}],
Nc = 1 - MN/sqrt(rho^2 + wn^2)) maximizando  min(beta)/max|beta|  (con el signo que da beta_inf > 0).
beta > 0 en todo el eje <=> objetivo > 0.
Uso: python optimiza_masa.py [generaciones]
"""
import os
import sys
import json
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution

AQUI = os.path.dirname(os.path.abspath(__file__))
L = 200.0
# parametros: M, frac_n (n = frac_n * 0.95 * 0.5/(M - 1/2)), MN/M, Nc(0), w, w2, lam, log a, wq
LIM = [(0.52, 2.0), (0.05, 1.0), (1.0, 3.0), (0.1, 0.9), (0.5, 5.0), (0.08, 0.8), (0.0, 0.5), (np.log(2), np.log(500)),
       (0.2, 3.0)]


def beta(p):
    M, fn, mrat, N0, w1, w2, lam, la, wq = p
    n = fn * 0.95 * 0.5 / (M - 0.5)
    MN = M * mrat
    wn = MN / (1 - N0)
    a = np.exp(la)

    def m(r):
        return M - (M - 0.5) * r ** (-n)

    def mr(r):
        return n * (M - 0.5) * r ** (-n - 1)

    def perf(x):
        z1, z2 = x / w1, x / w2
        e1, e2 = np.exp(-z1 ** 2 / 2), np.exp(-z2 ** 2 / 2)
        v = z1 * e1 + lam * z2 * e2
        v1 = (1 - z1 ** 2) * e1 / w1 + lam * (1 - z2 ** 2) * e2 / w2
        v2 = (z1 ** 3 - 3 * z1) * e1 / w1 ** 2 + lam * (z2 ** 3 - 3 * z2) * e2 / w2 ** 2
        g = np.exp(-x ** 2 / wq ** 2)
        A = -2 * (1 + a * g)
        A1 = -2 * a * g * (-2 * x / wq ** 2)
        s_ = np.sqrt(x ** 2 + wn ** 2)
        N = 1 - MN / s_
        N1 = MN * x / s_ ** 3
        return v, v1, v2, A / N, A1 / N - A * N1 / N ** 2

    def rhs(x, y):
        r, r1, bh, ch, bp, cp = y
        r2 = m(r) / r ** 2 - mr(r) / r
        K = 1 + r * r2 - r1 ** 2
        v, v1, v2, Q, Q1 = perf(x)
        u1 = 2 * v * v1 / r ** 2 - 2 * v ** 2 * r1 / r ** 3
        u2 = (2 * v1 ** 2 + 2 * v * v2) / r ** 2 - 8 * v * v1 * r1 / r ** 3 - 2 * v ** 2 * r2 / r ** 3 \
            + 6 * v ** 2 * r1 ** 2 / r ** 4
        S = (Q1 * r ** 4 * u1 + Q * 4 * r ** 3 * r1 * u1 + Q * r ** 4 * u2) / (4 * r ** 3)
        return [r1, r2, ch * r, -K * bh / r ** 3, cp * r, -K * bp / r ** 3 - S]

    s = solve_ivp(rhs, (0, L), [1, 0, 1, 0, 0, 0], dense_output=True, rtol=1e-9, atol=1e-12, max_step=0.05)
    YL = s.sol(L)
    cc = -YL[5] / YL[3]
    xs = np.linspace(0, L, 4001)
    Y = s.sol(xs)
    b = cc * Y[2] + Y[4]
    binf = cc * YL[2] + YL[4]
    return xs, b * np.sign(binf), dict(M=M, n=n, MN=MN, wn=wn, w=w1, w2=w2, lam=lam, a=a, wq=wq)


def objetivo(p):
    try:
        _, b, _ = beta(p)
        return -float(b.min() / np.max(np.abs(b)))
    except Exception:                                  # noqa: BLE001
        return 2.0


if __name__ == "__main__":
    gens = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    t0 = time.time()
    hist = []

    def cb(xk, convergence=None):
        hist.append(float(objetivo(xk)))
        print(f"[{time.time() - t0:6.0f} s] mejor objetivo (-min/max) = {hist[-1]:+.4f}", flush=True)

    res = differential_evolution(objetivo, LIM, maxiter=gens, popsize=8, seed=3, tol=1e-6, polish=False, callback=cb)
    xs, b, prm = beta(res.x)
    print("optimo:", prm, " min beta/max|beta| =", float(b.min() / np.max(np.abs(b))), flush=True)
    json.dump(dict(objetivo=float(-res.fun), parametros=prm, historia=hist),
              open(os.path.join(AQUI, "resultados", "optimiza_masa.json"), "w"), indent=1)
