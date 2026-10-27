# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): criba de gargantas con galileon cubico y flujo con cola de energia positiva.

Geometria por la masa de Misner-Sharp, m(r) = M - (M - r0/2)(r0/r)^n (r0 = 1): crece hacia fuera si M > 1/2;
abocinamiento r''(0) = 1/2 - m_r(1) > 0 si n (M - 1/2) < 1/2.  r(rho) de r'' = m/r^2 - m_r/r, r(0) = 1, r'(0) = 0.
Flujo v = V0 [z1 e^{-z1^2/2} + lam z2 e^{-z2^2/2}];  A4 = -[1 + a e^{-rho^2/wq^2}] (Q = 2 A4/Nc);
lapso del corte Nc = 1 - MN/sqrt(rho^2 + wn^2) (cola 1/rho, MN > M para que lejos la materia efectiva cumpla el NEC:
(Q r')' = -2 (r'/Nc)' > 0, la condicion cinetica del regimen sin flujo, informe §6);  beta de su EDO (par, acotada, beta_inf = 1 fija V0).
Simbolo principal en forma de perfiles (simbolo_perfiles.py): P = P2 c^2 + P1 c + P0, D = det M_aux < 0;  simbolo
efectivo de d proporcional a s = P/D = s2 c^2 + s1 c + s0 (normalizacion 2A om^2 - 2B k^2):
    cinetica en el corte tau positiva <=> s2 > 0;  gradiente positivo <=> s0 < 0;  hiperbolico <=> s1^2 - 4 s2 s0 >= 0.
En la zona lejana (v, Q' -> 0) el signo cinetico es el de beta' (r K)' (K = 1 + r r'' - r'^2): positivo con una cola de
masa creciente (r K = 3m - r m_rho).
Uso: python criba_masa.py [n]
"""
import os
import sys
import json
import time
import pickle
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
L = 300.0
pp = pickle.load(open(os.path.join(AQUI, "resultados", "simbolo_perfiles.pkl"), "rb"))
x = sp.Symbol("rho", real=True)
NOMF = ("r", "v", "Q", "Nc", "beta")
SIM = {}
for n in NOMF:
    f = sp.Function(n)(x)
    SIM[f] = sp.Symbol(f"{n}0")
    for k in range(1, 5):
        SIM[sp.Derivative(f, (x, k))] = sp.Symbol(f"{n}{k}")


def a_simbolos(e):
    e = e.xreplace({k: v for k, v in SIM.items() if isinstance(k, sp.Derivative)})
    return e.xreplace({k: v for k, v in SIM.items() if not isinstance(k, sp.Derivative)})


Dd = a_simbolos(pp["D"])
s_expr = [sp.together(a_simbolos(pp["P"][j][0] / pp["P"][j][1]) / Dd) for j in (2, 1, 0)]
ARG = sorted({s for e in s_expr for s in e.free_symbols}, key=str)
fS = [sp.lambdify(ARG, sp.simplify(e), "numpy") for e in s_expr]
print("argumentos del simbolo:", [str(a) for a in ARG], flush=True)


def masa(M, n):
    def m(r):
        return M - (M - 0.5) * r ** (-n)

    def mr(r):
        return n * (M - 0.5) * r ** (-n - 1)

    def mrr(r):
        return -n * (n + 1) * (M - 0.5) * r ** (-n - 2)
    return m, mr, mrr


def evalua(prm):
    m, mr, mrr = masa(prm["M"], prm["n"])
    w1, w2, lam, a, wq, MN, wn = (prm[k] for k in ("w", "w2", "lam", "a", "wq", "MN", "wn"))

    def perf(xx):
        z1, z2 = xx / w1, xx / w2
        e1, e2 = np.exp(-z1 ** 2 / 2), np.exp(-z2 ** 2 / 2)
        v = z1 * e1 + lam * z2 * e2
        v1 = (1 - z1 ** 2) * e1 / w1 + lam * (1 - z2 ** 2) * e2 / w2
        v2 = (z1 ** 3 - 3 * z1) * e1 / w1 ** 2 + lam * (z2 ** 3 - 3 * z2) * e2 / w2 ** 2
        v3 = (-z1 ** 4 + 6 * z1 ** 2 - 3) * e1 / w1 ** 3 + lam * (-z2 ** 4 + 6 * z2 ** 2 - 3) * e2 / w2 ** 3
        g = np.exp(-xx ** 2 / wq ** 2)
        A = -2 * (1 + a * g)                          # 2 A4
        A1 = -2 * a * g * (-2 * xx / wq ** 2)
        A2 = -2 * a * g * (4 * xx ** 2 / wq ** 4 - 2 / wq ** 2)
        s_ = np.sqrt(xx ** 2 + wn ** 2)
        N = 1 - MN / s_
        N1 = MN * xx / s_ ** 3
        N2 = MN * (1 / s_ ** 3 - 3 * xx ** 2 / s_ ** 5)
        Q = A / N
        Q1 = A1 / N - A * N1 / N ** 2
        Q2 = A2 / N - 2 * A1 * N1 / N ** 2 - A * N2 / N ** 2 + 2 * A * N1 ** 2 / N ** 3
        return v, v1, v2, v3, Q, Q1, Q2, N, N1, N2

    def geo(r, r1):
        r2 = m(r) / r ** 2 - mr(r) / r
        r3 = r1 * (-2 * m(r) / r ** 3 + mr(r) / r ** 2 + mr(r) / r ** 2 - mrr(r) / r)
        return r2, r3

    def S(xx, r, r1, r2, V0=1.0):
        v, v1, v2 = perf(xx)[:3]
        v, v1, v2 = V0 * v, V0 * v1, V0 * v2
        Q, Q1 = perf(xx)[4:6]
        u1 = 2 * v * v1 / r ** 2 - 2 * v ** 2 * r1 / r ** 3
        u2 = (2 * v1 ** 2 + 2 * v * v2) / r ** 2 - 8 * v * v1 * r1 / r ** 3 - 2 * v ** 2 * r2 / r ** 3 \
            + 6 * v ** 2 * r1 ** 2 / r ** 4
        return (Q1 * r ** 4 * u1 + Q * 4 * r ** 3 * r1 * u1 + Q * r ** 4 * u2) / (4 * r ** 3)

    def rhs(xx, y):
        r, r1, bh, ch, bp, cp = y                     # beta = (bh, ch = beta'/r) homogenea y (bp, cp) particular
        r2, _ = geo(r, r1)
        K = 1 + r * r2 - r1 ** 2
        return [r1, r2, ch * r, -K * bh / r ** 3, cp * r, -K * bp / r ** 3 - S(xx, r, r1, r2)]

    sol = solve_ivp(rhs, (0, L), [1.0, 0.0, 1.0, 0.0, 0.0, 0.0], dense_output=True, rtol=1e-10, atol=1e-13,
                    max_step=0.05)
    if not sol.success or np.any(sol.y[0] <= 0):
        return dict(prm, B=False, error="geometria")
    # modo creciente (beta ~ r^2: beta'/r -> 2) eliminado con beta'/r en rho = L (el acotado tiene beta'/r ~ rho^-3);
    # con masa, beta ~ r^2 lleva correcciones rho ln rho y un ajuste en (rho^2, 1, 1/rho) no sirve
    YL = sol.sol(L)
    cc = -YL[5] / YL[3]
    binf = cc * YL[2] + YL[4]
    xx = np.linspace(1e-3, L, 30001)
    Y = sol.sol(xx)
    bet = cc * Y[2] + Y[4]
    out = dict(prm, B=bool(binf > 0 and bet.min() > 0))
    if not out["B"]:
        return out
    V0 = 1 / np.sqrt(binf)
    xs = np.concatenate([np.geomspace(0.003, 1, 150), np.linspace(1.01, 60, 600)])
    Y = sol.sol(xs)
    r, r1 = Y[0], Y[1]
    r2, r3 = geo(r, r1)
    b = (cc * Y[2] + Y[4]) / binf
    b1 = (cc * Y[3] + Y[5]) * r / binf
    v, v1, v2, v3, Q, Q1, Q2, N, N1, N2 = perf(xs)
    val = dict(r0=r, r1=r1, r2=r2, r3=r3, v0=V0 * v, v1=V0 * v1, v2=V0 * v2, v3=V0 * v3, Q0=Q, Q1=Q1, Q2=Q2,
               Nc0=N, Nc1=N1, Nc2=N2, beta0=b, beta1=b1)
    s2, s1, s0 = (np.asarray(f(*[val[str(a_)] for a_ in ARG]), float) for f in fS)
    esc = np.abs(s2) + np.abs(s1) + np.abs(s0) + 1e-300
    disc = (s1 ** 2 - 4 * s2 * s0) / esc ** 2
    rK = r * (1 + r * r2 - r1 ** 2)
    QR = Q1 * r1 + Q * r2                             # (Q r')': cinetica del regimen sin flujo (informe, §6)
    lapso_ok = bool(np.all(N > np.abs(V0 * v)) and np.all(N > 0))
    out.update(V0=float(V0), beta0=float(b[0]), vmax=float(np.max(np.abs(V0 * v))), rpp0=float(r2[0]),
               H=bool(np.all(disc >= -1e-8)), cin_pos=float(np.mean(s2 > 0)), grad_pos=float(np.mean(s0 < 0)),
               cin_pos_todo=bool(np.all(s2 > 0)), lejos_cin=bool(np.all(s2[xs > 20] > 0)),
               rK_crece_lejos=bool(np.all(np.diff(rK[xs > 20]) > 0)), lapso_ok=lapso_ok,
               QRp_lejos=bool(np.all(QR[xs > 5] > 0)))
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    rng = np.random.default_rng(23)
    RES = []
    t0 = time.time()
    for i in range(n):
        M = float(rng.uniform(0.52, 2.0))
        nn = float(rng.uniform(0.2, min(4.0, 0.5 / (M - 0.5) * 0.95)))
        prm = dict(M=M, n=nn, MN=float(M * rng.uniform(1.0, 3.0)), w=float(np.exp(rng.uniform(np.log(0.8), np.log(5)))),
                   a=float(np.exp(rng.uniform(np.log(5), np.log(300)))),
                   wq=float(np.exp(rng.uniform(np.log(0.3), np.log(3)))),
                   wn=0.0,
                   w2=float(np.exp(rng.uniform(np.log(0.08), np.log(0.6)))), lam=float(rng.uniform(0, 0.3)))
        prm["wn"] = float(prm["MN"] / rng.uniform(0.3, 0.9))     # Nc(0) = 1 - MN/wn en (0.1, 0.7)
        try:
            c = evalua(prm)
        except Exception as ex_:                       # noqa: BLE001
            c = dict(prm, B=False, error=str(ex_))
        RES.append(c)
        if c.get("B"):
            print(f"[{time.time() - t0:6.0f} s] " + str({k_: (float(f"{v_:.4g}") if isinstance(v_, float) else v_)
                                                       for k_, v_ in c.items()}), flush=True)
    nB = sum(1 for c in RES if c.get("B"))
    print(f"{n} fondos: beta > 0 en {nB}; hiperbolicos {sum(1 for c in RES if c.get('H'))}; cinetica positiva en todo "
          f"rho {sum(1 for c in RES if c.get('cin_pos_todo'))}; cinetica positiva lejos "
          f"{sum(1 for c in RES if c.get('lejos_cin'))}; (Q r')' > 0 lejos {sum(1 for c in RES if c.get('QRp_lejos'))}; "
          f"todo a la vez {sum(1 for c in RES if c.get('H') and c.get('QRp_lejos') and c.get('lapso_ok'))}", flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", "criba_masa.json"), "w"), indent=1)
