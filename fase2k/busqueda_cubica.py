# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): busqueda de fondos con galileon cubico y flujo que pasen las condiciones necesarias:
  (B) beta = Nc G4 > 0 en todo el eje (con beta_inf = 1), A4 < 0 (Q < 0);
  (G) garganta (garganta_cubica.py, R0 = 1): velocidades reales y energia cinetica positiva del monopolo,
        Gam > 0  y  Lam^ < 0 ,
      Gam ~ Q0 V1 (R2 V3 - V1 (R4 + R2 - R2^2)) - 2 beta0 (1 + R2)(R4 + R2),
      Lam^ ~ -[Q0^2 V1^4 + 2 Q0 R2 V1^2 beta0 - 4 beta0^2 (1 + R2)^2].
Familia: r = sqrt(rho^2 + 1) (R2 = 1, R4 = -3);  v = V0 [z1 e^{-z1^2/2} + lam z2 e^{-z2^2/2}], z_i = rho/w_i;
Q = -2 [1 + a e^{-rho^2/wq^2}].  V0 queda fijado por beta_inf = 1.  (El lapso Nc es libre y no entra en (B), (G).)
Uso: python busqueda_cubica.py [n]
"""
import os
import sys
import json
import numpy as np
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
L = 300.0


def perfiles(w1, w2, lam, a, wq):
    def f(x):
        R = np.sqrt(x ** 2 + 1)
        R1, R2 = x / R, 1 / R ** 3
        v = v1 = v2 = 0.0
        for c, w in ((1.0, w1), (lam, w2)):
            z = x / w
            e = np.exp(-z ** 2 / 2)
            v = v + c * z * e
            v1 = v1 + c * (1 - z ** 2) * e / w
            v2 = v2 + c * (z ** 3 - 3 * z) * e / w ** 2
        g = np.exp(-x ** 2 / wq ** 2)
        Q = -2 * (1 + a * g)
        Q1 = -2 * a * g * (-2 * x / wq ** 2)
        return R, R1, R2, v, v1, v2, Q, Q1
    return f


def resolver(w1, w2, lam, a, wq):
    f = perfiles(w1, w2, lam, a, wq)

    def S(x):
        R, R1, R2, v, v1, v2, Q, Q1 = f(x)
        u1 = 2 * v * v1 / R ** 2 - 2 * v ** 2 * R1 / R ** 3
        u2 = (2 * v1 ** 2 + 2 * v * v2) / R ** 2 - 8 * v * v1 * R1 / R ** 3 - 2 * v ** 2 * R2 / R ** 3 \
            + 6 * v ** 2 * R1 ** 2 / R ** 4
        return (Q1 * R ** 4 * u1 + Q * 4 * R ** 3 * R1 * u1 + Q * R ** 4 * u2) / (4 * R ** 3)

    def rhs(x, y, fuente):
        R, R1, R2 = f(x)[:3]
        return [y[1] * R, -(1 + R * R2 - R1 ** 2) * y[0] / R ** 3 - (S(x) if fuente else 0)]

    kw = dict(dense_output=True, rtol=1e-9, atol=1e-12, max_step=0.05)
    sh = solve_ivp(lambda x, y: rhs(x, y, False), (0, L), [1.0, 0.0], **kw)
    sp_ = solve_ivp(lambda x, y: rhs(x, y, True), (0, L), [0.0, 0.0], **kw)
    xa = np.linspace(L / 2, L, 200)
    Mx = np.vstack([xa ** 2, np.ones_like(xa), 1 / xa]).T
    ah = np.linalg.lstsq(Mx, sh.sol(xa)[0], rcond=None)[0]
    ap = np.linalg.lstsq(Mx, sp_.sol(xa)[0], rcond=None)[0]
    c = -ap[0] / ah[0]
    binf = c * ah[1] + ap[1]
    x = np.linspace(0, L, 30001)
    beta = c * sh.sol(x)[0] + sp_.sol(x)[0]
    out = dict(w1=w1, w2=w2, lam=lam, a=a, wq=wq)
    if not (binf > 0 and beta.min() > 0):
        out.update(B=False)
        return out
    V0 = 1 / np.sqrt(binf)
    b0 = beta[0] / binf
    Q0 = -2 * (1 + a)
    V1 = V0 * (1 / w1 + lam / w2)
    V3 = V0 * (-3 / w1 ** 3 - 3 * lam / w2 ** 3)
    R2, R4 = 1.0, -3.0
    Gam = Q0 * V1 * (R2 * V3 - V1 * (R4 + R2 - R2 ** 2)) - 2 * b0 * (1 + R2) * (R4 + R2)
    Lam = -(Q0 ** 2 * V1 ** 4 + 2 * Q0 * R2 * V1 ** 2 * b0 - 4 * b0 ** 2 * (1 + R2) ** 2)
    out.update(B=True, V0=float(V0), beta0=float(b0), V1=float(V1), V3=float(V3), Gam=float(Gam), Lam=float(Lam),
               G=bool(Gam > 0 and Lam < 0), vmax=float(np.max(np.abs(V0 * f(x)[3]))))
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    rng = np.random.default_rng(7)
    RES = []
    for i in range(n):
        w1 = float(np.exp(rng.uniform(np.log(0.15), np.log(1.2))))
        w2 = float(np.exp(rng.uniform(np.log(1.0), np.log(6.0))))
        lam = float(rng.uniform(-1.0, 3.0))
        a = float(np.exp(rng.uniform(np.log(0.5), np.log(200.0))))
        wq = float(np.exp(rng.uniform(np.log(0.2), np.log(4.0))))
        try:
            c = resolver(w1, w2, lam, a, wq)
        except Exception as ex:                       # noqa: BLE001
            c = dict(w1=w1, w2=w2, lam=lam, a=a, wq=wq, B=False, error=str(ex))
        RES.append(c)
        if c.get("B"):
            print({k: (float(f"{x:.4g}") if isinstance(x, float) else x) for k, x in c.items()}, flush=True)
    nB = sum(1 for c in RES if c.get("B"))
    nG = sum(1 for c in RES if c.get("G"))
    print(f"{n} formas: beta > 0 en {nB}; ademas garganta sana (Gam > 0, Lam^ < 0) en {nG}", flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", "busqueda_cubica.json"), "w"), indent=1)
