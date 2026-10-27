# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): fondos con galileon cubico y flujo (G3 != 0, psi' != 0) con ondas gravitatorias sanas en el corte
(A4 < 0, G4 > 0).

Ecuacion de beta = Nc G4 (nogo_beta.py), con u = v^2/r^2:
    (beta'/r)' + K beta/r^3 = -(Q r^4 u')'/(4 r^3),    K = 1 + r r'' - r'^2,  Q = 2 A4/Nc .
Identidad integral (beta > 0, dividir por beta, integrar y completar el cuadrado):
    int (1/r)(beta'/beta + Q r^2 u'/(8 beta))^2 + int (1 + r'^2)/r^3 = -(3/4) int Q u' r'/beta
                                                                  + (1/64) int Q^2 r^3 u'^2/beta^2   (identidad_cubica.py),
cuyo ultimo termino es positivo sea cual sea el signo de Q: un flujo fuerte puede pagar la garganta.
Reconstruccion: r = sqrt(rho^2 + 1), v = V0 x e^{-x^2/2} (x = rho/w), Q = -2 [1 + a e^{-rho^2/wq^2}] (A4 < 0 con
Nc > 0 si a > -1), beta par acotada (unica: la ecuacion es inhomogenea), V0 fijado por beta_inf = 1 (G4 -> 1, A4 -> -1).
Se buscan formas con beta > 0 en todo el eje y se da max |v| (el lapso del corte, libre, debe superarlo).
Uso: python reconstruccion_cubica.py
"""
import os
import json
import numpy as np
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
L = 400.0


def perfiles(w, a, wq):
    def f(x):
        R = np.sqrt(x ** 2 + 1)
        R1, R2 = x / R, 1 / R ** 3
        z = x / w
        e = np.exp(-z ** 2 / 2)
        v = z * e
        v1 = (1 - z ** 2) * e / w
        v2 = (z ** 3 - 3 * z) * e / w ** 2
        g = np.exp(-x ** 2 / wq ** 2)
        Q = -2 * (1 + a * g)
        Q1 = -2 * a * g * (-2 * x / wq ** 2)
        return R, R1, R2, v, v1, v2, Q, Q1
    return f


def fuente(f, x):
    """S = (Q r^4 u')'/(4 r^3), u = v^2/r^2 (V0 = 1), con u'' exacta."""
    R, R1, R2, v, v1, v2, Q, Q1 = f(x)
    u1 = 2 * v * v1 / R ** 2 - 2 * v ** 2 * R1 / R ** 3
    u2 = (2 * v1 ** 2 + 2 * v * v2) / R ** 2 - 8 * v * v1 * R1 / R ** 3 - 2 * v ** 2 * R2 / R ** 3 \
        + 6 * v ** 2 * R1 ** 2 / R ** 4
    return (Q1 * R ** 4 * u1 + Q * 4 * R ** 3 * R1 * u1 + Q * R ** 4 * u2) / (4 * R ** 3)


def resolver(w, a, wq):
    f = perfiles(w, a, wq)

    def hom(x, y):
        R, R1, R2 = f(x)[:3]
        return [y[1] * R, -(1 + R * R2 - R1 ** 2) * y[0] / R ** 3]

    def inh(x, y):
        R, R1, R2 = f(x)[:3]
        return [y[1] * R, -(1 + R * R2 - R1 ** 2) * y[0] / R ** 3 - fuente(f, x)]

    kw = dict(dense_output=True, rtol=1e-10, atol=1e-13, max_step=0.05)
    sh = solve_ivp(hom, (0, L), [1.0, 0.0], **kw)
    sp_ = solve_ivp(inh, (0, L), [0.0, 0.0], **kw)
    xa = np.linspace(L / 2, L, 300)
    Mx = np.vstack([xa ** 2, np.ones_like(xa), 1 / xa]).T
    ah = np.linalg.lstsq(Mx, sh.sol(xa)[0], rcond=None)[0]
    ap = np.linalg.lstsq(Mx, sp_.sol(xa)[0], rcond=None)[0]
    c = -ap[0] / ah[0]
    x = np.linspace(0, L, 40001)
    beta = c * sh.sol(x)[0] + sp_.sol(x)[0]
    binf = c * ah[1] + ap[1]
    out = dict(w=w, a=a, wq=wq, beta_inf_V1=float(binf))
    if binf > 0 and beta.min() > 0:
        V0 = 1 / np.sqrt(binf)
        out.update(positiva=True, V0=float(V0), vmax=float(V0 * np.exp(-0.5)), beta0=float(beta[0] / binf),
                   beta_min=float(beta.min() / binf))
    else:
        out.update(positiva=False, beta_min_rel=float(beta.min() / abs(binf)) if binf != 0 else None)
    return out


if __name__ == "__main__":
    RES = []
    for w in (0.3, 0.5, 1.0, 2.0):
        for a in (-0.8, -0.5, 0.0, 2.0, 5.0, 20.0):
            for wq in (0.3, 1.0, 3.0):
                c = resolver(w, a, wq)
                RES.append(c)
                if c["positiva"]:
                    print({k: (float(f"{x:.4g}") if isinstance(x, float) else x) for k, x in c.items()}, flush=True)
    print(f"{sum(c['positiva'] for c in RES)} de {len(RES)} formas con beta > 0 (G4 > 0) y Q < 0 (A4 < 0)")
    json.dump(RES, open(os.path.join(AQUI, "resultados", "reconstruccion_cubica.json"), "w"), indent=1)
