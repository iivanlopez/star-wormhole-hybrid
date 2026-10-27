# -*- coding: utf-8 -*-
"""
Fase 2k: existencia de gargantas con phi = qt t (psi' = 0) -- reconstruccion.

Con v = 0 el fondo solo ve A2 = G2 y B4 = G4 (A4, es decir F4 en la GLPV, queda libre).  Con p = N A2, beta = N B4 y
la longitud propia rho, la lagrangiana estatica Lest = p r^2 + 2 beta (1 + r'^2) + 4 beta' r r' no depende de rho y
la ecuacion de h (fijada por el gauge) es la ligadura H = 0 de su primera integral:
    H = 2 beta (r'^2 - 1) + 4 beta' r r' - p r^2 = 0,
    E_r: p r = 2 (beta r')' + 2 beta'' r      ->     r^2 beta'' - r r' beta' + (1 + r r'' - r'^2) beta = 0,
y E_N se sigue de dH/drho = -(E_r r' + E_N N').  N no aparece: dada r(rho), beta es un modo cero de
    (r u')' + (1 + 2 r r'' - 2 r'^2) u / r = 0,      beta = r u,
que debe decaer (u ~ 1/rho, beta -> cte) en los dos extremos: hace falta afinar un parametro de forma de r:
    r = sqrt(rho^2 + r0^2) (1 + eps e^{-rho^2/w^2}).
El lapso N(rho) es libre; se toma N = exp(-M/sqrt(rho^2 + a^2)) (monotona en rho > 0: X = qt^2/2N^2 monotona).
G4 = beta/N, G2 = p/N.  Se comprueba:
  - beta > 0 (modo cero sin nodos) y G4 -> 1, G2 -> 0 lejos;
  - residuo de E_r y de E_N (forma de tercer orden) sobre la solucion;
  - la prueba dX = 0: d = r f con (r^3 f')' + r f = 0 (es D[d] = 0) y energia -2 int beta r^2 f'^2 < 0;
  - signos de Q = 2 A4/N y (Q r')' en la garganta con F4 = 0 (A4 = -G4 + 2X G4_X).
Uso: python reconstruccion_qt.py
"""
import os
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

AQUI = os.path.dirname(os.path.abspath(__file__))
L = 400.0


def geom(eps, r0=1.0, wf=1.5):
    def rr(x):
        g = np.exp(-x ** 2 / wf ** 2)
        s = np.sqrt(x ** 2 + r0 ** 2)
        s1, s2 = x / s, r0 ** 2 / s ** 3
        g1 = -2 * x / wf ** 2 * g
        g2 = (4 * x ** 2 / wf ** 4 - 2 / wf ** 2) * g
        f, f1, f2 = 1 + eps * g, eps * g1, eps * g2
        return s * f, s1 * f + s * f1, s2 * f + 2 * s1 * f1 + s * f2
    return rr


def beta_sol(rr):
    def f(x, y):
        b, b1 = y
        R, R1, R2 = rr(x)
        return [b1, (R * R1 * b1 - (1 + R * R2 - R1 ** 2) * b) / R ** 2]
    return solve_ivp(f, (0, L), [1.0, 0.0], dense_output=True, rtol=1e-12, atol=1e-14)


def coef_rho2(eps, **kw):
    s = beta_sol(geom(eps, **kw))
    x = np.linspace(L / 2, L, 300)
    Mx = np.vstack([x ** 2, np.ones_like(x), 1 / x]).T
    return np.linalg.lstsq(Mx, s.sol(x)[0], rcond=None)[0][0]


def caso(r0=1.0, wf=1.5, Mm=0.5, aa=1.0, qt=1.0, rango=(-0.9, 3.0)):
    es = np.linspace(*rango, 80)
    cs = [coef_rho2(e, r0=r0, wf=wf) for e in es]
    raices = [brentq(lambda e: coef_rho2(e, r0=r0, wf=wf), es[i], es[i + 1], xtol=1e-14)
              for i in range(len(es) - 1) if cs[i] * cs[i + 1] < 0]
    out = []
    for eps in raices:
        rr = geom(eps, r0, wf)
        s = beta_sol(rr)
        x = np.linspace(0, L, 80001)
        b, b1 = s.sol(x)
        b, b1 = b / b[-1], b1 / b[-1]
        R, R1, R2 = rr(x)
        b2 = (R * R1 * b1 - (1 + R * R2 - R1 ** 2) * b) / R ** 2
        p = (2 * b * (R1 ** 2 - 1) + 4 * b1 * R * R1) / R ** 2
        Er = p * R - 2 * (b1 * R1 + b * R2) - 2 * b2 * R
        # E_N en la forma p' r^2 + beta' (2 - 4 r r'' - 2 r'^2) = 0 (derivada numerica)
        EN = np.gradient(p, x) * R ** 2 + b1 * (2 - 4 * R * R2 - 2 * R1 ** 2)
        N = np.exp(-Mm / np.sqrt(x ** 2 + aa ** 2))
        N = N / N[-1]
        G4, G2 = b / N, p / N
        X = qt ** 2 / (2 * N ** 2)
        G4X = np.gradient(G4, x)[1:] / np.gradient(X, x)[1:]
        A4_0 = -G4[1] + 2 * X[1] * G4X[0]
        Q0 = 2 * A4_0 / N[0]
        # prueba dX = 0: d = r f, (r^3 f')' + r f = 0, f(0) = 1, f'(0) = 0
        sd = solve_ivp(lambda xx, y: [y[1] / rr(xx)[0] ** 3, -rr(xx)[0] * y[0]], (0, L), [1.0, 0.0],
                       dense_output=True, rtol=1e-12, atol=1e-14)
        f_, r3f1 = sd.sol(x)
        f1 = r3f1 / R ** 3
        E = 2 * np.trapezoid(-2 * b * R ** 2 * f1 ** 2, x)
        nodos = int(np.sum(np.diff(np.sign(b)) != 0))
        out.append(dict(eps=float(eps), r_garganta=float(R[0]), rpp_garganta=float(R2[0]), beta_min=float(b.min()),
                        nodos_beta=nodos, G4_min=float(G4.min()), G4_0=float(G4[0]), G2_0=float(G2[0]),
                        G2_inf=float(G2[-1]), res_Er=float(np.max(np.abs(Er[:-10]))),
                        res_EN=float(np.max(np.abs(EN[5:-5]))), A4_0=float(A4_0), Q0=float(Q0),
                        QR2_0=float(Q0 * R2[0]), E_prueba=float(E), f_L=float(f_[-1])))
    return out


if __name__ == "__main__":
    RES = []
    for (r0, wf) in ((1.0, 1.5), (1.0, 0.8), (1.0, 3.0), (2.0, 1.5)):
        for c in caso(r0=r0, wf=wf):
            c.update(r0=r0, w=wf)
            RES.append(c)
            print({k: (float(f"{v:.5g}") if isinstance(v, float) else v) for k, v in c.items()}, flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", "reconstruccion_qt.json"), "w"), indent=1)
