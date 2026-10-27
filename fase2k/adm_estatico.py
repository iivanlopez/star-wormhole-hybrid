# -*- coding: utf-8 -*-
"""
Fase 2k: el sector estatico en las variables ADM del corte phi = cte, en forma cerrada.

Con phi = qt tau (tau = t + psi/qt) y todo independiente de tau, en la longitud propia rho del corte:
    ds^2 = -(Nc^2 - v^2) dtau^2 + 2 v dtau drho + drho^2 + r^2 dOmega^2,      X = qt^2/(2 Nc^2),  y = sqrt(2X) = qt/Nc,
(Nc el lapso del corte, v = N^rho el desplazamiento).  La lagrangiana de la GLPV cuartica en el gauge unitario es
    L = A2 + A3 K + A4 (K^2 - K_ij K^ij) + B4 R3,
    A2 = G2,  A3 = y G3 - int_0^y G3 dy,  A4 = -G4 + 2 X G4_X + 4 X^2 F4,  B4 = G4
(derivacion en el informe, §3), y en el sector estatico l = 0 (K_rho^rho = -v'/Nc, K_th^th = -r' v/(r Nc)):
    Lest = Nc r^2 A2 - A3 (r^2 v)' + (2 A4/Nc) r' (v^2 r)' + Nc B4 (2 - 4 r r'' - 2 r'^2).
Comprobaciones:
  (1) Lest es la lagrangiana covariante reducida (fase2j/q_simbolico.py) evaluada en esa metrica salvo una derivada
      total: el operador de Euler-Lagrange de la diferencia se anula identicamente para Nc, v, r (teoria polinomica
      con coeficientes simbolicos);
  (2) GR (A4 = -1, B4 = 1): Schwarzschild en coordenadas de Painleve-Gullstrand (Nc = 1, v = sqrt(2M/r), r = rho) es
      solucion.
Uso: python adm_estatico.py
"""
import os
import sys
import json
import time
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "fase2j"))
import q_simbolico as S        # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


rho = S.l
tau = S.t
qt = sp.Symbol("qt", positive=True)
Nc, v, r = (sp.Function(n)(rho) for n in ("Nc", "v", "r"))
RES = {}


def euler_lagrange(L, campos):
    out = {}
    for f_ in campos:
        e = sp.euler_equations(L, [f_], [rho])[0].lhs
        out[str(f_.func)] = sp.simplify(sp.expand(e.xreplace({sp.sqrt(Nc ** 2): Nc})))
    return out


def lest(A2, A3, A4, B4):
    return (Nc * r ** 2 * A2 - A3 * sp.diff(r ** 2 * v, rho) + 2 * A4 / Nc * sp.diff(r, rho) * sp.diff(v ** 2 * r, rho)
            + Nc * B4 * (2 - 4 * r * sp.diff(r, rho, 2) - 2 * sp.diff(r, rho) ** 2))


if __name__ == "__main__":
    # ------------------------------------------------------------------ (1) frente a la lagrangiana covariante
    Xs = sp.Symbol("X", positive=True)
    ys = sp.Symbol("y", positive=True)
    k = sp.symbols("k1:4", real=True)
    c = sp.symbols("c1:3", real=True)
    b = sp.symbols("b1:3", real=True)
    f = sp.symbols("f0:2", real=True)
    G2 = k[0] * Xs + k[1] * Xs ** 2 + k[2] * Xs ** 3
    G3 = c[0] * Xs + c[1] * Xs ** 2
    G4 = 1 + b[0] * Xs + b[1] * Xs ** 2
    F4 = f[0] + f[1] * Xs

    def teoria(phi, X):
        sx = {Xs: X}
        return (G2.subs(sx), G3.subs(sx), G4.subs(sx), sp.diff(G4, Xs).subs(sx), F4.subs(sx), 0, 0, 0)

    g = S.metrica_esferica(Nc ** 2 - v ** 2, v, sp.Integer(1), r ** 2)
    gi, Gam, Ric, Rs = S.geometria(g)
    Lcov, esc = S.lagrangiana(g, gi, Gam, Ric, Rs, qt * tau, teoria)
    Lcov = sp.simplify(Lcov.subs(S.th, sp.pi / 3).xreplace({sp.sqrt(Nc ** 2): Nc}))
    assert sp.simplify(esc["X"] - qt ** 2 / (2 * Nc ** 2)) == 0
    log("lagrangiana covariante en la metrica estacionaria")
    G3y = G3.subs(Xs, ys ** 2 / 2)
    A3y = ys * G3y - sp.integrate(G3y, (ys, 0, ys))
    sy = {ys: qt / Nc}
    sX = {Xs: qt ** 2 / (2 * Nc ** 2)}
    A2 = G2.subs(sX)
    A3 = A3y.subs(sy)
    A4 = (-G4 + 2 * Xs * sp.diff(G4, Xs) + 4 * Xs ** 2 * F4).subs(sX)
    B4 = G4.subs(sX)
    dif = Lcov - lest(A2, A3, A4, B4)
    el = euler_lagrange(dif, [Nc, v, r])
    RES["EL_diferencia"] = {k_: str(x) for k_, x in el.items()}
    log(f"(1) Euler-Lagrange de L_cov - L_est: {el}")
    assert all(x == 0 for x in el.values())
    # ------------------------------------------------------------------ (2) Schwarzschild en Painleve-Gullstrand
    M = sp.Symbol("M", positive=True)
    egr = euler_lagrange(lest(0, 0, -1, 1), [Nc, v, r])
    pg = {Nc: 1, v: sp.sqrt(2 * M / rho), r: rho}
    res = {k_: sp.simplify(x.subs(pg).doit()) for k_, x in egr.items()}
    RES["PG"] = {k_: str(x) for k_, x in res.items()}
    log(f"(2) GR, Schwarzschild en Painleve-Gullstrand: {res}")
    assert all(x == 0 for x in res.values())
    json.dump(RES, open(os.path.join(AQUI, "resultados", "adm_estatico.json"), "w"), indent=1)
