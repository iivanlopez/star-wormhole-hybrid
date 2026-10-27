# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): la lagrangiana ADM dependiente del tiempo (gauge h_rhorho = 1, phi = qt tau) frente a la covariante.
    ds^2 = -(Nc^2 - v^2) dtau^2 + 2 v dtau drho + drho^2 + r^2 dOmega^2,   Nc, v, r funciones de (tau, rho),
    L_ADM = Nc r^2 [A2 + A3 K + A4 (4 K_r K_th + 2 K_th^2) + B4 R3],
    K_r = -v'/Nc,  K_th = (r_tau - v r')/(r Nc),  K = K_r + 2 K_th,  R3 = (2 - 4 r r'' - 2 r'^2)/r^2,
con A2 = G2, A3 = y G3 - int_0^y G3 dy, A4 = -G4 + 2X G4_X + 4X^2 F4, B4 = G4 (y = sqrt(2X) = qt/Nc).
Comprobacion: el operador de Euler-Lagrange (en tau y rho) de L_cov - L_ADM se anula identicamente para Nc, v, r
(teoria polinomica con coeficientes simbolicos).  Es la extension dependiente del tiempo y con desplazamiento de
adm_estatico.py (estatico) y adm_dinamico.py (v = 0).
Uso: python adm_tiempo_cov.py
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
rho, tau = S.l, S.t
qt = sp.Symbol("qt", positive=True)
Nc, v, r = (sp.Function(n)(tau, rho) for n in ("Nc", "v", "r"))
Xs = sp.Symbol("X", positive=True)
ys = sp.Symbol("y", positive=True)
k = sp.symbols("k1:3", real=True)
c = sp.symbols("c1:3", real=True)
b = sp.symbols("b1:3", real=True)
f = sp.symbols("f0:2", real=True)
G2 = k[0] * Xs + k[1] * Xs ** 2
G3 = c[0] * Xs + c[1] * Xs ** 2
G4 = 1 + b[0] * Xs + b[1] * Xs ** 2
F4 = f[0] + f[1] * Xs


def teoria(phi, X):
    sx = {Xs: X}
    return (G2.subs(sx), G3.subs(sx), G4.subs(sx), sp.diff(G4, Xs).subs(sx), F4.subs(sx), 0, 0, 0)


g = S.metrica_esferica(Nc ** 2 - v ** 2, v, sp.Integer(1), r ** 2)
gi, Gam, Ric, Rs = S.geometria(g)
Lcov, esc = S.lagrangiana(g, gi, Gam, Ric, Rs, qt * tau, teoria)
raiz = {sp.sqrt(Nc ** 2): Nc}
Lcov = sp.expand(Lcov.subs(S.th, sp.pi / 3).xreplace(raiz))
print(f"[{time.time() - T0:6.1f} s] lagrangiana covariante", flush=True)
G3y = G3.subs(Xs, ys ** 2 / 2)
A3 = (ys * G3y - sp.integrate(G3y, (ys, 0, ys))).subs(ys, qt / Nc)
sX = {Xs: qt ** 2 / (2 * Nc ** 2)}
A2 = G2.subs(sX)
A4 = (-G4 + 2 * Xs * sp.diff(G4, Xs) + 4 * Xs ** 2 * F4).subs(sX)
B4 = G4.subs(sX)
D = lambda e, x=rho, n=1: sp.diff(e, x, n)              # noqa: E731
Kr = -D(v) / Nc
Kth = (D(r, tau) - v * D(r)) / (r * Nc)
R3 = (2 - 4 * r * D(r, rho, 2) - 2 * D(r) ** 2) / r ** 2
Ladm = Nc * r ** 2 * (A2 + A3 * (Kr + 2 * Kth) + A4 * (4 * Kr * Kth + 2 * Kth ** 2) + B4 * R3)
dif = sp.expand(Lcov - Ladm)
RES = {}
for fn in (Nc, v, r):
    e = sp.euler_equations(dif, [fn], [tau, rho])
    val = 0 if not e else sp.simplify(sp.expand(e[0].lhs.xreplace(raiz)))
    RES[str(fn.func)] = str(val)
    print(f"[{time.time() - T0:6.1f} s] E-L de L_cov - L_ADM respecto de {fn.func}: {val}", flush=True)
    assert val == 0
json.dump(RES, open(os.path.join(AQUI, "resultados", "adm_tiempo_cov.json"), "w"), indent=1)
