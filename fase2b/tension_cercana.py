# -*- coding: utf-8 -*-
"""
Cascara con tension dependiente del escalar, Lambda_s(phi), en el limite de zona cercana (§13 de ley_analitica.py).

Cascara:  L_s = (1/8 pi) [ -Lambda_s(phi) + kappa F(phi) Rbar ] - n m(phi),  kappa = 4 zeta a,  Lambda_s(phi0) = -2a.
  - La tension aporta -Lambda_s a 8 pi S^t_t y a 8 pi S^th_th, y +2 Lambda_s' a la carga escalar sigma_phi = -16 pi dL_s/dphi.
  - Empalme del escalar del fondo: 2 rho beta = d[2] + 4 kappa F'/r0^2 - 2 Lambda_s'(phi0)  (el polvo compensa la carga
    de la tension; con Lambda_s' ~ a^2 es rho beta ~ -Lambda_s').
  - Perturbacion: delta(8 pi S^t_t) y delta(8 pi S^th_th) ganan -Lambda_s' q0; delta sigma_phi gana
    2 (Lambda_s'' + rho beta mu) q0 (tension y polvo entran juntos: solo cuenta la suma).
Escalas de zona cercana (r0 = 1):  Lambda_s'(phi0) = a^2 l1,  Lambda_s''(phi0) + rho beta mu = a s2.
Con l1 = s2 = 0 se recupera el polinomio de ley_analitica.py.
Salida: resultados/tension_cercana.json (polinomio P(s; zeta, u, F0, F2, l1, s2)).   Uso: python tension_cercana.py
"""
import sys
import json
import sympy as sp
import egb_F_gen as H

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

X = sp.symbols("X", nonnegative=True)
a, u, F2, F3, F0, zeta, W, f, c, l1, s2 = sp.symbols("a u F2 F3 F0 zeta W f c l1 s2", real=True)
s = sp.symbols("s", positive=True)
Z = 1 + X
L = sp.log(Z)
G0 = (1 + u**2 / 4) / (4 * u)
w = -a**2 * W


def orden_dominante(e, nmax=4, nmin=-3):
    e = sp.expand(e)
    for n in range(nmax, nmin - 1, -1):
        cc = sp.simplify(sp.limit(e / a**n, a, sp.oo))
        if cc != 0:
            return n, cc
    return None, sp.Integer(0)


# solucion de zona cercana (ley_analitica.py, §13.2): desplazamiento f + onda escalar c Z^-s
S = sp.Rational(1, 2) * (Z - (1 + 4 * F2 * u**2 * L) / Z)
A0 = f * (W - 1)
Q = u * A0 / ((W - 1) * Z) + c * Z**(-s)
V = (2 * u**2 * (Q + Z * sp.diff(Q, X)) - u * A0 * (4 * Z**2 + u**2) / Z) / (1 - W) - 16 * F2 * u**2 * L * Q
D = (V - (4 + u**2) * Q) / (8 * u)
y0 = dict(a1=a**2 * A0 / Z**2, d0=D / a, d1=sp.diff(D, X), q0=Q, q1=a * sp.diff(Q, X))
y0 = {k: sp.simplify(v.subs(X, 0)) for k, v in y0.items()}
bg0 = dict(r=sp.Integer(1), r1=sp.Integer(0), r2=sp.simplify(sp.diff(S, X).subs(X, 0)), P1=a, P2=-a**2, p1=u,
           p2=-a * u, e2=sp.Integer(1))
F1_0 = G0 / a
lam1 = a**2 * l1          # Lambda_s'(phi0)
sig2 = a * s2             # Lambda_s''(phi0) + rho beta mu


def capa(zeta):
    kappa, Lam = 4 * zeta * a, -2 * a
    M0 = sp.Matrix(H.est_principal(F1_0, bg0["r"], bg0["r1"], bg0["P1"], bg0["p1"]))
    dd = M0 * sp.Matrix([0, 2 * a, 2 * u])
    assert sp.simplify(dd[1] + Lam) == 0
    rho = -dd[0] + 2 * kappa * F0 - Lam
    rhobeta = (dd[2] + 4 * kappa * F1_0 - 2 * lam1) / 2
    return kappa, rho, rhobeta


def fuentes(zeta, d0, q0):
    kappa, rho, rhobeta = capa(zeta)
    Stt = 2 * rho * d0 - rhobeta * q0 + 2 * kappa * (F1_0 * q0 - 2 * F0 * d0) - lam1 * q0
    Sthth = -2 * kappa * w * (F0 * d0 + F1_0 * q0) - lam1 * q0
    dR = -4 * d0 - 4 * w * d0
    sphi = -2 * kappa * (F2 * q0 * 2 + F1_0 * dR) + 2 * rhobeta * (-2 * d0) + 2 * sig2 * q0
    return Stt, Sthth, sphi


rr, rr1, rr2, PP1, PP2, pp1, pp2, ee2 = (bg0[k] for k in ("r", "r1", "r2", "P1", "P2", "p1", "p2", "e2"))
yy = [y0[k] for k in ("a1", "d0", "d1", "q0", "q1")]
M0 = sp.Matrix(H.lin_principal(F1_0, rr, rr1, PP1, pp1))
Pp, pp = sp.symbols("Pp pp")
R = sp.Matrix(H.lin_resto(F1_0, F2, F3, w, rr, rr1, rr2, PP1, Pp, pp1, pp, ee2, *yy))
R0 = R.subs({Pp: 0, pp: 0})
delta = M0 * sp.Matrix([2 * y0["d1"], 2 * y0["a1"], 2 * y0["q1"]]) + (R.subs({Pp: 1, pp: 0}) - R0) * 2 * PP1 \
    + (R.subs({Pp: 0, pp: 1}) - R0) * 2 * pp1
lig0 = H.lin_ligaduras(F1_0, F2, F3, w, rr, rr1, rr2, PP1, PP2, pp1, pp2, ee2, *yy)
St = fuentes(zeta, y0["d0"], y0["q0"])
filas = dict(lig_t=lig0[0], lig_l=lig0[1], J_tt=delta[0] - St[0], J_thth=delta[1] - St[1], J_phi=delta[2] - St[2])
print("Empalmes y ligaduras en l = 0+ a orden dominante, con Lambda_s' = a^2 l1 y Lambda_s'' + rho beta mu = a s2:")
coef = {}
for k, e in filas.items():
    n, cc = orden_dominante(e)
    cc = sp.expand(sp.simplify(cc.subs(W, s**2)))
    cf, cs = sp.factor(cc.coeff(f)), sp.factor(cc.coeff(c))
    coef[k] = (cf, cs)
    print(f"   {k:7s} (orden a^{n}):  coef f = {cf};   coef c = {cs}")
tf, tc = coef["J_thth"]
pf, pc = coef["J_phi"]
Ps = sp.factor(sp.numer(sp.together(tf * pc - tc * pf)))
Ppol = sp.Poly(sp.expand(Ps), s)
print("\nPolinomio de grado", Ppol.degree(), "; coeficiente principal:", sp.factor(Ppol.LC()))
print("zeta = 0:", sp.factor(Ps.subs(zeta, 0)))
print("l1 = s2 = 0 (comprobacion frente a ley_analitica.json):")
ley = json.load(open("resultados/ley_analitica.json", encoding="utf-8"))
P0 = sp.sympify(ley["polinomio"], locals=dict(s=s, u=u, F0=F0, F2=F2, zeta=zeta))
print("   cociente =", sp.simplify(Ps.subs({l1: 0, s2: 0}) / P0))
json.dump(dict(polinomio=str(sp.expand(Ps)), variable="s", parametros=["u", "F0", "F2", "zeta", "l1", "s2"],
               escalas="Lambda_s'(phi0) = a^2 l1 r0;  Lambda_s''(phi0) + rho beta mu = a s2",
               zeta_0=str(sp.factor(Ps.subs(zeta, 0))), coeficiente_principal=str(sp.factor(Ppol.LC())),
               coef=dict((k, [str(v[0]), str(v[1])]) for k, v in coef.items())),
          open("resultados/tension_cercana.json", "w", encoding="utf-8"), indent=1)
print("\nGuardado: resultados/tension_cercana.json")
