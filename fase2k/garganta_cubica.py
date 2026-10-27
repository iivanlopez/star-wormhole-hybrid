# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): el simbolo principal del monopolo en una garganta con flujo (v = V1 rho + ..., A3_N != 0) en
funcion de los datos de Taylor del fondo en la garganta, usando las ecuaciones del fondo.

Simbolo (simbolo_garganta.py): S_eff ~ k^2 [Gam c^2 + Lam] / (N0^2 R2^2 V1^2), con
    Gam = N0^2 N2^2 kappa - N0^2 N2^2 R0^2 V1 a3_NN + V1^2 (8 N0 N2 R0 R2 a4_N - 12 N0 R2^2 a4 - 8 N2 R0 R2 a4),
    Lam = -4 N0^2 R2 [N0 N2^2 R0 beta_N^2 + V1^2 a4 (2 N0 N2 R0 b4_N - N0 R2 b4 + 2 N2 R0 b4)] / a4,
    kappa = R0^2 p_NN + (2 - 4 R0 R2) beta_NN   (p = N A2, beta = N B4);
velocidades reales y energia cinetica positiva (sin fantasma, sin inclinacion en la garganta) <=> Gam > 0, Lam > 0
... con la normalizacion S_eff = 2 A om^2 - 2 B k^2: A ~ Gam, B ~ -Lam.  (Se recalcula aqui con los datos de Taylor.)
Fondo en la garganta (rho = 0; pares N, r, Q, beta, p; impar v):
    H = 2 beta (r'^2 - 1) + 4 beta' r r' + Q r' (v^2 r)' - p r^2 = 0,
    (beta'/r)' + K beta/r^3 = -(Q r^4 (v^2/r^2)')'/(4 r^3),      r A3' = 2 v (Q r')'.
Uso: python garganta_cubica.py
"""
import os
import json
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
x = sp.Symbol("rho", real=True)
O = 6
R0, R2, R4 = sp.symbols("R0 R2 R4", real=True)
N0, N2, N4 = sp.symbols("N0 N2 N4", real=True)
Q0, Q2, Q4 = sp.symbols("Q0 Q2 Q4", real=True)
V1, V3, V5 = sp.symbols("V1 V3 V5", real=True)
b0, b2, b4_, p0, p2, p4 = sp.symbols("beta0 beta2 beta4 p0 p2 p4", real=True)
r = R0 + R2 * x ** 2 / 2 + R4 * x ** 4 / 24
N = N0 + N2 * x ** 2 / 2 + N4 * x ** 4 / 24
Q = Q0 + Q2 * x ** 2 / 2 + Q4 * x ** 4 / 24
v = V1 * x + V3 * x ** 3 / 6 + V5 * x ** 5 / 120
be = b0 + b2 * x ** 2 / 2 + b4_ * x ** 4 / 24
p = p0 + p2 * x ** 2 / 2 + p4 * x ** 4 / 24
D = lambda f: sp.diff(f, x)                         # noqa: E731


def serie(e, n=O):
    return sp.expand(sp.series(e, x, 0, n).removeO())


H = serie(2 * be * (D(r) ** 2 - 1) + 4 * D(be) * r * D(r) + Q * D(r) * D(v ** 2 * r) - p * r ** 2, 5)
K = 1 + r * D(D(r)) - D(r) ** 2
ecb = serie((D(D(be) / r) + K * be / r ** 3 + D(Q * r ** 4 * D(v ** 2 / r ** 2)) / (4 * r ** 3)) * r ** 3, 3)
sol = sp.solve([H.coeff(x, 0), H.coeff(x, 2), H.coeff(x, 4), ecb.coeff(x, 0), ecb.coeff(x, 2)],
               [p0, p2, p4, b2, b4_], dict=True)[0]
print("fondo en la garganta:", {k: sp.factor(v_) for k, v_ in sol.items() if k in (p0, b2)}, flush=True)
A3p = serie(2 * v * D(Q * D(r)) / r, 5)                 # A3' (impar)
A3 = sp.integrate(A3p, x)                              # salvo constante


def derivadas_N(f):
    """f(N) a partir de las series pares f(rho), N(rho): f_N(0) y f_NN(0)."""
    s = sp.Symbol("s")                                # s = rho^2
    fs = sp.expand(f.subs(x, sp.sqrt(s)))
    Ns = sp.expand(N.subs(x, sp.sqrt(s)))
    f1 = sp.diff(fs, s) / sp.diff(Ns, s)
    f2 = sp.diff(f1, s) / sp.diff(Ns, s)
    return sp.simplify(f1.subs(s, 0)), sp.simplify(f2.subs(s, 0))


b4f = be / N
a4f = Q * N / 2
a2f = p / N
bN, bNN = derivadas_N(be)
pN, pNN = derivadas_N(p)
b4N, b4NN = derivadas_N(b4f)
a4N, _ = derivadas_N(a4f)
a3N, a3NN = derivadas_N(A3)
kappa = R0 ** 2 * pNN + (2 - 4 * R0 * R2) * bNN
a4v = Q0 * N0 / 2
b4v = b0 / N0
Gam = (N0 ** 2 * N2 ** 2 * kappa - N0 ** 2 * N2 ** 2 * R0 ** 2 * V1 * a3NN
       + V1 ** 2 * (8 * N0 * N2 * R0 * R2 * a4N - 12 * N0 * R2 ** 2 * a4v - 8 * N2 * R0 * R2 * a4v))
Lamh = N0 * N2 ** 2 * R0 * bN ** 2 + V1 ** 2 * a4v * (2 * N0 * N2 * R0 * b4N - N0 * R2 * b4v + 2 * N2 * R0 * b4v)
Gs = sp.factor(sp.simplify(Gam.subs(sol)))
Ls = sp.factor(sp.simplify(Lamh.subs(sol)))
print("Gam (kinetico, con el fondo):", Gs, flush=True)
print("Lam^ (gradiente; B > 0 <=> Lam^/A4 > 0):", Ls, flush=True)
# comprobacion de a3_N con la ligadura usada en simbolo_garganta.py: a3_N N2 = 4 V1 a4 R2/(N0 R0)
print("a3_N - 4 V1 a4 R2/(N0 R0 N2) =", sp.simplify(a3N - 4 * V1 * a4v * R2 / (N0 * R0 * N2)), flush=True)
json.dump({"fondo": {str(k): str(sp.factor(v_)) for k, v_ in sol.items()}, "Gam": str(Gs), "Lam": str(Ls)},
          open(os.path.join(AQUI, "resultados", "garganta_cubica.json"), "w"), indent=1)
