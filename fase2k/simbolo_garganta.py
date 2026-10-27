# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): simbolo principal del monopolo en la garganta con desplazamiento de fondo (v != 0).

Lagrangiana ADM dependiente del tiempo (adm_tiempo_cov.py), gauge h = 1, perturbaciones Nc = N + nu, v = V + w,
r = R + d (funciones de t y rho).  Ecuaciones lineales de Euler-Lagrange; onda plana e^{i(k rho - om t)} con los
coeficientes congelados en la garganta: R' = V = N' = 0 y funciones de la teoria pares (a_i' = 0, a_i'' = a_iN N'').
det M(om, k) para (nu, w, d); con om = c k y k -> oo, el coeficiente dominante P(c) da las velocidades caracteristicas
radiales del monopolo en el corte phi = cte.  Caso v = 0 (V' = 0) como control: c^2 = 8 beta_N^2 (Q r'')/(kappa r Q^2)
(informe, §6).  Ademas: hacer V'(0) != 0 (flujo que atraviesa la garganta).
Uso: python simbolo_garganta.py
"""
import os
import json
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
t, rho, e = sp.symbols("t rho e", real=True)
G = lambda n: sp.Function(n)(rho)                       # noqa: E731
H = lambda n: sp.Function(n)(t, rho)                    # noqa: E731
N, V, R = G("N"), G("V"), G("R")
NOMS = ("a2", "a3", "a4", "b4")
fun = {n: [G(n), G(n + "N"), G(n + "NN")] for n in NOMS}
nu, w, d = H("nu"), H("w"), H("d")
Dr = lambda x, k_=1: sp.diff(x, rho, k_)                # noqa: E731
Dt = lambda x: sp.diff(x, t)                            # noqa: E731
Nc, v, r = N + e * nu, V + e * w, R + e * d
A = {n: fun[n][0] + fun[n][1] * e * nu + fun[n][2] * (e * nu) ** 2 / 2 for n in NOMS}
Kr = -Dr(v) / Nc
Kth = (Dt(r) - v * Dr(r)) / (r * Nc)
R3 = (2 - 4 * r * Dr(r, 2) - 2 * Dr(r) ** 2) / r ** 2
L = Nc * r ** 2 * (A["a2"] + A["a3"] * (Kr + 2 * Kth) + A["a4"] * (4 * Kr * Kth + 2 * Kth ** 2) + A["b4"] * R3)
L2 = sp.expand(sp.diff(L, e, 2).subs(e, 0) / 2)
EL = {str(f.func): sp.expand(sp.euler_equations(L2, [f], [t, rho])[0].lhs) for f in (nu, w, d)}

# onda plana y congelacion en la garganta
om, k = sp.symbols("omega k")
amp = {nu: sp.Symbol("Nh"), w: sp.Symbol("Wh"), d: sp.Symbol("Dh")}
fase = sp.exp(sp.I * (k * rho - om * t))
# valores en la garganta: funciones pares (R, N, a_i) e impar (V)
s = {}
R0, R2, R4 = sp.symbols("R0 R2 R4", positive=True)
N0, N2, N4 = sp.symbols("N0 N2 N4", real=True)
V1, V3 = sp.symbols("V1 V3", real=True)
s.update({R: R0, Dr(R): 0, Dr(R, 2): R2, Dr(R, 3): 0, Dr(R, 4): R4,
          N: N0, Dr(N): 0, Dr(N, 2): N2, Dr(N, 3): 0, Dr(N, 4): N4,
          V: 0, Dr(V): V1, Dr(V, 2): 0, Dr(V, 3): V3, Dr(V, 4): 0})
vals = {}
for n in NOMS:
    f0, f1, f2 = fun[n]
    a0, a1, a2_, a3_ = sp.symbols(f"{n}_0 {n}_N {n}_NN {n}_NNN", real=True)
    vals[n] = (a0, a1, a2_, a3_)
    # f0(rho) = A(N(rho)): f0' = A_N N', f0'' = A_NN N'^2 + A_N N''; f1 = A_N(N(rho)): f1' = A_NN N', f1'' = A_NNN N'^2 + A_NN N''
    s.update({f0: a0, Dr(f0): 0, Dr(f0, 2): a1 * N2, Dr(f0, 3): 0,
              f1: a1, Dr(f1): 0, Dr(f1, 2): a2_ * N2, Dr(f1, 3): 0,
              f2: a2_, Dr(f2): 0, Dr(f2, 2): a3_ * N2})


def congela(ex):
    ex = ex.subs({f: amp[f] * fase for f in (nu, w, d)}).doit()
    ex = sp.expand(sp.simplify(ex / fase))
    # sustituir primero las derivadas mas altas
    for n_ in range(4, -1, -1):
        ex = ex.subs({kk: vv for kk, vv in s.items() if isinstance(kk, sp.Derivative) and kk.derivative_count == n_})
    return sp.expand(ex.subs({kk: vv for kk, vv in s.items() if not isinstance(kk, sp.Derivative)}).subs(rho, 0))


M = sp.Matrix([[sp.expand(congela(EL[q]).coeff(amp[f])) for f in (nu, w, d)] for q in ("nu", "w", "d")])
det = sp.expand(M.det())
c = sp.Symbol("c")
P = sp.Poly(sp.expand(det.subs(om, c * k)), k)
grado = P.degree()
lider = sp.factor(P.coeff_monomial(k ** grado))
print("grado en k del determinante:", grado, flush=True)
print("coeficiente dominante P(c):", lider, flush=True)
# control v = 0
lider0 = sp.factor(lider.subs(V1, 0))
print("con V'(0) = 0:", lider0, flush=True)
# ecuacion efectiva para d: S_eff = det M / det M_aux (M_aux = bloque de nu y w); partes principales
Mx = M[:2, :2]
detx = sp.expand(Mx.det())
Px = sp.Poly(sp.expand(detx.subs(om, c * k)), k)
lid_x = sp.factor(Px.coeff_monomial(k ** Px.degree()))
print(f"det M_aux ~ k^{Px.degree()} * [{lid_x}]", flush=True)
print(f"S_eff ~ k^{grado - Px.degree()} * P(c)/D(c)", flush=True)
# ligadura de momento en la garganta: A3' = 2 V (Q R')'/R  ->  a3_N N2 = 2 V1 Q(0) R2/R0, Q = 2 a4/N
a3N, a4_0 = vals["a3"][1], vals["a4"][0]
lig = {a3N: 4 * V1 * a4_0 * R2 / (N0 * R0 * N2)}
P_l = sp.factor(lider.subs(lig))
D_l = sp.factor(lid_x.subs(lig))
print("P(c) con la ligadura de momento:", P_l, flush=True)
print("D(c) con la ligadura de momento:", D_l, flush=True)
Pc = sp.Poly(sp.expand(sp.numer(sp.together(P_l))), c)
cin, gra = sp.factor(Pc.coeff_monomial(c ** 2)), sp.factor(Pc.coeff_monomial(1))
print("P: coeficiente de c^2:", cin, flush=True)
print("P: termino sin c:", gra, flush=True)
json.dump({"grado": int(grado), "P": str(lider), "P_v0": str(lider0), "D": str(lid_x), "grado_D": int(Px.degree()),
           "P_ligadura": str(P_l), "D_ligadura": str(D_l), "cin": str(cin), "gra": str(gra),
           "M": [[str(x) for x in row] for row in M.tolist()]},
          open(os.path.join(AQUI, "resultados", "simbolo_garganta.json"), "w"), indent=1)
