# -*- coding: utf-8 -*-
"""
Fase 2n, parte 1: lemas exactos sobre el galileon cubico con flujo (fase2k, Parte II).  Todo simbolico (sympy).

  (L1) r^3 L[r'] = (r K)'  con  L[f] = (f'/r)' + K f / r^3,  K = 1 + r r'' - r'^2   (identidad para todo r(rho)).
  (L2) Con la masa de Misner-Sharp del corte, r'^2 = 1 - 2 m(r)/r:  (r K)' = -r' r^3 d/dr (m_r / r^2)
       = -4 pi r' r^3 de/dr, con e = m_r/(4 pi r^2) la densidad efectiva.
  (L3) Limite de flujo debil del coeficiente cinetico: s2 = -2 B/(r v^2 r''^2) con B -> -2 beta' (r K)' cuando
       (v, v', v'') = eps (a, b, c), eps -> 0  (a partir del simbolo exacto de fase2k/resultados/simbolo_perfiles.pkl).
       Se comprueba tambien que B es afin en Q' con coeficiente r^4 w' (r' w' - r'' w), w = v/r.
  (L4) Garganta: la EDO de beta da beta''(0) = -(beta0/R0^2) [1 + R0 R2 + x/2], x = Q0 R0^2 V1^2 / beta0, y la
       condicion de gradiente sano con A4(0) < 0 (Lambda_hat < 0) implica beta''(0) > R2 beta0 / (2 R0) > 0.
  (L5) Ninguna identidad lineal  int mu B = (sin beta)  con peso mu independiente de Q tiene integrando de signo
       definido: tras quitar Q' por partes queda Q * C con C conteniendo -2 mu r^4 r' w' w''.
Uso: python lemas.py
"""
import os
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
x = sp.Symbol("rho", real=True)
F = lambda n: sp.Function(n)(x)                      # noqa: E731
D = lambda f, n=1: sp.diff(f, x, n)                  # noqa: E731
res = {}

# (L1)
r = F("r")
K = 1 + r * D(r, 2) - D(r) ** 2
Lop = lambda f: D(D(f) / r) + K * f / r ** 3          # noqa: E731
res["L1"] = sp.simplify(r ** 3 * Lop(D(r)) - D(r * K)) == 0
print("(L1) r^3 L[r'] - (rK)' = 0:", res["L1"])

# (L2): r'^2 = 1 - 2m(r)/r  =>  r'' = m/r^2 - m_r/r ;  r''' = r' d/dr(r'')
R = sp.Symbol("R", positive=True)
m = sp.Function("m")(R)
rp = sp.sqrt(1 - 2 * m / R)
rpp = sp.diff(rp ** 2, R) / 2
rppp = rp * sp.diff(rpp, R)
rK_der = rp - rp ** 3 + R ** 2 * rppp                 # (rK)' = r' - r'^3 + r^2 r'''
obj = -rp * R ** 3 * sp.diff(sp.diff(m, R) / R ** 2, R)
res["L2"] = sp.simplify(rK_der - obj) == 0
KMS = sp.simplify(1 + R * rpp - rp ** 2)
res["L2_K"] = str(KMS)
print("(L2) (rK)' = -r' r^3 (m_r/r^2)_r:", res["L2"], "  K =", KMS)

# (L3): simbolo exacto de la fase 2k
pp = pickle.load(open(os.path.join(AQUI, "..", "fase2k", "resultados", "simbolo_perfiles.pkl"), "rb"))
SIM = {}
for n in ("r", "v", "Q", "Nc", "beta"):
    f = sp.Function(n)(x)
    SIM[f] = sp.Symbol(f"{n}0")
    for k in range(1, 5):
        SIM[sp.Derivative(f, (x, k))] = sp.Symbol(f"{n}{k}")


def a_simbolos(e):
    e = e.xreplace({k: v for k, v in SIM.items() if isinstance(k, sp.Derivative)})
    return e.xreplace({k: v for k, v in SIM.items() if not isinstance(k, sp.Derivative)})


s2 = sp.together(a_simbolos(pp["P"][2][0] / pp["P"][2][1]) / a_simbolos(pp["D"]))
r0, r1, r2, r3, v0, v1, v2, Q0, Q1, b1 = sp.symbols("r0 r1 r2 r3 v0 v1 v2 Q0 Q1 beta1")
B = sp.factor(-s2 * r0 * v0 ** 2 * r2 ** 2 / 2)      # B definido por s2 = -2 B / (r v^2 r''^2)
rKp = r1 - r1 ** 3 + r0 ** 2 * r3
w0, w1 = v0 / r0, (v1 * r0 - r1 * v0) / r0 ** 2
w2 = sp.diff(w1, v1) * v2 + sp.diff(w1, v0) * v1 + sp.diff(w1, r0) * r1 + sp.diff(w1, r1) * r2
W, W1 = r0 ** 2 * w1, sp.Symbol("W1")
# B de fase2k: -2 beta'(rK)' - Q v r''(r^2 W)'/r + Q' W (r' W - r v r'')
r2W_der = 2 * r0 * r1 * W + r0 ** 2 * (r0 ** 2 * w2 + 2 * r0 * r1 * w1)
B_2k = -2 * b1 * rKp - Q0 * v0 * r2 * r2W_der / r0 + Q1 * W * (r1 * W - r0 * v0 * r2)
res["B_igual_2k"] = sp.simplify(B - B_2k) == 0
eps, a_, b_, c_ = sp.symbols("epsilon a b c")
lim = sp.limit(B.subs({v0: eps * a_, v1: eps * b_, v2: eps * c_}), eps, 0)
res["L3_limite"] = str(sp.factor(lim))
res["L3"] = sp.simplify(lim + 2 * b1 * rKp) == 0
coefQ1 = sp.factor(sp.diff(B, Q1))
res["L3_coefQ1"] = sp.simplify(coefQ1 - r0 ** 4 * w1 * (r1 * w1 - r2 * w0)) == 0
print("(L3) B coincide con la forma de fase2k:", res["B_igual_2k"], "; limite flujo debil:", lim,
      "; coef. de Q' = r^4 w'(r'w' - r''w):", res["L3_coefQ1"])

# (L4): serie en la garganta
R0, R2, R4, V1, V3, Q0s, Q2, b0, b2 = sp.symbols("R0 R2 R4 V1 V3 Q0 Q2 beta0 beta2")
rs = R0 + R2 * x ** 2 / 2 + R4 * x ** 4 / 24
vs = V1 * x + V3 * x ** 3 / 6
Qs = Q0s + Q2 * x ** 2 / 2
bs = b0 + b2 * x ** 2 / 2
Ks = 1 + rs * sp.diff(rs, x, 2) - sp.diff(rs, x) ** 2
ode = sp.diff(sp.diff(bs, x) / rs, x) + Ks * bs / rs ** 3 + sp.diff(Qs * rs ** 4 * sp.diff(vs ** 2 / rs ** 2, x), x) / (4 * rs ** 3)
b2sol = sp.solve(sp.series(ode, x, 0, 1).removeO(), b2)[0]
xx = sp.Symbol("x")
res["L4_beta2"] = str(sp.factor(b2sol))
res["L4_forma"] = sp.simplify(b2sol + b0 / R0 ** 2 * (1 + R0 * R2 + xx / 2)).subs(xx, Q0s * R0 ** 2 * V1 ** 2 / b0) == 0
res["L4_forma"] = sp.simplify(b2sol + b0 / R0 ** 2 * (1 + R0 * R2 + Q0s * R0 ** 2 * V1 ** 2 / (2 * b0))) == 0
# Lambda_hat de fase2k §12 en la variable x:  [..] = beta0^2 (x^2 + 2 R0 R2 x - 4 (1 + R0 R2)^2)
corch = Q0s ** 2 * R0 ** 4 * V1 ** 4 + 2 * Q0s * R0 ** 3 * R2 * V1 ** 2 * b0 - 4 * b0 ** 2 * (1 + R0 * R2) ** 2
res["L4_corchete"] = sp.simplify(corch - b0 ** 2 * (xx ** 2 + 2 * R0 * R2 * xx - 4 * (1 + R0 * R2) ** 2)
                                 .subs(xx, Q0s * R0 ** 2 * V1 ** 2 / b0)) == 0
# x < x_- = -R0R2 - sqrt(R0^2R2^2 + 4(1+R0R2)^2)  =>  -x/2 - (1+R0R2) > R0R2/2 + [sqrt(..) - 2(1+R0R2)]/2 >= R0R2/2
a = sp.Symbol("a", positive=True)                     # a = R0 R2 > 0
cota = sp.sqrt(a ** 2 + 4 * (1 + a) ** 2) / 2 - (1 + a)
res["L4_cota"] = bool(sp.simplify((sp.sqrt(a ** 2 + 4 * (1 + a) ** 2)) ** 2 - (2 * (1 + a)) ** 2) == a ** 2)
print("(L4) beta''(0) =", b2sol, "; forma:", res["L4_forma"], "; corchete:", res["L4_corchete"],
      "; sqrt(a^2+4(1+a)^2) > 2(1+a):", res["L4_cota"])

# (L5): int mu B, quitando Q' por partes: coeficiente de w'' en C
mu, wf, Qf = F("mu"), F("w"), F("Q")
flujo = -Qf * wf * D(r, 2) * D(r ** 4 * D(wf)) + D(Qf) * r ** 4 * D(wf) * (D(r) * D(wf) - D(r, 2) * wf)
integ = mu * flujo
# mu Q' G = (mu Q G)' - Q (mu G)'
G = r ** 4 * D(wf) * (D(r) * D(wf) - D(r, 2) * wf)
C = sp.expand((integ - mu * D(Qf) * G - Qf * D(mu * G)) / Qf)
res["L5_sinQp"] = sp.simplify(sp.diff(C, D(Qf))) == 0
c_w2 = sp.factor(sp.expand(C).coeff(D(wf, 2)))
res["L5_coef_wpp"] = str(c_w2)
res["L5"] = sp.simplify(c_w2 + 2 * mu * r ** 4 * D(r) * D(wf)) == 0
print("(L5) coeficiente de w'' en C:", c_w2, " (= -2 mu r^4 r' w'):", res["L5"])

json.dump(res, open(os.path.join(AQUI, "resultados", "lemas.json"), "w"), indent=1)
assert all(v for k, v in res.items() if isinstance(v, bool)), res
print("todas las comprobaciones: OK")
