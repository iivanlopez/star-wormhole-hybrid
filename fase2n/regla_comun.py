# -*- coding: utf-8 -*-
"""
Fase 2n: la regla de signos comun de las gargantas sostenidas por un campo temporal (2k §6 y 2m §4.5).

Las dos acciones cineticas del monopolo se escriben  T = int [ -(2 G/N) d_t^2 + k ((G d_t / N)')^2 ]:
  2k (phi = q t, psi' = 0):  Q = 2 A4 / N  con  A4 = -G4 + 2 X G4X + 4 X^2 F4 = -Gcal,   T = int [Q d^2 + r ((Q d)')^2/(Q r')'],
                             luego k = 4 r / (Q r')';  en la garganta (Q r')' = Q r''.
  2m (A1 = 0, G3 = 0):       k = k2 = -4 r^2 N^4 / (U^2 ct),  ct(0) = N0 r0 Gcal r''(0) / X  (X > 0).
Se comprueba que en la garganta k(0) = -kappa0 / (Gcal(0) r''(0)) con kappa0 > 0 en los dos casos: los dos terminos
cineticos tienen el signo de -Gcal(0) cuando r''(0) > 0.
Uso: python regla_comun.py
"""
import os
import json
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
x = sp.Symbol("rho", real=True)
G4, G4X, F4, X = sp.symbols("G4 G4X F4 X", real=True)
Gcal = G4 - 2 * X * G4X - 4 * X ** 2 * F4
A4 = -G4 + 2 * X * G4X + 4 * X ** 2 * F4
res = {"A4_igual_menos_G": sp.simplify(A4 + Gcal) == 0}

# 2k: T con Q = 2 A4/N = -2 G/N  ->  forma comun
G, N, r, dd = (sp.Function(n)(x) for n in ("G", "N", "r", "d"))
Q = -2 * G / N
T2k = Q * dd ** 2 + r * sp.diff(Q * dd, x) ** 2 / sp.diff(Q * sp.diff(r, x), x)
k2k = 4 * r / sp.diff(Q * sp.diff(r, x), x)
Tcom = -2 * G / N * dd ** 2 + k2k * sp.diff(G * dd / N, x) ** 2
res["2k_forma_comun"] = sp.simplify(T2k - Tcom) == 0
# en la garganta: r' = 0  ->  (Q r')' = Q r''
r0, rpp, G0, N0, U0, kap = sp.symbols("r0 rpp G0 N0 U0 kappa", positive=True)
Gs = sp.Symbol("Gs", real=True)                   # Gcal(0) con signo libre
k2k_0 = 4 * r0 / ((-2 * Gs / N0) * rpp)
res["2k_k0"] = str(sp.simplify(k2k_0))
# 2m: k2 = -4 r^2 N^4/(U^2 ct), ct(0) = N0 r0 G rpp / X
Xp = sp.Symbol("Xp", positive=True)
k2m_0 = -4 * r0 ** 2 * N0 ** 4 / (U0 ** 2 * (N0 * r0 * Gs * rpp / Xp))
res["2m_k0"] = str(sp.simplify(k2m_0))
for nom, k0 in (("2k", k2k_0), ("2m", k2m_0)):
    kappa0 = sp.simplify(-k0 * Gs * rpp)
    res[f"{nom}_kappa0"] = str(kappa0)
    res[f"{nom}_kappa0_pos"] = bool(kappa0.is_positive)
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(AQUI, "resultados", "regla_comun.json"), "w"), indent=1)
assert all(v for v in res.values() if isinstance(v, bool))
print("regla comun: OK")
