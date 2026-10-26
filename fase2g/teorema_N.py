# -*- coding: utf-8 -*-
"""
Fase 2g, parte I: teorema del taquion del monopolo en la teoria cuartica mas alla de Horndeski (GLPV), valido con
cualquier intensidad de campo.  Usa la lagrangiana cuadratica y la reduccion del monopolo de la Fase 2f
(fase2f/resultados/bH_simbolico.pkl y monopolo_reduccion.pkl).

Gauge unitario (delta phi = 0).  La ligadura de momento es
    W = gam c + al1 d' + al0 d,     con   al0 = -Phi' al1   (T2 de la Fase 2f, identidad exacta),
de modo que  W = gam c + al1 N (d/N)'.  La funcion de prueba d = N (el lapso) cumple la ligadura con c = 0, y la accion
reducida se queda en la parte de d sola, L_rest(c = 0, d) = L_d, cuyos coeficientes solo contienen Gcal:
    p_d = -2 N Gcal,   u_d = -(2/r)(N r' Gcal)',   k_d = -2 G4/N      (T3 de la Fase 2f).
Identidad puntual (exacta, fuera de la capa de masas):
    -L_rest(c = 0, d = N)|_estatico = -2 N^3 Gcal (r'/r - Phi')^2 + D_l[ ... ],
y, como la derivada total no contribuye (es impar en la garganta y se anula en el infinito),
    E[N] = int (P N'^2 + U N^2) dl = -2 int N^3 Gcal (r'/r - Phi')^2 dl = -2 int N^5 Gcal [(r/N)']^2/r^2 dl  < 0
si Gcal > 0 (condicion de que el sector impar no tenga fantasmas, T4).  Con K, P > 0 y K N^2 integrable, el principio
variacional da un taquion:  omega_0^2 <= E[N] / int K N^2 < 0.
Comprobaciones:
  1. al0 + Phi' al1 = 0 (exacta) y, por tanto, W(c = 0, d = N) = 0 identicamente.
  2. L_rest(c = 0, d) tiene la forma k_d d_t^2 - pt d'^2 - ut d^2 - mt d d' (sin d d''), y p_d = pt, u_d = ut - mt'/2
     coinciden con -2 N Gcal y -(2/r)(N r' Gcal)'.
  3. La identidad puntual con la derivada total explicita  D_l[ mt N^2/2 - 2 N^3 r' Gcal/r ].
Escribe resultados/teorema_N.json.
Uso: python teorema_N.py
"""
import os
import sys
import json
import time
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
from chorros import Phi, rr, PP, Ee, eP, Xs, D       # noqa: E402

T0 = time.time()
RES = {}


def log(m):
    print(f"[{time.time() - T0:6.1f} s] {m}", flush=True)


def S(n):
    return sp.Symbol(n, real=True)


mr = pickle.load(open(os.path.join(F2F, "resultados", "monopolo_reduccion.pkl"), "rb"))
Xb = {Xs: -PP[1] ** 2 / 2}
N = sp.sqrt(Ee)

# ------------------------------------------------------------------ 1. ligadura de momento con d = N, c = 0
al1 = sp.expand(mr["al1"].xreplace(Xb))
al0 = sp.expand(mr["al0"].xreplace(Xb))
assert sp.simplify(al0 + Phi[1] * al1) == 0
dN = [N, N * Phi[1], N * (Phi[2] + Phi[1] ** 2)]                     # N, N', N''
W = sp.simplify(al1 * dN[1] + al0 * dN[0])
assert W == 0
RES["ligadura"] = "al0 = -Phi' al1 exacta; W(c = 0, d = N) = 0"
log("1. al0 = -Phi' al1 (exacta): la ligadura de momento se cumple con d = N, c = 0")

# ------------------------------------------------------------------ 2. L_rest(c = 0, d): coeficientes
L = sp.expand(mr["L_rest"].xreplace(Xb))
ceros_c = {S(n): 0 for n in ("c_00", "c_01", "c_10", "c_20")}
Ld = sp.expand(L.xreplace(ceros_c))
d00, d01, d02, d10, d20 = (S(n) for n in ("d_00", "d_01", "d_02", "d_10", "d_20"))
mon = sp.Poly(Ld, d00, d01, d02, d10, d20)
monos = {m: sp.expand(c) for m, c in zip(mon.monoms(), mon.coeffs())}
log("2. monomios de L_d: " + ", ".join(str(m) for m in monos))
permitidos = {(2, 0, 0, 0, 0), (0, 2, 0, 0, 0), (1, 1, 0, 0, 0), (0, 0, 0, 2, 0)}
assert set(monos) <= permitidos, set(monos) - permitidos
kd = monos.get((0, 0, 0, 2, 0), 0)
pt = -monos.get((0, 2, 0, 0, 0), 0)
ut = -monos.get((2, 0, 0, 0, 0), 0)
mt = -monos.get((1, 1, 0, 0, 0), 0)
G4, G401, F400 = S("G4_00"), S("G4_01"), S("F4_00")
p1 = PP[1]
Gc = G4 + G401 * p1 ** 2 - F400 * p1 ** 4                 # Gcal = G4 - 2X G4_X - 4X^2 F4 con X = -p1^2/2
p_d = pt
u_d = sp.expand(ut - D(mt) / 2)
assert sp.simplify(kd + 2 * G4 / N) == 0
assert sp.simplify(p_d + 2 * N * Gc) == 0
assert sp.simplify(u_d + 2 * D(N * rr[1] * Gc) / rr[0]) == 0
RES["L_d"] = "k_d = -2 G4/N, p_d = -2 N Gcal, u_d = -(2/r)(N r' Gcal)' (exactas)"
log("   k_d = -2 G4/N, p_d = -2 N Gcal, u_d = ut - mt'/2 = -(2/r)(N r' Gcal)'  (exactas)")

# ------------------------------------------------------------------ 3. identidad puntual
menosL = sp.expand(-Ld.xreplace({d00: dN[0], d01: dN[1], d02: dN[2], d10: 0, d20: 0}))
objetivo = -2 * N ** 3 * Gc * (rr[1] / rr[0] - Phi[1]) ** 2
td = D(mt * N ** 2 / 2 - 2 * N ** 3 * rr[1] * Gc / rr[0])
dif = sp.simplify(sp.expand(menosL - objetivo - td))
assert dif == 0, dif
RES["identidad"] = "-L_rest(c=0, d=N) = -2 N^3 Gcal (r'/r - Phi')^2 + D_l[mt N^2/2 - 2 N^3 r' Gcal/r] (exacta)"
log("3. -L_rest(c = 0, d = N) = -2 N^3 Gcal (r'/r - Phi')^2 + D_l[ mt N^2/2 - 2 N^3 r' Gcal / r ]  (exacta)")

# forma de mt (para la discusion de los terminos de borde): debe ser par en la garganta (impar x impar)
RES["mt"] = str(sp.factor(mt))
log(f"   mt = {sp.factor(mt)}")
json.dump(RES, open(os.path.join(AQUI, "resultados", "teorema_N.json"), "w"), indent=1)
log("comprobado; resultados/teorema_N.json")
