# -*- coding: utf-8 -*-
"""
Fase 2m: Proca generalizada cuartica en simetria esferica: lagrangiana reducida (SymPy).

  S = (1/16 pi) int d^4x sqrt(-g) L,
  L = G2(X, F, Y) + G3(X) nabla_m A^m + G4(X) R + G4_X [(nabla_m A^m)^2 - nabla_r A_s nabla^s A^r],
  X = -(1/2) A_m A^m,   F = -(1/4) F_mn F^mn,   Y = A^m A^n F_m^a F_na,   signatura (-,+,+,+)
(De Felice-Heisenberg-Kase-Mukohyama-Tsujikawa-Zhang; el termino c2 G4_X nabla_r A_s nabla^r A^s difiere del de arriba
en c2 G4_X (-2F) y se absorbe en G2).  Con A_m = d_m phi (F = 0) es la Horndeski cuartica con G2(X, 0, 0), G3^H = -G3,
G4^H = G4 (las convenciones de fase2f/bH_simbolico.py).

Metrica esferica general (sector l = 0) y vector, funciones de (t, l):
  ds^2 = -A dt^2 + 2 b dt dl + B dl^2 + C dOmega^2,   A_m = (A0, A1, 0, 0).
En l = 0, F_m^a F_na = -2F g^(2)_mn (solo F_tl != 0), asi que Y = 4 X F identicamente (se comprueba): G2(X, F, Y) entra
como G2(X, F, 4XF) y basta una funcion de (X, F).
"""
import os
import sys
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "fase2j"))
from q_simbolico import t, l, th, vp, COORD, metrica_esferica, geometria      # noqa: E402,F401


def escalares_A(g, gi, Gam, A):
    """X, F, Y, nabla.A y nabla_r A_s nabla^s A^r para el covector A (lista de 4 componentes)."""
    n = 4
    DA = sp.Matrix(n, n, lambda m, k: sp.diff(A[k], COORD[m]) - sum(Gam[s][m][k] * A[s] for s in range(n)))  # nabla_m A_k
    Fmn = sp.Matrix(n, n, lambda m, k: sp.diff(A[k], COORD[m]) - sp.diff(A[m], COORD[k]))
    up = [sum(gi[m, k] * A[k] for k in range(n)) for m in range(n)]
    X = -sum(up[m] * A[m] for m in range(n)) / 2
    Fup = gi * Fmn * gi
    F = -sum(Fmn[m, k] * Fup[m, k] for m in range(n) for k in range(n)) / 4
    FF = Fmn * gi * Fmn.T                                           # F_m^a F_na = F_mb g^ba F_na
    Y = sum(up[m] * up[k] * FF[m, k] for m in range(n) for k in range(n))
    div = sum(gi[m, k] * DA[m, k] for m in range(n) for k in range(n))
    M = DA * gi
    trT = sum((M * M)[a, a] for a in range(n))                      # nabla_r A_s g^sa nabla_a A_b g^br
    return dict(DA=DA, X=X, F=F, Y=Y, div=div, trT=trT)


def lagrangiana(g, gi, Gam, Rs, A, teoria):
    """Densidad sqrt(-g) L / sin(theta): teoria(X, F) -> (G2, G3, G4, G4X)."""
    e = escalares_A(g, gi, Gam, A)
    G2, G3, G4, G4X = teoria(e["X"], e["F"])
    L = G2 + G3 * e["div"] + G4 * Rs + G4X * (e["div"] ** 2 - e["trT"])
    sqrtg = sp.sqrt(g[0, 1] ** 2 - g[0, 0] * g[1, 1]) * g[2, 2]
    return sqrtg * L, e
