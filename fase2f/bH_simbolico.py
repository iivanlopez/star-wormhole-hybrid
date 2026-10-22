# -*- coding: utf-8 -*-
"""
Teoria cuartica mas alla de Horndeski (GLPV) en simetria esferica: lagrangiana reducida (SymPy).

  S = (1/16 pi) int d^4x sqrt(-g) L,
  L = G2(phi, X) - G3(phi, X) box phi + G4(phi, X) R + G4_X [(box phi)^2 - phi_mn phi^mn] + F4(phi, X) E4,
  E4 = eps^{mnr}_s eps^{m'n'r's} phi_m phi_m' phi_nn' phi_rr'
     = 2 X [(box phi)^2 - phi_mn^2] + 2 box phi phi^m phi_mn phi^n - 2 phi^m phi_mn phi^nl phi_l,
  X = -(1/2) g^{mn} phi_m phi_n,  signatura (-,+,+,+).
RG con escalar canonico: G4 = 1, G2 = X - V, G3 = F4 = 0 (las convenciones de las Fases 2b-2e).

Metrica esferica general (sector l = 0), funciones de (t, l):
  ds^2 = -A dt^2 + 2 b dt dl + B dl^2 + C dOmega^2,   phi(t, l).
La lagrangiana reducida es Lred = sqrt(A B + b^2) C L (la integral angular da un factor constante); por criticidad
simetrica (grupo compacto SO(3)) sus ecuaciones de Euler-Lagrange son las ecuaciones de campo del sector l = 0.

Las funciones de la teoria G2, G3, G4, F4 son funciones generales de (phi, X): generar_bH.py las evalua en X como
simbolo y recupera la dependencia en X con la regla de la cadena; en el fondo pasan a simbolos G_ij (i derivadas en
phi, j en X) evaluados en (phi_b, X_b).
"""
import itertools
import sympy as sp

t, l, th, vp = sp.symbols("t l theta varphi", real=True)
COORD = [t, l, th, vp]


def signo_perm(p):
    p = list(p)
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def geometria(g):
    """Inversa, Christoffel Gam[a][b][c] = Gamma^a_bc y escalar de Ricci para una metrica 4x4 general."""
    n = 4
    gi = sp.simplify(g.inv())
    dg = [[[sp.diff(g[a, b], COORD[c]) for c in range(n)] for b in range(n)] for a in range(n)]
    Gam = [[[sp.simplify(sum(gi[a, d] * (dg[d][b][c] + dg[d][c][b] - dg[b][c][d]) for d in range(n)) / 2)
             for c in range(n)] for b in range(n)] for a in range(n)]
    # Ricci R_bd = d_a Gam^a_bd - d_d Gam^a_ba + Gam^a_ae Gam^e_bd - Gam^a_de Gam^e_ba
    Ric = sp.zeros(n, n)
    for b_ in range(n):
        for d in range(b_, n):
            v = sum(sp.diff(Gam[a][b_][d], COORD[a]) - sp.diff(Gam[a][b_][a], COORD[d]) for a in range(n))
            v += sum(Gam[a][a][e] * Gam[e][b_][d] - Gam[a][d][e] * Gam[e][b_][a] for a in range(n) for e in range(n))
            Ric[b_, d] = Ric[d, b_] = v
    Rs = sum(gi[a, b_] * Ric[a, b_] for a in range(n) for b_ in range(n))
    return gi, Gam, Rs


def escalares_phi(g, gi, Gam, phi):
    """phi_m, phi_mn y los invariantes que usa la lagrangiana."""
    n = 4
    d1 = [sp.diff(phi, x) for x in COORD]
    d2 = sp.Matrix(n, n, lambda m, k: sp.diff(phi, COORD[m], COORD[k]) - sum(Gam[s][m][k] * d1[s] for s in range(n)))
    up = [sum(gi[m, k] * d1[k] for k in range(n)) for m in range(n)]                  # phi^m
    S = sp.Matrix(n, n, lambda a, k: sum(gi[a, m] * d2[m, k] for m in range(n)))       # phi^a_k
    X = -sum(up[m] * d1[m] for m in range(n)) / 2
    box = sum(S[a, a] for a in range(n))
    trS2 = sum(S[a, k] * S[k, a] for a in range(n) for k in range(n))                  # phi_mn phi^mn
    pSp = sum(up[m] * d2[m, k] * up[k] for m in range(n) for k in range(n))            # phi^m phi_mn phi^n
    pS2p = sum(up[m] * d2[m, k] * S[k, j] * up[j] for m in range(n) for k in range(n) for j in range(n))
    E4 = 2 * X * (box**2 - trS2) + 2 * box * pSp - 2 * pS2p
    return dict(d1=d1, d2=d2, X=X, box=box, trS2=trS2, E4=E4)


def E4_levi_civita(g, d1, d2):
    """E4 = eps^{mnr a} eps^{m'n'r' s} g_as phi_m phi_m' phi_nn' phi_rr' con la suma explicita (1/(-g) sgn sgn)."""
    detg = g.det()
    tot = 0
    perms = list(itertools.permutations(range(4)))
    for p1 in perms:
        m, n_, r_, a = p1
        if d1[m] == 0:
            continue
        s1 = signo_perm(p1)
        for p2 in perms:
            m2, n2, r2, s_ = p2
            if d1[m2] == 0 or g[a, s_] == 0 or d2[n_, n2] == 0 or d2[r_, r2] == 0:
                continue
            tot += s1 * signo_perm(p2) * g[a, s_] * d1[m] * d1[m2] * d2[n_, n2] * d2[r_, r2]
    return tot / (-detg)


def metrica_esferica(A, b, B, C):
    return sp.Matrix([[-A, b, 0, 0], [b, B, 0, 0], [0, 0, C, 0], [0, 0, 0, C * sp.sin(th) ** 2]])


def lagrangiana(g, gi, Gam, Rs, phi, teoria):
    """Densidad sqrt(-g) L / sin(theta) para la teoria dada: teoria(phi, X) -> (G2, G3, G4, G4X, F4)."""
    e = escalares_phi(g, gi, Gam, phi)
    G2, G3, G4, G4X, F4 = teoria(phi, e["X"])
    L = G2 - G3 * e["box"] + G4 * Rs + G4X * (e["box"] ** 2 - e["trS2"]) + F4 * e["E4"]
    # sqrt(-g)/sin(theta) = sqrt(A B + b^2) C para la metrica esferica (C = r^2 > 0: sin valor absoluto)
    sqrtg = sp.sqrt(g[0, 1] ** 2 - g[0, 0] * g[1, 1]) * g[2, 2]
    return sqrtg * L, e
