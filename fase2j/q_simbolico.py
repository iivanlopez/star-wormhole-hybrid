# -*- coding: utf-8 -*-
"""
Fase 2j: teoria quintica (Horndeski y mas alla, GLPV) en simetria esferica: lagrangiana reducida (SymPy).

  S = (1/16 pi) int d^4x sqrt(-g) L,
  L = G2 - G3 box phi + G4 R + G4_X [(box phi)^2 - phi_mn^2] + F4 E4
      + G5 G_mn phi^mn - (G5_X/6) [(box phi)^3 - 3 box phi phi_mn^2 + 2 phi_mn^3] + F5 E5,
  E4 = eps^{mnr}_s eps^{m'n'r's} phi_m phi_m' phi_nn' phi_rr',   E5 = eps^{mnrs} eps^{m'n'r's'} phi_m phi_m' phi_nn' phi_rr' phi_ss',
  X = -(1/2) g^{mn} phi_m phi_n,  signatura (-,+,+,+)  (las convenciones de la Fase 2f, fase2f/bH_simbolico.py).
Horndeski: F4 = F5 = 0.  La teoria cuartica de la Fase 2f es G5 = F5 = 0.  Un acoplo f(phi) G_GB (Gauss-Bonnet) equivale
a G5 = -4 f' ln|X| (mas terminos en G2, G3, G4 si f'' != 0); la teoria de Lu-Pang (4D EGB) es
  G2 = 8 alfa X^2,  G3 = -8 alfa X,  G4 = 1 + 4 alfa X,  G5 = -4 alfa ln|X|.
Metrica esferica general (sector l = 0), funciones de (t, l):  ds^2 = -A dt^2 + 2 b dt dl + B dl^2 + C dOmega^2,  phi(t, l).
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


def metrica_esferica(A, b, B, C):
    return sp.Matrix([[-A, b, 0, 0], [b, B, 0, 0], [0, 0, C, 0], [0, 0, 0, C * sp.sin(th) ** 2]])


def geometria(g):
    """Inversa, Christoffel Gam[a][b][c] = Gamma^a_bc, tensor de Ricci R_bd y escalar de Ricci."""
    n = 4
    gi = sp.simplify(g.inv())
    dg = [[[sp.diff(g[a, b], COORD[c]) for c in range(n)] for b in range(n)] for a in range(n)]
    Gam = [[[sp.simplify(sum(gi[a, d] * (dg[d][b][c] + dg[d][c][b] - dg[b][c][d]) for d in range(n)) / 2)
             for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.zeros(n, n)
    for b_ in range(n):
        for d in range(b_, n):
            v = sum(sp.diff(Gam[a][b_][d], COORD[a]) - sp.diff(Gam[a][b_][a], COORD[d]) for a in range(n))
            v += sum(Gam[a][a][e] * Gam[e][b_][d] - Gam[a][d][e] * Gam[e][b_][a] for a in range(n) for e in range(n))
            Ric[b_, d] = Ric[d, b_] = v
    Rs = sum(gi[a, b_] * Ric[a, b_] for a in range(n) for b_ in range(n))
    return gi, Gam, Ric, Rs


def escalares_phi(g, gi, Gam, Ric, Rs, phi):
    """phi_m, phi_mn y los invariantes de la lagrangiana quintica."""
    n = 4
    d1 = [sp.diff(phi, x) for x in COORD]
    d2 = sp.Matrix(n, n, lambda m, k: sp.diff(phi, COORD[m], COORD[k]) - sum(Gam[s][m][k] * d1[s] for s in range(n)))
    up = [sum(gi[m, k] * d1[k] for k in range(n)) for m in range(n)]                  # phi^m
    S = sp.Matrix(n, n, lambda a, k: sum(gi[a, m] * d2[m, k] for m in range(n)))       # phi^a_k
    X = -sum(up[m] * d1[m] for m in range(n)) / 2
    box = sum(S[a, a] for a in range(n))
    S2 = S * S
    trS2 = sum(S2[a, a] for a in range(n))
    trS3 = sum((S2 * S)[a, a] for a in range(n))
    pSp = sum(up[m] * d2[m, k] * up[k] for m in range(n) for k in range(n))            # phi^m phi_mn phi^n
    pS2p = sum(up[m] * d2[m, k] * S[k, j] * up[j] for m in range(n) for k in range(n) for j in range(n))
    E4 = 2 * X * (box ** 2 - trS2) + 2 * box * pSp - 2 * pS2p
    # G^a_k phi^k_a con G^a_k = R^a_k - delta R/2
    Rup = sp.Matrix(n, n, lambda a, k: sum(gi[a, m] * Ric[m, k] for m in range(n)))
    Gphi = sum(Rup[a, k] * S[k, a] for a in range(n) for k in range(n)) - Rs * box / 2
    return dict(d1=d1, d2=d2, X=X, box=box, trS2=trS2, trS3=trS3, E4=E4, Gphi=Gphi, pSp=pSp, pS2p=pS2p, S=S, up=up)


def E5_levi_civita(g, d1, d2):
    """E5 = eps^{mnrs} eps^{m'n'r's'} phi_m phi_m' phi_nn' phi_rr' phi_ss' (suma explicita, 1/(-g) sgn sgn)."""
    detg = g.det()
    tot = 0
    perms = list(itertools.permutations(range(4)))
    for p1 in perms:
        m, n_, r_, s_ = p1
        if d1[m] == 0:
            continue
        s1 = signo_perm(p1)
        for p2 in perms:
            m2, n2, r2, s2 = p2
            if d1[m2] == 0 or d2[n_, n2] == 0 or d2[r_, r2] == 0 or d2[s_, s2] == 0:
                continue
            tot += s1 * signo_perm(p2) * d1[m] * d1[m2] * d2[n_, n2] * d2[r_, r2] * d2[s_, s2]
    return tot / (-detg)


def lagrangiana(g, gi, Gam, Ric, Rs, phi, teoria, con_F5=False):
    """Densidad sqrt(-g) L / sin(theta): teoria(phi, X) -> (G2, G3, G4, G4X, F4, G5, G5X, F5)."""
    e = escalares_phi(g, gi, Gam, Ric, Rs, phi)
    G2, G3, G4, G4X, F4, G5, G5X, F5 = teoria(phi, e["X"])
    L = (G2 - G3 * e["box"] + G4 * Rs + G4X * (e["box"] ** 2 - e["trS2"]) + F4 * e["E4"]
         + G5 * e["Gphi"] - G5X / 6 * (e["box"] ** 3 - 3 * e["box"] * e["trS2"] + 2 * e["trS3"]))
    if con_F5:
        L += F5 * E5_levi_civita(g, e["d1"], e["d2"])
    sqrtg = sp.sqrt(g[0, 1] ** 2 - g[0, 0] * g[1, 1]) * g[2, 2]
    return sqrtg * L, e
