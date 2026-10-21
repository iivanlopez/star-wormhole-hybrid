# -*- coding: utf-8 -*-
"""
Ecuaciones de campo de Einstein-dilaton-Gauss-Bonnet (EdGB) para metricas esfericas dependientes de (t, l):

    S = (1/16 pi) int sqrt(-g) [ R - (1/2) (d phi)^2 + alpha F(phi) G ],     F = e^{-gamma phi}   (KKK 2011)
    G = R^2 - 4 R_mn R^mn + R_mnrs R^mnrs

    ds^2 = -e^{2A} dt^2 + e^{2B} dl^2 + R^2 dOmega^2,      A, B, R, phi funciones de (t, l)

Ecuaciones:
    E_mn = G_mn + Gam_mn - T_mn = 0,      T_mn = (1/2) d_m phi d_n phi - (1/4) g_mn (d phi)^2
    E_phi = box phi + alpha F'(phi) G = 0
con el tensor de Gauss-Bonnet (Doneva-Yazadjiev 2018, ec. 3), Psi_m = alpha F'(phi) d_m phi:
    Gam_mn = -R (nabla_m Psi_n + nabla_n Psi_m) - 4 (nabla.Psi) (R_mn - R g_mn / 2)
             + 4 R_ma nabla^a Psi_n + 4 R_na nabla^a Psi_m - 4 g_mn R^ab nabla_a Psi_b + 4 R^b_{man} nabla^a Psi_b
Identidad (invariancia difeomorfica) que fija signos y factores:  nabla^m E_mn = -(1/2) E_phi d_n phi.
"""
import sympy as sp


def tensores(g, X):
    """Christoffel, Riemann R^a_{bcd}, Ricci R_bd, escalar, para una metrica diagonal g en coordenadas X."""
    n = len(X)
    gi = sp.diag(*[1 / g[i, i] for i in range(n)])
    Gam = [[[sum(gi[a, s] * (sp.diff(g[s, b], X[c]) + sp.diff(g[s, c], X[b]) - sp.diff(g[b, c], X[s]))
                 for s in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]
    Riem = [[[[sp.diff(Gam[a][d][b], X[c]) - sp.diff(Gam[a][c][b], X[d])
               + sum(Gam[a][c][e] * Gam[e][d][b] - Gam[a][d][e] * Gam[e][c][b] for e in range(n))
               for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = sp.Matrix(n, n, lambda b, d: sum(Riem[a][b][a][d] for a in range(n)))
    Rs = sum(gi[b, b] * Ric[b, b] for b in range(n))
    return gi, Gam, Riem, Ric, Rs


def gauss_bonnet(g, gi, Riem, Ric, Rs, n=4):
    K = sum(Riem[a][b][c][d] ** 2 * g[a, a] * gi[b, b] * gi[c, c] * gi[d, d]
            for a in range(n) for b in range(n) for c in range(n) for d in range(n))
    RicSq = sum(Ric[a, b] ** 2 * gi[a, a] * gi[b, b] for a in range(n) for b in range(n))
    return Rs**2 - 4 * RicSq + K


def ecuaciones(A, B, R, phi, X, alpha, gamma, F=None, Fp=None):
    """Devuelve (E, Ephi, GB): E[m][n] covariante (4x4), Ephi escalar y el invariante de Gauss-Bonnet.
    Acoplamiento alpha F(phi) G con F = e^{-gamma phi} por defecto.  Para un acoplamiento general se puede pasar
    Fp = F'(phi) como funcion arbitraria de (t, l) (campo auxiliar chi): las ecuaciones solo dependen de F' y de sus
    primeras derivadas, de modo que luego se sustituye chi = F'(phi), d chi = F''(phi) d phi."""
    t, l, th, ph = X
    n = 4
    g = sp.diag(-sp.exp(2 * A), sp.exp(2 * B), R**2, R**2 * sp.sin(th) ** 2)
    gi, Gam, Riem, Ric, Rs = tensores(g, X)
    GB = gauss_bonnet(g, gi, Riem, Ric, Rs)
    if Fp is not None:
        pass
    elif F is None:
        F = sp.exp(-gamma * phi)
        Fp = -gamma * sp.exp(-gamma * phi)
    else:
        Fp = sp.diff(F, phi)
    dphi = [sp.diff(phi, x) for x in X]
    Psi = [alpha * Fp * dphi[m] for m in range(n)]
    nPsi = sp.Matrix(n, n, lambda m, k: sp.diff(Psi[k], X[m]) - sum(Gam[s][m][k] * Psi[s] for s in range(n)))
    divPsi = sum(gi[m, m] * nPsi[m, m] for m in range(n))
    upPsi = sp.Matrix(n, n, lambda a, k: gi[a, a] * nPsi[a, k])          # nabla^a Psi_k
    Ricup = sp.Matrix(n, n, lambda a, b: gi[a, a] * gi[b, b] * Ric[a, b])  # R^ab
    s_RnP = sum(Ricup[a, b] * nPsi[a, b] for a in range(n) for b in range(n))
    Ein = Ric - g * Rs / 2
    Gm = sp.zeros(n, n)
    for m in range(n):
        for k in range(m, n):
            v = (-Rs * (nPsi[m, k] + nPsi[k, m]) - 4 * divPsi * Ein[m, k]
                 + 4 * sum(Ric[m, a] * upPsi[a, k] + Ric[k, a] * upPsi[a, m] for a in range(n))
                 - 4 * g[m, k] * s_RnP
                 + 4 * sum(Riem[b][m][a][k] * upPsi[a, b] for a in range(n) for b in range(n)))
            Gm[m, k] = Gm[k, m] = v
    dphi2 = sum(gi[m, m] * dphi[m] ** 2 for m in range(n))
    T = sp.Matrix(n, n, lambda m, k: dphi[m] * dphi[k] / 2 - g[m, k] * dphi2 / 4)
    E = Ein + Gm - T
    sqrtg = sp.exp(A + B) * R**2 * sp.sin(th)
    box = sum(sp.diff(sqrtg * gi[m, m] * dphi[m], X[m]) for m in range(n)) / sqrtg
    Ephi = box + alpha * Fp * GB
    return E, Ephi, GB, (g, gi, Gam)


def divergencia(E, g, gi, Gam, X, nu):
    """nabla^m E_{m nu} para una metrica diagonal."""
    n = len(X)
    return sum(gi[m, m] * (sp.diff(E[m, nu], X[m])
                           - sum(Gam[s][m][m] * E[s, nu] + Gam[s][m][nu] * E[m, s] for s in range(n)))
               for m in range(n))
