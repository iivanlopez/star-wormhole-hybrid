# -*- coding: utf-8 -*-
"""
Modos radiales pares de los agujeros de KKK con la garganta libre (delta r(0) != 0) y la capa de la garganta dinamica.

Capa (accion compatible con las condiciones de empalme de KKK, PRD 85 044007):
    L_s = (1/8 pi) [ -Lambda_s + kappa F(phi) Rbar ] - n m(phi),      F = e^{-gamma phi}
    tension -> Lambda_s = -2 Phi'_0;  gravedad inducida con dilaton -> kappa = 2 alpha lambda_0 = 4 alpha Phi'_0 (lambda_0 = lambda_1, KKK);
    polvo con masa m(phi): rho = 8 pi n m = 2 Phi'_0 + 8 alpha F0 (Phi'_0 + gamma phi'_0)/r0^2,  beta = m'/m = phi'_0/rho,
    mu = m''/m' libre (acoplamiento exponencial m ~ e^{beta phi}: mu = beta).  n r^2 = cte (sin flujo: simetria Z2).
Empalmes: parte delta(l) de las ecuaciones linealizadas (E^t_t, E^th_th, E_phi) = (8 pi S^t_t, 8 pi S^th_th, sigma_phi),
con saltos [d'] = 2 d'(0+), [a'] = 2 a'(0+), [q'] = 2 q'(0+) (modos pares) y los saltos del fondo [Phi'] = 2 Phi'_0, [phi'] = 2 phi'_0.
Ligaduras E^l_t = E^l_l = 0 en 0+.  Bianchi en la capa hace redundante una de las 5 condiciones (se comprueba).
Infinito: q ~ e^{-k l}/r con k^2 = -w e^{-2 Phi_inf}; D(w) = coeficiente de la solucion creciente.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import egb_lineal_gen as G
from kkk import FondoKKK, segundas
from perturb_kkk import fondo_en, derivadas, ligaduras

PI8 = 8.0 * math.pi


class CapaKKK:
    """zeta escala el termino de gravedad inducida respecto de la eleccion de KKK (lambda_0 = lambda_1): kappa = zeta 4 alpha Phi'_0.
    Las tres condiciones de empalme del fondo fijan entonces Lambda_s, rho y beta (la carga dilatonica restante la lleva el polvo)."""

    def __init__(self, f, mu=None, zeta=1.0):
        self.f = f
        y0 = f.y0()
        self.r0, self.P1, self.p1 = y0[0], y0[3], y0[5]
        al, ga = f.alpha, f.gamma
        r0 = self.r0
        self.F0 = math.exp(-ga * y0[4])
        a0 = al * ga * self.F0
        self.zeta = zeta
        self.kappa = zeta * 4.0 * al * self.P1
        self.Lam = -2.0 * self.P1                                                     # empalme thth
        self.rho = 8 * a0 * self.p1 / r0**2 + 2 * self.kappa * self.F0 / r0**2 + 2 * self.P1   # empalme tt
        sphi = 2 * self.p1 + 16 * a0 * self.P1 / r0**2                                 # empalme del dilaton
        self.beta = (sphi - 4 * self.kappa * ga * self.F0 / r0**2) / (2 * self.rho)
        self.mu = self.beta if mu is None else mu

    def fuentes(self, d0, q0, w):
        """(8 pi dS^t_t, 8 pi dS^th_th, d sigma_phi) lineales en (d0, q0)."""
        r0, ka, F, ga, rho, be = self.r0, self.kappa, self.F0, self.f.gamma, self.rho, self.beta
        Fp, Fpp = -ga * F, ga * ga * F
        Stt = rho * (2 * d0 / r0 - be * q0) + 2 * ka * (Fp * q0 - 2 * F * d0 / r0) / r0**2
        Sthth = -2 * ka * w * (F * d0 / r0 + Fp * q0)
        dR = -4 * d0 / r0**3 - 4 * w * d0 / r0
        sphi = -2 * ka * (Fpp * q0 * 2 / r0**2 + Fp * dR) + 2 * rho * be * (-2 * d0 / r0 + self.mu * q0)
        return Stt, Sthth, sphi

    def fondo_ok(self):
        """Comprueba los empalmes del fondo: parte delta de las ecuaciones = fuentes de la capa estatica."""
        f = self.f
        bg = fondo_en(f, 0.0)
        M = np.array(G.principal(f.alpha, f.gamma, *bg))
        salto = np.array([0.0, 2 * self.P1, 2 * self.p1])
        delta = M @ salto
        Stt = -(self.rho - 2 * self.kappa * self.F0 / self.r0**2 + self.Lam)
        Sthth = -self.Lam
        sphi = -2 * self.kappa * (-f.gamma * self.F0) * 2 / self.r0**2 + 2 * self.rho * self.beta
        return delta, np.array([Stt, Sthth, sphi])


def matriz_garganta(f, capa, w):
    """Filas: (E^l_t, E^l_l, empalme tt, empalme thth, empalme phi); columnas: u = (a', d, d', q, q') en 0+."""
    al, ga = f.alpha, f.gamma
    bg = list(fondo_en(f, 0.0))
    M = np.array(G.principal(al, ga, *bg))
    K = np.zeros((5, 5), dtype=complex if isinstance(w, complex) else float)
    for k in range(5):
        u = np.zeros(5)
        u[k] = 1.0
        a1, d0, d1, q0, q1 = u
        lig = G.ligaduras(al, ga, w, *bg, 0.0, a1, d0, d1, q0, q1)
        # parte delta del resto: lineal en P2 y p2
        bgP, bgp, bg0 = list(bg), list(bg), list(bg)
        bgP[4], bgP[7] = 1.0, 0.0
        bgp[4], bgp[7] = 0.0, 1.0
        bg0[4], bg0[7] = 0.0, 0.0
        R0 = np.array(G.resto(al, ga, w, *bg0, 0.0, a1, d0, d1, q0, q1))
        dP = np.array(G.resto(al, ga, w, *bgP, 0.0, a1, d0, d1, q0, q1)) - R0
        dp = np.array(G.resto(al, ga, w, *bgp, 0.0, a1, d0, d1, q0, q1)) - R0
        delta = M @ np.array([2 * d1, 2 * a1, 2 * q1]) + dP * 2 * capa.P1 + dp * 2 * capa.p1
        fu = np.array(capa.fuentes(d0, q0, w))
        K[:, k] = [lig[0], lig[1], *(delta - fu)]
    return K


def datos_garganta(f, capa, w):
    """Solucion con d(0) = 1 de (E^l_t, E^l_l, empalme thth, empalme phi); devuelve u y el residuo del empalme tt."""
    K = matriz_garganta(f, capa, w)
    filas = [0, 1, 3, 4]
    A = K[np.ix_(filas, [0, 2, 3, 4])]
    b = -K[filas, 1]
    x = np.linalg.solve(A, b)
    u = np.array([x[0], 1.0, x[1], x[2], x[3]])
    res_tt = K[2] @ u / (np.abs(K[2]) @ np.abs(u))
    return u, res_tt, K


def vector_garganta(f, capa, w, ref=None):
    """Vector nulo unitario de (E^l_t, E^l_l, empalme thth, empalme phi): no tiene polos (la normalizacion d(0) = 1
    los tiene donde la solucion admisible tiene d(0) = 0).  El signo se fija por continuidad con ref."""
    K = matriz_garganta(f, capa, w)
    v = np.linalg.svd(K[[0, 1, 3, 4]])[2][-1]
    if ref is not None and v @ ref < 0:
        v = -v
    res_tt = K[2] @ v / (np.abs(K[2]) @ np.abs(v))
    return v, res_tt


def barrido(f, capa, ws, **kw):
    """D(w) sobre una malla ordenada, con el signo del dato inicial seguido por continuidad."""
    ref, out = None, []
    for w in ws:
        v, res = vector_garganta(f, capa, w, ref)
        ref = v
        out.append((D(f, capa, w, u=v, **kw), res))
    return out


def D(f, capa, w, l_extra=30.0, rtol=1e-10, devolver=False, u=None):
    """Coeficiente de la solucion creciente de q en el infinito (w < 0, reloj de la garganta)."""
    if u is None:
        u, res_tt = vector_garganta(f, capa, w)
    else:
        res_tt = None
    kap = math.sqrt(max(-w * math.exp(-2 * f.Phi_inf), 1e-30))
    l1 = min(10.0 * f.r0 + l_extra / kap, 0.9 * f.l_end)
    if devolver:
        sol = solve_ivp(lambda l, y: derivadas(f, l, y, w), (0.0, l1), u, method="DOP853", rtol=rtol, atol=1e-14,
                        dense_output=True)
        y_fin = sol.y[:, -1]
    else:
        # renormalizacion por tramos (factor > 0: no cambia el signo de D) para |w| grande
        def grande(l, y):
            return np.max(np.abs(y)) - 1e80
        grande.terminal = True
        y_fin, l, escala = np.array(u, dtype=float), 0.0, 0.0
        y_fin = y_fin / np.max(np.abs(y_fin))
        while l < l1 * (1 - 1e-14):
            s_ = solve_ivp(lambda x, y: derivadas(f, x, y, w), (l, l1), y_fin, method="DOP853", rtol=rtol,
                           atol=1e-14 * np.max(np.abs(y_fin)), events=[grande])
            y_fin, l = s_.y[:, -1], s_.t[-1]
            y_fin = y_fin / np.max(np.abs(y_fin))
    a1, d0, d1, q0, q1 = y_fin
    r, s = f.sol.sol(l1)[:2]
    val = (q1 + (kap + s / r) * q0) * r * math.exp(-kap * l1)
    if devolver:
        return val, sol, res_tt
    return val
