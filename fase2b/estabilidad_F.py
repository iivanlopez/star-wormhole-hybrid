# -*- coding: utf-8 -*-
"""
Modos radiales pares de agujeros EsGB con acoplamiento general F(phi) (fondo_F.FondoF), con la garganta (o el ecuador)
libre y la capa de l = 0 dinamica.  Generaliza estabilidad_kkk.py:
  capa  L_s = (1/8 pi) [ -Lambda_s(phi) + kappa F(phi) Rbar ] - n m(phi),   kappa = zeta 4 Phi'_0
        (zeta = 1: lambda0 = lambda1, la eleccion de KKK y de Antoniou et al.);  F es el acoplamiento completo.
        Lambda_s(phi0) = -2 Phi'_0;  Lambda_s'(phi0) = dLam y Lambda_s''(phi0) = d2Lam (por defecto 0: tension constante).
  empalmes = parte delta(l) de (E^t_t, E^th_th, E_phi);  ligaduras E^l_t = E^l_l = 0 en 0+.
  infinito: q ~ e^{-k l}/r, k^2 = -w e^{-2 Phi_inf}.  w con el reloj de l = 0 (Phi(0) = 0).
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
import egb_F_gen as H
from fondo_F import segundas


def fondo_en(f, l):
    y = f.sol.sol(l)
    x, _ = segundas(y, f.Fd)
    F0, F1, F2, F3 = f.Fd(y[4])
    return (F1, F2, F3), (y[0], y[1], x[0], y[3], x[1], y[5], x[2], math.exp(-2.0 * y[2]))


def derivadas(f, l, y, w):
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = fondo_en(f, l)
    a1, d0, d1, q0, q1 = y
    M = np.array(H.lin_principal(F1, r, r1, P1, p1))
    b = -np.array(H.lin_resto(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1))
    d2, a2, q2 = np.linalg.solve(M, b)
    return [a2, d1, d2, q1, q2]


def ligaduras(f, l, y, w):
    (F1, F2, F3), bg = fondo_en(f, l)
    return H.lin_ligaduras(F1, F2, F3, w, *bg, *y)


class CapaF:
    def __init__(self, f, zeta=1.0, mu=None, dLam=0.0, d2Lam=0.0):
        self.f, self.zeta = f, zeta
        y0 = f.y0()
        self.r0, self.P1, self.p1 = y0[0], y0[3], y0[5]
        self.F0, self.F1, self.F2, _ = f.Fd(y0[4])
        r0 = self.r0
        self.kappa = zeta * 4.0 * self.P1
        self.Lam = -2.0 * self.P1
        self.dLam, self.d2Lam = dLam, d2Lam
        d = self.delta_fondo()
        # empalmes del fondo: d = (8 pi S^t_t, 8 pi S^th_th, sigma_phi) de la capa estatica; la tension aporta
        # -Lambda_s a 8 pi S^t_t y 8 pi S^th_th y +2 Lambda_s' a sigma_phi = -16 pi dL_s/dphi (el polvo compensa esta carga)
        self.rho = -d[0] + 2 * self.kappa * self.F0 / r0**2 - self.Lam
        self.beta = (d[2] + 4 * self.kappa * self.F1 / r0**2 - 2 * dLam) / (2 * self.rho)
        self.mu = self.beta if mu is None else mu
        assert abs(d[1] + self.Lam) < 1e-9 * max(1.0, abs(self.Lam))

    def delta_fondo(self):
        (F1, F2, F3), bg = fondo_en(self.f, 0.0)
        M = np.array(H.lin_principal(F1, bg[0], bg[1], bg[3], bg[5]))
        return M @ np.array([0.0, 2 * self.P1, 2 * self.p1])

    def fuentes(self, d0, q0, w):
        r0, ka, F0, F1, F2, rho, be = self.r0, self.kappa, self.F0, self.F1, self.F2, self.rho, self.beta
        Stt = rho * (2 * d0 / r0 - be * q0) + 2 * ka * (F1 * q0 - 2 * F0 * d0 / r0) / r0**2 - self.dLam * q0
        Sthth = -2 * ka * w * (F0 * d0 / r0 + F1 * q0) - self.dLam * q0
        dR = -4 * d0 / r0**3 - 4 * w * d0 / r0
        sphi = -2 * ka * (F2 * q0 * 2 / r0**2 + F1 * dR) + 2 * rho * be * (-2 * d0 / r0 + self.mu * q0) \
            + 2 * self.d2Lam * q0
        return Stt, Sthth, sphi


def matriz_garganta(f, capa, w):
    (F1, F2, F3), bg = fondo_en(f, 0.0)
    r, r1, r2, P1, P2, p1, p2, e2 = bg
    M = np.array(H.lin_principal(F1, r, r1, P1, p1))
    K = np.zeros((5, 5), dtype=complex if isinstance(w, complex) else float)
    for k in range(5):
        u = np.zeros(5)
        u[k] = 1.0
        a1, d0, d1, q0, q1 = u
        lig = H.lin_ligaduras(F1, F2, F3, w, *bg, a1, d0, d1, q0, q1)
        R0 = np.array(H.lin_resto(F1, F2, F3, w, r, r1, r2, P1, 0.0, p1, 0.0, e2, a1, d0, d1, q0, q1))
        dP = np.array(H.lin_resto(F1, F2, F3, w, r, r1, r2, P1, 1.0, p1, 0.0, e2, a1, d0, d1, q0, q1)) - R0
        dp = np.array(H.lin_resto(F1, F2, F3, w, r, r1, r2, P1, 0.0, p1, 1.0, e2, a1, d0, d1, q0, q1)) - R0
        delta = M @ np.array([2 * d1, 2 * a1, 2 * q1]) + dP * 2 * capa.P1 + dp * 2 * capa.p1
        K[:, k] = [lig[0], lig[1], *(delta - np.array(capa.fuentes(d0, q0, w)))]
    return K


def vector_garganta(f, capa, w, ref=None):
    K = matriz_garganta(f, capa, w)
    v = np.linalg.svd(K[[0, 1, 3, 4]])[2][-1]
    if ref is not None and np.real(np.vdot(ref, v)) < 0:
        v = -v
    res_tt = abs(K[2] @ v) / (np.abs(K[2]) @ np.abs(v))
    return v, res_tt


def vector_menores(f, capa, w):
    K = matriz_garganta(f, capa, w)[[0, 1, 3, 4]]
    return np.array([(-1) ** k * np.linalg.det(np.delete(K, k, axis=1)) for k in range(5)])


def _l1(f, kap, l_extra):
    return min(10.0 * f.r0 + l_extra / max(abs(kap), 1e-3), 0.9 * f.l_end)


def D(f, capa, w, u=None, l_extra=30.0, rtol=1e-10):
    """Coeficiente de la solucion creciente en el infinito (w real < 0), integracion renormalizada."""
    if u is None:
        u = vector_garganta(f, capa, w)[0]
    kap = math.sqrt(max(-w * math.exp(-2 * f.Phi_inf), 1e-30))
    l1 = _l1(f, kap, l_extra)

    def grande(x, y):
        return np.max(np.abs(y)) - 1e80
    grande.terminal = True
    y, l = np.array(u, dtype=float) / np.max(np.abs(u)), 0.0
    while l < l1 * (1 - 1e-14):
        s = solve_ivp(lambda x, yy: derivadas(f, x, yy, w), (l, l1), y, method="DOP853", rtol=rtol,
                      atol=1e-14 * np.max(np.abs(y)), events=[grande])
        y, l = s.y[:, -1], s.t[-1]
        y = y / np.max(np.abs(y))
    a1, d0, d1, q0, q1 = y
    r, s_ = f.sol.sol(l1)[:2]
    return (q1 + (kap + s_ / r) * q0) * r * math.exp(-kap * l1)


def Dc(f, capa, w, l_extra=30.0, rtol=1e-10):
    """D analitica (complejo, dato inicial por menores, rama principal de k)."""
    w = complex(w)
    u = vector_menores(f, capa, w).astype(complex)
    kap = np.sqrt(-w * math.exp(-2 * f.Phi_inf))
    if kap.real < 0:
        kap = -kap
    l1 = _l1(f, kap, l_extra)
    s = solve_ivp(lambda x, yy: derivadas(f, x, yy, w), (0.0, l1), u, method="DOP853", rtol=rtol, atol=1e-14)
    a1, d0, d1, q0, q1 = s.y[:, -1]
    r, s_ = f.sol.sol(l1)[:2]
    return (q1 + (kap + s_ / r) * q0) * r * np.exp(-kap * l1)


def raices(f, capa, ws, refinar=True):
    """Raices reales de D en la malla ordenada ws (signo del dato inicial seguido por continuidad)."""
    from scipy.optimize import brentq
    ref, vals, res = None, [], 0.0
    for w in ws:
        v, r_ = vector_garganta(f, capa, w, ref)
        ref = v
        vals.append(D(f, capa, w, u=v))
        res = max(res, r_)
    out = []
    for i in range(len(ws) - 1):
        if vals[i] * vals[i + 1] < 0:
            if not refinar:
                out.append(0.5 * (ws[i] + ws[i + 1]))
                continue
            r0 = vector_garganta(f, capa, ws[i])[0]
            try:
                out.append(brentq(lambda x: D(f, capa, x, u=vector_garganta(f, capa, x, r0)[0]), ws[i], ws[i + 1],
                                  xtol=1e-9, rtol=1e-10))
            except ValueError:
                out.append(float("nan"))
    return out, res


def secante(f, capa, w0, it=40):
    w0 = complex(w0)
    w1 = w0 * (1 + 1e-3) + 1e-2j
    D0, D1 = Dc(f, capa, w0), Dc(f, capa, w1)
    paso = 0.0
    for _ in range(it):
        if D1 == D0:
            break
        paso = -D1 * (w1 - w0) / (D1 - D0)
        if abs(paso) > 0.3 * abs(w1):
            paso *= 0.3 * abs(w1) / abs(paso)
        w0, D0 = w1, D1
        w1 = w1 + paso
        D1 = Dc(f, capa, w1)
        if abs(paso) < 1e-10 * abs(w1):
            break
    return w1, abs(paso), abs(D1)


def velocidades(f, l):
    """Cuadrados de las velocidades caracteristicas locales (limite de alta frecuencia): autovalores finitos lambda = w/k^2
    del problema M x = lambda N x, con M la parte principal en l y N los coeficientes de w (d, q).  lambda > 0 en todo
    punto <=> sistema hiperbolico; lambda < 0 en algun punto <=> problema de valores iniciales mal planteado."""
    import scipy.linalg as sl
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = fondo_en(f, l)
    M = np.array(H.lin_principal(F1, r, r1, P1, p1))
    N = np.zeros((3, 3))
    base = np.array(H.lin_resto(F1, F2, F3, 0.0, r, r1, r2, P1, P2, p1, p2, e2, 0, 0, 0, 0, 0))
    for j, (d0, q0) in ((0, (1.0, 0.0)), (2, (0.0, 1.0))):
        N[:, j] = np.array(H.lin_resto(F1, F2, F3, 1.0, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0)) \
            - np.array(H.lin_resto(F1, F2, F3, 0.0, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0))
    lam = sl.eig(M, N, right=False)
    # el lapso no tiene termino en w: su autovalor es infinito y a veces sale numericamente como +-1e19
    return sorted(x.real for x in lam if np.isfinite(x) and abs(x) < 1e12 and abs(x.imag) < 1e-9 * max(1, abs(x)))


def raices_menores(f, capa, ws, refinar=True):
    """Raices reales de D en la malla ws con el dato inicial por menores 4x4 (analitico, sin ambiguedad de signo) e
    integracion renormalizada (factores > 0).  Robusto frente a giros rapidos del vector nulo, que hacen fallar el
    seguimiento de signo de raices()."""
    from scipy.optimize import brentq

    def g(w):
        return D(f, capa, w, u=np.real(vector_menores(f, capa, float(w))))
    vals = [g(w) for w in ws]
    out = []
    for i in range(len(ws) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            if not refinar:
                out.append(0.5 * (ws[i] + ws[i + 1]))
                continue
            try:
                out.append(brentq(g, ws[i], ws[i + 1], xtol=1e-9, rtol=1e-10))
            except ValueError:
                out.append(float("nan"))
    return out
