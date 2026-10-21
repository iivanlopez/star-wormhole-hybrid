# -*- coding: utf-8 -*-
"""
Perturbaciones radiales de los agujeros de KKK en longitud propia (g_ll = 1, shift nulo):
    A = Phi + a(l) e^{i omega t},   R = r + d(l) e^{i omega t},   phi = p + q(l) e^{i omega t},   w = omega^2
(w medido con el tiempo t normalizado con Phi(0) = 0: reloj de la garganta).
Sistema en el bulto (l > 0):  y = (a', d, d', q, q')  ->  (d'', a'', q'') del sistema principal 3x3 (tt, thth, dilaton);
ligaduras E^l_t y E^l_l (sin derivadas segundas) que se imponen en l = 0+ y se propagan por Bianchi.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from kkk import segundas
import egb_lineal_gen as G


def fondo_en(f, l):
    """bg = (r, r1, r2, P1, P2, p, p1, p2, e2) del fondo f (FondoKKK resuelto) en l."""
    y = f.sol.sol(l)
    r2, P2, p2, _ = segundas(y, f.alpha, f.gamma)
    return (y[0], y[1], r2, y[3], P2, y[4], y[5], p2, math.exp(-2.0 * y[2]))


def derivadas(f, l, y, w):
    bg = fondo_en(f, l)
    a1, d0, d1, q0, q1 = y
    M = np.array(G.principal(f.alpha, f.gamma, *bg))
    b = -np.array(G.resto(f.alpha, f.gamma, w, *bg, 0.0, a1, d0, d1, q0, q1))
    d2, a2, q2 = np.linalg.solve(M, b)
    return [a2, d1, d2, q1, q2]


def ligaduras(f, l, y, w):
    a1, d0, d1, q0, q1 = y
    return G.ligaduras(f.alpha, f.gamma, w, *fondo_en(f, l), 0.0, a1, d0, d1, q0, q1)


def integrar(f, y0, w, l1, rtol=1e-10, densa=False):
    return solve_ivp(lambda l, y: derivadas(f, l, y, w), (0.0, l1), y0, method="DOP853", rtol=rtol, atol=1e-14,
                     dense_output=densa)
