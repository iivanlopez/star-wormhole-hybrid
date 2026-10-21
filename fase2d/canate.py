# -*- coding: utf-8 -*-
"""
Agujero de gusano liso de tipo Ellis en EsGB + Maxwell (Canate, Sultana, Kazanas, arXiv:1907.09463, PRD 100, 064007),
en las convenciones de fase2b (accion (1/16 pi) [R - (dphi)^2/2 + F(phi) G], longitud propia l, N = e^Phi):

  ds^2 = -dt^2 + dl^2 + (l^2 + Q^2) dOmega^2          (Ellis sin masa, ultraestatico: Phi = 0)
  phi(l) = (2 Qs/Q) arctan(l/Q)                       (su ec. 32 con rho = sqrt(r^2 - Q^2) = l; impar en l)
  F(phi) = f(phi) de su ec. (33)                       (par en phi; F'(0) = 0 en la garganta)
  campo electrico radial con carga conservada: 8 pi T^mu_nu = (Qe^2/r^4) diag(-1, -1, 1, 1), Qe^2 = Q^2 + Qs^2
  (normalizacion comprobada en verificar_canate.py: la ligadura ll en la garganta la fija).

Maxwell en simetria esferica: la carga se conserva y F_{tr} F^{tr} = -Qe^2/r^4 para cualquier metrica esferica, asi
que el tensor de energia (componentes mixtas) solo depende de r: las ecuaciones de fase2b reciben +-Qe^2/r^4 y las
perturbadas -+4 Qe^2 d/r^5, exactamente (no hay grados de libertad electromagneticos radiales).

Contenido:
  acoplamiento(Q, Qs)        phi -> (F, F', F'', F''') (SymPy, lambdify)
  segundas_Q, est_ligadura_Q fondo con Maxwell
  FondoCanate                fondo exacto para l de cualquier signo; .sol.sol(l), .Fd, .ecuaciones (gancho para
                             fase2c/estabilidad_lisa: derivadas, ligaduras y fondo_en con Maxwell)
  derivadas_Q, ligaduras_Q, fondo_en_Q, velocidades_Q
"""
import os
import sys
import math
from types import SimpleNamespace
import numpy as np
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(AQUI, "..", "fase2b"))
sys.path.append(os.path.join(AQUI, "..", "fase2c"))
import egb_F_gen as H                      # noqa: E402


# ------------------------------------------------------------------ acoplamiento de Canate et al. (ec. 33)
def _f_simbolica():
    x, Q, Qs = sp.symbols("phi Q Q_s", real=True)
    k = Q / (2 * Qs)
    c, s = sp.cos(k * x), sp.sin(k * x)
    f = (Q**2 + Qs**2) * (c - c**3 * sp.log(c**4) - (1 + 2 * c**2) * (Q / Qs) * x * s) / (12 * c**3)
    d = [f]
    for _ in range(3):
        d.append(sp.simplify(sp.diff(d[-1], x)))
    return x, Q, Qs, d


_X, _Q, _QS, _FD = _f_simbolica()
_FD_NUM = sp.lambdify((_X, _Q, _QS), _FD, modules="numpy")


def acoplamiento(Q, Qs):
    def Fd(phi):
        return tuple(float(v) for v in _FD_NUM(phi, Q, Qs))
    return Fd


def comprobar_derivadas_publicadas():
    """f' y f'' de la ec. (34) del articulo frente a las derivadas simbolicas de (33)."""
    x, Q, Qs = _X, _Q, _QS
    k = Q / (2 * Qs)
    f1p = -(Q**2 + Qs**2) * Q**2 * x / (8 * Qs**2) / sp.cos(k * x) ** 4
    f2p = -(Q**2 + Qs**2) * Q**2 / (8 * Qs**3) * (2 * Q * x * sp.sin(k * x) + Qs * sp.cos(k * x)) / sp.cos(k * x) ** 5
    return sp.simplify(_FD[1] - f1p), sp.simplify(_FD[2] - f2p)


# ------------------------------------------------------------------ ecuaciones con Maxwell
def segundas_Q(y, Fd, Qe2):
    """(r'', Phi'', phi'') y det M del fondo con el campo electrico (y = (r, r', Phi, Phi', phi, phi', ...))."""
    r, s, Phi, q, p, u = y[:6]
    F0, F1, F2, F3 = Fd(p)
    M = np.array(H.est_principal(F1, r, s, q, u))
    b = np.array(H.est_resto(F1, F2, r, s, q, u))
    b[0] += Qe2 / r**4                 # E^t_t - 8 pi T^t_t,  8 pi T^t_t = -Qe^2/r^4
    b[1] -= Qe2 / r**4                 # E^th_th - 8 pi T^th_th, 8 pi T^th_th = +Qe^2/r^4
    return np.linalg.solve(M, -b), np.linalg.det(M)


def est_ligadura_Q(y, Fd, Qe2):
    r, s, Phi, q, p, u = y[:6]
    return H.est_ligadura(Fd(p)[1], r, s, q, u) + Qe2 / r**4


def fondo_en_Q(f, l):
    y = f.sol.sol(l)
    x, _ = segundas_Q(y, f.Fd, f.Qe2)
    F0, F1, F2, F3 = f.Fd(y[4])
    return (F1, F2, F3), (y[0], y[1], x[0], y[3], x[1], y[5], x[2], math.exp(-2.0 * y[2]))


def derivadas_Q(f, l, y, w):
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = fondo_en_Q(f, l)
    a1, d0, d1, q0, q1 = y
    M = np.array(H.lin_principal(F1, r, r1, P1, p1))
    b = np.array(H.lin_resto(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1))
    em = 4.0 * f.Qe2 * d0 / r**5
    b = b + np.array([-em, em, 0.0])  # delta(8 pi T^t_t) = +4 Qe^2 d/r^5, delta(8 pi T^th_th) = -4 Qe^2 d/r^5
    d2, a2, q2 = np.linalg.solve(M, -b)
    return [a2, d1, d2, q1, q2]


def ligaduras_Q(f, l, y, w):
    (F1, F2, F3), bg = fondo_en_Q(f, l)
    lt, ll = H.lin_ligaduras(F1, F2, F3, w, *bg, *y)
    return [lt, ll - 4.0 * f.Qe2 * y[1] / bg[0] ** 5]   # E^l_l: delta(8 pi T^l_l) = +4 Qe^2 d/r^5


def velocidades_Q(f, l):
    """Como estabilidad_F.velocidades, con las derivadas segundas del fondo de las ecuaciones con Maxwell
    (Maxwell no aporta terminos en w ni derivadas segundas: solo cambia el fondo)."""
    import scipy.linalg as sl
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = fondo_en_Q(f, l)
    M = np.array(H.lin_principal(F1, r, r1, P1, p1))
    N = np.zeros((3, 3))
    for j, (d0, q0) in ((0, (1.0, 0.0)), (2, (0.0, 1.0))):
        N[:, j] = np.array(H.lin_resto(F1, F2, F3, 1.0, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0)) \
            - np.array(H.lin_resto(F1, F2, F3, 0.0, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0))
    lam = sl.eig(M, N, right=False)
    return sorted(x.real for x in lam if np.isfinite(x) and abs(x) < 1e12 and abs(x.imag) < 1e-9 * max(1, abs(x)))


ECUACIONES_Q = SimpleNamespace(derivadas=derivadas_Q, ligaduras=ligaduras_Q, fondo_en=fondo_en_Q,
                               velocidades=velocidades_Q)


# ------------------------------------------------------------------ fondo exacto
class _Exacta:
    def __init__(self, Q, Qs):
        self.Q, self.Qs = Q, Qs

    def sol(self, l):
        Q, Qs = self.Q, self.Qs
        r = math.sqrt(l * l + Q * Q)
        return np.array([r, l / r, 0.0, 0.0, 2 * Qs / Q * math.atan(l / Q), 2 * Qs / r**2, math.asinh(l / Q)])


class FondoCanate:
    """Fondo exacto (l de cualquier signo).  Atributos usados por fase2c/estabilidad_lisa: sol.sol, Fd, ecuaciones,
    rpp0, phi1, Phi_inf; der/izq con Phi_inf y l_end para las colas."""

    def __init__(self, Qs, Q=1.0, l_max=1e6):
        self.Q, self.Qs = Q, Qs
        self.Qe2 = Q * Q + Qs * Qs
        self.Fd = acoplamiento(Q, Qs)
        self.sol = _Exacta(Q, Qs)
        self.ecuaciones = ECUACIONES_Q
        self.r0 = Q
        self.a = 0.0
        self.rpp0 = 1.0 / Q
        self.phi1 = 2 * Qs / Q**2
        self.Phi_inf = 0.0
        self.der = SimpleNamespace(Phi_inf=0.0, l_end=l_max)
        self.izq = SimpleNamespace(Phi_inf=0.0, l_end=l_max)
        self.l_end = l_max

    def exactas_segundas(self, l):
        Q, Qs = self.Q, self.Qs
        r = math.sqrt(l * l + Q * Q)
        return np.array([Q * Q / r**3, 0.0, -4 * Qs * l / r**4])
