# -*- coding: utf-8 -*-
"""
Agujeros de gusano de Kanti-Kleihaus-Kunz (EdGB, PRL 107 271101; PRD 85 044007) en longitud propia:

    ds^2 = -e^{2 Phi(l)} dt^2 + dl^2 + r(l)^2 dOmega^2,    S = (1/16pi) int sqrt(-g) [R - (dphi)^2/2 + alpha e^{-gamma phi} G]

Ecuaciones estaticas (derivadas y verificadas con egb_simbolico.py / verificar_egb.py), con a = alpha gamma e^{-gamma phi}:
  tt :  (2r + 8a phi' r') r'' + 4a (r'^2 - 1) phi''                    = -[r^2 phi'^2/4 + r'^2 - 1 + 4a gamma phi'^2 (1 - r'^2)]
  thth: (1 + 4a Phi' phi') r'' + 4a Phi' r' phi'' + (r + 4a phi' r') Phi'' = -[r Phi'^2 + r phi'^2/4 + Phi' r' - 4a gamma Phi' phi'^2 r' + 4a Phi'^2 phi' r']
  dil:  -16a Phi' r' r'' + r^2 phi'' + 8a (1 - r'^2) Phi''               = -[8a Phi'^2 (1 - r'^2) + r^2 Phi' phi' + 2 r r' phi']
  ll (ligadura): (r'^2 - 1)/r^2 + 2 Phi' r'/r - phi'^2/4 + 4a Phi' phi' (3r'^2 - 1)/r^2 = 0

Garganta (l = 0+, lado derecho de la solucion simetrica): r = r0, r' = 0, Phi = 0, phi = phi0, phi' = phi1 < 0 libre;
la ligadura da Phi'(0) = -(1 + r0^2 phi1^2/4) / (4 a0 phi1).  La solucion Z2 tiene un salto en phi' y Phi' en l = 0,
sostenido por una capa de materia (KKK, sec. III).
Parametros de KKK (phi_inf = 0): f0 = 1/(r0 r''(0)) y alpha/r0^2; por la simetria phi -> phi + c, l -> e^{-gamma c/2} l
basta integrar con r0 = 1, phi0 = 0 y despues alpha_KKK = alpha e^{-gamma phi_inf}.
"""
import math
from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp


def segundas(y, al, ga):
    """Resuelve el sistema lineal 3x3 para (r'', Phi'', phi''). Devuelve (r'', Phi'', phi'', det)."""
    r, s, Phi, q, p, u = y              # s = r', q = Phi', u = phi'
    a = al * ga * math.exp(-ga * p)
    M = np.array([[2 * r + 8 * a * u * s, 0.0, 4 * a * (s * s - 1)],
                  [1 + 4 * a * q * u, r + 4 * a * u * s, 4 * a * q * s],
                  [-16 * a * q * s, 8 * a * (1 - s * s), r * r]])
    b = -np.array([r * r * u * u / 4 + s * s - 1 + 4 * a * ga * u * u * (1 - s * s),
                   r * q * q + r * u * u / 4 + q * s - 4 * a * ga * q * u * u * s + 4 * a * q * q * u * s,
                   8 * a * q * q * (1 - s * s) + r * r * q * u + 2 * r * s * u])
    det = np.linalg.det(M)
    x = np.linalg.solve(M, b)
    return x[0], x[1], x[2], det


def ligadura(y, al, ga):
    r, s, Phi, q, p, u = y
    a = al * ga * math.exp(-ga * p)
    return (s * s - 1) / r**2 + 2 * q * s / r - u * u / 4 + 4 * a * q * u * (3 * s * s - 1) / r**2


@dataclass
class FondoKKK:
    alpha: float                # acoplamiento (con r0 = 1, phi0 = 0)
    phi1: float                 # phi'(0+) < 0
    gamma: float = 1.0
    r0: float = 1.0
    l_max: float = 2.0e4
    rtol: float = 1e-11
    info: dict = field(default_factory=dict)

    def y0(self):
        a0 = self.alpha * self.gamma
        q0 = -(1 + self.r0**2 * self.phi1**2 / 4) / (4 * a0 * self.phi1)
        return np.array([self.r0, 0.0, 0.0, q0, 0.0, self.phi1])

    def rhs(self, l, y):
        rpp, Ppp, ppp, _ = segundas(y, self.alpha, self.gamma)
        return [y[1], rpp, y[3], Ppp, y[5], ppp]

    def resolver(self):
        al, ga = self.alpha, self.gamma
        y0 = self.y0()
        rpp0, Ppp0, ppp0, det0 = segundas(y0, al, ga)
        self.rpp0 = rpp0
        self.f0 = 1.0 / (self.r0 * rpp0) if rpp0 > 0 else float("inf")

        def ev_det(l, y):
            return segundas(y, al, ga)[3] / (abs(det0) + 1e-300)
        ev_det.terminal = True

        def ev_rmin(l, y):                     # r' < 0 lejos de la garganta: segundo cuello / ecuador
            return y[1] + 1e-12 if l > 1e-6 else 1.0
        ev_rmin.terminal = True
        sol = solve_ivp(self.rhs, (0.0, self.l_max), y0, method="DOP853", rtol=self.rtol, atol=1e-14,
                        events=[ev_det, ev_rmin], dense_output=True)
        self.sol = sol
        y = sol.y[:, -1]
        self.l_end = sol.t[-1]
        self.regular = sol.status == 0 and sol.t[-1] >= self.l_max * (1 - 1e-12)
        r, s, Phi, q, p, u = y
        # asintotica: r' -> 1 - M/r ..., phi -> phi_inf - D/r, Phi -> Phi_inf - M/r
        self.M = r * (1 - s * s) / 2                       # masa de Misner-Sharp (unidades r0 = 1, phi0 = 0)
        self.D = r * r * u                                 # carga escalar: phi ~ phi_inf - D/r  ->  phi' ~ D/r^2
        self.phi_inf = p + self.D / r
        self.Phi_inf = Phi + q * r                          # Phi ~ Phi_inf - M/r -> Phi' ~ M/r^2
        self.ligadura_max = max(abs(ligadura(sol.y[:, i], al, ga)) for i in range(0, sol.y.shape[1], max(1, sol.y.shape[1] // 200)))
        # parametros de KKK (phi_inf = 0, r0 reescalado)
        lam = math.exp(ga * self.phi_inf / 2)
        self.alpha_kkk = al * math.exp(-ga * self.phi_inf)     # alpha / r0'^2 con r0' = lam r0
        self.M_kkk = self.M * lam / (lam * self.r0)            # M / r0' (adimensional)
        self.D_kkk = self.D * lam / (lam * self.r0)
        return self

    def capa(self):
        """Capa de polvo en la garganta (KKK, PRD 85 044007): rho r0 y lambda0 r0 = lambda1 r0 = 2 Phi'_0 r0.
        rho r0 es invariante bajo phi -> phi + c, l -> e^{-gamma c/2} l, asi que se calcula con r0 = 1, phi0 = 0:
            rho = 2 Phi'_0 + 8 (alpha e^{-gamma phi0}/r0^2) (Phi'_0 + gamma phi'_0)."""
        y0 = self.y0()
        q0, u0 = y0[3], y0[5]
        return 2 * q0 + 8 * self.alpha * (q0 + self.gamma * u0), 2 * q0
