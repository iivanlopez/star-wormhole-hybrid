# -*- coding: utf-8 -*-
"""
Estrella de neutrones ordinaria (sin agujero de gusano), codigo independiente para validar:
  - TOV en coordenada areal con la EoS de DFKK (tabla "Without a wormhole" de R2).
  - Oscilaciones radiales en la forma de primer orden de Chanmugam (1977) / Kokkotas & Ruoff (2001),
    variables xi = Delta r / r y Delta P.  Criterio de control: omega_0^2 = 0 en el maximo de M(rho_c).
Unidades G = c = 1, km.
"""
import math
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from hibrido import PI4, PI8, MSUN_KM


class EstrellaNS:
    def __init__(self, eos, rho_bc, theta_stop=1e-7):
        self.eos, self.rho_bc = eos, rho_bc
        self.sig = eos.sigma(rho_bc)
        self.c1 = self.sig * (eos.n + 1.0)
        self.theta_stop = theta_stop
        self._tov()

    def theta(self, Phi):
        return ((1.0 + self.c1) * math.exp(-Phi) - 1.0) / self.c1

    def bg(self, r, y):
        m, Phi = y[0], y[1]
        th = self.theta(Phi)
        eps, P, _, _ = self.eos.termo(max(th, 0.0), self.rho_bc, self.sig)
        dPhi = (m + PI4 * r**3 * P) / (r * (r - 2.0 * m))
        return m, Phi, eps, P, dPhi

    def _tov(self):
        e_c, P_c, _, _ = self.eos.termo(1.0, self.rho_bc, self.sig)
        r1 = 1e-5
        y1 = [PI4 / 3.0 * e_c * r1**3, 2.0 * math.pi / 3.0 * (e_c + 3.0 * P_c) * r1**2, r1]

        def rhs(r, y):
            m, Phi, eps, P, dPhi = self.bg(r, y)
            return [PI4 * r * r * eps, dPhi, 1.0 / math.sqrt(1.0 - 2.0 * m / r)]

        def sup(r, y):
            return self.theta(y[1]) - self.theta_stop
        sup.terminal, sup.direction = True, -1
        sol = solve_ivp(rhs, (r1, 1e3), y1, method="DOP853", rtol=1e-12, atol=1e-16,
                        events=[sup], dense_output=True)
        self.sol, self.r1 = sol.sol, r1
        self.R = sol.t[-1]
        self.M, Phi_R, self.R_prop = sol.y[0, -1], sol.y[1, -1], sol.y[2, -1]
        self.dPhi0 = Phi_R - 0.5 * math.log(1.0 - 2.0 * self.M / self.R)   # Phi_fis = Phi - dPhi0
        self.Gam = self.eos.Gam

    def resumen(self):
        return dict(R=self.R, R_prop=self.R_prop, M=self.M / MSUN_KM, r_g=2.0 * self.M)

    def F(self, w2):
        """Delta P en la superficie para omega^2 = w2 (medido en el infinito)."""
        e_c, P_c, _, _ = self.eos.termo(1.0, self.rho_bc, self.sig)
        y0 = [1.0, -3.0 * self.Gam * P_c]

        def rhs(r, y):
            xi, DP = y
            m, Phi, eps, P, dPhi = self.bg(r, self.sol(r))
            hP = eps + P
            Pr = -hP * dPhi
            el = 1.0 / (1.0 - 2.0 * m / r)
            enu = math.exp(2.0 * (Phi - self.dPhi0))
            dxi = -(3.0 * xi + DP / (self.Gam * P)) / r - Pr * xi / hP
            dDP = (xi * (w2 * el / enu * hP * r - 4.0 * Pr)
                   + xi * (Pr * Pr * r / hP - PI8 * el * hP * P * r)
                   + DP * (Pr / hP - PI4 * hP * r * el))
            return [dxi, dDP]
        sol = solve_ivp(rhs, (self.r1, self.R), y0, method="DOP853", rtol=1e-10, atol=[1e-14, 1e-20])
        return sol.y[1, -1]

    def modo_fundamental(self, w_lo, w_hi, n=40):
        import numpy as np
        ws = np.linspace(w_lo, w_hi, n)
        Fs = [self.F(w) for w in ws]
        for i in range(n - 1):
            if Fs[i] * Fs[i + 1] < 0:
                return brentq(self.F, ws[i], ws[i + 1], xtol=1e-16, rtol=1e-11)
        return None
