# -*- coding: utf-8 -*-
"""
Agujeros de gusano de Einstein-escalar-Gauss-Bonnet con acoplamiento general F(phi) (Antoniou et al. 2019,
arXiv:1904.13091), en longitud propia:  ds^2 = -e^{2 Phi} dt^2 + dl^2 + r^2 dOmega^2.
Ecuaciones: egb_F_gen.py (generar_F.py), validadas contra el caso de KKK.
Centro l = 0 (lado l > 0 de la solucion Z2): r = r0, r' = 0, Phi = 0, phi = phi0, phi' = phi1;
la ligadura ll da Phi'(0) = (1 + r0^2 phi1^2/4) / (4 F'(phi0) phi1).
r''(0) > 0: garganta (una garganta);  r''(0) < 0: ecuador, con una garganta lisa en l_t > 0 (doble garganta).
Comparacion con Antoniou: coordenada eta con ds_esp^2 = e^{f1} [deta^2 + (eta^2 + eta0^2) dOmega^2] y f1(inf) = 0:
eta0 = lim 2 r exp(-int_0^l dl'/r); su parametro es alpha/eta0^2.
"""
import math
from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import egb_F_gen as H


def acoplamiento(tipo, alpha, **kw):
    """Devuelve una funcion phi -> (F, F', F'', F''')."""
    if tipo == "exp":
        g = kw.get("gamma", 1.0)
        return lambda x: tuple(alpha * (-g) ** k * math.exp(-g * x) for k in range(4))
    if tipo == "cuadratico":
        return lambda x: (alpha * x * x, 2 * alpha * x, 2 * alpha, 0.0)
    if tipo == "gauss":
        g = kw.get("gamma", 1.0)

        def F(x):
            e = math.exp(-g * x * x)
            return (alpha * e, alpha * e * (-2 * g * x), alpha * e * (4 * g * g * x * x - 2 * g),
                    alpha * e * (-8 * g**3 * x**3 + 12 * g * g * x))
        return F
    raise ValueError(tipo)


def segundas(y, Fd):
    r, s, Phi, q, p, u = y[:6]
    F0, F1, F2, F3 = Fd(p)
    M = np.array(H.est_principal(F1, r, s, q, u))
    b = -np.array(H.est_resto(F1, F2, r, s, q, u))
    return np.linalg.solve(M, b), np.linalg.det(M)


@dataclass
class FondoF:
    tipo: str
    alpha: float
    phi0: float
    phi1: float
    r0: float = 1.0
    l_max: float = 2.0e4
    rtol: float = 1e-11
    kw: dict = field(default_factory=dict)

    def __post_init__(self):
        self.Fd = acoplamiento(self.tipo, self.alpha, **self.kw)

    def y0(self):
        F1 = self.Fd(self.phi0)[1]
        q0 = (1 + self.r0**2 * self.phi1**2 / 4) / (4 * F1 * self.phi1)
        return np.array([self.r0, 0.0, 0.0, q0, self.phi0, self.phi1, 0.0])   # ultima: int dl / r

    def rhs(self, l, y):
        x, _ = segundas(y, self.Fd)
        return [y[1], x[0], y[3], x[1], y[5], x[2], 1.0 / y[0]]

    def resolver(self):
        y0 = self.y0()
        x0, det0 = segundas(y0, self.Fd)
        self.rpp0 = x0[0]
        self.ecuador = x0[0] < 0
        self.f0 = 1.0 / (self.r0 * x0[0])

        def ev_det(l, y):
            return segundas(y, self.Fd)[1] / abs(det0)
        ev_det.terminal = True

        def ev_r(l, y):                         # r -> 0: singular
            return y[0] - 1e-3 * self.r0
        ev_r.terminal = True

        def ev_rp(l, y):                        # cambios de signo de r' (garganta lisa o nuevo ecuador)
            return y[1] if l > 1e-9 else (1.0 if not self.ecuador else -1.0)
        ev_rp.terminal = False

        sol = solve_ivp(self.rhs, (0.0, self.l_max), y0, method="DOP853", rtol=self.rtol, atol=1e-14,
                        events=[ev_det, ev_r, ev_rp], dense_output=True)
        self.sol = sol
        self.l_end = sol.t[-1]
        cruces = sol.t_events[2]
        self.l_garganta = cruces[0] if (self.ecuador and len(cruces) >= 1) else (0.0 if not self.ecuador else None)
        n_cruces_esperados = 1 if self.ecuador else 0
        self.regular = sol.status == 0 and len(cruces) == n_cruces_esperados
        r, s, Phi, q, p, u, I = sol.y[:, -1]
        self.M = r * (1 - s * s) / 2
        self.D = r * r * u
        self.phi_inf = p + self.D / r
        self.Phi_inf = Phi + q * r
        self.eta0 = 2 * r * math.exp(-I) * math.exp(-self.M / r * math.log(r)) if r > 0 else float("nan")
        self.eta0 = 2 * r * math.exp(-I)
        self.alpha_eta = self.alpha / self.eta0**2
        self.nodos = int(np.sum(np.diff(np.sign(sol.y[4] - self.phi_inf)) != 0))
        if self.l_garganta:
            self.r_garganta = float(sol.sol(self.l_garganta)[0])
        return self

    def capa(self, zeta=1.0):
        """Capa de polvo en l = 0 con lambda0 = lambda1 (Antoniou, KKK) escalada por zeta:
        rho = 2 Phi'_0 + 2 kappa F0/r0^2 - 8 F'(phi0) phi'_0/r0^2,  kappa = zeta 4 Phi'_0."""
        y0 = self.y0()
        q0, u0 = y0[3], y0[5]
        F0, F1 = self.Fd(self.phi0)[:2]
        kappa = zeta * 4 * q0
        return 2 * q0 + 2 * kappa * F0 / self.r0**2 - 8 * F1 * u0 / self.r0**2


def disparo_phi_inf(tipo, alpha, phi0, bracket, objetivo=0.0, **kw):
    """phi1 tal que phi_inf = objetivo (con phi0 y alpha fijos, r0 = 1)."""
    def g(u):
        f = FondoF(tipo, alpha, phi0, u, l_max=4e3, rtol=1e-10, kw=kw).resolver()
        if not f.regular:
            return float("nan")
        return f.phi_inf - objetivo
    u = brentq(g, *bracket, xtol=1e-13, rtol=1e-12)
    return FondoF(tipo, alpha, phi0, u, kw=kw).resolver()
