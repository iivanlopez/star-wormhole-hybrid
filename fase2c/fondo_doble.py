# -*- coding: utf-8 -*-
"""
Agujeros de gusano EsGB LISOS (sin cascara) y en general asimetricos: el sistema estatico de fase2b/fondo_F.py
(egb_F_gen.est_principal / est_resto, funcion segundas) integrado desde la garganta l = 0 hacia los dos lados.

Datos lisos en l = 0 (r0 = 1):  r = 1, r' = 0, Phi = 0, phi = phi0, phi' = phi1,  Phi' = (1 + phi1^2/4)/(4 F'(phi0) phi1)
(ligadura ll).  Las ecuaciones son invariantes bajo l -> -l (todas las ecuaciones de fondo son pares en el numero de
derivadas), asi que el lado l < 0 es el lado "derecho" de los datos (phi0, -phi1) [Phi'(0) -> -Phi'(0)], leido con
l~ = -l.  Las dos mitades empalman con todas sus derivadas: la solucion es analitica a traves de la garganta.

Clases:
  Lado   una mitad (l~ >= 0) con clasificacion del final: 'AF' (asintoticamente plana), 'horizonte' (e^Phi -> 0 a l
         finito: Phi' -> -inf), 'r0' (r -> 0), 'det' (parte principal singular), 'phi' (phi' diverge), 'fallo'.
  Doble  las dos mitades pegadas; .sol.sol(l) para l de cualquier signo, compatible con fase2b/estabilidad_F
         (fondo_en, derivadas, ligaduras, velocidades).
"""
import os
import sys
import math
from dataclasses import dataclass, field
import numpy as np
from scipy.integrate import solve_ivp

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fase2b"))
from fondo_F import acoplamiento as acoplamiento_F, segundas      # noqa: E402
import egb_F_gen as H                                              # noqa: E402

PHIP_MAX = 1e7      # |Phi'| o |phi'| por encima de esto: fin singular (horizonte o divergencia del escalar)
DET_MIN = 1e-4      # fallo del integrador con |det M| < DET_MIN |det M(0)|: degeneracion de la parte principal


def acoplamiento(tipo, alpha, **kw):
    """fondo_F.acoplamiento mas el acoplamiento lineal F = alpha phi (sGB con simetria de traslacion)."""
    if tipo == "lineal":
        return lambda x: (alpha * x, alpha, 0.0, 0.0)
    return acoplamiento_F(tipo, alpha, **kw)


@dataclass
class Lado:
    tipo: str
    alpha: float
    phi0: float
    phi1: float
    l_max: float = 2.0e3
    rtol: float = 1e-10
    kw: dict = field(default_factory=dict)

    def __post_init__(self):
        self.Fd = acoplamiento(self.tipo, self.alpha, **self.kw)

    def y0(self):
        F1 = self.Fd(self.phi0)[1]
        q0 = (1 + self.phi1**2 / 4) / (4 * F1 * self.phi1)
        return np.array([1.0, 0.0, 0.0, q0, self.phi0, self.phi1, 0.0])

    def rhs(self, l, y):
        x, _ = segundas(y, self.Fd)
        return [y[1], x[0], y[3], x[1], y[5], x[2], 1.0 / y[0]]

    def resolver(self):
        y0 = self.y0()
        x0, det0 = segundas(y0, self.Fd)
        self.a = y0[3]                  # Phi'(0) en la direccion de este lado
        self.rpp0 = x0[0]               # r''(0) > 0: garganta; < 0: ecuador

        def ev_det(l, y):
            return segundas(y, self.Fd)[1] / abs(det0)
        ev_det.terminal = True

        def ev_r(l, y):
            return y[0] - 1e-3
        ev_r.terminal = True

        def ev_P(l, y):
            return PHIP_MAX - max(abs(y[3]), abs(y[5]))
        ev_P.terminal = True

        def ev_rp(l, y):
            return y[1] if l > 1e-9 else (1.0 if self.rpp0 > 0 else -1.0)

        sol = solve_ivp(self.rhs, (0.0, self.l_max), y0, method="DOP853", rtol=self.rtol, atol=1e-14,
                        events=[ev_det, ev_r, ev_P, ev_rp], dense_output=True)
        self.sol = sol
        self.l_end = sol.t[-1]
        self.cruces = list(sol.t_events[3])
        r, s, Phi, q, p, u, I = sol.y[:, -1]
        self.det_fin = segundas(sol.y[:, -1], self.Fd)[1] / det0
        if sol.status == 0:
            # AF: r' > 0 y masa de Misner-Sharp m = r (1 - r'^2)/2 convergida (m(l_max) frente a m(l_max/2)) y << r;
            # excluye las geometrias conicas (r' -> cte < 1, m proporcional a r)
            rh, sh = sol.sol(0.5 * self.l_end)[:2]
            m, mh = r * (1 - s * s) / 2, rh * (1 - sh * sh) / 2
            ok = s > 0 and abs(m) < 0.05 * r and abs(m - mh) <= 0.05 * max(1.0, abs(m))
            self.fin = "AF" if ok else "noAF"
        elif len(sol.t_events[0]):
            self.fin = "det"
        elif len(sol.t_events[1]):
            self.fin = "r0"
        elif len(sol.t_events[2]):
            self.fin = "horizonte" if abs(q) >= abs(u) and q < 0 else ("phi" if abs(u) > abs(q) else "Phi+")
        elif q < -1e3:
            self.fin = "horizonte"
        elif abs(self.det_fin) < DET_MIN:
            self.fin = "det"            # det M -> 0 sin cambiar de signo: el integrador se detiene antes del evento
        else:
            self.fin = "fallo"
        self.regular = self.fin == "AF"
        # magnitudes asintoticas (validas si AF): Phi ~ Phi_inf - M/r, phi ~ phi_inf - D/r  (D = r^2 phi')
        self.M = r * (1 - s * s) / 2
        self.M_Phi = q * r * r * (1 - 2 * self.M / r) ** 0.5 if r > 2 * self.M else float("nan")
        self.D = r * r * u
        self.phi_inf = p + self.D / r
        self.Phi_inf = Phi + q * r
        self.nodos = int(np.sum(np.diff(np.sign(sol.y[4] - self.phi_inf)) != 0))
        return self


class _Densa:
    """Salida densa de las dos mitades en la coordenada l de la solucion completa."""
    SIGNO = np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0])     # r, r', Phi, Phi', phi, phi', int dl/r

    def __init__(self, der, izq):
        self.der, self.izq = der, izq

    def sol(self, l):
        if l >= 0.0:
            return self.der.sol.sol(l)
        return self.izq.sol.sol(-l) * self.SIGNO


class Doble:
    """Agujero liso completo: l > 0 con (phi0, phi1), l < 0 con (phi0, -phi1) reflejado."""

    def __init__(self, tipo, alpha, phi0, phi1, l_max=2.0e3, rtol=1e-10, kw=None):
        kw = kw or {}
        self.tipo, self.alpha, self.phi0, self.phi1, self.kw = tipo, alpha, phi0, phi1, kw
        self.der = Lado(tipo, alpha, phi0, phi1, l_max, rtol, kw).resolver()
        self.izq = Lado(tipo, alpha, phi0, -phi1, l_max, rtol, kw).resolver()
        self.Fd = self.der.Fd
        self.sol = _Densa(self.der, self.izq)
        self.r0 = 1.0
        self.a = self.der.a                       # Phi'(0) (el mismo a los dos lados en la coordenada l)
        self.rpp0 = self.der.rpp0
        self.regular = self.der.regular and self.izq.regular
        self.l_end = min(self.der.l_end, self.izq.l_end)

    def resumen(self):
        d, i = self.der, self.izq
        return dict(tipo=self.tipo, alpha=self.alpha, phi0=self.phi0, phi1=self.phi1, a=self.a, rpp0=self.rpp0,
                    fin_der=d.fin, fin_izq=i.fin, l_der=d.l_end, l_izq=i.l_end,
                    M_der=d.M, M_izq=i.M, D_der=d.D, D_izq=i.D, phi_inf_der=d.phi_inf, phi_inf_izq=i.phi_inf,
                    Phi_inf_der=d.Phi_inf, Phi_inf_izq=i.Phi_inf, cruces_der=len(d.cruces), cruces_izq=len(i.cruces),
                    nodos_der=d.nodos, nodos_izq=i.nodos)
