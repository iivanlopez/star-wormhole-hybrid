# -*- coding: utf-8 -*-
"""
Fase 2j: problema reducido del monopolo en el marco de la semilla (agujero negro de Lu-Pang en r > r0), para el sector par
de las gargantas de Bakopoulos-Charmousis-Kanti.

Con la accion reducida S = int [K d_t^2 - P d_lbar^2 - U d^2] dlbar (reduccion_q.py) y dlbar = dr/sqrt(h):
    -(P_r y')' + U_r y = w K_r y   en r > r0,   P_r = P sqrt(h),  U_r = U/sqrt(h),  K_r = K/sqrt(h),   w = omega^2,
con condicion natural y'(r0) = 0: las funciones pares de la longitud propia l de la garganta son funciones regulares de
r - r0 ~ l^2, sin condicion en r0, y la condicion natural del problema variacional es P y'(r0) = 0.
Los modos inestables son los w < 0; se cuentan con el angulo de Prufer (como fase2f/espectro.py).
"""
import math
import numpy as np
from scipy.integrate import solve_ivp


class SLr:
    """-(P y')' + U y = w K y en [r0, rmax] con funciones P, U, K de r (vectorizadas)."""

    def __init__(self, r0, rmax, fP, fU, fK):
        self.r0, self.rmax = r0, rmax
        self.fP, self.fU, self.fK = fP, fU, fK

    def prufer(self, w, rmax=None, rtol=1e-10):
        rmax = self.rmax if rmax is None else rmax
        th0 = math.pi / 2                          # y = 1, P y' = 0 (Neumann)

        def rhs(s, th):                            # s = ln(r - r0 + 1)
            r = self.r0 + math.exp(s) - 1
            c, sn = math.cos(th[0]), math.sin(th[0])
            return [math.exp(s) * (c * c / self.fP(r) + (w * self.fK(r) - self.fU(r)) * sn * sn)]
        sol = solve_ivp(rhs, (0.0, math.log(rmax - self.r0 + 1)), [th0], method="DOP853", rtol=rtol, atol=1e-12)
        return th0, sol.y[0][-1]

    def conteo(self, w, rmax=None):
        th0, thL = self.prufer(w, rmax)
        return int(math.floor(thL / math.pi)) - int(math.floor(th0 / math.pi))

    def autovalores_negativos(self, n=3, wmin=-1.0, rmax=None, tol=1e-9):
        while self.conteo(wmin, rmax) > 0:
            wmin *= 4
        out = []
        for k in range(n):
            a, b = wmin, -1e-14
            if self.conteo(b, rmax) <= k:
                break
            while b - a > tol * max(1e-6, abs(a)):
                m = (a + b) / 2
                if self.conteo(m, rmax) > k:
                    b = m
                else:
                    a = m
            out.append((a + b) / 2)
        return out
