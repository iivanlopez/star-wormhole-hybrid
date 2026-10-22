# -*- coding: utf-8 -*-
"""
Espectro lineal de los sectores l = 0 (monopolo, accion reducida de monopolo_reduccion.py) e impar (l >= 2, accion
de chi) sobre una garganta reconstruida (reconstruccion.Garganta).

Ambos sectores dan un problema de Sturm-Liouville   -(P y')' + U y = w K y,   l in (0, inf),   w = omega^2,
con paridad en la garganta (fondo simetrico):  par y'(0) = 0,  impar y(0) = 0.
Con K > 0 y P > 0 (sin fantasmas ni inestabilidad de gradiente: problema bien planteado) el espectro es real y
los modos inestables -- taquiones, lentos o rapidos -- son exactamente los autovalores w < 0.  Su numero es el de
ceros en (0, inf) de la solucion de w = 0 regular en la garganta (teorema de oscilacion de Sturm; el continuo empieza
en w = 0).  Los ceros se cuentan con el angulo de Prufer  y = rho sin(th), P y' = rho cos(th),
    th' = cos^2(th)/P + (w K - U) sin^2(th)      (th crece al cruzar th = m pi porque P > 0),
y mas alla de la malla con la forma de Liouville  x = int sqrt(K/P) dl,  Y = (K P)^{1/4} y,  Y'' = V Y:  para w = 0
la solucion es Y ~ A + B x (monopolo) y hay un cero mas si Y y Y_x tienen signos opuestos en el borde.
Velocidad radial fisica: c^2 = P/(K N^2) (la de la luz es 1).
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline


def malla(lmin=1e-3, L=400.0, n=6000):
    a = np.geomspace(lmin, 1.0, n // 4)
    b = np.geomspace(1.0, L, 3 * n // 4 + 1)[1:]
    return np.concatenate([a, b])


def coeficientes(g, sector, ls, LL=2, trozo=300):
    """K, P, U, N en la malla (por trozos para limitar la memoria de los chorros)."""
    K, P, U, N = [], [], [], []
    for i in range(0, len(ls), trozo):
        l = ls[i:i + trozo]
        o = g.monopolo(l) if sector == "monopolo" else g.impar(l, LL)
        K.append(o["K"].deriv(0))
        P.append(o["P"].deriv(0))
        U.append(o["U"].deriv(0))
        N.append(o["V"]["N"].deriv(0))
    return tuple(np.concatenate(x) for x in (K, P, U, N))


class SL:
    def __init__(self, ls, K, P, U):
        self.ls, self.K, self.P, self.U = ls, K, P, U
        x = np.log(ls)
        self.fK = CubicSpline(x, np.log(np.abs(K)) if np.all(K > 0) else K)
        self.fP = CubicSpline(x, np.log(np.abs(P)) if np.all(P > 0) else P)
        self.fU = CubicSpline(x, U)
        self.logK = np.all(K > 0)
        self.logP = np.all(P > 0)

    def _k(self, l):
        v = self.fK(np.log(l))
        return np.exp(v) if self.logK else v

    def _p(self, l):
        v = self.fP(np.log(l))
        return np.exp(v) if self.logP else v

    def prufer(self, w, paridad, L=None, rtol=1e-11):
        """Angulo de Prufer de la solucion regular en la garganta desde lmin hasta L."""
        l0 = self.ls[0]
        L = self.ls[-1] if L is None else L
        P0, U0, K0 = self._p(l0), self.fU(np.log(l0)), self._k(l0)
        if paridad == "par":
            c2 = (U0 - w * K0) / (2 * P0)
            y, yp = 1 + c2 * l0 ** 2, 2 * c2 * l0
        else:
            c3 = (U0 - w * K0) / (6 * P0)
            y, yp = l0 + c3 * l0 ** 3, 1 + 3 * c3 * l0 ** 2
        th0 = math.atan2(y, P0 * yp)

        def rhs(s, th):                       # variable s = log l
            l = math.exp(s)
            c, sn = math.cos(th[0]), math.sin(th[0])
            return [l * (c * c / self._p(l) + (w * self._k(l) - float(self.fU(s))) * sn * sn)]
        sol = solve_ivp(rhs, (math.log(l0), math.log(L)), [th0], method="DOP853", rtol=rtol, atol=1e-12)
        return th0, sol.y[0][-1], sol

    def ceros(self, paridad, L=None):
        """Numero de ceros en (0, L) de la solucion de w = 0 y diagnostico del borde (forma de Liouville)."""
        th0, thL, _ = self.prufer(0.0, paridad, L)
        n = int(math.floor(thL / math.pi)) - int(math.floor(th0 / math.pi))
        L = self.ls[-1] if L is None else L
        # Liouville en L: y'/y = cot(th)/P;  Y_x/Y = sqrt(P/K) [ (KP)^{1/4}'/(KP)^{1/4} + y'/y ]
        h = 1e-4 * L
        kp = lambda l: (self._k(l) * self._p(l)) ** 0.25                                  # noqa: E731
        dlog = (math.log(kp(L + h)) - math.log(kp(L - h))) / (2 * h)
        ypy = math.cos(thL) / math.sin(thL) / self._p(L)
        Yx_Y = math.sqrt(self._p(L) / self._k(L)) * (dlog + ypy)
        return n, Yx_Y

    def conteo(self, w, paridad, L=None):
        """Numero de ceros en (0, L) de la solucion regular con energia w; para w < 0 y L mucho mayor que la longitud
        de decaimiento es el numero de autovalores por debajo de w (oscilacion de Sturm)."""
        th0, thL, _ = self.prufer(w, paridad, L)
        return int(math.floor(thL / math.pi)) - int(math.floor(th0 / math.pi))

    def autovalores_negativos(self, paridad, n, wmin=-1.0, L=None, tol=1e-10):
        """Los n autovalores w < 0 (de menor a mayor) por biseccion en conteo(w)."""
        while self.conteo(wmin, paridad, L) > 0:
            wmin *= 4
        out = []
        for k in range(n):
            a, b = wmin, 0.0                    # conteo(a) <= k < conteo(b-)
            b = -1e-14
            if self.conteo(b, paridad, L) <= k:
                break
            while b - a > tol * max(1.0, abs(a)):
                m = (a + b) / 2
                if self.conteo(m, paridad, L) > k:
                    b = m
                else:
                    a = m
            out.append((a + b) / 2)
        return out
