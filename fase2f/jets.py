# -*- coding: utf-8 -*-
"""
Aritmetica de chorros (series de Taylor truncadas) vectorizada en numpy, para derivar exactamente (salvo redondeo)
composiciones de funciones a lo largo del fondo.

  Jet(c)   c[k] = f^(k)(l0)/k!,  forma (N+1, npts);  +, -, *, /, ** (exponente real), neg, comparacion de orden
  var(l0, N)               chorro de la variable l en los puntos l0
  const(x, N, npts)        constante
  Jet.d()                  derivada (pierde un orden)
  Jet.integral(c0)         primitiva con valor c0 en l0 (gana un orden; se trunca a N)
  Jet.deriv(k)             f^(k)(l0)
  exp, log, sqrt           funciones elementales
Los chorros se pueden pasar a funciones generadas por sympy.lambdify (operadores + modulos {'sqrt': sqrt, ...}).
"""
import math
import numpy as np


class Jet:
    __array_priority__ = 1000

    def __init__(self, c):
        self.c = np.asarray(c, dtype=float)

    @property
    def N(self):
        return self.c.shape[0] - 1

    # ---------------------------------------------------------------- utilidades
    @staticmethod
    def _como(x, ref):
        if isinstance(x, Jet):
            return x
        c = np.zeros_like(ref.c)
        c[0] = x
        return Jet(c)

    def deriv(self, k=0):
        return self.c[k] * math.factorial(k)

    def d(self):
        N = self.N
        c = np.zeros_like(self.c)
        for k in range(N):
            c[k] = (k + 1) * self.c[k + 1]
        c[N] = np.nan                      # orden perdido
        return Jet(c)

    def integral(self, c0):
        N = self.N
        c = np.zeros_like(self.c)
        c[0] = c0
        for k in range(1, N + 1):
            c[k] = self.c[k - 1] / k
        return Jet(c)

    # ---------------------------------------------------------------- aritmetica
    def __add__(self, o):
        o = Jet._como(o, self)
        return Jet(self.c + o.c)

    __radd__ = __add__

    def __sub__(self, o):
        o = Jet._como(o, self)
        return Jet(self.c - o.c)

    def __rsub__(self, o):
        o = Jet._como(o, self)
        return Jet(o.c - self.c)

    def __neg__(self):
        return Jet(-self.c)

    def __pos__(self):
        return self

    def __mul__(self, o):
        if not isinstance(o, Jet):
            return Jet(self.c * o)
        a, b = self.c, o.c
        N = min(a.shape[0], b.shape[0]) - 1
        c = np.zeros((N + 1,) + a.shape[1:])
        for k in range(N + 1):
            c[k] = np.sum(a[:k + 1] * b[k::-1], axis=0)
        return Jet(c)

    __rmul__ = __mul__

    def __truediv__(self, o):
        if not isinstance(o, Jet):
            return Jet(self.c / o)
        a, b = self.c, o.c
        N = min(a.shape[0], b.shape[0]) - 1
        c = np.zeros((N + 1,) + a.shape[1:])
        for k in range(N + 1):
            s = a[k] - (np.sum(b[1:k + 1] * c[k - 1::-1][:k], axis=0) if k else 0.0)
            c[k] = s / b[0]
        return Jet(c)

    def __rtruediv__(self, o):
        return Jet._como(o, self) / self

    def __pow__(self, p):
        if isinstance(p, Jet):
            return exp(log(self) * p)
        if float(p).is_integer() and p >= 0:
            p = int(p)
            out = Jet._como(1.0, self)
            base = self
            while p:
                if p & 1:
                    out = out * base
                base = base * base
                p >>= 1
            return out
        if float(p).is_integer() and p < 0:
            return 1.0 / (self ** (-int(p)))
        x = self.c
        N = self.N
        y = np.zeros_like(x)
        y[0] = x[0] ** p
        for k in range(1, N + 1):
            s = np.zeros_like(x[0])
            for i in range(1, k + 1):
                s = s + (p * i - (k - i)) * x[i] * y[k - i]
            y[k] = s / (k * x[0])
        return Jet(y)

    def __rpow__(self, b):
        return exp(self * math.log(b))


def var(l0, N):
    l0 = np.atleast_1d(np.asarray(l0, dtype=float))
    c = np.zeros((N + 1, l0.size))
    c[0] = l0
    if N >= 1:
        c[1] = 1.0
    return Jet(c)


def const(x, N, npts):
    c = np.zeros((N + 1, npts))
    c[0] = x
    return Jet(c)


def exp(x):
    if not isinstance(x, Jet):
        return np.exp(x)
    a = x.c
    N = x.N
    y = np.zeros_like(a)
    y[0] = np.exp(a[0])
    for k in range(1, N + 1):
        y[k] = sum(i * a[i] * y[k - i] for i in range(1, k + 1)) / k
    return Jet(y)


def log(x):
    if not isinstance(x, Jet):
        return np.log(x)
    a = x.c
    N = x.N
    y = np.zeros_like(a)
    y[0] = np.log(a[0])
    for k in range(1, N + 1):
        y[k] = (a[k] - sum(i * y[i] * a[k - i] for i in range(1, k)) / k) / a[0]
    return Jet(y)


def sqrt(x):
    if not isinstance(x, Jet):
        return np.sqrt(x)
    return x ** 0.5


MODULOS = [{"sqrt": sqrt, "exp": exp, "log": log}]
