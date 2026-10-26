# -*- coding: utf-8 -*-
"""
Fase 2i: simbolos y utilidades exactas del solido elastico sobre los chorros de la Fase 2f (chorros.py).
  ff_k        d^k f0/dl^k (coordenada radial material del fondo)
  Fs(i,j,k)   d^i_{X_r} d^j_{X_t} d^k_f F en el fondo (X_r = ff1^2, X_t = ff0^2/r^2)
  D(e, amps)  derivada total en l (chorros.D + ff_k -> ff_(k+1) + regla de la cadena de Fs)
  ValoresS    valores racionales exactos (rellena G_ij y Fs_ijk al azar)
  chorro_capa_solido(fe, semilla, orden)
              chorro en la capa de masas con solido: G2_00 de E_B = 0 y (Phi_k, r_k, P_k, ff_k), k = 2..orden, de
              E_A = E_C = E_p = E_f = 0 y sus derivadas
"""
import re
import sympy as sp
import chorros as CH
from chorros import Phi, rr, PP, Ee, eP, Xs

ff = sp.symbols("ff0:12", real=True)
FPAT = re.compile(r"^Fs_(\d)(\d)(\d)$")


def Fs(i, j, k):
    return sp.Symbol(f"Fs_{i}{j}{k}", real=True)


Xrb = ff[1] ** 2
Xtb = ff[0] ** 2 / rr[0] ** 2
dXr = 2 * ff[1] * ff[2]
dXt = 2 * ff[0] * ff[1] / rr[0] ** 2 - 2 * ff[0] ** 2 * rr[1] / rr[0] ** 3
rho = -Fs(0, 0, 0)
p_r = Fs(0, 0, 0) - 2 * Xrb * Fs(1, 0, 0)
p_t = Fs(0, 0, 0) - Xtb * Fs(0, 1, 0)


def D(e, amps=()):
    """Derivada total en l: la de chorros.D mas ff_k -> ff_(k+1) y la regla de la cadena para Fs_ijk."""
    out = CH.D(e, amps)
    fs = e.free_symbols
    for k in range(len(ff) - 1):
        if ff[k] in fs:
            out += sp.diff(e, ff[k]) * ff[k + 1]
    for s in fs:
        m = FPAT.match(s.name)
        if m:
            i, j, k = (int(x) for x in m.groups())
            out += sp.diff(e, s) * (Fs(i + 1, j, k) * dXr + Fs(i, j + 1, k) * dXt + Fs(i, j, k + 1) * ff[1])
    return sp.expand(out)


class ValoresS(CH.Valores):
    """Como chorros.Valores, rellenando tambien los Fs_ijk al azar."""

    def rellena_G(self, e, salvo=()):
        super().rellena_G(e, salvo)
        for s in sorted(e.free_symbols, key=str):
            if FPAT.match(s.name) and s not in self and s not in salvo:
                self[s] = self.r()


def chorro_capa_solido(fe, semilla, orden=3):
    v = ValoresS(semilla)
    v[Ee] = sp.Rational(121, 100)
    v[eP] = sp.Rational(11, 10)
    v[Phi[1]] = v.r()
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[rr[1]] = v.r()
    v[PP[0]] = v.r()
    v[PP[1]] = v.r()
    v[Xs] = -v[PP[1]] ** 2 / 2
    v[ff[0]] = sp.Rational(v.rnd.randint(300, 1500), 1000)
    v[ff[1]] = sp.Rational(v.rnd.randint(300, 1500), 1000)
    g200 = CH.gsym("G2", 0, 0)
    eB = v.ev(fe["B"], salvo=(g200,))
    c1 = eB.coeff(g200)
    v[g200] = -sp.expand(eB - c1 * g200) / c1
    eqs = [fe["A"], fe["C"], fe["p"], fe["f"]]
    for k in range(2, orden + 1):
        inc = [Phi[k], rr[k], PP[k], ff[k]]
        ev = [v.ev(e) for e in eqs]
        Mk = sp.Matrix([[x.coeff(u) for u in inc] for x in ev])
        bk = sp.Matrix([sp.expand(x - sum(x.coeff(u) * u for u in inc)) for x in ev])
        assert all(not x.free_symbols for x in bk) and Mk.det() != 0, [x.free_symbols for x in bk]
        sol = Mk.LUsolve(-bk)
        for u, s in zip(inc, sol):
            s = sp.nsimplify(s) if not s.is_Rational else s
            assert s.is_Rational, s
            v[u] = s
        if k < orden:
            eqs = [D(e) for e in eqs]
    assert v.ev(fe["B"]) == 0
    return v


def chorro_equilibrio(fe_f, semilla, orden=4):
    """Chorro con solo la ecuacion del solido E_f = 0 impuesta (metrica arbitraria): ff_k de D^(k-2) E_f = 0."""
    v = ValoresS(semilla)
    v[Ee] = sp.Rational(121, 100)
    v[eP] = sp.Rational(11, 10)
    for k in range(1, 9):
        v[Phi[k]] = v.r()
        v[rr[k]] = v.r()
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[ff[0]] = sp.Rational(v.rnd.randint(300, 1500), 1000)
    v[ff[1]] = sp.Rational(v.rnd.randint(300, 1500), 1000)
    e = fe_f
    for k in range(2, orden + 1):
        x = v.ev(e)
        c1 = x.coeff(ff[k])
        assert c1 != 0
        v[ff[k]] = -sp.expand(x - c1 * ff[k]) / c1
        e = D(e)
    return v
