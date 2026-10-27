# -*- coding: utf-8 -*-
"""
Fase 2k: utilidades exactas sobre chorros del fondo estatico con escalar lineal en el tiempo, phi = qt t + psi(l)
(como fase2j/chorros_q.py).  Simbolos G_0j (j derivadas en X) de G2, G3, G4, F4, G5.
  D(e, amps)   derivada total en l (fondo, X_b, G_0j por la regla de la cadena, listas de amplitudes)
  Valores      valores racionales exactos; los G_0j se rellenan al azar al aparecer (o se anulan los de 'cero')
  chorro_en_capa(fe, J, semilla, orden, cero)
               chorro en la capa de masas: G2_00 de E_B = 0 y (Phi_k, r_k, P_k), k = 2..orden, de E_A = E_C = J = 0
               (J = 0 da E_b = E_p = 0; la identidad de Noether hace E_B constante).
"""
import re
import random
import sympy as sp

Phi = sp.symbols("Phi0:9", real=True)
rr = sp.symbols("r0:9", real=True)
PP = sp.symbols("P0:10", real=True)
Ee = sp.Symbol("E", positive=True)
Xs = sp.Symbol("Xs", real=True)
qt = sp.Symbol("qt", real=True)
XB = qt ** 2 / (2 * Ee) - PP[1] ** 2 / 2
DXB = -qt ** 2 * Phi[1] / Ee - PP[1] * PP[2]
GPAT = re.compile(r"^(G2|G3|G4|F4|G5)_0(\d)$")


def gsym(nombre, j):
    return sp.Symbol(f"{nombre}_0{j}", real=True)


def es_G(s):
    return GPAT.match(s.name) is not None


def D(e, amps=()):
    e = sp.sympify(e)
    fs = e.free_symbols
    out = 0
    for seq in (Phi, rr, PP) + tuple(amps):
        for k in range(len(seq) - 1):
            if seq[k] in fs:
                out += sp.diff(e, seq[k]) * seq[k + 1]
    if Ee in fs:
        out += sp.diff(e, Ee) * 2 * Phi[1] * Ee
    if Xs in fs:
        out += sp.diff(e, Xs) * DXB
    for s in fs:
        m = GPAT.match(s.name)
        if m:
            out += sp.diff(e, s) * gsym(m.group(1), int(m.group(2)) + 1) * DXB
    return sp.expand(out)


class Valores(dict):
    def __init__(self, semilla, cero=()):
        super().__init__()
        self.rnd = random.Random(semilla)
        self.cero = tuple(cero)

    def r(self, a=-900, b=900):
        return sp.Rational(self.rnd.randint(a, b), 997)

    def rellena_G(self, e, salvo=()):
        for s in sorted(e.free_symbols, key=str):
            if es_G(s) and s not in self and s not in salvo:
                self[s] = 0 if GPAT.match(s.name).group(1) in self.cero else self.r()

    def ev(self, e, salvo=()):
        e = sp.sympify(e)
        self.rellena_G(e, salvo)
        return sp.expand(e.xreplace(self))


_CACHE = {}


def chorro_en_capa(fe, J, semilla, orden=3, cero=("G5",), fijos=None):
    """Chorro exacto de la capa de masas.  El sistema (E_A, E_B, J) es degenerado en sus derivadas mas altas (una
    combinacion es una ligadura de orden menor), asi que se usa G2 como incognita (reconstruccion local):
      nivel 2: (r_2, P_2, G2_00); nivel k >= 3: (r_k, P_k, G2_0(k-1)) de (D^{k-2} J, D^{k-2} E_A, D^{k-2} E_B)
      (E_B ya contiene G2_01, que se fija al azar en el nivel 2);  Phi_k (k >= 2) al azar.
    E_p = -D J y, por la identidad de Noether, E_C = 0 (r' != 0).  Se llega hasta r, P de orden orden + 1."""
    v = Valores(semilla, cero)
    v[Ee] = sp.Rational(121, 100)
    v[qt] = sp.Rational(v.rnd.randint(200, 900), 1000)
    v[Phi[1]] = v.r()
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[rr[1]] = v.r()
    v[PP[1]] = v.r()
    for k in range(2, orden + 2):
        v[Phi[k]] = v.r()
    if fijos:
        v.update(fijos)
    v[Xs] = v[qt] ** 2 / (2 * v[Ee]) - v[PP[1]] ** 2 / 2
    clave = (id(fe), id(J))
    lista = _CACHE.setdefault(clave, [[J, fe["A"], fe["B"]]])
    for k in range(2, orden + 2):
        while len(lista) <= k - 2:
            lista.append([D(e) for e in lista[-1]])
        eqs = lista[k - 2]
        inc = [rr[k], PP[k], gsym("G2", 0 if k == 2 else k - 1)]
        ev = [v.ev(e, salvo=(inc[2],)) for e in eqs]
        Mk = sp.Matrix([[x.coeff(u) for u in inc] for x in ev])
        bk = sp.Matrix([sp.expand(x - sum(x.coeff(u) * u for u in inc)) for x in ev])
        assert all(not x.free_symbols for x in bk), [x.free_symbols for x in bk]
        assert all(not x.free_symbols for x in Mk), Mk
        assert Mk.det() != 0, (k, Mk)
        sol = Mk.LUsolve(-bk)
        for u, s_ in zip(inc, sol):
            v[u] = sp.Rational(s_)
    return v
