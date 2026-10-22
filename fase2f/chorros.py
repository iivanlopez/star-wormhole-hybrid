# -*- coding: utf-8 -*-
"""
Utilidades exactas sobre chorros del fondo estatico (Fase 2f), compartidas por las comprobaciones y las reducciones.

  simbolos   Phi_k, r_k, P_k (k-esima derivada en l de Phi, r, phi), E = e^{2 Phi0}, eP = e^{Phi0}, Xs = X_b y G_ij
             (i derivadas en phi, j en X de G2, G3, G4, F4 en (phi_b, X_b)); mismos nombres que generar_bH.py
  D(e, amps) derivada total en l (fondo, G_ij por la regla de la cadena y listas de amplitudes de perturbaciones)
  Valores    valores racionales exactos de un chorro; los G_ij se rellenan al azar al aparecer
  chorro_en_capa(fe, semilla, orden)
             chorro en la capa de masas: G2_00 de la ligadura E_B = 0 y (Phi_k, r_k, P_k), k = 2..orden, de
             E_A = E_C = E_p = 0 y sus derivadas en l (aritmetica racional exacta)
"""
import re
import random
import sympy as sp

Phi = sp.symbols("Phi0:9", real=True)
rr = sp.symbols("r0:9", real=True)
PP = sp.symbols("P0:10", real=True)
Ee = sp.Symbol("E", positive=True)
eP = sp.Symbol("eP", positive=True)
Xs = sp.Symbol("Xs", real=True)
GPAT = re.compile(r"^(G2|G3|G4|F4)_(\d)(\d)$")


def gsym(nombre, i, j):
    return sp.Symbol(f"{nombre}_{i}{j}", real=True)


def es_G(s):
    return GPAT.match(s.name) is not None


def D(e, amps=()):
    """Derivada total en l: fondo (Phi_k, r_k, P_k -> k+1), E -> 2 Phi1 E, eP -> Phi1 eP, X_b -> -P1 P2,
    G_ij -> G_(i+1)j P1 - G_i(j+1) P1 P2, y las amplitudes de 'amps' (listas ordenadas por derivada)."""
    fs = e.free_symbols
    assert Phi[0] not in fs
    out = 0
    for seq in (Phi, rr, PP) + tuple(amps):
        for k in range(len(seq) - 1):
            if seq[k] in fs:
                out += sp.diff(e, seq[k]) * seq[k + 1]
    if Ee in fs:
        out += sp.diff(e, Ee) * 2 * Phi[1] * Ee
    if eP in fs:
        out += sp.diff(e, eP) * Phi[1] * eP
    if Xs in fs:
        out += sp.diff(e, Xs) * (-PP[1] * PP[2])
    for s in fs:
        m = GPAT.match(s.name)
        if m:
            nm, i, j = m.group(1), int(m.group(2)), int(m.group(3))
            out += sp.diff(e, s) * (gsym(nm, i + 1, j) * PP[1] - gsym(nm, i, j + 1) * PP[1] * PP[2])
    return sp.expand(out)


class Valores(dict):
    """Valores racionales exactos de un chorro; los simbolos G_ij se rellenan al azar al aparecer."""

    def __init__(self, semilla):
        super().__init__()
        self.rnd = random.Random(semilla)

    def r(self, a=-900, b=900):
        return sp.Rational(self.rnd.randint(a, b), 997)

    def rellena_G(self, e, salvo=()):
        for s in sorted(e.free_symbols, key=str):
            if es_G(s) and s not in self and s not in salvo:
                self[s] = self.r()

    def ev(self, e, salvo=()):
        self.rellena_G(e, salvo)
        return sp.expand(e.xreplace(self))


def chorro_en_capa(fe, semilla, orden=3):
    """Chorro de fondo en la capa de masas hasta 'orden' (ecuaciones del fondo fe y sus derivadas), exacto."""
    v = Valores(semilla)
    v[Ee] = sp.Rational(121, 100)
    v[eP] = sp.Rational(11, 10)
    v[Phi[1]] = v.r()
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[rr[1]] = v.r()
    v[PP[0]] = v.r()
    v[PP[1]] = v.r()
    v[Xs] = -v[PP[1]] ** 2 / 2
    g200 = gsym("G2", 0, 0)
    eB = v.ev(fe["B"], salvo=(g200,))
    c1 = eB.coeff(g200)
    assert c1 != 0 and sp.expand(eB - c1 * g200).free_symbols == set()
    v[g200] = -sp.expand(eB - c1 * g200) / c1
    eqs = [fe["A"], fe["C"], fe["p"]]
    for k in range(2, orden + 1):
        inc = [Phi[k], rr[k], PP[k]]
        ev = [v.ev(e) for e in eqs]
        Mk = sp.Matrix([[x.coeff(u) for u in inc] for x in ev])
        bk = sp.Matrix([sp.expand(x - sum(x.coeff(u) * u for u in inc)) for x in ev])
        assert all(not x.free_symbols for x in bk) and Mk.det() != 0
        sol = Mk.LUsolve(-bk)
        for u, s in zip(inc, sol):
            v[u] = sp.nsimplify(s)
        if k < orden:
            eqs = [D(e) for e in eqs]
    return v
