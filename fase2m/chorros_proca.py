# -*- coding: utf-8 -*-
"""
Fase 2m: utilidades exactas sobre chorros del fondo estatico de la Proca generalizada (como fase2j/chorros_q.py).
  D(e, amps)    derivada total en l (fondo, G_ij por la regla de la cadena con X_b y F_b, listas de amplitudes)
  Valores       valores racionales exactos; los G_ij se rellenan al azar al aparecer
  chorro_en_capa(fe, semilla, orden, caso)
     caso "general": A0, A1 != 0.  G2_00 de E_B = 0; en cada orden k, (N_k, r_k, U_k, G2_k0) de
                     D^(k-2){E_A, E_C, E_u} = 0 y D^(k-1) E_w = 0  (V_(k-1) libre: un chorro de alguna teoria).
     caso "R":       A0 = 0 (vector radial).  E_u = 0 identicamente; (N_k, r_k, G2_k0) de E_A, E_C, E_w.
     caso "E":       A1 = 0 (vector temporal), G3 = 0.  E_w = 0 identicamente; (N_k, r_k, U_k) de E_A, E_C, E_u.
  Por la identidad de Noether radial, E_B = 0 en el orden 0 basta.
"""
import re
import random
import sympy as sp

NMAXB = 9
NN = sp.symbols(f"N0:{NMAXB}", positive=True)
rr = sp.symbols(f"r0:{NMAXB}", real=True)
UU = sp.symbols(f"U0:{NMAXB}", real=True)
VV = sp.symbols(f"V0:{NMAXB}", real=True)
GPAT = re.compile(r"^(G2|G3|G4)_(\d+)(\d)$")
Xb = UU[0] ** 2 / (2 * NN[0] ** 2) - VV[0] ** 2 / 2
Fb = UU[1] ** 2 / (2 * NN[0] ** 2)


def gsym(nombre, i, j):
    return sp.Symbol(f"{nombre}_{i}{j}", real=True)


def es_G(s):
    return GPAT.match(s.name) is not None


def _D0(e):
    out = 0
    fs = e.free_symbols
    for seq in (NN, rr, UU, VV):
        for k in range(NMAXB - 1):
            if seq[k] in fs:
                out += sp.diff(e, seq[k]) * seq[k + 1]
    return out


dXb, dFb = sp.expand(_D0(Xb)), sp.expand(_D0(Fb))


def D(e, amps=()):
    fs = e.free_symbols
    out = _D0(e)
    for seq in amps:
        for k in range(len(seq) - 1):
            if seq[k] in fs:
                out += sp.diff(e, seq[k]) * seq[k + 1]
    for s in fs:
        m = GPAT.match(s.name)
        if m:
            nm, i, j = m.group(1), int(m.group(2)), int(m.group(3))
            out += sp.diff(e, s) * (gsym(nm, i + 1, j) * dXb + (gsym(nm, i, j + 1) * dFb if nm == "G2" else 0))
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
        self.rellena_G(e, salvo)
        return sp.expand(e.xreplace(self))


def chorro_en_capa(fe, semilla, orden=6, caso="general", fijos=None):
    cero = ("G3",) if caso == "E" else ()
    v = Valores(semilla, cero)
    if fijos:
        v.update(fijos)
    v[NN[0]] = sp.Rational(v.rnd.randint(500, 1500), 1000)
    v[NN[1]] = v.r()
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[rr[1]] = v.r()
    for k in range(NMAXB):
        if caso == "R":
            v[UU[k]] = 0
        if caso == "E":
            v[VV[k]] = 0
    if caso != "R":
        v[UU[0]] = v.r()
        v[UU[1]] = v.r()
    if caso != "E":
        v[VV[0]] = v.r()
    g200, g210 = gsym("G2", 0, 0), gsym("G2", 1, 0)
    con_w = caso != "E"
    if con_w:                                       # E_w no contiene G2_00: primero G2_10, luego G2_00 de E_B
        ew = v.ev(fe["w"], salvo=(g210,))
        c1 = ew.coeff(g210)
        assert c1 != 0 and not sp.expand(ew - c1 * g210).free_symbols, ew.free_symbols
        v[g210] = -sp.expand(ew - c1 * g210) / c1
    eB = v.ev(fe["B"], salvo=(g200,))
    c1 = eB.coeff(g200)
    assert c1 != 0 and not sp.expand(eB - c1 * g200).free_symbols
    v[g200] = -sp.expand(eB - c1 * g200) / c1
    eqs = [fe["A"], fe["C"]] + ([fe["u"]] if caso != "R" else [])
    ew = D(fe["w"]) if con_w else None
    for k in range(2, orden + 1):
        if caso != "E":
            v.setdefault(VV[k - 1], v.r())
        inc = [NN[k], rr[k]] + ([UU[k]] if caso != "R" else []) + ([gsym("G2", k, 0)] if con_w else [])
        lista = eqs + ([ew] if con_w else [])
        ev = [v.ev(e, salvo=inc) for e in lista]
        Mk = sp.Matrix([[x.coeff(u) for u in inc] for x in ev])
        bk = sp.Matrix([sp.expand(x - sum(x.coeff(u) * u for u in inc)) for x in ev])
        assert all(not x.free_symbols for x in bk), [x.free_symbols for x in bk]
        assert all(not x.free_symbols for x in Mk)
        assert Mk.det() != 0, (k, Mk)
        sol = Mk.LUsolve(-bk)
        for u, s in zip(inc, sol):
            v[u] = s
        if k < orden:
            eqs = [D(e) for e in eqs]
            if con_w:
                ew = D(ew)
    return v
