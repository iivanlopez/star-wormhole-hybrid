# -*- coding: utf-8 -*-
"""
Exporta las ecuaciones de resultados/bH_simbolico.pkl (generar_bH.py) a codigo numerico: escribe bH_gen.py con
  SIMB_G                       nombres de los simbolos G_ij (i derivadas en phi, j en X) que usan las ecuaciones
  fondo(G, Phi, r, P)          [E_A, E_B, E_C, E_p] del fondo estatico
  lineal(G, Phi, r, P, w)      coeficientes de las ecuaciones lineales l = 0 en el gauge beta = c = 0:
                               dict ec -> {(v, k): coef} con v en (a, d, q) y k el orden de derivada en l;
                               ec en (a, d, q, be, c); la de be se devuelve dividida por i omega
  G = dict nombre -> valor (deben estar todos los de SIMB_G); Phi = (Phi0, Phi1, Phi2), r = (r0, r1, r2),
  P = (P0, P1, P2): derivadas en l del fondo (no aparecen ordenes mayores); w = omega^2 con el t de la metrica.
  E = e^{2 Phi0} se calcula de Phi0: los terminos en w llevan w/E = omega^2 e^{-2 Phi} (frecuencia local).
Uso: python exportar_bH.py
"""
import os
import sys
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
d = pickle.load(open(os.path.join(AQUI, "resultados", "bH_simbolico.pkl"), "rb"))
fe, lin, GS, Phi, rr, PP, Ee, Xs, iw, amp = (d[k] for k in ("fondo_eq", "lineales", "GS", "Phi", "rr", "PP", "Ee", "Xs",
                                                            "iw", "amp"))
w = sp.Symbol("w", real=True)
gauge = {amp[(v, k)]: 0 for v in ("be", "c") for k in range(6)}
nombresG = sorted({str(s) for s in GS.values()})
simbG = [sp.Symbol(n, real=True) for n in nombresG]
bg = list(Phi) + list(rr) + list(PP)
sub_fijo = {Xs: -PP[1] ** 2 / 2}


def limpio(e):
    return sp.expand(e.subs(sub_fijo))


fondo_list = [limpio(fe[n]) for n in ("A", "B", "C", "p")]
coef = {}
for ec in ("a", "d", "q", "be", "c"):
    e = limpio(lin[ec].subs(gauge))
    pol = sp.Poly(e, iw)
    pot = sorted(set(m[0] for m in pol.monoms()))
    if ec == "be":
        assert all(p % 2 == 1 for p in pot), pot          # E_be es impar en i omega
        e = sp.expand(e / iw)
    else:
        assert all(p % 2 == 0 for p in pot), (ec, pot)
    e = sp.expand(e.subs(iw ** 2, -w))
    assert not e.has(iw)
    for v in ("a", "d", "q"):
        for k in range(6):
            c = e.coeff(amp[(v, k)])
            if c != 0:
                coef[(ec, v, k)] = c
    resto = sp.expand(e - sum(coef[(ec, v, k)] * amp[(v, k)] for v in ("a", "d", "q") for k in range(6)
                              if (ec, v, k) in coef))
    assert resto == 0, (ec, resto)
claves = sorted(coef)
# solo los simbolos G que quedan en las expresiones finales (los demas se cancelan)
usados = set().union(*[e.free_symbols for e in fondo_list], *[c.free_symbols for c in coef.values()])
nombresG = [n for n in nombresG if sp.Symbol(n, real=True) in usados]
simbG = [sp.Symbol(n, real=True) for n in nombresG]
orden = max(k for k in range(len(Phi)) if any(Phi[k] in e.free_symbols or rr[k] in e.free_symbols or
                                              PP[k] in e.free_symbols for e in fondo_list + list(coef.values())))
assert orden == 2, orden
bg = list(Phi[:3]) + list(rr[:3]) + list(PP[:3])
args = simbG + bg + [Ee, w]
lineas = ['# -*- coding: utf-8 -*-', '"""Generado por exportar_bH.py: no editar a mano."""', "import math", "from math import sqrt", "",
          f"SIMB_G = {nombresG!r}", f"CLAVES = {claves!r}", ""]
f_fondo = sp.lambdify(args, fondo_list, modules="math", cse=True)
f_lin = sp.lambdify(args, [coef[k] for k in claves], modules="math", cse=True)
import inspect                                        # noqa: E402
src_f = inspect.getsource(f_fondo).replace("def _lambdifygenerated(", "def _fondo(")
src_l = inspect.getsource(f_lin).replace("def _lambdifygenerated(", "def _lineal(")
lineas += [src_f, "", src_l, ""]
lineas += ["def fondo(G, Phi, r, P):",
           "    return _fondo(*[G[n] for n in SIMB_G], *Phi[:3], *r[:3], *P[:3], math.exp(2 * Phi[0]), 0.0)", "",
           "def lineal(G, Phi, r, P, w):",
           "    vals = _lineal(*[G[n] for n in SIMB_G], *Phi[:3], *r[:3], *P[:3], math.exp(2 * Phi[0]), w)",
           "    out = {}",
           "    for (ec, v, k), c in zip(CLAVES, vals):",
           "        out.setdefault(ec, {})[(v, k)] = c",
           "    return out", ""]
open(os.path.join(AQUI, "bH_gen.py"), "w", encoding="utf-8").write("\n".join(lineas))
print(f"escrito bH_gen.py: {len(nombresG)} simbolos G, {len(claves)} coeficientes lineales")
for ec in ("a", "d", "q", "be", "c"):
    ks = [(v, k) for (e_, v, k) in claves if e_ == ec]
    print(f"  E_{ec}: " + ", ".join(f"{v}{k}" for v, k in ks))
