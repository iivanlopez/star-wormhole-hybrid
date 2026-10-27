# -*- coding: utf-8 -*-
"""
Fase 2k: estructura de la accion cuadratica l = 0 en el gauge unitario temporal (dphi = 0) con phi = qt t + psi(l).
Lista los pares de chorros de (a, beta, c, d) con coeficiente no nulo en la capa de masas (2 chorros exactos) y separa
las partes par e impar en el numero de derivadas temporales.
Uso: python estructura_t.py [G5]
"""
import os
import sys
import pickle
import itertools
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_t import chorro_en_capa   # noqa: E402

CON_G5 = len(sys.argv) > 1 and sys.argv[1] == "G5"
NOMBRE = "t5" if CON_G5 else "t"
d = pickle.load(open(os.path.join(AQUI, "resultados", f"{NOMBRE}_simbolico.pkl"), "rb"))
fe, J, L2 = d["fondo_eq"], d["J"], d["L2"]
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
pj = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in ("a", "be", "c", "d", "q") for o in ORD}
L2u = sp.expand(L2.xreplace({pj[("q", o)]: 0 for o in ORD}))
chs = [chorro_en_capa(fe, J, 700 + s, orden=4, cero=() if CON_G5 else ("G5",)) for s in range(2)]
P = sp.Poly(L2u, *[pj[(v, o)] for v in ("a", "be", "c", "d") for o in ORD])
filas = []
for mon, coef in zip(P.monoms(), P.coeffs()):
    gens = [g for g, e in zip(P.gens, mon) for _ in range(e)]
    vals = [ch.ev(coef) for ch in chs]
    nt = sum(int(str(g).split("_")[1][0]) for g in gens)
    filas.append((" * ".join(str(g) for g in gens), nt, vals[0] != 0 or vals[1] != 0, len(sp.Add.make_args(coef))))
for f in sorted(filas, key=lambda f: (f[1], f[0])):
    print(f"{f[0]:16s} nt={f[1]}  {'NO NULO' if f[2] else 'nulo en capa'}  ({f[3]} terminos)")
