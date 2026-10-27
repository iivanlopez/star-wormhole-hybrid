# -*- coding: utf-8 -*-
"""
Fase 2j: por que K = 0 en el agujero negro de Lu-Pang.  Comprueba, en chorros con las funciones de la teoria de Lu-Pang
(G2 = 8 al X^2, G3 = -8 al X, G4 = 1 + 4 al X, G5 = -4 al ln|X|) y datos iniciales (r0, r', phi', Phi') al azar,
  1. si el sistema del fondo E_A = E_C = E_p = 0 determina (Phi'', r'', phi'') (rango de la matriz de esas segundas
     derivadas); en la 4D-EGB regularizada la ecuacion del escalar es una combinacion de las de la metrica;
  2. si K = 0 en cualquier solucion de esa teoria (no solo en el agujero negro), resolviendo E_B para Phi'.
Uso: python degeneracion_lp.py
"""
import os
import sys
import json
import random
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_q import Phi, rr, PP, Ee, Xs, D        # noqa: E402

q = pickle.load(open(os.path.join(AQUI, "resultados", "q_simbolico.pkl"), "rb"))
R = pickle.load(open(os.path.join(AQUI, "resultados", "reduccion_q.pkl"), "rb"))
fe = q["fondo_eq"]
al = sp.Rational(1, 2)
GLP = {"G2": 8 * al * Xs ** 2, "G3": -8 * al * Xs, "G4": 1 + 4 * al * Xs, "F4": sp.Integer(0),
       "G5": -4 * al * sp.log(-Xs), "F5": sp.Integer(0)}


def teoria(expr):
    """Sustituye los G_ij por las derivadas de las funciones de Lu-Pang en Xs (simbolo) y X = -P1^2/2."""
    rep = {}
    for s in expr.free_symbols:
        n = s.name
        if "_" in n and n.split("_")[0] in GLP and len(n.split("_")[1]) == 2:
            nom, ij = n.split("_")
            i, j = int(ij[0]), int(ij[1])
            rep[s] = 0 if i > 0 else sp.diff(GLP[nom], Xs, j)
    return expr.xreplace(rep)


RES = {}
fl = {k: sp.expand(teoria(v)) for k, v in fe.items()}
# 1. matriz de las segundas derivadas en (E_A, E_C, E_p)
inc = [Phi[2], rr[2], PP[2]]
Mk = sp.Matrix([[sp.diff(fl[e], u) for u in inc] for e in ("A", "C", "p")])
rnd = random.Random(7)
rangos = []
for _ in range(4):
    val = {Ee: sp.Rational(121, 100), rr[0]: sp.Rational(rnd.randint(1500, 4000), 1000),
           rr[1]: sp.Rational(rnd.randint(300, 900), 1000), PP[1]: -sp.Rational(rnd.randint(50, 300), 1000),
           Phi[1]: sp.Rational(rnd.randint(50, 400), 1000)}
    val[Xs] = -val[PP[1]] ** 2 / 2
    Mv = Mk.xreplace(val).applyfunc(lambda x: sp.nsimplify(sp.N(x, 40), rational=False))
    Mn = sp.Matrix(Mv).evalf(40)
    rangos.append(int(Mn.rank(iszerofunc=lambda x: abs(x) < 1e-30)))
RES["rango_segundas"] = rangos
print("rango de d(E_A, E_C, E_p)/d(Phi'', r'', phi'') con las funciones de Lu-Pang (4 puntos al azar):", rangos)
# 2. soluciones genericas de la teoria de Lu-Pang: E_B = 0 para Phi', y (Phi_k, r_k, phi_k) de E_A = E_C = E_p = 0
g_, a1_, a0_ = R["extra"]
sub = {}
for simb, expr in ((g_, R["gam"]), (a1_, R["al1"]), (a0_, R["al0"])):
    e = expr
    for k in range(4):
        sub[simb[k]] = e
        e = D(e)
Kexp = sp.expand(teoria(R["K"].xreplace(sub)))
Pexp = sp.expand(teoria(R["P"].xreplace(sub)))
eqs0 = [fl["A"], fl["C"], fl["p"]]
sols = []
for _ in range(3):
    val = {Ee: sp.Rational(121, 100), rr[0]: sp.Rational(rnd.randint(1500, 4000), 1000),
           rr[1]: sp.Rational(rnd.randint(300, 900), 1000), PP[1]: -sp.Rational(rnd.randint(50, 300), 1000)}
    val[Xs] = -val[PP[1]] ** 2 / 2
    eB = sp.expand(fl["B"].xreplace(val))
    raices = [x for x in sp.Poly(sp.nsimplify(eB, rational=True), Phi[1]).nroots(n=40) if abs(sp.im(x)) < 1e-30]
    if not raices:
        continue
    val[Phi[1]] = sp.re(raices[0])
    eqs = eqs0
    ok = True
    for k in range(2, 6):
        inc = [Phi[k], rr[k], PP[k]]
        ev = [sp.expand(e.xreplace(val).subs(Xs, val[Xs])) for e in eqs]
        Mk_ = sp.Matrix([[x.coeff(u) for u in inc] for x in ev]).evalf(40)
        bk = sp.Matrix([sp.expand(x - sum(x.coeff(u) * u for u in inc)) for x in ev]).evalf(40)
        if abs(Mk_.det()) < 1e-30:
            ok = False
            break
        so = Mk_.LUsolve(-bk)
        for u, s in zip(inc, so):
            val[u] = s
        if k < 5:
            eqs = [D(e) for e in eqs]
    if not ok:
        sols.append(dict(Phi1=float(val[Phi[1]]), degenerado=True))
        print("solucion generica de Lu-Pang: sistema de segundas derivadas degenerado en la ligadura")
        continue
    Kv = float(sp.N(Kexp.xreplace(val).subs(Xs, val[Xs]), 30))
    Pv = float(sp.N(Pexp.xreplace(val).subs(Xs, val[Xs]), 30))
    sols.append(dict(r=float(val[rr[0]]), r1=float(val[rr[1]]), P1=float(val[PP[1]]), Phi1=float(val[Phi[1]]), K=Kv, P=Pv))
    print(f"solucion generica de Lu-Pang (r {float(val[rr[0]]):.3f}, r' {float(val[rr[1]]):.3f}, phi' {float(val[PP[1]]):.3f}, "
          f"Phi' {float(val[Phi[1]]):.4f}): K = {Kv:+.4e}, P = {Pv:+.4e}")
RES["soluciones_genericas"] = sols
json.dump(RES, open(os.path.join(AQUI, "resultados", "degeneracion_lp.json"), "w"), indent=1)
