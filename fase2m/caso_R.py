# -*- coding: utf-8 -*-
"""
Fase 2m: vector radial (caso R: A0 = 0, A = A1(l) dl).  El fondo es el de Horndeski con phi' = A1 (verificar_proca.py).
Gauge radial w = delta A1 = 0 (el analogo del gauge unitario); quedan a, beta, c, d y al = delta A0.
Comprueba (exacto):
  1. -L_rest(c = 0, d = N f) estatico (a = beta = al = 0) = -2 N^3 Gcal (f' - sig f)^2, Gcal = G4 + A1^2 G4_X  (fuera de
     la capa de masas; con G2(X, F) y G3 cualesquiera).
  2. La prueba d = N, c = a = beta = al = 0 cumple las ligaduras: E_a = E_beta = E_al = 0 en chorros de la capa de masas.
  3. al no se acopla a (a, beta) salvo derivadas totales:  a al_t + a_t al  (mismo coeficiente) es d_t(a al), y los
     terminos mu beta al + nu (beta al' + beta' al) son -(nu' - mu) beta al con nu' = mu en la capa de masas.
     Por tanto  L2 = -F W + L_rest(c, d) + L_al,  L_al = k_F al'^2 + m_al al^2 + n al al' + al J[c_t, d_t]:
     al solo se acopla a velocidades.  Su eliminacion no cambia la energia potencial (la del vacio de Horndeski) y
     resta -J M^-1 J/4 a la cinetica.
Uso: python caso_R.py
"""
import os
import sys
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import chorros_proca as CP                              # noqa: E402
from chorros_proca import NN, rr, UU, VV, gsym          # noqa: E402

d = pickle.load(open(os.path.join(AQUI, "resultados", "proca_simbolico.pkl"), "rb"))
L2, pj, fe, lin, amp, iw = (d[k] for k in ("L2", "pj", "fondo_eq", "lineales", "amp", "iw"))
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
cU = {UU[k]: 0 for k in range(CP.NMAXB)}
LR = sp.expand(L2.xreplace(cU).xreplace({pj[("w", o)]: 0 for o in ORD}))
CHK = {}

# 1. identidad de la familia c = 0
f = sp.symbols("f0:4", real=True)


def Df(e):
    out = CP.D(e)
    for k in range(3):
        out += sp.diff(e, f[k]) * f[k + 1]
    return sp.expand(out)


dN = [NN[0] * f[0]]
for _ in range(2):
    dN.append(Df(dN[-1]))
est = {pj[(v, o)]: 0 for v in ("a", "be", "c", "al") for o in ORD}
est.update({pj[("d", o)]: 0 for o in ORD if o[0] > 0})
est.update({pj[("d", (0, 0))]: dN[0], pj[("d", (0, 1))]: dN[1], pj[("d", (0, 2))]: dN[2]})
menosL = sp.expand(-LR.xreplace(est))
Gc = gsym("G4", 0, 0) + VV[0] ** 2 * gsym("G4", 1, 0)
sig = rr[1] / rr[0] - NN[1] / NN[0]
dif = sp.simplify(menosL + 2 * NN[0] ** 3 * Gc * (f[1] - sig * f[0]) ** 2)
CHK["identidad_c0"] = str(dif)
assert dif == 0, dif
print("1. -L_rest(0, N f) = -2 N^3 Gcal (f' - sig f)^2 con Gcal = G4 + A1^2 G4_X: identidad exacta (puntual)")

# 2. ligaduras en la prueba d = N
sub = {amp[(v, k)]: 0 for v in ("a", "be", "c", "al", "w") for k in range(6)}
e = NN[0]
for k in range(6):
    sub[amp[("d", k)]] = e
    e = CP._D0(e)
res = []
for sem in range(3):
    v = CP.chorro_en_capa(fe, 300 + sem, orden=3, caso="R")
    res.append({n: float(v.ev(sp.expand(lin[n].xreplace(cU).xreplace(sub).xreplace({iw: 0}))))
                for n in ("a", "be", "al")})
CHK["ligaduras_prueba"] = res
print("2. E_a, E_beta, E_al en la prueba d = N (3 chorros de la capa de masas):", res)
assert all(abs(x) == 0 for r_ in res for x in r_.values())

# 3. acoplos de al con (a, beta)
A0, A10 = pj[("a", (0, 0))], pj[("a", (1, 0))]
B0, B1 = pj[("be", (0, 0))], pj[("be", (0, 1))]
Al0, Al1, Al10 = pj[("al", (0, 0))], pj[("al", (0, 1))], pj[("al", (1, 0))]
k1, k2 = LR.coeff(A0).coeff(Al10), LR.coeff(A10).coeff(Al0)
mu = LR.coeff(B0).coeff(Al0)
nu1, nu2 = LR.coeff(B0).coeff(Al1), LR.coeff(B1).coeff(Al0)
otros = [t for t in sp.Add.make_args(LR) if (t.has(A0) or t.has(A10) or t.has(B0) or t.has(B1)
                                             or any(t.has(pj[(x, o)]) for x in ("a", "be") for o in ORD[3:])
                                             or t.has(pj[("a", (0, 1))]) or t.has(pj[("be", (1, 0))]))
         and any(t.has(pj[("al", o)]) for o in ORD)]
assert sp.expand(k1 - k2) == 0 and sp.expand(nu1 - nu2) == 0
PJ = set(pj.values())
monos = {sp.Mul(*sorted([s for s in t.free_symbols if s in PJ], key=str)) for t in otros}
assert monos == {A0 * Al10, A10 * Al0, B0 * Al0, B0 * Al1, B1 * Al0}, monos   # a al_t, a_t al, beta al, beta al', beta' al
nm = []
for sem in range(3):
    v = CP.chorro_en_capa(fe, 400 + sem, orden=4, caso="R")
    nm.append(float(v.ev(sp.expand(CP.D(nu1) - mu))))
CHK["nu_prima_menos_mu"] = nm
print("3. a al_t + a_t al con el mismo coeficiente (d_t(a al)); nu' - mu en la capa de masas:", nm)
assert all(x == 0 for x in nm)
cin = {str(t.as_coeff_Mul()[1]): str(sp.factor_terms(t.as_coeff_Mul()[0])) for t in [sp.Add(*[
    t for t in sp.Add.make_args(LR) if any(t.has(pj[("al", o)]) for o in ORD) and not
    any(t.has(pj[(x, o)]) for x in ("a", "be") for o in ORD)])]}
LA = sp.Add(*[t for t in sp.Add.make_args(LR) if any(t.has(pj[("al", o)]) for o in ORD)
              and not any(t.has(pj[(x, o)]) for x in ("a", "be") for o in ORD)])
pol = sp.Poly(LA, *[pj[(x, o)] for x in ("c", "d", "al") for o in ORD])
CHK["L_al"] = {str(m): str(sp.factor_terms(c)) for m, c in zip(pol.monoms(), pol.coeffs())}
print("   L_al (sin a, beta):")
gens = [pj[(x, o)] for x in ("c", "d", "al") for o in ORD]
for m, c in zip(pol.monoms(), pol.coeffs()):
    print("     ", "*".join(str(gens[i]) + ("^2" if e_ == 2 else "") for i, e_ in enumerate(m) if e_), "->",
          sp.factor_terms(c))
json.dump(CHK, open(os.path.join(AQUI, "resultados", "caso_R.json"), "w"), indent=1)
