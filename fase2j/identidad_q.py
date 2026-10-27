# -*- coding: utf-8 -*-
"""
Fase 2j: la energia de la familia c = 0 (d = N f) en la teoria quintica de Horndeski, en forma canonica.

-L_rest(0, N f)|_estatico (reduccion_q.py) es una forma cuadratica en (f, f', f'').  Integrando por partes,
    m02 f f'' = -m02 f'^2 - m02' f f' + (m02 f f')',   m12 f' f'' = -(m12'/2) f'^2 + (m12 f'^2/2)',
    b f f' = -(b'/2) f^2 + (b f^2/2)',
queda la forma canonica (unica) A f'^2 + C f^2.  En la teoria cuartica (Fase 2g) A = -2 N^3 Gcal y
C = -2 N^3 Gcal sig^2 - (2 N^3 Gcal sig)' (la identidad puntual -2 N^3 Gcal (f' - sig f)^2), y con f = 1 la energia es
E[N] = int C = -2 int N^3 Gcal sig^2 < 0.  Aqui se calculan A y C con G5 y se comparan con los coeficientes del sector
impar de Horndeski (Kobayashi-Motohashi-Suyama, arXiv:1202.4893), en longitud propia y con X = -phi'^2/2:
    G_KMS/2 = G4 - 2X G4_X + X (Phi' phi' G5_X + G5_phi),   H_KMS/2 = G4 - 2X G4_X + X ((r'/r) phi' G5_X + G5_phi),
    F_KMS/2 = G4 + (1/2) phi' X' G5_X - X G5_phi.
Escribe resultados/identidad_q.json.
Uso: python identidad_q.py
"""
import os
import sys
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_q import Phi, rr, PP, Ee, Xs, D, chorro_en_capa   # noqa: E402


def S(n):
    return sp.Symbol(n, real=True)


R = pickle.load(open(os.path.join(AQUI, "resultados", "reduccion_q.pkl"), "rb"))
q = pickle.load(open(os.path.join(AQUI, "resultados", "q_simbolico.pkl"), "rb"))
horn = {s: 0 for e in q["fondo_eq"].values() for s in e.free_symbols if s.name.startswith(("F4_", "F5_"))}
fe = {k: sp.expand(v.xreplace(horn)) for k, v in q["fondo_eq"].items()}
f = sp.symbols("f0:5", real=True)
mL = sp.expand(R["menosL"])
N = sp.sqrt(Ee)
p1, p2 = PP[1], PP[2]
Xb = -p1 ** 2 / 2
G4, G4X, G5, G5X, G5p = S("G4_00"), S("G4_01"), S("G5_00"), S("G5_01"), S("G5_10")
sig = rr[1] / rr[0] - Phi[1]
RES = {}

# coeficientes m_ij de -L_rest(0, N f) en (f, f', f'')
m = {}
for i in range(3):
    for j in range(i, 3):
        c = mL.coeff(f[i], 1).coeff(f[j], 1) if i != j else mL.coeff(f[i], 2)
        m[(i, j)] = sp.expand(c)
resto = sp.expand(mL - sum(m[k] * f[k[0]] * f[k[1]] for k in m))
assert resto == 0
assert m[(2, 2)] == 0, "aparece f''^2"
A = sp.expand(m[(1, 1)] - m[(0, 2)] - D(m[(1, 2)]) / 2)
b = sp.expand(m[(0, 1)] - D(m[(0, 2)]))
C = sp.expand(m[(0, 0)] - D(b) / 2)
GKMS = G4 - 2 * Xb * G4X + Xb * (Phi[1] * p1 * G5X + G5p)
HKMS = G4 - 2 * Xb * G4X + Xb * (rr[1] / rr[0] * p1 * G5X + G5p)
FKMS = G4 + sp.Rational(1, 2) * p1 * (-p1 * p2) * G5X - Xb * G5p
Gtil = sp.expand(R["Gtil"])
print("Gtil - H_KMS/2 =", sp.simplify(Gtil - HKMS))
for nombre, cand in (("H_KMS/2", HKMS), ("G_KMS/2", GKMS), ("F_KMS/2", FKMS)):
    print(f"A + 2 N^3 ({nombre}) =", sp.factor(sp.simplify(A + 2 * N ** 3 * cand)))
dC = sp.factor(sp.simplify(C + 2 * N ** 3 * GKMS * sig ** 2 + D(2 * N ** 3 * GKMS * sig)))
print("C - [-2 N^3 G_KMS sig^2 - (2 N^3 G_KMS sig)'] =", dC)
cero5 = {s_: 0 for s_ in C.free_symbols if s_.name.startswith('G5_')}
print('   con G5 = 0:', sp.simplify(dC.xreplace(cero5)))
dC2 = sp.factor(sp.simplify(C + 2 * N ** 3 * HKMS * sig ** 2 + D(2 * N ** 3 * HKMS * sig)))
print("C - [-2 N^3 H_KMS sig^2 - (2 N^3 H_KMS sig)'] =", dC2)
RES.update(A=str(sp.factor(A)), dC_GKMS=str(dC), dC_HKMS=str(dC2))
# dependencia en G5 sin derivar (la teoria es invariante bajo G5 -> G5 + cte salvo derivadas totales)
print("dA/dG5_00 =", sp.simplify(sp.diff(A, G5)), "   dC/dG5_00 =", sp.factor(sp.simplify(sp.diff(C, G5))))
RES["dC_dG5"] = str(sp.factor(sp.simplify(sp.diff(C, G5))))
pickle.dump(dict(A=A, C=C, b=b, m=m, GKMS=GKMS, HKMS=HKMS), open(os.path.join(AQUI, "resultados", "identidad_q.pkl"), "wb"))
json.dump(RES, open(os.path.join(AQUI, "resultados", "identidad_q.json"), "w"), indent=1)
