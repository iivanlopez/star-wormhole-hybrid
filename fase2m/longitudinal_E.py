# -*- coding: utf-8 -*-
"""
Fase 2m: el sector longitudinal en la garganta (caso E).  Tras eliminar F (L_eff = -W^2/4kap + L_rest):
  - at = delta A0 - U a es auxiliar par (sin at_t^2 ni at'^2): masa estatica m_at (energia_E.py);
  - wt = delta A1 + U beta/N^2 es auxiliar impar: el termino wt_t^2 de Maxwell se cancela con el de -W^2/4kap, y queda
    c_t wt^2 con c_t = c_ww - (e_ww)'/2  (e_ww: coeficiente de wt wt').
En el limite local (desacoplo, coeficientes congelados) el par (at, wt) con la parte de Maxwell g da
    w^2 = -(c_t/m_at) k^2 - 2 c_t N^2 / (g r^2 N ...)   ->  sano si m_at > 0 y c_t < 0.
Calcula exactamente m_at(0) y c_t(0) en una garganta simetrica (E_B y las segundas derivadas en la capa de masas) y
comprueba la cancelacion de wt_t^2 y de at'^2.
Uso: python longitudinal_E.py
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
en = pickle.load(open(os.path.join(AQUI, "resultados", "energia_E.pkl"), "rb"))
re = pickle.load(open(os.path.join(AQUI, "resultados", "reduccion_E.pkl"), "rb"))
fe = d["fondo_eq"]
q, Lr, kap, Wj = re["q"], re["Lr"], en["kap"], en["Wj"]
cv = {VV[k]: 0 for k in range(CP.NMAXB)}
Fj = NN[0] ** 2 * q[("a", (0, 1))] + q[("be", (1, 0))]
Leff = sp.expand((-Wj ** 2 / (4 * kap) + Lr).xreplace(cv))
CHK = {}
wt10, wt0, wt1, at1 = q[("wt", (1, 0))], q[("wt", (0, 0))], q[("wt", (0, 1))], q[("at", (0, 1))]
CHK["coef_wt_t^2"] = str(sp.simplify(Leff.coeff(wt10, 2)))
CHK["coef_at'^2"] = str(sp.simplify(Leff.coeff(at1, 2)))
assert CHK["coef_wt_t^2"] == "0" and CHK["coef_at'^2"] == "0"
print("wt_t^2 y at'^2 se cancelan en L_eff: wt y at son auxiliares")
cww = sp.simplify(Leff.coeff(wt0, 2).xreplace({s: 0 for s in q.values()}))
eww = sp.simplify(Leff.coeff(wt0).coeff(wt1))
ct = sp.expand(cww - CP.D(sp.expand(eww)) / 2)
th = {VV[k]: 0 for k in range(9)}
th.update({NN[1]: 0, rr[1]: 0, UU[1]: 0})
G3z = lambda e: e.xreplace({s: 0 for s in e.free_symbols if s.name.startswith("G3")})  # noqa: E731
E = {n: sp.expand(G3z(fe[n].xreplace(th))) for n in fe}
sol = sp.solve([E["A"], E["C"], E["u"]], [NN[2], rr[2], UU[2]], dict=True)[0]
g200 = sp.solve(E["B"], gsym("G2", 0, 0))[0]


def garg(e):
    return sp.factor(sp.simplify(G3z(e).xreplace(th).xreplace(sol).xreplace({gsym("G2", 0, 0): g200})))


X, r0, N0 = sp.symbols("X r0 N0", positive=True)
mt0 = garg(en["mt"])
ct0 = garg(ct)
rpp = garg(sol[rr[2]])
aX = {UU[0]: N0 * sp.sqrt(2 * X), NN[0]: N0, rr[0]: r0}
G2X, G2XX, G2F, G4, G4X, G4XX = (gsym(*a) for a in (("G2", 1, 0), ("G2", 2, 0), ("G2", 0, 1), ("G4", 0, 0),
                                                     ("G4", 1, 0), ("G4", 2, 0)))
Gc = G4 - 2 * X * G4X
Th = G2X * r0 ** 2 + 2 * G4X
for nom, e in (("r''(0)", rpp), ("m_at(0)", mt0), ("c_t(0)", ct0)):
    e = sp.factor(sp.simplify(e.xreplace(aX)))
    CHK[nom] = str(e)
    print(nom, "=", e)
# formas compactas
f_r = -2 * X * Th / (4 * r0 * Gc)
f_m = (G4 * ((G2X + 2 * X * G2XX) * r0 ** 2 + 2 * (G4X + 2 * X * G4XX)) + 4 * X ** 2 * r0 ** 2 *
            (G2X * G4XX - G2XX * G4X)) / (2 * N0 * Gc)
print("r''(0) = -2X Theta/(4 r0 Gcal):", sp.simplify(rpp.xreplace(aX) - f_r) == 0)
print("m_at(0) = [G4((G2X+2X G2XX) r0^2 + 2(G4X+2X G4XX)) + 4X^2 r0^2 (G2X G4XX - G2XX G4X)]/(2 N0 Gcal):",
      sp.simplify(mt0.xreplace(aX) - f_m) == 0)
ct_G4lin = sp.factor(sp.simplify(ct0.xreplace(aX).xreplace({G4XX: 0, gsym("G4", 3, 0): 0})))
CHK["c_t(0)_G4XX=0"] = str(ct_G4lin)
print("c_t(0) con G4_XX = 0:", ct_G4lin)
# c_t(0) en funcion de Theta (G4_XX = 0)
Ths = sp.Symbol("Theta")
ct_Th = sp.factor(sp.simplify(ct_G4lin.xreplace({G2X: (Ths - 2 * G4X) / r0 ** 2})))
CHK["c_t(0)_Theta"] = str(ct_Th)
print("c_t(0) (G4_XX = 0) con G2_X = (Theta - 2 G4_X)/r0^2:", ct_Th)
json.dump(CHK, open(os.path.join(AQUI, "resultados", "longitudinal_E.json"), "w"), indent=1)
