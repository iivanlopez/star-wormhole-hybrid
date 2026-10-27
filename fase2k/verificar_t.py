# -*- coding: utf-8 -*-
"""
Fase 2k: comprobaciones de las ecuaciones del fondo y de la accion cuadratica con phi = qt t + psi(l).
  1. Solucion "stealth" de Babichev-Charmousis (arXiv:1312.3204, JHEP 08 (2014) 106): en la teoria
        zeta R + beta G^{mn} phi_m phi_n   (equivale a G4 = zeta + beta X, G2 = G3 = F4 = 0),
     Schwarzschild h = 1 - mu/r con phi = qt t + psi, psi_r = +- qt sqrt(mu/r)/h (X = qt^2/2 constante) cumple
     E_A = E_B = E_C = J = 0 (40 cifras; en longitud propia psi' = sqrt(h) psi_r).
  2. Modos gauge exactos en chorros de la capa de masas (aritmetica racional exacta):
       xi^t = zeta(l) e^{i w t}:  a = i w zeta, beta = -e^{2Phi} zeta', dphi = qt zeta;
       xi^l = xi(l) e^{i w t}:    a = Phi' xi, beta = i w xi, c = xi', d = r' xi, dphi = psi' xi;
     las cinco ecuaciones lineales se anulan.
Uso: python verificar_t.py [G5]
"""
import os
import sys
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_t import Phi, rr, PP, Ee, Xs, qt, D, chorro_en_capa, gsym   # noqa: E402

CON_G5 = len(sys.argv) > 1 and sys.argv[1] == "G5"
NOMBRE = "t5" if CON_G5 else "t"
d = pickle.load(open(os.path.join(AQUI, "resultados", f"{NOMBRE}_simbolico.pkl"), "rb"))
fe, J, lin, iw, amp = d["fondo_eq"], d["J"], d["lineales"], d["iw"], d["amp"]
RES = {}

# ------------------------------------------------------------------ 1. stealth Schwarzschild
mu, rv, ze, be, q0 = sp.symbols("mu rv zeta beta q0", positive=True)
h = 1 - mu / rv
sq = sp.sqrt(h)
dl = lambda e: sp.diff(e, rv) * sq                      # noqa: E731
res_st = {}
for signo in (1, -1):
    ePhi = [sp.log(h) / 2]
    er = [rv]
    ep = [None, signo * q0 * sp.sqrt(mu / rv) / sq]
    for k in range(1, 6):
        ePhi.append(dl(ePhi[-1]))
        er.append(dl(er[-1]))
    for k in range(2, 7):
        ep.append(dl(ep[-1]))
    for (m_, r_, z_, b_, qq) in ((1, sp.Rational(27, 10), 1, sp.Rational(3, 10), sp.Rational(7, 10)),
                                 (2, sp.Rational(53, 10), sp.Rational(3, 2), -sp.Rational(11, 10), sp.Rational(13, 10))):
        pt = {mu: m_, rv: r_, q0: qq}
        num = {Ee: sp.N(h.subs(pt), 50), rr[0]: sp.N(r_, 50), qt: qq}
        num.update({Phi[k]: sp.N(ePhi[k].subs(pt), 50) for k in range(1, 6)})
        num.update({rr[k]: sp.N(er[k].subs(pt), 50) for k in range(1, 6)})
        num.update({PP[k]: sp.N(ep[k].subs(pt), 50) for k in range(1, 7)})
        Xn = sp.N(qq ** 2 / (2 * num[Ee]) - num[PP[1]] ** 2 / 2, 50)
        assert abs(Xn - qq ** 2 / 2) < 1e-40
        num[Xs] = Xn
        for (nom, i, j), s in d["GS"].items():
            val = 0
            if nom == "G4":
                val = (z_ + b_ * Xn) if j == 0 else (b_ if j == 1 else 0)
            num[s] = val
        for n, e in list(fe.items()) + [("J", J)]:
            res_st[n] = max(res_st.get(n, 0.0), abs(float(sp.N(e.xreplace(num), 40))))
RES["stealth_residuos"] = res_st
print("1. stealth Schwarzschild (Babichev-Charmousis), 2 signos x 2 puntos:", {k: f"{v:.1e}" for k, v in res_st.items()},
      flush=True)
assert max(res_st.values()) < 1e-35

# control: con G4 = zeta + beta X y la misma metrica, psi_r con otra amplitud no es solucion
num2 = dict(num)
num2[PP[1]] = num[PP[1]] * sp.Rational(11, 10)
num2[PP[2]] = num[PP[2]] * sp.Rational(11, 10)
num2[PP[3]] = num[PP[3]] * sp.Rational(11, 10)
num2[Xs] = qt.subs(num) ** 2 / (2 * num[Ee]) - num2[PP[1]] ** 2 / 2
for (nom, i, j), s in d["GS"].items():
    num2[s] = 0
    if nom == "G4":
        num2[s] = (z_ + b_ * num2[Xs]) if j == 0 else (b_ if j == 1 else 0)
ctrl = max(abs(float(sp.N(e.xreplace(num2), 20))) for e in list(fe.values()) + [J])
RES["stealth_control"] = ctrl
print(f"   control (psi' x 1.1): residuo maximo {ctrl:.2e} (no nulo)", flush=True)

# ------------------------------------------------------------------ 2. modos gauge en la capa de masas
z = sp.symbols("z0:9", real=True)
x = sp.symbols("x0:9", real=True)
NM = 5


def amps_de(base):
    out = {}
    for v, e in base.items():
        for k in range(NM + 1):
            out[amp[(v, k)]] = e
            e = D(e, (z, x))
    return out


gauge = {
    "xi_t": {"a": iw * z[0], "be": -Ee * z[1], "c": 0, "d": 0, "q": qt * z[0]},
    "xi_l": {"a": Phi[1] * x[0], "be": iw * x[0], "c": x[1], "d": rr[1] * x[0], "q": PP[1] * x[0]},
}
RES["gauge"] = {}
for sem in range(2):
    v = chorro_en_capa(fe, J, 500 + sem, orden=7, cero=() if CON_G5 else ("G5",))
    for k in range(9):
        v[z[k]] = v.r()
        v[x[k]] = v.r()
    v[iw] = v.r()
    for nom, base in gauge.items():
        sub = amps_de(base)
        for e in lin:
            val = v.ev(lin[e].xreplace(sub))
            RES["gauge"][f"{nom}_E{e}_{sem}"] = str(val)
            assert val == 0, (nom, e, val)
print("2. modos gauge xi^t, xi^l: las 5 ecuaciones lineales se anulan exactamente en 2 chorros de la capa de masas",
      flush=True)
json.dump(RES, open(os.path.join(AQUI, "resultados", f"verificar_{NOMBRE}.json"), "w"), indent=1, default=str)
