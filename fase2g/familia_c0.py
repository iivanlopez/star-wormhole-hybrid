# -*- coding: utf-8 -*-
"""
Fase 2g, parte II: la familia de pruebas con c = 0 (delta g_ll = 0) en el sector l = 0 con fluido.

Con el fluido, la ligadura de momento es W = gam c + al1 N (d/N)' = -2 kap_b xi (xi: desplazamiento del fluido).  Para
d = N h(l) y c = 0 basta mover el fluido con xi = -al1 N h'/(2 kap_b) = 4 N Gcal h'/(r (rho~ + p~)) (dentro de la
estrella).  En esta familia NO aparece ninguna funcion libre de la teoria salvo Gcal y G4:
    T[h] = int [ kap_b xi_t^2 + k_d d_t^2 ]          = int [ 8 N Gcal^2 h_t'^2/(rho~ + p~) - 2 G4 N h_t^2 ]
    E[h] = int [ Q^2/(4 kap_a) - L_rest(0, N h) ]     = int [ c_s^2 Q^2/(2 N r^2 (rho~ + p~)) - 2 N^3 Gcal (h' - sig h)^2 ]
    Q    = (4 N^2 r Gcal h')' + 2 N^2 r (rho~ + p~) h,      sig = r'/r - Phi' = (ln(r/N))'.
Fuera de la estrella (sin fluido) W = 0 exige h = cte; con h = 1 se recupera la prueba d = N, xi = 0.
Este script comprueba simbolicamente (exacto, fuera de la capa de masas) que
    -L_rest(c = 0, d = N f)|_estatico = -2 N^3 Gcal (f' - sig f)^2      (sin derivada total)
y que la parte cinetica de L_rest(c = 0, d) es k_d d_t^2 con k_d = -2 G4/N.
Uso: python familia_c0.py
"""
import os
import sys
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
from chorros import Phi, rr, PP, Ee, Xs, D       # noqa: E402

mr = pickle.load(open(os.path.join(F2F, "resultados", "monopolo_reduccion.pkl"), "rb"))
Xb = {Xs: -PP[1] ** 2 / 2}
N = sp.sqrt(Ee)


def S(n):
    return sp.Symbol(n, real=True)


L = sp.expand(mr["L_rest"].xreplace(Xb))
Ld = sp.expand(L.xreplace({S(n): 0 for n in ("c_00", "c_01", "c_10", "c_20")}))
f = sp.symbols("f0:4", real=True)                      # f, f', f'', f'''


def Df(e):
    out = D(e)
    for k in range(3):
        out += sp.diff(e, f[k]) * f[k + 1]
    return sp.expand(out)


d0 = N * f[0]
d1 = Df(d0)
d2 = Df(d1)
menosL = sp.expand(-Ld.xreplace({S("d_00"): d0, S("d_01"): d1, S("d_02"): d2, S("d_10"): 0, S("d_20"): 0}))
G4, G401, F400 = S("G4_00"), S("G4_01"), S("F4_00")
Gc = G4 + G401 * PP[1] ** 2 - F400 * PP[1] ** 4
sig = rr[1] / rr[0] - Phi[1]
dif = sp.simplify(menosL + 2 * N ** 3 * Gc * (f[1] - sig * f[0]) ** 2)
assert dif == 0, dif
cin = sp.expand(Ld.coeff(S("d_10"), 2))
assert sp.simplify(cin + 2 * G4 / N) == 0
assert sp.expand(Ld - cin * S("d_10") ** 2).has(S("d_10")) is False
print("-L_rest(c = 0, d = N f) = -2 N^3 Gcal (f' - sig f)^2  (exacta, puntual) y cinetica k_d = -2 G4/N:  comprobado")
json.dump({"familia_c0": "-L_rest(0, N f) = -2 N^3 Gcal (f' - sig f)^2 exacta; k_d = -2 G4/N"},
          open(os.path.join(AQUI, "resultados", "familia_c0.json"), "w"), indent=1)
