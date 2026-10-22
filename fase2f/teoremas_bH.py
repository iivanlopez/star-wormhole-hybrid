# -*- coding: utf-8 -*-
"""
Identidades exactas de la Fase 2f (teoria cuartica mas alla de Horndeski, fondo estatico esferico con lapso
N = e^Phi, longitud propia l, X_b = -phi'^2/2), comprobadas en SymPy sobre chorros arbitrarios (fuera de la capa de
masas salvo que se indique):

  Notacion:  Gcal = G4 - 2X G4_X - 4X^2 F4,   H = G4 - 2X G4_X,
             Gamma1 = G4 - 4X G4_X - 4X^2 G4_XX - 20 X^2 F4 - 8 X^3 F4_X,
             Gamma2 = phi'(phi'^2 G3_X + 2 G4_phi - 2 phi'^2 G4_phiX),
             E_A, E_B, E_C, E_p: ecuaciones del fondo (generar_bH.py).
  T1 (tercera ecuacion del fondo):   -N^2 E_A + (r^2/2) E_C = [N r Gcal (r' - r Phi')]' - N G4.
  T2 (monopolo, coeficientes de la ligadura y de la accion reducida, monopolo_reduccion.py):
        gam = r (4 r' Gamma1 + r Gamma2)/N,   al1 = -4 r Gcal/N,   al0 = -Phi' al1,
        K = (2/N^2) [ (N xi)' - N G4 ],   xi = -al1 r H/gam = 4 r Gcal H/(4 r' Gamma1 + r Gamma2);
     con T1, en la capa de masas:  N^2 K/2 = W',   W = N [ xi - r Gcal (r' - r Phi') ].
  T3 (potencial y gradiente del monopolo):  con A1 = al1/gam, A0 = al0/gam = -Phi' A1 y los coeficientes de L_rest
        (E_c^rest = k_c c + e0 d + e1 d' + e_tt d_tt,  L_d = k_d d_t^2 - p_d d'^2 - u_d d^2):
        U = u_d - (k_c/2) A0^2 + A0 e0 + (1/2)(k_c A1 A0 - A1 e0 - A0 e1)',   P = p_d + A1 e1 - (k_c/2) A1^2,
        u_d = -(2/r) (N r' Gcal)'.
  T4 (sector impar, L = alfa Q^2 + beta h0^2 + gamma h1^2):
        alfa = N_l Gcal/N,  beta = N_l [lam_l G4/N + 2 E_A]/r^2,  gamma = N_l [-lam_l N Gcal - 2 E_B]/r^2,
        N_l = l(l+1)/(2l+1),  lam_l = (l-1)(l+2);  (l = 2, 3).
Consecuencias (fase2f_resultados.md):  int N^2 K dl = 2 [W] = 6 (M+ + M-) (masas de los dos lados, asintotica RG);
        int U dl = -2 int N Gcal r'^2/r^2 dl + int Phi'^2 (P - p_d - A1 e1) dl - int Phi' A1 e0 dl.
Uso: python teoremas_bH.py      (escribe resultados/teoremas_bH.json)
"""
import os
import sys
import json
import time
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros import Phi, rr, PP, Ee, eP, Xs, D, chorro_en_capa       # noqa: E402

T0 = time.time()
RES = {}


def log(m):
    print(f"[{time.time() - T0:6.1f} s] {m}", flush=True)


d = pickle.load(open(os.path.join(AQUI, "resultados", "bH_simbolico.pkl"), "rb"))
fe, lin, amp, iw = d["fondo_eq"], d["lineales"], d["amp"], d["iw"]
mr = pickle.load(open(os.path.join(AQUI, "resultados", "monopolo_reduccion.pkl"), "rb"))
g_, a1_, a0_ = mr["extra"]
Xb = {Xs: -PP[1] ** 2 / 2}


def S(n):
    return sp.Symbol(n, real=True)


N = sp.sqrt(Ee)
G4, G401, G402, F400, F401 = (S(n) for n in ("G4_00", "G4_01", "G4_02", "F4_00", "F4_01"))
p1 = PP[1]
Gc = G4 + G401 * p1 ** 2 - F400 * p1 ** 4
Hc = G4 + G401 * p1 ** 2
Gam1 = G4 + 2 * p1 ** 2 * G401 - p1 ** 4 * G402 - 5 * p1 ** 4 * F400 + p1 ** 6 * F401
Gam2 = p1 * (p1 ** 2 * S("G3_01") + 2 * S("G4_10") - 2 * p1 ** 2 * S("G4_11"))
EA, EB, EC = (sp.expand(fe[k].xreplace(Xb)) for k in ("A", "B", "C"))

# ------------------------------------------------------------------ T1
t1 = sp.expand(-Ee * EA + rr[0] ** 2 / 2 * EC - (D(N * rr[0] * Gc * (rr[1] - rr[0] * Phi[1])) - N * G4))
RES["T1"] = str(t1)
assert t1 == 0
log("T1: -N^2 E_A + (r^2/2) E_C = [N r Gcal (r' - r Phi')]' - N G4  (identidad exacta)")

# ------------------------------------------------------------------ T2
gam_e, al1_e, al0_e = (sp.expand(mr[k].xreplace(Xb)) for k in ("gam", "al1", "al0"))
assert sp.simplify(gam_e - rr[0] * (4 * rr[1] * Gam1 + rr[0] * Gam2) / N) == 0
assert sp.simplify(al1_e + 4 * rr[0] * Gc / N) == 0
assert sp.simplify(al0_e + Phi[1] * al1_e) == 0
F = a1_[0] * rr[0] * Hc / g_[0]                                    # = -xi con los chorros abstractos
K_cand = 2 * (Phi[1] * F - D(F, mr["extra"]) + 2 * a0_[0] * rr[0] * Hc / g_[0] - G4) / N
dK = sp.expand(mr["K"].xreplace(Xb) - K_cand)
assert dK == 0
# con al0 = -Phi' al1:  K = (2/N^2) [ (N xi)' - N G4 ],  xi = -F
xi = -F
K_T2 = 2 / N ** 2 * (D(N * xi, mr["extra"]) - N * G4)
dK2 = sp.simplify((K_cand - K_T2).subs(a0_[0], -Phi[1] * a1_[0]))
assert dK2 == 0
RES["T2"] = "gam, al1, al0 y K = (2/N^2)[(N xi)' - N G4] exactos"
log("T2: gam = r(4 r' Gamma1 + r Gamma2)/N, al1 = -4 r Gcal/N, al0 = -Phi' al1, K = (2/N^2)[(N xi)' - N G4] "
    "(identidades exactas)")
# en la capa de masas: N^2 K/2 = W' (con los chorros explicitos)
W = N * (4 * rr[0] * Gc * Hc / (4 * rr[1] * Gam1 + rr[0] * Gam2) - rr[0] * Gc * (rr[1] - rr[0] * Phi[1]))
xi_e = 4 * rr[0] * Gc * Hc / (4 * rr[1] * Gam1 + rr[0] * Gam2)
K_exp = 2 / N ** 2 * (D(N * xi_e) - N * G4)
for sem in range(3):
    v = chorro_en_capa(fe, 600 + sem, orden=3)
    x = v.ev(sp.together(Ee * K_exp / 2 - D(W)))
    assert x == 0, x
RES["T2_capa"] = "N^2 K/2 = W' exacta en 3 chorros en la capa de masas"
log("T2 en la capa de masas: N^2 K/2 = W',  W = N [xi - r Gcal (r' - r Phi')]  (exacta, 3 chorros)")

# ------------------------------------------------------------------ T3
z = {amp[(v_, k)]: 0 for v_ in ("a", "be", "c", "q") for k in range(6)}
Ed = sp.expand(lin["d"].xreplace(z))
u_d = sp.expand(-Ed.coeff(amp[("d", 0)]).subs(iw, 0) / 2).xreplace(Xb)
p_d = sp.expand(Ed.coeff(amp[("d", 2)]).subs(iw, 0) / 2).xreplace(Xb)
zc = {amp[(v_, k)]: 0 for v_ in ("a", "be", "q") for k in range(6)}
Ec = sp.expand(lin["c"].xreplace(zc))
e0 = sp.expand(Ec.coeff(amp[("d", 0)]).subs(iw, 0)).xreplace(Xb)
e1 = sp.expand(Ec.coeff(amp[("d", 1)])).xreplace(Xb)
k_c = sp.expand(mr["k_c"].xreplace(Xb))
assert sp.expand(u_d + 2 * D(N * rr[1] * Gc) / rr[0]) == 0
log("T3: u_d = -(2/r)(N r' Gcal)'  (identidad exacta)")
A1 = a1_[0] / g_[0]
A0 = a0_[0] / g_[0]
U_T3 = u_d - k_c / 2 * A0 ** 2 + A0 * e0 + sp.Rational(1, 2) * D(k_c * A1 * A0 - A1 * e0 - A0 * e1, mr["extra"])
P_T3 = p_d + A1 * e1 - k_c / 2 * A1 ** 2
dP = sp.simplify(mr["P"].xreplace(Xb) - P_T3)
dU = sp.simplify(sp.together(mr["U"].xreplace(Xb) - U_T3))
RES["T3_P"] = str(dP)
RES["T3_U"] = str(dU)
assert dP == 0 and dU == 0
log("T3: U = u_d - (k_c/2) A0^2 + A0 e0 + (1/2)(k_c A1 A0 - A1 e0 - A0 e1)',  P = p_d + A1 e1 - (k_c/2) A1^2  "
    "(identidades exactas con los chorros abstractos de gam, al1, al0)")

# ------------------------------------------------------------------ T4
import simbolico_num as SN                                      # noqa: E402
for LL in (2, 3):
    Nl = sp.Rational(LL * (LL + 1), 2 * LL + 1)
    lam = (LL - 1) * (LL + 2)
    al, be, ga = (SN.IMPAR_EXPR[LL][k].subs(eP, N) for k in ("alfa", "beta", "gamma"))
    assert sp.simplify(al - Nl * Gc / N) == 0
    assert sp.simplify(be - Nl * (lam * G4 / N + 2 * EA) / rr[0] ** 2) == 0
    assert sp.simplify(ga - Nl * (-lam * N * Gc - 2 * EB) / rr[0] ** 2) == 0
RES["T4"] = "alfa = N_l Gcal/N, beta = N_l[lam G4/N + 2 E_A]/r^2, gamma = N_l[-lam N Gcal - 2 E_B]/r^2 (l = 2, 3)"
log("T4: sector impar alfa = N_l Gcal/N, beta = N_l [lam_l G4/N + 2 E_A]/r^2, gamma = N_l [-lam_l N Gcal - 2 E_B]/r^2 "
    "(exactas, l = 2, 3)")
json.dump(RES, open(os.path.join(AQUI, "resultados", "teoremas_bH.json"), "w"), indent=1)
log("todas las identidades comprobadas; resultados/teoremas_bH.json")
