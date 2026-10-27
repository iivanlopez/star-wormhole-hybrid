# -*- coding: utf-8 -*-
"""
Fase 2j: las gargantas de Bakopoulos-Charmousis-Kanti (arXiv:2111.09857) vistas en el marco de la semilla.

La garganta es la imagen disforme g = gbar - D(Xbar) dphi dphi del exterior r > r0 del agujero negro de Lu-Pang
(Horndeski: G2 = 8 al X^2, G3 = -8 al X, G4 = 1 + 4 al X, G5 = -4 al ln|X|),
    hbar = 1 + (r^2/2al)(1 - sqrt(1 + 8 al M/r^3)),   dphi/dlbar = (sqrt(hbar) - 1)/r   (longitud propia de la semilla),
con la garganta en hbar(r0) = (1 - lam)^2.  Con phi estatico la transformacion solo cambia g_ll: el lapso es el mismo y,
con c = 0, delta phi = 0 (y beta = 0), delta X = 0 a primer y segundo orden.  La prueba d = N f, c = 0 del agujero de
gusano es la prueba d = N f, cbar = 0 de la semilla en r > r0 (dos copias) y su energia es la de la semilla con la
densidad tal cual (sin integrar por partes: en r0 las dos copias no se pegan de forma lisa en el marco de la semilla):
    E[N] = 2 int_{r0}^inf m00 dlbar,     m00 = -L_rest(0, N)|_{f = 1}.
Por identidad_q.py, m00 = C + b'/2 con C = -2 N^3 Gcal sig^2 - (2 N^3 Gcal sig)', Gcal = G4 - 2X G4_X + X Phi' phi' G5_X,
de modo que
    E[N] = 2 [ -2 int_{r0}^inf N^3 Gcal sig^2 dlbar + 2 N^3 Gcal sig (r0) - b(r0)/2 ]      (volumen + borde).
Ademas: K del problema reducido (sin fantasmas: K > 0) y P (sin inestabilidad de gradiente: P > 0) en r > r0.
Uso: python bck_semilla.py
"""
import os
import sys
import json
import math
import pickle
import numpy as np
import sympy as sp
from scipy.optimize import brentq
from scipy.integrate import quad

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_q import Phi, rr, PP, Ee, Xs, D        # noqa: E402

R = pickle.load(open(os.path.join(AQUI, "resultados", "reduccion_q.pkl"), "rb"))
I = pickle.load(open(os.path.join(AQUI, "resultados", "identidad_q.pkl"), "rb"))
al_, M_, rv = sp.symbols("alpha M rv", positive=True)
h = 1 + rv ** 2 / (2 * al_) * (1 - sp.sqrt(1 + 8 * al_ * M_ / rv ** 3))
sq = sp.sqrt(h)


def dl(e):
    return sp.diff(e, rv) * sq


ePhi = [None, sp.diff(h, rv) / (2 * sq)]
er = [rv, sq]
ep = [None, (sq - 1) / rv]
for k in range(2, 8):
    ePhi.append(dl(ePhi[-1]))
    er.append(dl(er[-1]))
    ep.append(dl(ep[-1]))
GLP = {"G2": 8 * al_ * Xs ** 2, "G3": -8 * al_ * Xs, "G4": 1 + 4 * al_ * Xs, "F4": sp.Integer(0),
       "G5": -4 * al_ * sp.log(-Xs), "F5": sp.Integer(0)}


def a_fondo(expr):
    """Sustituye en una expresion de los chorros el fondo de Lu-Pang (funcion de rv, alpha, M)."""
    rep = {}
    for s in expr.free_symbols:
        n = s.name
        if "_" in n and n.split("_")[0] in GLP and len(n.split("_")[1]) == 2:
            nom, ij = n.split("_")
            i, j = int(ij[0]), int(ij[1])
            rep[s] = 0 if i > 0 else sp.diff(GLP[nom], Xs, j)
    e = expr.xreplace(rep)
    e = e.xreplace({Phi[k]: ePhi[k] for k in range(1, 8)}).xreplace({rr[k]: er[k] for k in range(1, 8)})
    e = e.xreplace({PP[k]: ep[k] for k in range(1, 8)}).xreplace({rr[0]: rv, Ee: h})
    return e.subs(Xs, -ep[1] ** 2 / 2)


fs = sp.symbols("f0:5", real=True)
m00 = a_fondo(R["menosL"].xreplace({fs[0]: 1, fs[1]: 0, fs[2]: 0, fs[3]: 0, fs[4]: 0}))
N = sp.sqrt(Ee)
sig = rr[1] / rr[0] - Phi[1]
Gc = a_fondo(I["GKMS"])
vol = a_fondo(-2 * N ** 3 * I["GKMS"] * sig ** 2)
borde = a_fondo(2 * N ** 3 * I["GKMS"] * sig - I["b"] / 2)
# coeficientes del problema reducido (con las derivadas de gam, al1, al0 en el fondo)
g_, a1_, a0_ = R["extra"]
sub = {}
for simb, expr in ((g_, R["gam"]), (a1_, R["al1"]), (a0_, R["al0"])):
    e = a_fondo(expr)
    for k in range(4):
        sub[simb[k]] = e
        e = dl(e)
Kf = a_fondo(R["K"]).xreplace(sub)
Pf = a_fondo(R["P"]).xreplace(sub)
f_m00 = sp.lambdify((rv, al_, M_), m00 / sq, "mpmath")      # por unidad de r: dlbar = dr/sqrt(h)
f_vol = sp.lambdify((rv, al_, M_), vol / sq, "mpmath")
f_borde = sp.lambdify((rv, al_, M_), borde, "mpmath")
f_G = sp.lambdify((rv, al_, M_), Gc, "mpmath")
f_K = sp.lambdify((rv, al_, M_), Kf, "mpmath")
f_P = sp.lambdify((rv, al_, M_), Pf, "mpmath")
f_h = sp.lambdify((rv, al_, M_), h, "mpmath")
f_sig = sp.lambdify((rv, al_, M_), a_fondo(sig), "mpmath")


def r0_de(al, M, lam):
    rp = M + math.sqrt(max(M * M - al, 0.0))
    return brentq(lambda r: float(f_h(r, al, M)) - (1 - lam) ** 2, rp * (1 + 1e-12) + 1e-12, 1e6)


RES = []
for al, M, lam in ((0.01, 1.0, 0.5), (0.1, 1.0, 0.5), (0.5, 1.0, 0.5), (0.9, 1.0, 0.5), (0.5, 1.0, 0.2),
                   (0.5, 1.0, 0.8), (0.5, 1.0, 0.95), (0.2, 1.0, 0.05), (0.99, 1.0, 0.3)):
    r0 = r0_de(al, M, lam)
    r0_bck = (M + math.sqrt(max(M * M - al * lam ** 3 * (2 - lam) ** 3, 0))) / (lam * (2 - lam))
    E = 2 * float(quad(lambda r: float(f_m00(r, al, M)), r0, np.inf, limit=500)[0])
    V = 2 * float(quad(lambda r: float(f_vol(r, al, M)), r0, np.inf, limit=500)[0])
    B = 2 * float(f_borde(r0, al, M))
    rs = np.concatenate([[r0 * (1 + 1e-9)], np.geomspace(r0 * (1 + 1e-6), 400 * r0, 600)])
    G = np.array([float(f_G(r, al, M)) for r in rs])
    K = np.array([float(f_K(r, al, M)) for r in rs])
    P = np.array([float(f_P(r, al, M)) for r in rs])
    s0 = float(f_sig(r0, al, M))
    RES.append(dict(alpha=al, M=M, lam=lam, r0=r0, r0_BCK=r0_bck, E_N=E, volumen=V, borde=B, sig_r0=s0,
                    G_min=float(G.min()), K_min=float(K.min()), K_max=float(K.max()), P_min=float(P.min()),
                    r_photon=float(rs[np.argmin(np.abs([float(f_sig(r, al, M)) for r in rs]))])))
    print(f"al {al} M {M} lam {lam}: r0 {r0:.5f} (BCK {r0_bck:.5f}); E[N] = {E:+.5e} = volumen {V:+.5e} + borde {B:+.5e} "
          f"(dif {E - V - B:+.1e}); sig(r0) {s0:+.4f}; Gcal en [{G.min():.4f}, {G.max():.4f}]; K en [{K.min():+.3e}, "
          f"{K.max():+.3e}]; min P {P.min():+.3e}", flush=True)
json.dump(RES, open(os.path.join(AQUI, "resultados", "bck_semilla.json"), "w"), indent=1)
