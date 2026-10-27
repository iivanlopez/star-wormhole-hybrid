# -*- coding: utf-8 -*-
"""
Fase 2k: el monopolo en las variables ADM del corte phi = cte (gauge unitario temporal, dphi = 0).

Con tau = t + psi/qt (phi = qt tau) el fondo es estacionario con desplazamiento:
    g_tt^tau = -e^{2Phi},  N_l = g_{tau l} = E s,  h_ll = 1 - E s^2,  N_tau^2 = E/(1 - E s^2),   s = psi'/qt.
Perturbaciones (en las variables estaticas a, beta, c, d con dA = 2E a, db = beta, dB = 2c, dC = 2 r d):
    beta~ = dN_l = beta + 2 E s a,      h = dh_ll/2 = c - s beta - s^2 E a,
    n = d(N_tau^2)/2 = E a + E s beta~/(1 - E s^2) - E^2 s^2 h/(1 - E s^2)^2 .
Ecuaciones de Euler-Lagrange en las variables nuevas (cambio puntual u = M v:  E_v = M^T E_u).
Comprobaciones en chorros exactos de la capa de masas:
  (1) E_n no contiene derivadas de n ni de beta~ (espaciales ni temporales): el lapso es auxiliar algebraico;
  (2) E_beta~ no contiene derivadas de beta~ ni de n: el desplazamiento tambien es algebraico en l = 0;
  (3) que derivadas de h y d aparecen.
Uso: python adm_t.py [G5]
"""
import os
import sys
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_t import Phi, rr, PP, Ee, Xs, qt, D, chorro_en_capa   # noqa: E402

CON_G5 = len(sys.argv) > 1 and sys.argv[1] == "G5"
NOMBRE = "t5" if CON_G5 else "t"
dat = pickle.load(open(os.path.join(AQUI, "resultados", f"{NOMBRE}_simbolico.pkl"), "rb"))
fe, J, lin, iw, amp = dat["fondo_eq"], dat["J"], dat["lineales"], dat["iw"], dat["amp"]
NM = 5
n_ = sp.symbols(f"n0:{NM + 3}", real=True)
B_ = sp.symbols(f"B0:{NM + 3}", real=True)
h_ = sp.symbols(f"h0:{NM + 3}", real=True)
d_ = sp.symbols(f"dd0:{NM + 3}", real=True)
AMPS = (n_, B_, h_, d_)
s = PP[1] / qt
H = 1 - Ee * s ** 2
a0 = (n_[0] - Ee * s * B_[0] / H + Ee ** 2 * s ** 2 * h_[0] / H ** 2) / Ee
be0 = B_[0] - 2 * Ee * s * a0
c0 = h_[0] + s * be0 + s ** 2 * Ee * a0
assert sp.simplify(c0 - s * be0 - s ** 2 * Ee * a0 - h_[0]) == 0
base = {"a": a0, "be": be0, "c": c0, "d": d_[0], "q": sp.Integer(0)}
sub = {}
for v, e in base.items():
    e = sp.together(e)
    for k in range(NM + 1):
        sub[amp[(v, k)]] = e
        e = D(e, AMPS)
var = {"n": n_[0], "B": B_[0], "h": h_[0]}


def ecuaciones_en(ch):
    """E_n, E_B, E_h, E_d evaluadas en el chorro (polinomios en amplitudes e iw con coeficientes racionales)."""
    subn = {k: ch.ev(e) for k, e in sub.items()}
    Eu = {v: sp.expand(ch.ev(lin[v]).xreplace(subn)) for v in ("a", "be", "c", "d")}
    out = {k: sp.expand(sum(Eu[u] * ch.ev(sp.diff(base[u], x)) for u in ("a", "be", "c"))) for k, x in var.items()}
    out["d"] = Eu["d"]
    return out


def coeficientes(x):
    out = {}
    for fam, nom in ((n_, "n"), (B_, "B"), (h_, "h"), (d_, "d")):
        for k, sym in enumerate(fam):
            c = sp.expand(x.coeff(sym))
            for p in range(4):
                cp = c.coeff(iw, p) if p > 0 else c.subs(iw, 0)
                if cp != 0:
                    out[f"{nom}{k}*iw^{p}"] = cp
    return out


if __name__ == "__main__":
    RES = {}
    for sem in range(2):
        ch = chorro_en_capa(fe, J, 900 + sem, orden=7, cero=() if CON_G5 else ("G5",))
        Enew = ecuaciones_en(ch)
        for k in ("n", "B", "h", "d"):
            cs = coeficientes(Enew[k])
            RES.setdefault(k, []).append(sorted(cs))
            print(f"chorro {sem}: E_{k}: {sorted(cs)}", flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", f"adm_{NOMBRE}.json"), "w"), indent=1)
