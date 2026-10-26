# -*- coding: utf-8 -*-
"""
Fase 2g, parte IV: la garganta cargada (reconstruccion_carga.GargantaCargada) cerca de la garganta.

Con carga, en una garganta simetrica (gam(0) = 0, H(0) = 0, k_c(0) = 0 por la regularidad de P) el coeficiente e0 de la
ecuacion de c no se anula en la garganta:
    e0_tot(0) = 2 N0 r0 S(0) + 2 q2 N0/r0^3 = 4 q2 N0/r0^3     (S(0) = q2/r0^4 por las ligaduras con carga),
de modo que la ecuacion de c en l = 0 se reduce a e0(0) d(0) = 0 (d = delta r: el radio de la garganta no puede variar)
y el potencial del problema reducido es singular, U ~ -C/l^2 con C/P(0) >> 1/4 (caida al centro).  Este script lo
muestra para la garganta de referencia de la Fase 2f con q2 = 0, 0.5, 1, 1.5, 2 y da la energia de la prueba d = N con carga,
    E[N] = int [ -2 N^3 Gcal sig^2 + 4 q2 N^3/r^4 ] dl  (formula; la prueba no es admisible si e0(0) != 0).
Escribe resultados/carga_garganta.json.
Uso: python carga_garganta.py
"""
import os
import sys
import json
import pickle
import numpy as np
import sympy as sp
from scipy.integrate import simpson

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
sys.path.insert(0, AQUI)
import simbolico_num as SN             # noqa: E402
import reconstruccion_carga as RC      # noqa: E402

BASE = dict(n=2, b=4.0, frac=0.9, wdb=0.5, M=2.0, a=2.0, wK=10.0, wPK=1.5, cm2=0.8)
cf = pickle.load(open(os.path.join(AQUI, "resultados", "coeficientes_fluido.pkl"), "rb"))
_P0 = sp.Symbol("Pf0", real=True)
E0 = SN.compilar(SN.con_Xb(cf["e0"]).subs(_P0, 0))
KC = SN.compilar(SN.con_Xb(cf["k_c"]).subs(_P0, 0))


def construir(q2):
    x = BASE
    n = x["n"]
    G40 = 2 * n / (2 * n - 1)
    return RC.GargantaCargada(q2=q2, phi1=1.0, b=x["b"], n=n, d0=x["frac"] * G40, wd=x["wdb"] * x["b"], M=x["M"],
                              a=x["a"], Lmax=4000.0, wK=x["wK"], wP=x["wPK"] * x["wK"], cm2=x["cm2"])


if __name__ == "__main__":
    out = []
    for q2 in (0.0, 0.5, 1.0, 1.5, 2.0):
        try:
            g = construir(q2)
        except AssertionError:
            print(f"q2 = {q2}: el fondo no existe (la carga cierra la garganta: r -> 0)", flush=True)
            out.append(dict(q2=q2, fondo="no existe"))
            continue
        l = np.array([1e-3, 3e-3, 1e-2, 3e-2])
        V = g.V(l)
        N, r = V["N"].deriv(0), V["r0"].deriv(0)
        e0 = E0(V).deriv(0) + 2 * q2 * N / r ** 3
        kc = KC(V).deriv(0) + q2 * N / r ** 2
        m = g._mono(V)
        U, P = m["U"].deriv(0), m["P"].deriv(0)
        ls = np.concatenate([np.linspace(1e-4, 30, 6001), np.geomspace(30, 3000, 2001)[1:]])
        B = g._fondo_basico(ls)
        Nl, rl, r1l = B["N"].deriv(0), B["r0"].deriv(0), B["r1"].deriv(0)
        Gc = g.perfiles(ls)["Gc"]
        sig = r1l / rl - B["Phi"].deriv(1)
        EN = simpson(-2 * Nl ** 3 * Gc * sig ** 2 + 4 * q2 * Nl ** 3 / rl ** 4, x=ls)
        r_ = dict(q2=q2, e0_tot=[float(v) for v in e0], e0_teo=float(4 * q2 * g.N0 / g.r0 ** 3), k_c=[float(v) for v in kc],
                  U_l2=[float(v) for v in U * l ** 2], P=[float(v) for v in P], EN_formula=float(EN))
        out.append(r_)
        print(f"q2 = {q2}: e0(l->0) = {e0[0]:+.5f} (4 q2 N0/r0^3 = {r_['e0_teo']:+.5f}), k_c(l->0) = {kc[0]:+.2e}, "
              f"U l^2 = {np.round(U * l ** 2, 4)}, P(0) = {P[0]:.4f};  E[N] (formula) = {EN:+.4f}", flush=True)
    json.dump(out, open(os.path.join(AQUI, "resultados", "carga_garganta.json"), "w"), indent=1)
