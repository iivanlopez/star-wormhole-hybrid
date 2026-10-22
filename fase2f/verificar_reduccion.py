# -*- coding: utf-8 -*-
"""
Comprobacion de extremo a extremo de la accion reducida del monopolo (monopolo_reduccion.py) sobre un fondo
reconstruido (reconstruccion.Garganta, en la capa de masas):

  1. se integra el sistema lineal completo l = 0 en el gauge beta = c = 0 (bH_gen.lineal: E_a, E_d, E_q dan
     (a'', d'', q''); E_beta y E_c son ligaduras que se imponen en el dato inicial) para un w dado;
  2. se forma la variable del gauge unitario  Psi = d - (r'/phi') q  (invariante bajo el modo gauge zeta);
  3. se comprueba que Psi cumple  -(P Psi')' + U Psi = w K Psi  y que las ligaduras se propagan.
Tambien se comprueba la identidad N^2 K/2 = W' con W = N [xi - r Gcal (r' - r Phi')] y los limites de W (masas).
Uso: python verificar_reduccion.py
"""
import os
import sys
import json
import math
import numpy as np
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import bH_gen                     # noqa: E402
import busqueda as B              # noqa: E402

X = dict(n=2, b=4.0, frac=0.9, wdb=0.5, M=2.0, a=2.0, wK=10.0, wPK=1.5, cm2=0.8)
g = B.construir(X)


def datos(l, w):
    V = g.V(np.array([l]))
    G = {k: float(V[k].c[0][0]) for k in bH_gen.SIMB_G}
    Ph = [float(V["Phi"].deriv(0)[0]), float(V["Phi1"].deriv(0)[0]), float(V["Phi2"].deriv(0)[0])]
    r_ = [float(V[f"r{k}"].deriv(0)[0]) for k in range(3)]
    P_ = [float(V[f"P{k}"].deriv(0)[0]) for k in range(3)]
    return bH_gen.lineal(G, Ph, r_, P_, w)


def ec(cs, ecu, vals):
    return sum(c * vals.get(k, 0.0) for k, c in cs.get(ecu, {}).items())


def rhs(l, y, w):
    a1, d0, d1, q0, q1 = y
    cs = datos(l, w)
    base = {("a", 0): 0.0, ("a", 1): a1, ("d", 0): d0, ("d", 1): d1, ("q", 0): q0, ("q", 1): q1}
    M = np.zeros((3, 3))
    b = np.zeros(3)
    for i, e in enumerate(("a", "d", "q")):
        b[i] = ec(cs, e, base)
        for j, v in enumerate(("a", "d", "q")):
            M[i, j] = cs.get(e, {}).get((v, 2), 0.0)
    a2, d2, q2 = np.linalg.solve(M, -b)
    return [a2, d1, d2, q1, q2]


def ligaduras(l, y, w):
    a1, d0, d1, q0, q1 = y
    cs = datos(l, w)
    vals = {("a", 0): 0.0, ("a", 1): a1, ("d", 0): d0, ("d", 1): d1, ("q", 0): q0, ("q", 1): q1}
    return ec(cs, "be", vals), ec(cs, "c", vals)


def dato_inicial(l, w, d0, q0, q1):
    cs = datos(l, w)
    be = cs["be"]
    d1 = -(be.get(("d", 0), 0) * d0 + be.get(("q", 0), 0) * q0 + be.get(("q", 1), 0) * q1) / be[("d", 1)]
    c = cs["c"]
    a1 = -(c.get(("d", 0), 0) * d0 + c.get(("d", 1), 0) * d1 + c.get(("q", 0), 0) * q0
           + c.get(("q", 1), 0) * q1) / c[("a", 1)]
    return [a1, d0, d1, q0, q1]


RES = {}
l0, l1 = 1.0, 4.0
for w in (-0.3, 0.0, 0.7):
    y0 = dato_inicial(l0, w, 0.3, -0.2, 0.5)
    sol = solve_ivp(rhs, (l0, l1), y0, args=(w,), method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True)
    ls = np.linspace(l0 + 0.05, l1 - 0.05, 9)
    m = g.monopolo(ls)
    V = m["V"]
    s = V["r1"] / V["P1"]                                # r'/phi'
    res, esc, lig = [], [], []
    for i, l in enumerate(ls):
        y = sol.sol(l)
        a2, _, d2, _, q2 = rhs(l, y, w)
        a1, d0, d1, q0, q1 = y
        s0, s1, s2 = (s.deriv(k)[i] for k in range(3))
        Psi = d0 - s0 * q0
        Psi1 = d1 - s1 * q0 - s0 * q1
        Psi2 = d2 - s2 * q0 - 2 * s1 * q1 - s0 * q2
        K, P, U = (m[k].deriv(0)[i] for k in ("K", "P", "U"))
        P1 = m["P"].deriv(1)[i]
        r_ = -P1 * Psi1 - P * Psi2 + U * Psi - w * K * Psi
        res.append(r_)
        esc.append(abs(P1 * Psi1) + abs(P * Psi2) + abs(U * Psi) + abs(w * K * Psi))
        lig.append(max(abs(x) for x in ligaduras(l, y, w)))
    rel = max(abs(a) / b for a, b in zip(res, esc))
    RES[f"w={w}"] = dict(residuo_SL_relativo=rel, ligaduras_max=max(lig))
    print(f"w = {w:+.2f}: residuo relativo de la ecuacion reducida en Psi = d - (r'/phi') q: {rel:.2e};  "
          f"ligaduras E_beta, E_c propagadas: max {max(lig):.1e}", flush=True)
    assert rel < 1e-6

# identidad del termino cinetico: N^2 K/2 = W'
ls = np.array([0.5, 1.0, 2.0, 5.0, 20.0])
m = g.monopolo(ls)
V = m["V"]
xi = 4 * V["r0"] ** 2 * V["Gc"] * V["H"] / (V["N"] * m["gam"])                                  # xi = -al1 r H/gam
W = V["N"] * (xi - V["r0"] * V["Gc"] * (V["r1"] - V["r0"] * V["Phi1"]))
err = np.max(np.abs((V["E"] * m["K"] / 2 - W.d()).deriv(0)) / np.abs((V["E"] * m["K"]).deriv(0)))
RES["identidad_W"] = float(err)
print(f"N^2 K/2 - W' (relativo, 5 puntos): {err:.1e}")
L = 3000.0
Vb = g.V(np.array([L]))
mb = g.monopolo(np.array([L]))
r, r1, N = (float(Vb[k].deriv(0)[0]) for k in ("r0", "r1", "N"))
Phi1 = float(Vb["Phi1"].deriv(0)[0])
xiL = float((4 * Vb["r0"] ** 2 * Vb["Gc"] * Vb["H"] / (Vb["N"] * mb["gam"])).deriv(0)[0])      # xi = -al1 r H/gam
WL = N * (xiL - r * float(Vb["Gc"].deriv(0)[0]) * (r1 - r * Phi1))
m_ms = r * (1 - r1 ** 2) / 2
m_k = r * r * Phi1 * N
RES["W_inf"] = dict(W=WL, masa_MS=m_ms, masa_komar=m_k, prediccion=2 * m_ms + m_k)
print(f"W(L = {L:.0f}) = {WL:.5f};  2 m_MS + M_Komar = {2 * m_ms + m_k:.5f}  (m_MS = {m_ms:.5f}, M_K = {m_k:.5f})")
json.dump(RES, open(os.path.join(AQUI, "resultados", "verificar_reduccion.json"), "w"), indent=1)
