# -*- coding: utf-8 -*-
"""
Comprobaciones de la Fase 2f (teoria cuartica mas alla de Horndeski sobre fondos estaticos esfericos).

  1. bH_gen.py (exportar_bH.py) reproduce las expresiones simbolicas de resultados/bH_simbolico.pkl.
  2. Identidades de Noether linealizadas del sector l = 0 (reparametrizaciones zeta(t, l) de l y tau(t, l) de t),
     en aritmetica racional exacta, para un chorro de fondo arbitrario FUERA de la capa de masas y perturbaciones
     arbitrarias (E_A, E_C, E_p del fondo; E_a, E_be, E_c, E_d, E_q las lineales en el gauge beta = c = 0):
       zeta:  Phi' E_a + r' E_d + phi' E_q - i w E_be - D_l E_c + 2 E E_A (2 Phi' a + a') + 2 E_C (r' d + r d')
              + E_p q' = 0
       tau:   -i w (E_a + 2 E E_A a - 2 r E_C d - E_p q) + D_l (E E_be) = 0
  3. Modos gauge sobre chorros EN la capa de masas (G_ij racionales aleatorios, G2 fijado por la ligadura E_B = 0,
     derivadas segundas y terceras del fondo despejadas de E_A = E_C = E_p = 0 y de sus derivadas en l):
       zeta = e^{i omega t}:  a = Phi' - w int e^{-2 Phi},  d = r',  q = phi';      tau:  a = constante.
  4. Limite RG + escalar + V (G4 = 1, G2 = X - V): las cinco ecuaciones lineales se anulan sobre los datos del
     sistema de fase2b/fase2e (E^t_t, E^th_th, E_phi con las ligaduras E^l_t, E^l_l) en estados aleatorios.
  5. Sector impar (impar_simbolico.py, impar_ecuaciones.py):
     RG en Schwarzschild (longitud propia): Kerr lineal y rotacion rigida (l = 1), solucion estatica
     h0 = r^2 (r - 2M) (l = 2), y controles que no son solucion;
     teoria completa: modo gauge de l = 1 (rotacion diferencial: H0 = w r^2 Lambda, H1 = -r^2 Lambda', con
     h1 -> i omega H1) sobre chorros en la capa de masas, para Lambda(l) arbitraria;
     dependencia en l: L2(l) / N_l es afin en (l - 1)(l + 2) con N_l = l (l + 1) / (2 l + 1) (l = 1, 2, 3).
Uso: python verificar_bH.py        (escribe resultados/verificar_bH.json)
"""
import os
import sys
import json
import math
import time
import pickle
import numpy as np
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
sys.path.append(os.path.join(AQUI, "..", "fase2e"))
T0 = time.time()
RES = {}


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


d = pickle.load(open(os.path.join(AQUI, "resultados", "bH_simbolico.pkl"), "rb"))
fe, lin, iw, amp = (d[k] for k in ("fondo_eq", "lineales", "iw", "amp"))
from chorros import Phi, rr, PP, Ee, eP, Xs, GPAT, D, Valores        # noqa: E402
from chorros import chorro_en_capa as _chorro                         # noqa: E402
assert Ee == d["Ee"] and Xs == d["Xs"]
gauge = {amp[(v, k)]: 0 for v in ("be", "c") for k in range(6)}
LIN = {v: sp.expand(lin[v].subs(gauge)) for v in ("a", "be", "c", "d", "q")}
AMPS = [[amp[(v, k)] for k in range(6)] for v in ("a", "d", "q")]


def chorro_en_capa(semilla, orden=3):
    return _chorro(fe, semilla, orden)


# ------------------------------------------------------------------ 1. bH_gen frente al pickle
import bH_gen                                              # noqa: E402
rng = np.random.default_rng(1)
err1 = 0.0
w_s = sp.Symbol("w", real=True)
for _ in range(4):
    G = {n: float(rng.uniform(-1, 1)) for n in bH_gen.SIMB_G}
    Ph = [float(x) for x in rng.uniform(-1, 1, 3)]
    r_ = [float(rng.uniform(0.5, 2)), float(rng.uniform(-1, 1)), float(rng.uniform(-1, 1))]
    P_ = [float(x) for x in rng.uniform(-1, 1, 3)]
    wv = float(rng.uniform(-2, 2))
    val = {sp.Symbol(n, real=True): sp.Float(x, 30) for n, x in G.items()}
    val.update({Phi[k]: Ph[k] for k in range(3)})
    val.update({rr[k]: r_[k] for k in range(3)})
    val.update({PP[k]: P_[k] for k in range(3)})
    val[Ee] = math.exp(2 * Ph[0])
    val[Xs] = -P_[1] ** 2 / 2
    num = bH_gen.fondo(G, Ph, r_, P_)
    for x, n in zip(num, ("A", "B", "C", "p")):
        ref = float(fe[n].xreplace(val))
        err1 = max(err1, abs(x - ref) / (1 + abs(ref)))
    coef = bH_gen.lineal(G, Ph, r_, P_, wv)
    for ec in ("a", "d", "q", "be", "c"):
        e = LIN[ec] / iw if ec == "be" else LIN[ec]
        e = sp.expand(e).subs(iw ** 2, -w_s)
        for vv in ("a", "d", "q"):
            for k in range(3):
                ref = float(sp.expand(e.coeff(amp[(vv, k)])).xreplace({**val, w_s: wv}))
                x = coef.get(ec, {}).get((vv, k), 0.0)
                err1 = max(err1, abs(x - ref) / (1 + abs(ref)))
RES["bH_gen_vs_simbolico"] = err1
log(f"1. bH_gen.py frente a las expresiones simbolicas (4 puntos, fondo y 35 coeficientes): {err1:.1e}")
assert err1 < 1e-12

# ------------------------------------------------------------------ 2. identidades de Noether linealizadas
E_A, E_C, E_p = fe["A"], fe["C"], fe["p"]
a, dd, q = AMPS
DEc = D(LIN["c"], AMPS)
DEbe = D(Ee * LIN["be"], AMPS)
id_zeta = (Phi[1] * LIN["a"] + rr[1] * LIN["d"] + PP[1] * LIN["q"] - iw * LIN["be"] - DEc
           + 2 * Ee * E_A * (2 * Phi[1] * a[0] + a[1]) + 2 * E_C * (rr[1] * dd[0] + rr[0] * dd[1]) + E_p * q[1])
id_tau = -iw * (LIN["a"] + 2 * Ee * E_A * a[0] - 2 * rr[0] * E_C * dd[0] - E_p * q[0]) + DEbe
for sem in range(3):
    v = Valores(100 + sem)
    v[Ee] = sp.Rational(121, 100)
    for k in range(3):
        v[Phi[k + 1]] = v.r()
        v[rr[k]] = v.r() if k else sp.Rational(v.rnd.randint(600, 2000), 1000)
        v[PP[k]] = v.r()
        for s in (a[k], dd[k], q[k]):
            v[s] = v.r()
    v[Xs] = -v[PP[1]] ** 2 / 2
    for nombre, e in (("zeta", id_zeta), ("tau", id_tau)):
        x = v.ev(e)
        assert x == 0, (nombre, sem, x)
RES["noether_lineal"] = "exacta (racional) en 3 chorros fuera de la capa de masas, polinomio en i omega"
log("2. identidades de Noether linealizadas (zeta y tau): residuo exactamente 0 en 3 chorros fuera de la capa "
    "de masas")

# ------------------------------------------------------------------ 3. modos gauge en la capa de masas
res3 = {}
for sem in range(3):
    v = chorro_en_capa(200 + sem, orden=3)
    # comprobacion de paso: la derivada de la ligadura se anula (identidad de Noether del fondo)
    assert v.ev(D(fe["B"])) == 0
    E_ = v[Ee]
    modo_z = {a[0]: v[Phi[1]], a[1]: v[Phi[2]] + iw ** 2 / E_, a[2]: v[Phi[3]] - 2 * iw ** 2 * v[Phi[1]] / E_,
              dd[0]: v[rr[1]], dd[1]: v[rr[2]], dd[2]: v[rr[3]], q[0]: v[PP[1]], q[1]: v[PP[2]], q[2]: v[PP[3]]}
    modo_t = {s: 0 for seq in AMPS for s in seq}
    modo_t[a[0]] = 1
    for ec in ("a", "be", "c", "d", "q"):
        for nombre, modo in (("zeta", modo_z), ("tau", modo_t)):
            mm = {s: 0 for seq in AMPS for s in seq}
            mm.update(modo)
            x = v.ev(LIN[ec].xreplace(mm))
            res3[f"{nombre}_{ec}_{sem}"] = str(x)
            assert x == 0, (nombre, ec, sem, x)
RES["modos_gauge_l0"] = "exactos en 3 chorros en la capa de masas (5 ecuaciones, todo omega)"
log("3. modos gauge zeta (a = Phi' - w int e^{-2Phi}, d = r', q = phi') y tau (a = cte) en la capa de masas: "
    "las 5 ecuaciones lineales se anulan exactamente (3 chorros, todo omega)")

# ------------------------------------------------------------------ 4. limite RG + escalar + V frente a fase2b
import potencial as PV                                     # noqa: E402
H = PV.H
rg_G = {n: 0.0 for n in bH_gen.SIMB_G}
err4 = 0.0
rng = np.random.default_rng(7)
for _ in range(20):
    Ph0, Ph1 = rng.uniform(-1, 1, 2)
    r0_, r1_ = rng.uniform(0.6, 2), rng.uniform(-0.9, 0.9)
    p0_, p1_ = rng.uniform(-1, 1, 2)
    V1, V2 = rng.uniform(-1, 1, 2)
    y = np.array([r0_, r1_, Ph0, Ph1, p0_, p1_])
    # V0 de la ligadura l l (fondo en la capa de masas): est_ligadura + V0/2 = 0
    V0 = -2 * H.est_ligadura(0.0, r0_, r1_, Ph1, p1_)
    Vd = (lambda V0=V0, V1=V1, V2=V2, p0_=p0_: (lambda x: (V0 + V1 * (x - p0_) + V2 * (x - p0_) ** 2 / 2,
                                                          V1 + V2 * (x - p0_), V2, 0.0)))()
    Fd0 = PV.acoplamiento("cuadratico", 0.0)
    (r2_, Ph2, p2_), _ = PV.segundas_V(y, Fd0, Vd)
    assert abs(PV.est_ligadura_V(y, Fd0, Vd)) < 1e-12
    e2 = math.exp(-2 * Ph0)
    wv = float(rng.uniform(-2, 2))
    bg = (r0_, r1_, r2_, Ph1, Ph2, p1_, p2_, e2)
    # datos (d0, q0, q1) al azar; (a1, d1) de las ligaduras E^l_t = E^l_l = 0 (lineales en a1, d1)

    def lig(a1, d1, d0, q0, q1):
        lt, ll = H.lin_ligaduras(0.0, 0.0, 0.0, wv, *bg, a1, d0, d1, q0, q1)
        return np.array([lt, ll + V1 * q0 / 2])
    d0_, q0_, q1_ = rng.uniform(-1, 1, 3)
    b0 = lig(0.0, 0.0, d0_, q0_, q1_)
    Mlig = np.column_stack([lig(1.0, 0.0, d0_, q0_, q1_) - b0, lig(0.0, 1.0, d0_, q0_, q1_) - b0])
    a1_, d1_ = np.linalg.solve(Mlig, -b0)
    M = np.array(H.lin_principal(0.0, r0_, r1_, Ph1, p1_))
    bb = np.array(H.lin_resto(0.0, 0.0, 0.0, wv, *bg, a1_, d0_, d1_, q0_, q1_)) + np.array([V1 * q0_ / 2, V1 * q0_ / 2,
                                                                                         -V2 * q0_])
    d2_, a2_, q2_ = np.linalg.solve(M, -bb)
    G = dict(rg_G)
    G.update({"G4_00": 1.0, "G2_00": -p1_ ** 2 / 2 - V0, "G2_01": 1.0, "G2_10": -V1, "G2_20": -V2})
    coef = bH_gen.lineal(G, [Ph0, Ph1, Ph2], [r0_, r1_, r2_], [p0_, p1_, p2_], wv)
    a0_ = float(rng.uniform(-1, 1))           # a sin derivar: su coeficiente debe anularse en la capa de masas
    vals = {("a", 0): a0_, ("a", 1): a1_, ("a", 2): a2_, ("d", 0): d0_, ("d", 1): d1_, ("d", 2): d2_,
            ("q", 0): q0_, ("q", 1): q1_, ("q", 2): q2_}
    for ec, cs in coef.items():
        tot = sum(c * vals[k] for k, c in cs.items())
        esc = sum(abs(c * vals[k]) for k, c in cs.items())
        err4 = max(err4, abs(tot) / esc)
RES["limite_RG_lineal"] = err4
log(f"4. limite RG + escalar + V: las 5 ecuaciones lineales sobre datos de fase2b/fase2e (20 estados en la capa "
    f"de masas): residuo relativo {err4:.1e}")
assert err4 < 1e-10

# ------------------------------------------------------------------ 5. sector impar
IMP = {}
for LL in (1, 2, 3):
    fn = os.path.join(AQUI, "resultados", f"impar_ec_L{LL}.pkl")
    if os.path.exists(fn):
        IMP[LL] = pickle.load(open(fn, "rb"))
R_, M_, J_ = sp.symbols("R M J", positive=True)
f_ = 1 - 2 * M_ / R_


def Dl_s(e):                                               # d/dl = r' d/dr en Schwarzschild, r' = sqrt(f)
    return sp.sqrt(f_) * sp.diff(e, R_)


def impar_RG(LL, H0e, H1e, wv):
    di = IMP[LL]
    ecs, amp_i, w_i = di["ecs"], di["amp"], di["w"]
    rk = [R_]
    Pk = [sp.log(f_) / 2]
    for _ in range(7):
        rk.append(sp.simplify(Dl_s(rk[-1])))
        Pk.append(sp.simplify(Dl_s(Pk[-1])))
    bg = {rr[k]: rk[k] for k in range(8)}
    bg.update({Phi[k]: Pk[k] for k in range(1, 8)})
    bg[eP] = sp.sqrt(f_)
    bg.update({PP[k]: 0 for k in range(9)})
    bg[Xs] = 0
    out = {}
    for v_ in ("h0", "h1"):
        e = ecs[v_].subs(w_i, wv)
        rep = {}
        for s in e.free_symbols:
            if GPAT.match(s.name):
                rep[s] = 1 if s.name == "G4_00" else 0
        sub = {}
        h = H0e
        for k in range(7):
            sub[amp_i[("h0", k)]] = h
            h = Dl_s(h)
        h = H1e
        for k in range(7):
            sub[amp_i[("h1", k)]] = h
            h = Dl_s(h)
        out[v_] = sp.simplify(e.xreplace(rep).xreplace(bg).xreplace(sub))
    return out


res5 = {}
if 1 in IMP:
    x = impar_RG(1, 2 * J_ / R_, 0, 0)
    res5["RG_Kerr_l1"] = str(x["h0"])
    assert x["h0"] == 0
    for wv in (0, sp.Rational(3, 7)):
        x = impar_RG(1, R_ ** 2, 0, wv)
        res5[f"RG_rotacion_rigida_l1_w{wv}"] = [str(x["h0"]), str(x["h1"])]
        assert x["h0"] == 0 and x["h1"] == 0
    x = impar_RG(1, 1 / R_ ** 2, 0, 0)
    res5["RG_control_l1"] = str(x["h0"])
    assert x["h0"] != 0
    log("5a. impar RG l = 1: Kerr lineal y rotacion rigida (w = 0, 3/7) son solucion; control h0 = 1/r^2 no")
if 2 in IMP:
    x = impar_RG(2, R_ ** 2 * (R_ - 2 * M_), 0, 0)
    res5["RG_estatica_l2"] = str(x["h0"])
    assert x["h0"] == 0
    x = impar_RG(2, R_ ** 3, 0, 0)
    res5["RG_control_l2"] = str(x["h0"])
    assert x["h0"] != 0
    log(f"5b. impar RG l = 2: h0 = r^2 (r - 2M) estatica es solucion; control h0 = r^3 da {x['h0']}")
if 1 in IMP:
    di = IMP[1]
    ecs, amp_i, w_i = di["ecs"], di["amp"], di["w"]
    Lam = sp.symbols("Lam0:8", real=True)
    nmax = max(k for (v_, k), s in amp_i.items() if any(s in ecs[u].free_symbols for u in ecs))
    for sem in range(3):
        v = chorro_en_capa(300 + sem, orden=3)
        for k in range(8):
            v[Lam[k]] = v.r()
        h0e = w_i * rr[0] ** 2 * Lam[0]
        h1e = -rr[0] ** 2 * Lam[1]
        sub = {}
        for k in range(nmax + 1):
            sub[amp_i[("h0", k)]] = h0e
            sub[amp_i[("h1", k)]] = h1e
            h0e, h1e = D(h0e, [Lam]), D(h1e, [Lam])
        for v_ in ("h0", "h1"):
            x = sp.expand(v.ev(ecs[v_].xreplace(sub)))
            assert x == 0, (sem, v_, x)
    res5["gauge_l1_teoria_completa"] = "exacto en 3 chorros en la capa de masas, Lambda(l) y omega arbitrarios"
    log("5c. impar l = 1, teoria completa: el modo gauge (rotacion diferencial) anula E_h0 y E_h1 exactamente "
        "(3 chorros en la capa de masas, Lambda y omega arbitrarios)")
if all(k in IMP for k in (1, 2, 3)):
    # L2 integrada: N_l [alfa + beta (l - 1)(l + 2)] con N_l = l (l + 1)/(2 l + 1)  (mismos simbolos en los tres)
    Ls = {}
    for LL in (1, 2, 3):
        Ls[LL] = pickle.load(open(os.path.join(AQUI, "resultados", f"impar_L{LL}.pkl"), "rb"))["L2"]
    N = {LL: sp.Rational(LL * (LL + 1), 2 * LL + 1) for LL in (1, 2, 3)}
    lam = {LL: (LL - 1) * (LL + 2) for LL in (1, 2, 3)}
    alfa = sp.expand(Ls[1] / N[1])
    beta = sp.expand((Ls[2] / N[2] - alfa) / lam[2])
    resto = sp.expand(Ls[3] / N[3] - alfa - lam[3] * beta)
    res5["dependencia_l"] = str(resto)
    assert resto == 0
    log("5d. dependencia en l de la lagrangiana impar: L2(3) = N_3 [alfa + 10 beta] con alfa, beta de l = 1, 2 "
        "(exacto)")
RES["impar"] = res5
json.dump(RES, open(os.path.join(AQUI, "resultados", "verificar_bH.json"), "w"), indent=1, default=str)
log("todas las comprobaciones superadas; resultados/verificar_bH.json")
