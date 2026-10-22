# -*- coding: utf-8 -*-
"""
Resultados numericos de la Fase 2f sobre gargantas reconstruidas (reconstruccion.Garganta):

  1. Ultraestaticas (M = 0): fondo en la capa de masas; K cambia de signo con int_0^inf K = W(inf) - W(0) = 0
     (fantasma obligado) e int U = -2 int Gcal r'^2/r^2 < 0 (taquion obligado: funcion de prueba d = const).
  2. Garganta masiva de referencia (M = 2, a = 2): K, P > 0 (hasta l = L y en la cola asintotica, cola_K),
     velocidades, identidad int N^2 K = 2 W y masas,
     descomposicion de int U, espectro completo del monopolo (taquiones rapidos y lentos) y del sector impar.
  3. Campo debil: sucesion M = 2 ... 0.05 (con phi1 ~ sqrt(M) para que K > 0 sea posible): w0 ~ -1/M^p, numero de
     taquiones, cota variacional int U/int K.
  4. Vias de escape: campo fuerte y superluminalidad (c^2 hasta 30).
Escribe resultados/gargantas_bH.json y resultados/perfiles_bH.npz (para figuras_fase2f.py).
Uso: python gargantas_bH.py
"""
import os
import sys
import json
import math
import time
import numpy as np
from scipy.integrate import simpson
from concurrent.futures import ProcessPoolExecutor

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import reconstruccion as R          # noqa: E402
import espectro as E                # noqa: E402
import busqueda as B                # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


BASE = dict(n=2, b=4.0, frac=0.9, wdb=0.5, M=2.0, a=2.0, wK=10.0, wPK=1.5, cm2=0.8)
ULTRA = [dict(b=1.0, n=2, d0=2.0, wd=1.0), dict(b=4.0, n=2, d0=1.2, wd=2.0), dict(b=2.0, n=1, d0=1.8, wd=2 ** 0.5)]


def malla_integral(L=3000.0):
    return np.concatenate([np.linspace(1e-3, 50, 20001), np.geomspace(50, L, 4001)[1:]])


def perfiles(g, ls, sector="monopolo", trozo=2000):
    claves = ("K", "P", "U")
    acc = {k: [] for k in claves + ("r", "r1", "N", "Phi1", "Gc", "H", "G4", "P1", "W")}
    for i in range(0, len(ls), trozo):
        o = g.monopolo(ls[i:i + trozo]) if sector == "monopolo" else g.impar(ls[i:i + trozo])
        V = o["V"]
        for k in claves:
            acc[k].append(o[k].deriv(0))
        for k, n in (("r", "r0"), ("r1", "r1"), ("N", "N"), ("Phi1", "Phi1"), ("Gc", "Gc"), ("H", "H"), ("G4", "G4"),
                     ("P1", "P1")):
            acc[k].append(V[n].deriv(0))
        if sector == "monopolo":
            xi = V["r0"] * V["Gc"] * V["H"] / (V["r1"] * (V["N"] * o["gam"] / (4 * V["r0"] * V["r1"])))
            W = V["N"] * (xi - V["r0"] * V["Gc"] * (V["r1"] - V["r0"] * V["Phi1"]))
            acc["W"].append(W.deriv(0))
    return {k: np.concatenate(v) for k, v in acc.items() if v}


def ultra(p):
    g = R.Garganta(phi1=1.0, Lmax=4000.0, **p)
    ls = malla_integral()
    f = perfiles(g, ls)
    res = g.residuos_fondo(np.array([1e-3, 0.5, 2.0, 10.0, 100.0]))
    IK = simpson(f["K"], x=ls)
    IKabs = simpson(np.abs(f["K"]), x=ls)
    IU = simpson(f["U"], x=ls)
    I1 = simpson(-2 * f["Gc"] * f["r1"] ** 2 / f["r"] ** 2, x=ls)
    return dict(p=p, res_fondo=float(np.max(np.abs(res))), minK=float(f["K"].min()), maxK=float(f["K"].max()),
                int_K=float(IK), int_absK=float(IKabs), W_L=float(f["W"][-1]), int_U=float(IU),
                menos2_int_Gr2=float(I1))


def masiva(x):
    g = B.construir(x)
    ls = malla_integral()
    f = perfiles(g, ls)
    IU = simpson(f["U"], x=ls)
    I1 = simpson(-2 * f["N"] * f["Gc"] * f["r1"] ** 2 / f["r"] ** 2, x=ls)
    I2 = simpson(f["Phi1"] ** 2 * f["P"], x=ls)
    IN2K = simpson(f["N"] ** 2 * f["K"], x=ls)
    L = ls[-1]
    mMS = f["r"][-1] * (1 - f["r1"][-1] ** 2) / 2
    mK = f["r"][-1] ** 2 * f["Phi1"][-1] * f["N"][-1]
    res = g.residuos_fondo(np.array([1e-3, 0.5, 2.0, 10.0, 100.0]))
    lm = E.malla(1e-3, 400.0, 2000)
    K, P, U, N = E.coeficientes(g, "monopolo", lm)
    sl = E.SL(lm, K, P, U)
    out = dict(x=x, res_fondo=float(np.max(np.abs(res))), minK=float(min(K.min(), f["K"].min())),
               minP=float(min(P.min(), f["P"].min())),                                   # malla y hasta l = L
               maxc2=float(np.max(P / (K * N ** 2))), N0=float(g.N0),
               int_N2K=float(IN2K), dos_W_L=float(2 * f["W"][-1]), masa_MS_L=float(mMS), masa_Komar_L=float(mK),
               int_U=float(IU), I1=float(I1), I2_Phi2P=float(I2), resto=float(IU - I1 - I2), int_K=float(simpson(f["K"], x=ls)),
               cola=g.cola_K())
    for par in ("par", "impar"):
        nz = sl.ceros(par)[0]
        out[f"w_{par}"] = [float(w) for w in sl.autovalores_negativos(par, nz, tol=1e-8)] if nz else []
    for LL in (2, 3):
        Ko, Po, Uo, No = E.coeficientes(g, "impar", lm, LL=LL)
        o = dict(minK=float(Ko.min()), minP=float(Po.min()), maxc2=float(np.max(Po / (Ko * No ** 2))))
        if Ko.min() > 0 and Po.min() > 0:
            so = E.SL(lm, Ko, Po, Uo)
            o["nodos"] = [so.ceros("par")[0], so.ceros("impar")[0]]
        out[f"impar_l{LL}"] = o
    return out, g, ls, f


def debil(M):
    x = dict(BASE, M=M, phi1=math.sqrt(M / 2))
    r = B.evalua(x, malla=E.malla(1e-3, 400.0, 1500))
    g = B.construir(x)
    ls = malla_integral()
    f = perfiles(g, ls)
    r["int_U_mitad"] = float(simpson(f["U"], x=ls))
    r["int_K_mitad"] = float(simpson(f["K"], x=ls))
    return r


def escape(args):
    M, a, wK, cm2 = args
    x = dict(BASE, M=M, a=a, wK=wK, cm2=cm2)
    return B.evalua(x, malla=E.malla(1e-3, 400.0, 900))


if __name__ == "__main__":
    RES = {}
    with ProcessPoolExecutor(max_workers=3) as ex:
        RES["ultraestaticas"] = list(ex.map(ultra, ULTRA))
    for u in RES["ultraestaticas"]:
        log(f"ultraestatica {u['p']}: |E_u| <= {u['res_fondo']:.1e}; K en [{u['minK']:+.3f}, {u['maxK']:+.3f}], "
            f"int_0^L K = {u['int_K']:+.2e} (int |K| = {u['int_absK']:.3f}), W(L) = {u['W_L']:+.2e}; "
            f"int_0^L U = {u['int_U']:+.4f} = -2 int Gcal r'^2/r^2 = {u['menos2_int_Gr2']:+.4f}")
    out, g, ls, f = masiva(BASE)
    RES["masiva"] = out
    log(f"masiva de referencia: |E_u| <= {out['res_fondo']:.1e}; N0 = {out['N0']:.3f}; min K = {out['minK']:+.2e}, "
        f"min P = {out['minP']:+.2e}, max c^2 = {out['maxc2']:.4f}")
    log(f"   int_0^L N^2 K = {out['int_N2K']:.5f} = 2 W(L) = {out['dos_W_L']:.5f};  m_MS(L) = {out['masa_MS_L']:.4f}, "
        f"M_Komar(L) = {out['masa_Komar_L']:.4f}  (M = {BASE['M']})")
    log(f"   int_0^L U = {out['int_U']:+.4f} = I1 {out['I1']:+.4f} + int Phi'^2 P {out['I2_Phi2P']:+.4f} + resto "
        f"{out['resto']:+.4f};  int K = {out['int_K']:.3f}")
    c = out["cola"]
    log(f"   cola (l > L): (r' Phi' - r'') l^3 = {c['A']:.2f} L^2 {c['alfa']:+.2f} L {c['beta']:+.2f} (C1 = {c['C1']:.4f}; "
        f"dispersion de beta {c['dispersion_beta']:.1e}); minimo para L >= ln L: {c['ymin']:+.3f}")
    log(f"   monopolo: w_par = {out['w_par']},  w_impar = {out['w_impar']}")
    for LL in (2, 3):
        log(f"   impar l = {LL}: {out[f'impar_l{LL}']}")
    # perfiles para las figuras
    lsf = np.concatenate([np.geomspace(1e-3, 1, 200), np.linspace(1.0, 40.0, 800)[1:]])
    ff = perfiles(g, lsf, trozo=1000)
    fo = perfiles(g, lsf, sector="impar", trozo=1000)
    gu = R.Garganta(phi1=1.0, Lmax=4000.0, **ULTRA[1])
    fu = perfiles(gu, lsf, trozo=1000)
    np.savez(os.path.join(AQUI, "resultados", "perfiles_bH.npz"), l=lsf,
             **{f"m_{k}": v for k, v in ff.items()}, **{f"o_{k}": v for k, v in fo.items()},
             **{f"u_{k}": v for k, v in fu.items()})
    Ms = [2.0, 1.0, 0.5, 0.25, 0.1, 0.05]
    with ProcessPoolExecutor(max_workers=6) as ex:
        RES["campo_debil"] = list(ex.map(debil, Ms))
    for M, r in zip(Ms, RES["campo_debil"]):
        if "n_par" in r:
            w0 = min(r["w_par"] + r["w_impar"] + [0.0])
            log(f"campo debil M = {M}: N0 = {r['N0']:.3f}; taquiones ({r['n_par']}, {r['n_impar']}); w0 = {w0:+.4f} "
                f"(M w0 = {M * w0:+.3f}); int U/int K (cota) = {r['int_U_mitad'] / r['int_K_mitad']:+.4f}")
        else:
            log(f"campo debil M = {M}: N0 = {r['N0']:.3f}; min K = {r['minK']:+.2e}, min P = {r['minP']:+.2e}, "
                f"minimo de la cola {r['cola_ymin']:+.2e} (fantasma)")
    combos = [(M, a, wK, cm2) for (M, a) in ((3.0, 2.0), (4.0, 3.0), (8.0, 4.0)) for wK in (6.0,) for cm2 in (1.0, 4.0, 30.0)]
    with ProcessPoolExecutor(max_workers=8) as ex:
        esc = list(ex.map(escape, combos))
    RES["escape"] = [dict(M=c[0], a=c[1], wK=c[2], cm2=c[3], **{k: v for k, v in r.items() if k != "x"}) for c, r in zip(combos, esc)]
    for c, r in zip(combos, esc):
        if "n_par" in r:
            log(f"escape M = {c[0]}, a = {c[1]}, cm2 = {c[3]}: N0 = {r['N0']:.3f}, max c^2 = {r['maxc2']:.2f}, "
                f"taquiones ({r['n_par']}, {r['n_impar']}): {r['w_par']} / {r['w_impar']}")
        else:
            log(f"escape M = {c[0]}, a = {c[1]}, cm2 = {c[3]}: no valido (min K = {r.get('minK')}, min P = "
                f"{r.get('minP')}, minimo de la cola {r.get('cola_ymin')})")
    json.dump(RES, open(os.path.join(AQUI, "resultados", "gargantas_bH.json"), "w"), indent=1, default=float)
    log("guardado resultados/gargantas_bH.json y resultados/perfiles_bH.npz")
