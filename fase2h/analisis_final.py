# -*- coding: utf-8 -*-
"""
Fase 2h: datos del informe (resultados/analisis_final.json y resultados/perfiles_fase2h.npz).
  A. Exterior (fondo de prueba de fondo_nu): asintotica de la teoria natural (U l^4, U_0 l^4, P l^2, K l^2) y autovalor
     de Dirichlet del exterior frente a theta (presupuesto de K movido fuera) y cm2 (P = cm2 N^2 K).
  B. Candidato rejilla_pozo (cumple la criba de la Fase 2g con la superficie fuera de la esfera de fotones): perfiles,
     criba con superficie fija y libre (y escala del muelle de superficie), prueba sin compresion, espectro completo
     frente a la rigidez Lam de las direcciones c != 0 (completo y familia c = 0, superficie fija y libre), modo
     inestable y velocidades caracteristicas.
  C. Semillas de la Fase 2g: criba con superficie libre y prueba sin compresion.
Uso: python analisis_final.py
"""
import os
import sys
import json
import time
import warnings
import numpy as np
from scipy.linalg import eigh

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2f"), os.path.join(AQUI, "..", "fase2g"), AQUI):
    sys.path.insert(0, p_)
import criba as C                     # noqa: E402
import fondo_nu as FN                 # noqa: E402
import teoria_nu as TN                # noqa: E402
import espectro2 as E2                # noqa: E402
import diseno_U as DU                 # noqa: E402
import velocidades as VE              # noqa: E402

warnings.filterwarnings("ignore")
T0 = time.time()
RES = {}
PERF = {}


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


PRUEBA = dict(wt=1.0, nu0=2.0, b=3.0, n=2, d0=0.5, wd=1.0, lam=1.0, lT=30.0, A=0.05, Ns=1.25, phi1=0.5)

# ------------------------------------------------------------------ A. exterior
F = FN.FondoNu(**PRUEBA)
T = TN.TeoriaNu(F, theta=0.0, b_in=0.05)
la = np.geomspace(1.2 * F.ls, 3500.0, 200)
p = DU.perfil_vacio(T, la)
Q2 = F.kappa * F.phi1 ** 2 * F.b ** 4
PERF.update(A_l=la, A_Ul4=p["U"] * la ** 4, A_U0l4=p["U0"] * la ** 4, A_Pl2=p["P"] * la ** 2, A_Kl2=p["K"] * la ** 2)
RES["A"] = dict(Q2=Q2, ls=F.ls, M=F.M_adm, U_sobre_P_l2_lejos=float(p["U"][-1] * la[-1] ** 2 / p["P"][-1]))
log(f"A: Q^2 = kappa phi1^2 b^4 = {Q2:.3f};  U l^2/P en l = 3500: {RES['A']['U_sobre_P_l2_lejos']:+.4f} (escalar sin masa: -2)")
le = np.concatenate([np.linspace(F.ls * (1 + 1e-6), 5 * F.ls, 700), np.geomspace(5 * F.ls, 1500, 700)[1:]])
RES["A"]["dirichlet"] = []
for th in (0.0, 0.25, 0.5, 0.75, 0.9):
    for cm2 in (1.0, 3.0):
        Tt = TN.TeoriaNu(F, theta=th, cm2=cm2, b_in=0.05)
        pt = DU.perfil_vacio(Tt, le)
        w, _ = DU.exterior_dirichlet(le, pt["K"], pt["P"], pt["U"], k=2)
        RES["A"]["dirichlet"].append(dict(theta=th, cm2=cm2, w0=float(w[0])))
        log(f"A: exterior Dirichlet theta = {th}, cm2 = {cm2}: w0 = {w[0]:+.3e}")

# ------------------------------------------------------------------ B. candidato rejilla_pozo
x = json.load(open(os.path.join(AQUI, "resultados", "rejilla_pozo.json")))["x"]
F = FN.FondoNu(**x)
ls = F.ls
lp = np.concatenate([np.linspace(1e-3, 1.6 * ls, 1500)])
f = F.campos(lp)
PERF.update(B_l=lp, B_N=f["N"], B_r=f["r"], B_r1=f["r1"], B_rp=f["rp"], B_Gc=f["Gc"], B_G4=f["G4"], B_sig=f["sig"],
            B_ls=np.array([ls]))
fs = F.campos(np.array([ls * (1 - 1e-9)]))
RES["B"] = dict(x=x, kappa=F.kappa, M_adm=F.M_adm, N0=F.N0, ls=ls, r_s=float(fs["r"][0]),
                rs_sobre_M=float(fs["r"][0] / F.M_adm), r_sig_ls=F.margen_exterior())
RES["B"]["criba_fija"] = C.criba_normalizada(F, n=160)
RES["B"]["criba_libre"] = C.criba_normalizada(F, n=160, libre=True)
RES["B"]["sin_compresion"] = C.prueba_sin_compresion(F)
log(f"B: criba fija {RES['B']['criba_fija']['margen']:+.4f}, libre {RES['B']['criba_libre']['margen']:+.4f}; "
    f"E_Q0 = {RES['B']['sin_compresion']['E_Q0']:+.3e}")
# escala del muelle de superficie
Es = C.energia_superficie(F)
nodos = C.malla_estrella(ls, 160)
fm = C.formas(F, nodos, que=("E", "M", "R"))
K = len(nodos)
Eo = F.E_out()
mu = []
for fac in (0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1e3, 1e4):
    Em, Rm = fm["E"].copy(), fm["R"].copy()
    Em[2 * (K - 1), 2 * (K - 1)] += Eo
    Rm[2 * (K - 1), 2 * (K - 1)] += abs(Eo)
    Em[2 * K - 1, 2 * K - 1] += fac * Es
    lib = [i for i in range(2 * K) if i != 1]
    w = eigh(Em[np.ix_(lib, lib)], Rm[np.ix_(lib, lib)], eigvals_only=True, subset_by_index=[0, 0])[0]
    mu.append(dict(factor=fac, Epar=float(w)))
RES["B"]["muelle"] = mu
log("B: escala del muelle: " + ", ".join(f"x{m['factor']:g}: {m['Epar']:+.4f}" for m in mu))
# prueba sin compresion: perfil y(l)
from scipy.integrate import solve_ivp         # noqa: E402


def rhs(l, z):
    g = F.campos(np.array([l]))
    return [z[1] / (4 * g["N"][0] ** 2 * g["r"][0] * g["Gc"][0]), -2 * g["N"][0] ** 2 * g["r"][0] * g["rp"][0] * z[0]]


sq = solve_ivp(rhs, (0, ls * (1 - 1e-9)), [1.0, 0.0], rtol=1e-10, atol=1e-13, dense_output=True)
ly = np.linspace(0, ls * (1 - 1e-9), 600)
PERF.update(B_ly=ly, B_yQ0=sq.sol(ly)[0])
# espectro frente a Lam
lc = ls
mallas = (240, 1500.0, 250)
RES["B"]["espectro"] = []
for Lam in (0.0, 0.01, 0.1, 1.0, 10.0, 100.0):
    Tt = TN.TeoriaNu(F, theta=0.0, eps=0.05, cm2=1.0, b_in=0.05, Lam=Lam)
    a_, b_ = E2.mallas(lc, mallas[0], mallas[1], mallas[2])
    cf = E2.Coefs(Tt, np.concatenate([a_, b_[1:]]))
    for libre in (False, True):
        for c0 in (False, True):
            for par in ("par", "impar"):
                o = E2.espectro(Tt, par, lc, *mallas, k=3, coefs=cf, libre=libre, c_cero=c0)
                RES["B"]["espectro"].append(dict(Lam=Lam, libre=libre, c_cero=c0, par=par, w=o.get("w"),
                                                 Tmin=o["Tmin"]))
    sel = [e for e in RES["B"]["espectro"] if e["Lam"] == Lam]
    log(f"B: Lam = {Lam}: " + "; ".join(f"{'libre' if e['libre'] else 'fija'}{' c=0' if e['c_cero'] else ''} {e['par']} "
                                        f"{e['w'][0]:+.2e}" for e in sel))
    if Lam == 10.0:
        # modo inestable con superficie libre (paridad par)
        E, Tm, nodos2, _, Tr, info = E2.ensamblar(Tt, "par", lc, *mallas, coefs=cf, libre=True, devolver_Tr=True)
        w, vec = eigh(E, Tm, subset_by_index=[0, 0])
        u = Tr @ vec[:, 0]
        h = u[0:info["nh"]:2]
        PERF.update(B_lmodo=nodos2, B_hmodo=h / h[0])
        lv = np.concatenate([np.geomspace(1e-2, ls, 300)[:-1], [0.999 * ls]])
        for Lv, Tv in ((0.0, TN.TeoriaNu(F, theta=0.0, eps=0.05, cm2=1.0, b_in=0.05, Lam=0.0)), (10.0, Tt)):
            v = VE.velocidades(Tv, lv)
            PERF.update(**{f"B_lv": lv, f"B_c2r_{int(Lv)}": v["c2_rapida"], f"B_c2l_{int(Lv)}": v["c2_lenta"],
                           f"B_detK_{int(Lv)}": v["detK_rel"], f"B_cphi2_{int(Lv)}": v["c2_vacio"], "B_cs2": v["cs2"]})
            RES["B"][f"vel_Lam{int(Lv)}"] = dict(c2_rapida_min=float(np.nanmin(v["c2_rapida"])),
                                                 c2_rapida_l1=float(np.interp(1.0, lv, v["c2_rapida"])),
                                                 c2_rapida_l01=float(np.interp(0.1, lv, v["c2_rapida"])),
                                                 c2_lenta_max=float(np.nanmax(v["c2_lenta"])),
                                                 cphi2_max=float(np.nanmax(v["c2_vacio"])),
                                                 detK_min=float(np.nanmin(v["detK_rel"])))
            log(f"B: velocidades (Lam = {Lv}): " + ", ".join(f"{k} {val:.3g}" for k, val in RES["B"][f"vel_Lam{int(Lv)}"].items()))

# ------------------------------------------------------------------ C. semillas de la Fase 2g
RES["C"] = {}
for sem in (1, 2):
    xs = json.load(open(os.path.join(AQUI, "..", "fase2g", "resultados", f"criba_optimiza_{sem}.json")))["x"]
    for etiqueta, xx in ((f"semilla{sem}", xs), (f"semilla{sem}_m1", dict(xs, m=1.0))):
        G = C.FondoFluido(**xx)
        RES["C"][etiqueta] = dict(m=xx["m"], criba_fija=C.criba_normalizada(G, n=160),
                                  criba_libre=C.criba_normalizada(G, n=160, libre=True),
                                  sin_compresion=C.prueba_sin_compresion(G))
        r = RES["C"][etiqueta]
        log(f"C: {etiqueta}: fija {r['criba_fija']['margen']:+.4f}, libre Epar {r['criba_libre']['Epar']:+.4f}, "
            f"E_Q0 = {r['sin_compresion']['E_Q0']:+.3e} (E[1] = {r['sin_compresion']['E1']:+.3e})")

json.dump(RES, open(os.path.join(AQUI, "resultados", "analisis_final.json"), "w"), indent=1, default=float)
np.savez(os.path.join(AQUI, "resultados", "perfiles_fase2h.npz"), **PERF)
log("escrito resultados/analisis_final.json y perfiles_fase2h.npz")
