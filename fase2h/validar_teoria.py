# -*- coding: utf-8 -*-
"""
Fase 2h: validacion de la reconstruccion (teoria.Teoria) sobre la mejor configuracion de la Fase 2g (semilla 2) con la
zona lejana canonica (fondo.Fondo):
  1. ecuaciones del fondo con fluido (E_A, E_B, E_C, E_p de fase2g/resultados/fluido_simbolico.pkl, compiladas aparte)
     evaluadas sobre los chorros de la teoria reconstruida: deben anularse;
  2. K de la formula T2 (simbolico_num.MONO["K"], con gam, al1 de la teoria) frente a la K de diseno;
  3. gam frente a r (4 r' Gamma1)/N con la Gamma1 de diseno;
  4. continuidad del fondo en lm (chorros interiores frente a lejanos);
  5. K > 0 en todas partes y K >= N (rho~+p~) xi^2/(2 Gcal^2)/N^2 dentro de la estrella (ausencia de fantasmas).
Uso: python validar_teoria.py
"""
import os
import sys
import json
import time
import pickle
import numpy as np
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2f"), os.path.join(AQUI, "..", "fase2g"), AQUI):
    sys.path.insert(0, p_)
import jets as J                      # noqa: E402
import simbolico_num as SN            # noqa: E402
import reconstruccion as R            # noqa: E402
import fondo as FH                    # noqa: E402
import teoria as TE                   # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:6.1f} s] {m}", flush=True)


fl = pickle.load(open(os.path.join(AQUI, "..", "fase2g", "resultados", "fluido_simbolico.pkl"), "rb"))
FONDO = {k: SN.compilar(SN.con_Xb(fl["fondo_eq"][k])) for k in ("A", "B", "C", "p")}

if __name__ == "__main__":
    x = json.load(open(os.path.join(AQUI, "..", "fase2g", "resultados", "criba_optimiza_2.json")))["x"]
    F = FH.Fondo(**x, kappa=1.0, phi1=1.0, Lmax=4000.0)
    log(f"fondo: ok {F.ok}, ls = {F.ls:.3f}, lm = {F.lm:.1f}, Delta = {F.Delta:.1f}, beta_nu = {F.beta_nu:.2f}, "
        f"M_adm = {F.M_adm:.4f}, N_inf = {F.N_inf:.5f}")
    T = TE.Teoria(F)
    log(f"teoria: ventana de K [{T.lw1:.1f}, {T.lw2:.1f}], de P [{T.wP1:.1f}, {T.wP2:.1f}], Lz = {T.Lz:.1f}, B = {T.B}, "
        f"exceso sin B = {T.exceso_sin_B:+.4f}")
    if T.B is None:
        sys.exit(0)
    l = np.array([0.05, 0.5, 2.0, 10.0, 20.0, 27.0, 40.0, 100.0, 170.0, 190.0, 300.0, 500.0, 900.0, 2000.0])
    V = T.V(l)
    res = {k: np.abs(f(V).c[0]) for k, f in FONDO.items()}
    escala = {k: np.abs(V["G4"].c[0]) * 1.0 for k in FONDO}
    for i, li in enumerate(l):
        log(f"  l = {li:7.2f}: |E_A| = {res['A'][i]:.1e}, |E_B| = {res['B'][i]:.1e}, |E_C| = {res['C'][i]:.1e}, "
            f"|E_p| = {res['p'][i]:.1e}")
    W = R.Garganta._W(V)
    Kf = SN.MONO["K"](W).c[0]
    Kd = V["K"].c[0]
    gam = SN.MONO["gam"](V).c[0]
    G1 = V["G1"].c[0]
    r, r1, N = V["r0"].c[0], V["r1"].c[0], V["N"].c[0]
    log("  K (T2) / K (diseno) - 1 = " + np.array2string((Kf / Kd - 1), precision=2))
    log("  gam / (4 r r' Gamma1/N) - 1 = " + np.array2string(gam / (4 * r * r1 * G1 / N) - 1, precision=2))
    # continuidad en lm
    e = 1e-6
    Vi = T.fondo(np.array([F.lm - e, F.lm + e]))
    for k in ("N", "r"):
        log(f"  continuidad en lm de {k}: " + ", ".join(f"d^{j}: {Vi[k].deriv(j)[0]:+.6e} / {Vi[k].deriv(j)[1]:+.6e}"
                                                         for j in range(4)))
    # K > 0 y ausencia de fantasmas con fluido
    lt = np.concatenate([np.linspace(1e-3, F.ls, 600), np.geomspace(F.ls, 3000, 600)[1:]])
    Ks, cond = [], []
    for i in range(0, lt.size, 300):
        Vt = T.V(lt[i:i + 300])
        Kt = Vt["K"].c[0]
        xi = Vt["xi"].c[0]
        rp = Vt["rp"].c[0]
        Nt = Vt["N"].c[0]
        Gct = Vt["Gc"].c[0]
        Ks.append(Kt)
        cond.append(Nt * Nt * Kt - Nt * rp * xi * xi / (2 * Gct * Gct))
    Ks, cond = np.concatenate(Ks), np.concatenate(cond)
    log(f"  min K = {Ks.min():+.3e} (en l = {lt[np.argmin(Ks)]:.2f});  min [N^2 K - N (rho~+p~) xi^2/(2 Gcal^2)] = "
        f"{cond.min():+.3e} (en l = {lt[np.argmin(cond)]:.2f})")
