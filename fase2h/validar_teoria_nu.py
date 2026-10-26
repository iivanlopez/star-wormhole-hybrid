# -*- coding: utf-8 -*-
"""
Fase 2h: validacion de la reconstruccion teoria_nu.TeoriaNu sobre un fondo fondo_nu.FondoNu:
  1. ecuaciones del fondo con fluido (E_A, E_B, E_C, E_p) sobre los chorros de la teoria: deben anularse;
  2. K de la formula T2 frente a la K de diseno, y gam frente a 4 r r' Gamma1/N;
  3. chorros del fondo frente a la solucion del ODE (N, N', r, r');
  4. K > 0 en todas partes y ausencia de fantasmas con fluido dentro de la estrella;
  5. perfil exterior: U, U_0, D, sig, K, K_nat (diseno_U.perfil_vacio).
Uso: python validar_teoria_nu.py [json de optimiza_nu]   (sin argumento: el fondo de prueba de fondo_nu)
"""
import os
import sys
import json
import time
import pickle
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2f"), os.path.join(AQUI, "..", "fase2g"), AQUI):
    sys.path.insert(0, p_)
import simbolico_num as SN            # noqa: E402
import reconstruccion as R            # noqa: E402
import fondo_nu as FN                 # noqa: E402
import teoria_nu as TN                # noqa: E402
import diseno_U as DU                 # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:6.1f} s] {m}", flush=True)


fl = pickle.load(open(os.path.join(AQUI, "..", "fase2g", "resultados", "fluido_simbolico.pkl"), "rb"))
FONDO = {k: SN.compilar(SN.con_Xb(fl["fondo_eq"][k])) for k in ("A", "B", "C", "p")}

if __name__ == "__main__":
    if len(sys.argv) > 1:
        d = json.load(open(sys.argv[1]))
        F = FN.FondoNu(**d["x"])
    else:
        F = FN.FondoNu(wt=1.0, nu0=2.0, b=3.0, n=2, d0=0.5, wd=1.0, lam=1.0, lT=30.0, A=0.05, Ns=1.25, phi1=0.5)
    log(f"fondo: ok {F.ok}, kappa = {F.kappa:.6f}, ls = {F.ls:.4f}, M_adm = {F.M_adm:.4f}, N0 = {F.N0:.4f}, "
        f"r sig(ls) = {F.margen_exterior():+.4f}")
    T = TN.TeoriaNu(F)
    log(f"teoria: l1 = {T.l1:.2f}, l2 = {T.l2:.2f}, LD = [{T.LDa:.1f}, {T.LD:.1f}], B = {T.B}, "
        f"exceso sin B = {T.exceso_sin_B:+.4f}, Z_obj = {T.Z_obj:.4f}, min K_nat exterior = {T.Kn_min_ext:+.3e}")
    if T.B is None:
        sys.exit(0)
    ls = F.ls
    l = np.array([0.05, 0.5, 0.3 * ls, 0.7 * ls, 0.99 * ls, 1.01 * ls, 1.1 * ls, 1.5 * ls, 3 * ls, 6 * ls, 10 * ls, 20 * ls,
                  1000.0, 2500.0])
    V = T.V(l)
    res = {k: np.abs(f(V).c[0]) for k, f in FONDO.items()}
    for i, li in enumerate(l):
        log(f"  l = {li:8.3f}: |E_A| = {res['A'][i]:.1e}, |E_B| = {res['B'][i]:.1e}, |E_C| = {res['C'][i]:.1e}, "
            f"|E_p| = {res['p'][i]:.1e}")
    W = R.Garganta._W(V)
    Kf = SN.MONO["K"](W).c[0]
    Kd = V["K"].c[0]
    gam = SN.MONO["gam"](V).c[0]
    r, r1, N = V["r0"].c[0], V["r1"].c[0], V["N"].c[0]
    log("  K (T2) / K (diseno) - 1 = " + np.array2string((Kf / Kd - 1), precision=2))
    log("  gam / (4 r r' Gamma1/N) - 1 = " + np.array2string(gam / (4 * r * r1 * V["G1"].c[0] / N) - 1, precision=2))
    Nu, Npu, rv, rpv = F.sol.sol(l)
    err = max(np.max(np.abs(N - Nu / F._Ninf)), np.max(np.abs(V["N"].c[1] - Npu / F._Ninf)),
              np.max(np.abs(r / rv - 1)), np.max(np.abs(r1 - rpv)))
    log(f"  chorros frente al ODE: max error {err:.1e}")
    lt = np.concatenate([np.linspace(1e-3, ls, 600), np.geomspace(ls, 3000, 600)[1:]])
    Ks, cond = [], []
    for i in range(0, lt.size, 300):
        Vt = T.V(lt[i:i + 300])
        Kt, xi, rp, Nt, Gct = (Vt[k].c[0] for k in ("K", "xi", "rp", "N", "Gc"))
        Ks.append(Kt)
        cond.append(Nt * Nt * Kt - Nt * rp * xi * xi / (2 * Gct * Gct))
    Ks, cond = np.concatenate(Ks), np.concatenate(cond)
    log(f"  min K = {Ks.min():+.3e} (en l = {lt[np.argmin(Ks)]:.2f});  min [N^2 K - N (rho~+p~) xi^2/(2 Gcal^2)] = "
        f"{cond.min():+.3e} (en l = {lt[np.argmin(cond)]:.2f})")
    le = np.concatenate([np.linspace(ls * (1 + 1e-6), 4 * ls, 300), np.geomspace(4 * ls, 0.95 * T.LD, 100)[1:]])
    p = DU.perfil_vacio(T, le)
    log(f"  exterior: max err A1 e1 = {p['err_a1e1'].max():.1e}, A1 e0 = {p['err_a1e0'].max():.1e}")
    for li in (1.001, 1.02, 1.05, 1.1, 1.2, 1.5, 2, 3, 4, 6):
        i = np.argmin(abs(le - li * ls))
        log(f"    l = {le[i]:8.3f}: sig = {p['sig'][i]:+.3e}  D = {p['D'][i]:+.3e}  D' = {p['Dp'][i]:+.3e}  "
            f"K = {p['K'][i]:+.3e}  K_nat = {p['Kn'][i]:+.3e}  P = {p['P'][i]:+.3e}  P_nat = {p['Pnat'][i]:+.3e}  "
            f"U = {p['U'][i]:+.3e}  U_0 = {p['U0'][i]:+.3e}")
