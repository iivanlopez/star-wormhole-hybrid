# -*- coding: utf-8 -*-
"""
Fase 2h: analisis completo de un candidato (fondo_nu.FondoNu) con el diseno de exterior natural (teoria_nu, theta = 0).
  1. fondo: kappa, M_adm, N0, ls, r_s/M, r sig(ls); criba de la Fase 2g (n = 160) y sector impar (criba.impar);
  2. teoria: B, ecuaciones del fondo y K de T2 sobre los chorros, K > 0 y cota de fantasmas con fluido;
  3. exterior: Sturm-Liouville con Dirichlet en la superficie (teoria natural);
  4. espectro completo de dos campos (espectro2 con superficie libre), dos mallas y dos cajas;
  5. velocidades caracteristicas radiales (velocidades.py).
Uso: python candidato_nu.py <json con "x"> [eps] [cm2] [b_in] [Lam]  -> resultados/candidato_<nombre>.json
"""
import os
import sys
import json
import time
import pickle
import warnings
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2f"), os.path.join(AQUI, "..", "fase2g"), AQUI):
    sys.path.insert(0, p_)
import simbolico_num as SN            # noqa: E402
import reconstruccion as R            # noqa: E402
import criba as C                     # noqa: E402
import fondo_nu as FN                 # noqa: E402
import teoria_nu as TN                # noqa: E402
import espectro2 as E2                # noqa: E402
import diseno_U as DU                 # noqa: E402
import velocidades as VE              # noqa: E402

warnings.filterwarnings("ignore")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


def analizar(x, eps=0.05, cm2=1.0, b_in=0.0, mallas=((120, 1500.0, 250), (200, 1500.0, 400), (120, 3000.0, 300)),
             nombre="candidato", Lam=0.0):
    res = dict(x=x, eps=eps, cm2=cm2, b_in=b_in, Lam=Lam)
    F = FN.FondoNu(**x)
    if not F.ok:
        log(f"fondo invalido: {F.motivo}")
        res["motivo"] = F.motivo
        return res
    f_s = F.campos(np.array([F.ls * (1 - 1e-9)]))
    rs = float(f_s["r"][0])
    res["fondo"] = dict(kappa=F.kappa, M_adm=F.M_adm, M_komar=F.M_komar, m_MS=F.m_MS, N0=F.N0, ls=F.ls, r_s=rs,
                        rs_sobre_M=rs / F.M_adm, r_sig_ls=F.margen_exterior(), rho_p0=float(F.campos(np.array([0.0]))["rp"][0]))
    log("fondo: " + ", ".join(f"{k} = {v:.5g}" for k, v in res["fondo"].items()))
    cr = C.criba_normalizada(F, n=160)
    res["criba"] = cr
    log("criba (superficie fija): " + ", ".join(f"{k} {v:+.4f}" for k, v in cr.items()))
    res["criba_libre"] = C.criba_normalizada(F, n=160, libre=True)
    res["sin_compresion"] = C.prueba_sin_compresion(F)
    log(f"criba (superficie libre): margen {res['criba_libre']['margen']:+.4f}; prueba sin compresion E_Q0 = "
        f"{res['sin_compresion']['E_Q0']:+.3e}")
    im = C.impar(F)
    res["impar"] = im
    log(f"sector impar: {im}")
    T = TN.TeoriaNu(F, theta=0.0, eps=eps, cm2=cm2, b_in=b_in, Lam=Lam)
    res["B"], res["exceso_sin_B"] = T.B, T.exceso_sin_B
    log(f"teoria: B = {T.B}, exceso sin B = {T.exceso_sin_B:+.4f}, min K_nat exterior = {T.Kn_min_ext:+.3e}")
    if T.B is None:
        return res
    # comprobaciones de la reconstruccion
    fl = pickle.load(open(os.path.join(AQUI, "..", "fase2g", "resultados", "fluido_simbolico.pkl"), "rb"))
    FONDO = {k: SN.compilar(SN.con_Xb(fl["fondo_eq"][k])) for k in ("A", "B", "C", "p")}
    ls = F.ls
    lch = np.array([0.05, 0.3 * ls, 0.7 * ls, 0.99 * ls, 1.05 * ls, 2 * ls, 10 * ls, 50 * ls])
    V = T.V(lch)
    res["err_fondo"] = float(max(np.max(np.abs(f(V).c[0])) for f in FONDO.values()))
    Kf = SN.MONO["K"](R.Garganta._W(V)).c[0]
    res["err_K_T2"] = float(np.max(np.abs(Kf / V["K"].c[0] - 1)))
    lt = np.concatenate([np.linspace(1e-3, ls, 800), np.geomspace(ls, 3000, 600)[1:]])
    Ks, cond = [], []
    for i in range(0, lt.size, 350):
        Vt = T.V(lt[i:i + 350])
        Kt, xi, rp, Nt, Gct = (Vt[k].c[0] for k in ("K", "xi", "rp", "N", "Gc"))
        Ks.append(Kt)
        cond.append(Nt * Nt * Kt - Nt * rp * xi * xi / (2 * Gct * Gct))
    Ks, cond = np.concatenate(Ks), np.concatenate(cond)
    res["min_K"], res["min_cota_fantasma"] = float(Ks.min()), float(cond.min())
    log(f"comprobaciones: ecuaciones del fondo {res['err_fondo']:.1e}, K(T2)/K - 1 = {res['err_K_T2']:.1e}, "
        f"min K = {Ks.min():+.3e}, min [N^2 K - cota] = {cond.min():+.3e}")
    # exterior
    le = np.concatenate([np.linspace(ls * (1 + 1e-6), 5 * ls, 700), np.geomspace(5 * ls, 1500, 700)[1:]])
    p = DU.perfil_vacio(T, le)
    w, _ = DU.exterior_dirichlet(le, p["K"], p["P"], p["U"])
    res["exterior_w"] = [float(v) for v in w]
    res["exterior_minU0"] = float(p["U0"].min())
    log(f"exterior (Dirichlet en ls, L = 1500): w = {np.array2string(w, precision=3)}; min U_0 = {p['U0'].min():+.2e}")
    # espectro completo
    lc = T.lc()
    res["espectro"] = []
    for n_in, L, n_out in mallas:
        a, b = E2.mallas(lc, n_in, L, n_out)
        cf = E2.Coefs(T, np.concatenate([a, b[1:]]))
        for libre in (True, False):
            for par in ("par", "impar"):
                o = E2.espectro(T, par, lc, n_in, L, n_out, k=6, coefs=cf, libre=libre)
                o.update(par=par, n_in=n_in, L=L, n_out=n_out, libre=libre)
                res["espectro"].append(o)
                log(f"espectro {par}, superficie {'libre' if libre else 'fija'} (n_in {n_in}, L {L:.0f}, n_out {n_out}): "
                    f"Tmin = {o['Tmin']:.2e}, "
                    + (f"w = {np.array2string(np.array(o['w']), precision=4)}" if "w" in o else f"fantasmas {o['n_T_neg']}"))
    # velocidades
    lv = np.concatenate([np.geomspace(1e-2, ls, 300)[:-1], [0.999 * ls]])
    v = VE.velocidades(T, lv)
    res["vel"] = dict(c2_rapida_max=float(np.nanmax(v["c2_rapida"])), c2_rapida_min=float(np.nanmin(v["c2_rapida"])),
                      c2_lenta_max=float(np.nanmax(v["c2_lenta"])), detK_min=float(np.nanmin(v["detK_rel"])),
                      c2_rapida_l1=float(np.interp(1.0, lv, v["c2_rapida"])),
                      c2_rapida_media=float(np.interp(0.5 * ls, lv, v["c2_rapida"])))
    log("velocidades (dentro): " + ", ".join(f"{k} {val:.3g}" for k, val in res["vel"].items()))
    json.dump(res, open(os.path.join(AQUI, "resultados", f"{nombre}.json"), "w"), indent=1, default=float)
    return res


if __name__ == "__main__":
    d = json.load(open(sys.argv[1]))
    eps = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    cm2 = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    b_in = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    Lam = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
    nombre = "candidato_" + os.path.splitext(os.path.basename(sys.argv[1]))[0]
    analizar(d["x"], eps, cm2, b_in, nombre=nombre, Lam=Lam,
             mallas=((240, 1500.0, 250), (400, 1500.0, 400), (240, 3000.0, 300)))
