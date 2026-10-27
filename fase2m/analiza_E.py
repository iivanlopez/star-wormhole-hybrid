# -*- coding: utf-8 -*-
"""
Fase 2m: refina las gargantas del barrido (garganta_E.py), comprueba que son gargantas de verdad (r''(0) > 0,
U -> 0, N -> N_inf > 0, r creciente) y evalua sobre ellas:
  - condiciones: G4 > 0, Gcal = G4 - 2X G4_X > 0, g = G2_F + 2F G2_FF > 0 (aqui 1), m_at > 0 (masa estatica de at,
    ~ termino cinetico del modo longitudinal, at ~ d_t pi) y c_ww < 0 (coeficiente de wt^2, ~ -gradiente, wt ~ pi');
  - la energia de la prueba d = N:  E[N] = int [vol + ele + J_N^2/(4 m_at)],
      vol = -2 N^3 G4 sig^2,   ele = (2N/g)(2 G4_X N (U/N)' - G2_F U')^2   (con la cola -2 N_inf^3/L mas alla de L).
Uso: python analiza_E.py   (lee resultados/garganta_E_barrido.json; escribe resultados/analiza_E.json, perfiles_E.npz)
"""
import os
import sys
import json
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import garganta_E as G                                  # noqa: E402


def refina(p, a, b, it=60):
    sa = G.tiro(a, p)
    for _ in range(it):
        m = 0.5 * (a + b)
        sm = G.tiro(m, p)
        if sm == sa:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def perfil(X0, p):
    s = G.integra(X0, p, L=60.0, denso=True)
    U = s.y[4]
    # corte: donde |U| es minimo antes de que el modo creciente domine
    i = int(np.argmin(np.abs(U) + 1e-300 * 0)) if s.t_events[1].size == 0 else len(s.t) - 1
    i = max(i, 10)
    Lc = s.t[i]
    l = np.linspace(0, Lc, 6001)
    y = s.sol(l)
    return l, y, Lc


def diagnostico(X0, p):
    l, y, Lc = perfil(X0, p)
    N, N1, r, r1, U, U1 = y
    d2 = np.array([G.rhs(0, y[:, k], p) for k in range(len(l))]).T
    N2, r2, U2 = d2[1], d2[3], d2[5]
    N3, r3, U3 = np.gradient(N2, l), np.gradient(r2, l), np.gradient(U2, l)
    Ninf = N[-1]
    args = (N, N1, r, r1, U, U1, N2, r2, U2, N3, r3, U3) + tuple(p)
    vol, ele, JN, mt, Gc, G4, cww = [np.broadcast_to(np.asarray(x, dtype=float), l.shape) for x in G.f_dens(*args)]
    corr = JN ** 2 / (4 * mt)
    w = np.full_like(l, l[1] - l[0])
    w[0] = w[-1] = w[0] / 2
    cola = -2 * Ninf ** 3 / r[-1]                         # int_L^inf -2 N^3 / r^2 dl  con r ~ l
    E = dict(vol=float((vol * w).sum() + cola), ele=float((ele * w).sum()), corr_at=float((corr * w).sum()))
    E["total"] = E["vol"] + E["ele"] + E["corr_at"]
    # normalizacion N_inf = 1: las tres densidades escalan como N^3
    E = {k: v / Ninf ** 3 for k, v in E.items()}
    rpp0 = G.rhs(0, y[:, 0], p)[3]
    return dict(p=list(p), X0=X0, r0=float(r[0]), Lc=float(Lc), N0_sobre_Ninf=float(N[0] / Ninf), rpp0=float(rpp0),
                r_monotona=bool(np.all(r1[1:] > 0)), U_final=float(U[-1] / U[0]), minN=float(N.min() / Ninf),
                min_Gcal=float(Gc.min()), min_G4=float(G4.min()), min_mt=float(mt.min()), max_mt=float(mt.max()),
                mt0=float(mt[0]), max_cww=float(cww.max()), cww0=float(cww[0]),
                mt_cambia_signo=bool(mt.min() < 0 < mt.max()), E=E), (l, y, vol, ele, corr, mt)


if __name__ == "__main__":
    br = json.load(open(os.path.join(AQUI, "resultados", "garganta_E_barrido.json")))
    res, perf = [], {}
    for k, o in enumerate(br):
        p = tuple(o["p"])
        X0 = refina(p, o["X0a"], o["X0b"])
        try:
            dg, pr = diagnostico(X0, p)
        except Exception as ex:                          # noqa: BLE001
            print("fallo", p, X0, ex)
            continue
        res.append(dg)
        perf[f"c{k}"] = np.vstack([pr[0], pr[1], pr[2], pr[3], pr[4], pr[5]])
        print(json.dumps(dg), flush=True)
    json.dump(res, open(os.path.join(AQUI, "resultados", "analiza_E.json"), "w"), indent=1)
    np.savez_compressed(os.path.join(AQUI, "resultados", "perfiles_E.npz"),
                        **{k: v[:, ::10].astype(np.float32) for k, v in perf.items()})   # 1 de cada 10 puntos (float32: los perfiles casi horizonte desbordan, no se usan)
