# -*- coding: utf-8 -*-
"""
Fase 2m (continuacion): la esfera ct = 0 y el fantasma en las gargantas electricas regulares (analiza_E.json).

Con la reduccion completa (reduccion_completa_E.py) la energia cinetica del monopolo es
    T = int [ -(2 Gcal/N) d_t^2 + k2 ((Gcal d_t/N)')^2 ],     k2 = -4 r^2 N^4/(U^2 ct).
Sobre cada garganta: ct(l), mt(l), k2; la region ct > 0 alrededor de la garganta (donde T < 0 para toda d con soporte
en ella), la primera esfera ct = 0 (k2 diverge y cambia de signo) y las esferas mt = 0 (el potencial
J^2/(4 mt) diverge).  Comprueba ct(0) = -N0 Theta/2 (longitudinal_E.py) y da T para una perturbacion de prueba
d_t = (1 - (l/s)^2)^3 en l < s = 0.95 l_c (soporte compacto en ct > 0).
Uso: python esfera_E.py   (escribe resultados/esfera_E.json)
"""
import os
import sys
import json
import pickle
import numpy as np
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import garganta_E as G                                  # noqa: E402
from chorros_proca import NN, rr, UU                    # noqa: E402

rc = pickle.load(open(os.path.join(AQUI, "resultados", "reduccion_completa_E.pkl"), "rb"))
ct = G.T(rc["ct"])
mt = G.T(rc["mt"])
Gc = G.T(rc["Gc"])
f = sp.lambdify(G.ARGS2, [ct, mt, Gc], "numpy")
# comprobacion del modo gauge residual d = r': coeficientes de L_s (cuadratica en d, d', d'') y de J
dj = [sp.Symbol(f"d_0{k}", real=True) for k in range(3)]
Ls = G.T(rc["Ls"])
Jt = G.T(rc["Jtot"])
cLs = [[sp.expand(sp.diff(Ls, dj[i], dj[j])) for j in range(3)] for i in range(3)]   # Ls = 1/2 d^T C d
jk = [sp.expand(Jt.coeff(dj[k])) for k in range(3)]
f_c = sp.lambdify(G.ARGS2, [cLs[i][j] for i in range(3) for j in range(3)] + jk, "numpy")


def ceros(l, y):
    s = np.sign(y)
    i = np.where(s[1:] * s[:-1] < 0)[0]
    return [float(l[k] - y[k] * (l[k + 1] - l[k]) / (y[k + 1] - y[k])) for k in i]


if __name__ == "__main__":
    import analiza_E as A
    res = []
    for o in json.load(open(os.path.join(AQUI, "resultados", "analiza_E.json"))):
        if not (o["r_monotona"] and o["rpp0"] > 0 and abs(o["U_final"]) < 1e-5 and o["N0_sobre_Ninf"] < 10):
            continue
        p = tuple(o["p"])
        l, y, Lc = A.perfil(o["X0"], p)
        N, N1, r, r1, U, U1 = y
        d2 = np.array([G.rhs(0, y[:, k], p) for k in range(len(l))]).T
        N2, r2, U2 = d2[1], d2[3], d2[5]
        z = np.zeros_like(l)
        cts, mts, Gcs = [np.broadcast_to(np.asarray(x, float), l.shape) for x in
                         f(N, N1, r, r1, U, U1, N2, r2, U2, z, z, z, *p)]
        b, m2, l2, l3 = p
        X0 = o["X0"]
        Th = (m2 + 2 * l2 * X0 + 3 * l3 * X0 ** 2) * o["r0"] ** 2 + 2 * b
        zc, zm = ceros(l, cts), ceros(l, mts)
        lc = zc[0] if zc else np.inf
        # prueba localizada en la region ct > 0
        s = 0.95 * lc if np.isfinite(lc) else 0.5 * l[-1]   # soporte compacto dentro de ct > 0
        dt = np.where(l < s, (1 - (l / s) ** 2) ** 3, 0.0)
        k2 = -4 * r ** 2 * N ** 4 / (U ** 2 * cts)
        Qd = Gcs * dt / N
        dQ = np.gradient(Qd, l)
        dens = -2 * Gcs / N * dt ** 2 + k2 * dQ ** 2
        w = np.gradient(l)
        Tbump = float(2 * (dens * w).sum())            # dos lados de la garganta
        resto = float(2 * (np.maximum(dens, 0) * w).sum())   # parte positiva (debe ser 0)
        # K r' y V r' por diferencias finitas en la zona sin la esfera mt = 0
        g_ = lambda x: np.gradient(x, l)                     # noqa: E731
        N3, r3, U3 = g_(N2), g_(r2), g_(U2)
        cc = [np.broadcast_to(np.asarray(x, float), l.shape) for x in
              f_c(N, N1, r, r1, U, U1, N2, r2, U2, N3, r3, U3, *p)]
        C = np.array(cc[:9]).reshape(3, 3, -1)
        j0, j1, j2 = cc[9:]
        x = [r1, r2, g_(r2)]                                  # d = r' y sus derivadas
        pk = [sum(C[k, j] * x[j] for j in range(3)) for k in range(3)]
        Es = pk[0] - g_(pk[1]) + g_(g_(pk[2]))                # dLs/dd
        ph = (j0 * x[0] + j1 * x[1] + j2 * x[2]) / (2 * mts)
        Jd = j0 * ph - g_(j1 * ph) + g_(g_(j2 * ph))
        Vr = -Es + Jd
        Qr = Gcs * r1 / N
        Kr = -4 * Gcs / N * r1 - 2 * Gcs / N * g_(k2 * g_(Qr))
        lm = zm[0] if zm else l[-1]
        zona = (l > 0.05 * lc) & (l < 0.8 * min(lc, lm))
        relK = float(np.max(np.abs(Kr[zona])) / np.max(np.abs(4 * Gcs / N * r1)[zona]))
        # control: la misma medida con d = N (no es un modo gauge) debe dar O(1)
        xN = [N, N1, N2]
        pN = [sum(C[k, j] * xN[j] for j in range(3)) for k in range(3)]
        EsN = pN[0] - g_(pN[1]) + g_(g_(pN[2]))
        phN = (j0 * xN[0] + j1 * xN[1] + j2 * xN[2]) / (2 * mts)
        JdN = j0 * phN - g_(j1 * phN) + g_(g_(j2 * phN))
        KN = -4 * Gcs / N * N - 2 * Gcs / N * g_(k2 * g_(Gcs))
        ctrlK = float(np.max(np.abs(KN[zona])) / np.max(np.abs(4 * Gcs)[zona]))
        ctrlV = float(np.max(np.abs((-EsN + JdN)[zona])) / np.max((np.abs(EsN) + np.abs(JdN))[zona]))
        relV = float(np.max(np.abs(Vr[zona])) / np.max((np.abs(Es) + np.abs(Jd))[zona]))
        fila = dict(p=list(p), gauge_Kr_rel=relK, gauge_Vr_rel=relV, control_KN_rel=ctrlK, control_VN_rel=ctrlV, X0=X0, r0=o["r0"], ct0=float(cts[0]), ct0_formula=float(-N[0] * Th / 2),
                    l_ct0=lc, r_ct0_sobre_r0=float(np.interp(lc, l, r) / r[0]) if np.isfinite(lc) else None,
                    ceros_ct=zc, ceros_mt=zm, mt0=float(mts[0]), min_Gcal=float(Gcs.min()),
                    T_prueba=Tbump, T_prueba_cola=resto, s=s)
        res.append(fila)
        print(json.dumps(fila), flush=True)
    json.dump(res, open(os.path.join(AQUI, "resultados", "esfera_E.json"), "w"), indent=1)
