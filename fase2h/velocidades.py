# -*- coding: utf-8 -*-
"""
Fase 2h: velocidades caracteristicas radiales del monopolo con fluido (simbolo principal del sistema de dos campos).

En las variables (xi, d) de la reduccion de la Fase 2g (xi: desplazamiento del fluido, W = -2 kap_b xi, d = N h):
    L = kap_b xi_t^2 + K d_t^2 + 2 X xi_t d_t - Q^2/(4 kap_a) - E_resto(c*, d),
    c* = (W - al1 d' - al0 d)/gam,  Q = (N^2 W)' + M,  M = -N (rho~+p~)(r^2 c* + 2 r d),
con X = kap_b e_tt/gam (de -e_tt N c_t h_t) y K la de la Fase 2f.  Las derivadas de mayor orden dan
    cinetica:  [[kap_b, X], [X, K]]
    gradiente: [[N^4 kap_b^2/kap_a, -N^3 kap_b (rho~+p~) r^2 A1/(2 kap_a)],
                [  .             , (N (rho~+p~) r^2 A1)^2/(4 kap_a) + P]],     A1 = al1/gam,
P = p_d + A1 e1 - (k_c/2) A1^2.  Las velocidades coordenadas v^2 son los autovalores generalizados y las propias
c^2 = v^2/N^2.  Sin acoplo (A1 = X = 0) son c_s^2 (fluido) y P/(N^2 K) (escalar).  Fuera del fluido c^2 = P/(N^2 K).
det(cinetica) > 0 es la condicion de ausencia de fantasmas con fluido (equivale a N^2 K >= N (rho~+p~) xi^2/(2 Gcal^2)).
"""
import numpy as np


def velocidades(T, l):
    l = np.atleast_1d(np.asarray(l, dtype=float))
    c, _ = T.coeficientes(l)
    N, r, rp, cs2 = c["N"], c["r"], c["rp"], c["cs2"]
    gam, al1, kc, e1, ett, K = c["gam"], c["al1"], c["kc"], c["e1"], c["ett"], c["K"]
    A1 = al1 / gam
    P = -2 * N * c["Gc"] + A1 * e1 - 0.5 * kc * A1 * A1
    out = dict(l=l, c2_vacio=P / (N * N * K), P=P, K=K)
    dentro = rp > 0
    c1 = np.full(l.size, np.nan)
    c2 = np.full(l.size, np.nan)
    detK = np.full(l.size, np.nan)
    if np.any(dentro):
        i = dentro
        kb = r[i] ** 2 * rp[i] / (2 * N[i])
        ka = N[i] * r[i] ** 2 * rp[i] / (2 * cs2[i])
        X = kb * ett[i] / gam[i]
        B = N[i] * rp[i] * r[i] ** 2 * A1[i]
        g11 = N[i] ** 4 * kb ** 2 / ka
        g12 = -N[i] ** 3 * kb * rp[i] * r[i] ** 2 * A1[i] / (2 * ka)
        g22 = B * B / (4 * ka) + P[i]
        k11, k12, k22 = kb, X, K[i]
        # det(G - v^2 Kc) = 0:  a v^4 - b v^2 + c = 0
        a = k11 * k22 - k12 * k12
        b = g11 * k22 + g22 * k11 - 2 * g12 * k12
        cc = g11 * g22 - g12 * g12
        disc = np.sqrt(np.maximum(b * b - 4 * a * cc, 0.0))
        v1 = (b - disc) / (2 * a)
        v2 = (b + disc) / (2 * a)
        c1[i] = v1 / N[i] ** 2
        c2[i] = v2 / N[i] ** 2
        detK[i] = a / (k11 * k22)
    out.update(c2_lenta=c1, c2_rapida=c2, detK_rel=detK, cs2=cs2)
    return out
