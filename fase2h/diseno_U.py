# -*- coding: utf-8 -*-
"""
Fase 2h: el potencial del monopolo en la region sin fluido (entre la estrella y la zona lejana) y su diseno.

Con G3 = 0 y G4 = G4(X) (Gamma2 = 0), en la region de vacio y con Z = N xi (xi = r Gcal H/(r' Gamma1)), las identidades T3
de la Fase 2f dan  A1 e1 = 4 N Gcal (1 + r Phi'/r')  y  A1 e0 = -4 N Gcal (r'/r + Phi') + 4 Z/r^2,  de modo que
    U = U_fondo - (2/N^2) (N^2 Z/r^2)' + (1/N) (N Phi' P)'                                         (*)
con U_fondo dependiente solo del fondo (N, r, Gcal).  Las dos funciones libres que quedan (Z, es decir K, y P) entran
solo como derivadas totales con pesos N^2 y N: pueden redistribuir U, no subirlo en conjunto.
Con (*) la condicion U >= 0 en [la, lb] se integra:
    N^2 Z/r^2 (l) >= N^2 Z/r^2 (lb) - int_l^lb [ N^2 U_fondo/2 + (N/2)(N Phi' P)' ] dl'      (cota inferior de Z)
y K > 0 (Z' > N G4) con Z(lb) fijo da la cota superior Z(l) <= Z(lb) - int_l^lb N G4.
Este modulo calcula U_fondo (comprobando (*) con dos disenos distintos) y las dos cotas.
"""
import os
import sys
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2f"), os.path.join(AQUI, "..", "fase2g"), AQUI):
    sys.path.insert(0, p_)
import completo as CO                 # noqa: E402
import jets as J                      # noqa: E402


def perfil_vacio(T, ls, trozo=400):
    """K, P, U (Sturm-Liouville, T3) y U_fondo segun (*) en la region sin fluido."""
    acc = {}
    for i in range(0, ls.size, trozo):
        l = ls[i:i + trozo]
        V = T.V(l)
        N, r, r1, Gc = V["N"], V["r0"], V["r1"], V["Gc"]
        gam, al1 = CO.COEF["gam"](V), CO.COEF["al1"](V)
        kc, e0, e1 = CO.COEF["k_c"](V), CO.COEF["e0"](V), CO.COEF["e1"](V)
        A1 = al1 / gam
        A0 = -V["Phi1"] * A1
        u_d = -(2 / r) * (N * r1 * Gc).d()
        p_d = -2 * N * Gc
        P = p_d + A1 * e1 - (kc / 2) * A1 * A1
        U = u_d - (kc / 2) * A0 * A0 + A0 * e0 + 0.5 * (kc * A1 * A0 - A1 * e0 - A0 * e1).d()
        Z = V["Zj"]
        Uz = -(2 / (N * N)) * (N * N * Z / (r * r)).d()
        Up = (1 / N) * (N * V["Phi1"] * P).d()
        Ub = U - Uz - Up
        # comprobacion de A1 e1 y A1 e0 (Gamma2 = 0, vacio)
        a1e1 = 4 * N * Gc * (1 + r * V["Phi1"] / r1)
        a1e0 = -4 * N * Gc * (r1 / r + V["Phi1"]) + 4 * Z / (r * r)
        d = dict(l=l, N=N.c[0], Phi1=V["Phi1"].c[0], r=r.c[0], r1=r1.c[0], G4=V["G4"].c[0], Gc=Gc.c[0],
                 Z=Z.c[0], K=V["K"].c[0], P=P.c[0], U=U.c[0], Ub=Ub.c[0], Uz=Uz.c[0], Up=Up.c[0],
                 err_a1e1=np.abs((A1 * e1 - a1e1).c[0]) / np.abs(a1e1.c[0]),
                 err_a1e0=np.abs((A1 * e0 - a1e0).c[0]) / (np.abs(a1e0.c[0]) + 1e-30), rp=V["rp"].c[0])
        if "Zn" in V:
            # U_0 = U_fondo + U_z(Z_nat): el U de la teoria natural sin el termino de P (y = N es supersolucion si
            # U_0 + (2/r^2)(2 sig D - D') >= 0, con cualquier P > 0)
            Zn = V["Zn"]
            U0 = Ub - (2 / (N * N)) * (N * N * Zn / (r * r)).d()
            d.update(U0=U0.c[0], Zn=Zn.c[0], Kn=V["Kn"].c[0], D=(Z - Zn).c[0], Dp=(Z - Zn).d().c[0],
                     sig=(r1 / r - V["Phi1"]).c[0], Pnat=V["P_nat"].c[0])
        for k, v in d.items():
            acc.setdefault(k, []).append(np.asarray(v))
    return {k: np.concatenate(v) for k, v in acc.items()}


def cotas(p, la_idx, lb_idx):
    """Cota inferior (U >= 0, sin el termino de P) y superior (K >= 0) de Z entre ls[la_idx] y ls[lb_idx]."""
    l = p["l"]
    N, r, G4 = p["N"], p["r"], p["G4"]
    sl = slice(la_idx, lb_idx + 1)
    ll = l[sl]
    W_b = (N * N * p["Z"] / r ** 2)[lb_idx]
    f = (N * N * p["Ub"] / 2)[sl]
    # int_l^lb f
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(ll))])
    int_l_lb = cum[-1] - cum
    Zlb = (r ** 2 / (N * N))[sl] * (W_b - int_l_lb)
    g = (N * G4)[sl]
    cumg = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(ll))])
    Zub = p["Z"][lb_idx] - (cumg[-1] - cumg)
    return ll, Zlb, Zub


def exterior_dirichlet(l, K, P, U, k=4):
    """Autovalores minimos de -(P d')' + U d = w K d en [l0, l1] con d = 0 en los extremos (elementos lineales,
    masa concentrada): el problema exterior desacoplado de la estrella."""
    from scipy.linalg import eigh
    h = np.diff(l)
    lm = 0.5 * (l[1:] + l[:-1])
    Pm, Um = np.interp(lm, l, P), np.interp(lm, l, U)
    n = l.size
    A = np.zeros((n, n))
    M = np.zeros((n, n))
    for e in range(n - 1):
        A[e:e + 2, e:e + 2] += Pm[e] / h[e] * np.array([[1, -1], [-1, 1]]) \
            + Um[e] * h[e] / 6 * np.array([[2, 1], [1, 2]])
        M[e, e] += 0.5 * h[e] * K[e]
        M[e + 1, e + 1] += 0.5 * h[e] * K[e + 1]
    w, v = eigh(A[1:-1, 1:-1], M[1:-1, 1:-1], subset_by_index=[0, k - 1])
    return w, v
