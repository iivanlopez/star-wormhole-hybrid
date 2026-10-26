# -*- coding: utf-8 -*-
"""
Fase 2i: criba de la familia c = 0 (sin funciones libres de la teoria) para materia elastica en el monopolo de la GLPV
cuartica, con la energia del solido en forma lagrangiana (comovil; reduccion_solido.py, parte C).

Familia c = 0: d = N h, delta g_ll = 0 y el material se mueve para cumplir la ligadura de momento,
    xi = 4 N Gcal h'/(r w),   w = rho~ + p~_r          (en el vacio W = 0 exige h' = 0: h constante)
    T[h] = int [ 8 N Gcal^2 h'^2/w - 2 G4 N h^2 ]
    E[h] = int [ -2 N^3 Gcal (h' - sig h)^2 - L_lag(xi, xi', d, d') ]
con L_lag = r^2 N [ G0 (z2 + z1^2/2) + G_r e_r2 + G_t e_t2 + (G_r e_r1 + G_t e_t1) z1 + (G_rr e_r1^2 + 2 G_rt e_r1 e_t1
                     + G_tt e_t1^2)/2 ],
    e_r1 = xi',  e_r2 = -xi'^2/2,  e_t1 = (r' xi + d)/r,
    e_t2 = [(r'^2 + r r'') xi^2 + 2 r' d xi + 2 r d' xi]/(2 r^2) - (r' xi + d)^2/r^2,
    z1 = e_r1 + 2 e_t1 + Phi' xi,  z2 = e_r2 + 2 e_t2 + Phi'' xi^2/2.
La integral es sobre el dominio material: las superficies libres no llevan terminos de superficie.
Modulos (deformaciones logaritmicas; e = energia por elemento = e_v(eps_r + 2 eps_t) + e_s(eps_r - eps_t)):
    S = e_v'' [[1, 2], [2, 4]] + e_s'' [[1, -1], [-1, 1]],   e_v'' = w c_L^2 - p~_r - (4/3) w c_T^2,  e_s'' = (4/3) w c_T^2,
    c_L^2 = (S_rr + p~_r)/w (velocidad longitudinal radial), c_T^2 = mu/w (cizalla);
    G0 = -rho~, G_r = w, G_t = 2 (rho~ + p~_t), G_rr = -S_rr - rho~ - 2 p~_r, G_rt = -S_rt - 2 rho~ - 2 p~_r - 2 p~_t,
    G_tt = -S_tt - 4 rho~ - 8 p~_t.   Fluido: c_T = 0 (S = (w c_s^2 - p~)[[1, 2], [2, 4]]).
Dominios: hueco de vacio [0, la] (h = h(la) constante; impar: h(la) = 0), material [la, lb] (elementos de Hermite
cubicos en h), exterior (h = h(lb), energia E_out h(lb)^2).  Con la = 0 la materia llega a la garganta (par: h'(0) = 0).
"""
import numpy as np
import sympy as sp
from scipy.linalg import eigh

GX, GW = np.polynomial.legendre.leggauss(8)


# ------------------------------------------------------------------ forma lagrangiana: matriz 4x4 en v = (xi, xi', d, d')
def _construir_forma():
    xi, xip, dd, ddp = sp.symbols("xi xip dd ddp", real=True)
    r, r1, r2, P1, P2, G0, Gr, Gt, Grr, Grt, Gtt = sp.symbols("r r1 r2 P1 P2 G0 Gr Gt Grr Grt Gtt", real=True)
    er1 = xip
    er2 = -xip ** 2 / 2
    et1 = (r1 * xi + dd) / r
    et2 = ((r1 ** 2 + r * r2) * xi ** 2 + 2 * r1 * dd * xi + 2 * r * ddp * xi) / (2 * r ** 2) - (r1 * xi + dd) ** 2 / r ** 2
    z1 = er1 + 2 * et1 + P1 * xi
    z2 = er2 + 2 * et2 + P2 * xi ** 2 / 2
    L = sp.expand(G0 * (z2 + z1 ** 2 / 2) + Gr * er2 + Gt * et2 + (Gr * er1 + Gt * et1) * z1
                  + (Grr * er1 ** 2 + 2 * Grt * er1 * et1 + Gtt * et1 ** 2) / 2)
    v = [xi, xip, dd, ddp]
    M = [[sp.diff(L, a, b) / 2 for b in v] for a in v]
    args = (r, r1, r2, P1, P2, G0, Gr, Gt, Grr, Grt, Gtt)
    return [[sp.lambdify(args, M[i][j], "numpy") for j in range(4)] for i in range(4)]


_FORMA = _construir_forma()


def forma_lag(p):
    """Matriz M (4, 4, npts) con L_lag/(r^2 N) = v^T M v, v = (xi, xi', d, d')."""
    args = (p["r"], p["r1"], p["r2"], p["Phi1"], p["Phi2"], p["G0"], p["Gr"], p["Gt"], p["Grr"], p["Grt"], p["Gtt"])
    n = np.size(p["r"])
    M = np.empty((4, 4, n))
    for i in range(4):
        for j in range(4):
            M[i, j] = np.broadcast_to(_FORMA[i][j](*args), (n,))
    return M


def modulos(rho, pr, pt, cL2, cT2):
    """(G0, G_r, G_t, G_rr, G_rt, G_tt) y (e_v'', e_s'') del solido isotropo de velocidades c_L, c_T."""
    w = rho + pr
    Es = 4.0 / 3.0 * cT2 * w
    Srr = w * cL2 - pr
    Ev = Srr - Es
    Srt = 2 * Ev - Es
    Stt = 4 * Ev + Es
    return dict(G0=-rho, Gr=w, Gt=2 * (rho + pt), Grr=-Srr - rho - 2 * pr, Grt=-Srt - 2 * rho - 2 * pr - 2 * pt,
                Gtt=-Stt - 4 * rho - 8 * pt, Ev=Ev, Es=Es)


# ------------------------------------------------------------------ elementos de Hermite
def hermite(t, H):
    v = np.array([2 * t ** 3 - 3 * t ** 2 + 1, H * (t ** 3 - 2 * t ** 2 + t), -2 * t ** 3 + 3 * t ** 2, H * (t ** 3 - t ** 2)])
    d1 = np.array([6 * t ** 2 - 6 * t, H * (3 * t ** 2 - 4 * t + 1), -6 * t ** 2 + 6 * t, H * (3 * t ** 2 - 2 * t)]) / H
    d2 = np.array([12 * t - 6, H * (6 * t - 4), -12 * t + 6, H * (6 * t - 2)]) / H ** 2
    return v, d1, d2


def formas(P, nodos):
    """Matrices globales (h, h' por nodo) de E y T en el material (sin hueco ni exterior) y las de referencia
    positivas R (de E: int 2 N^3 Gcal [(h' - sig h)^2 + h^2/r^2]) y RT (de T, con signos +)."""
    K = len(nodos)
    nd = 2 * K
    out = {k: np.zeros((nd, nd)) for k in ("E", "T", "R", "RT", "M")}
    for e in range(K - 1):
        x0, x1 = nodos[e], nodos[e + 1]
        H = x1 - x0
        t = 0.5 * (GX + 1)
        w = 0.5 * GW * H
        p = P.perfil(x0 + t * H)
        v, d1, d2 = hermite(t, H)
        idx = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        N, Gc, sig, r = p["N"], p["Gc"], p["sig"], p["r"]
        q = 4 * N * Gc / (r * p["w"])
        qp = q * (p["Phi1"] + p["Gcp"] / Gc - p["r1"] / r - p["wp"] / p["w"])
        # v = (xi, xi', d, d') en funcion de (h, h', h'') en los puntos: filas de la transformacion por funcion base
        xi_b = q * d1
        xip_b = qp * d1 + q * d2
        d_b = N * v
        dp_b = N * p["Phi1"] * v + N * d1
        V = np.stack([xi_b, xip_b, d_b, dp_b])                  # (4, 4 bases, npts)
        Mq = forma_lag(p)                                       # (4, 4, npts)
        pre = -(r ** 2) * N * w
        Lm = np.einsum("ian,ijn,jbn,n->ab", V, Mq, V, pre)      # -int L_lag
        D = d1 - sig * v
        B = 2 * N ** 3 * Gc
        out["E"][np.ix_(idx, idx)] += Lm - (D * B * w) @ D.T
        out["R"][np.ix_(idx, idx)] += (D * B * w) @ D.T + ((v / r) * B * w) @ (v / r).T
        wt = 8 * N * Gc ** 2 / p["w"]
        out["T"][np.ix_(idx, idx)] += (d1 * wt * w) @ d1.T - (v * 2 * p["G4"] * N * w) @ v.T
        out["RT"][np.ix_(idx, idx)] += (d1 * wt * w) @ d1.T + (v * 2 * p["G4"] * N * w) @ v.T
        out["M"][np.ix_(idx, idx)] += (v * w) @ v.T
    return out


def criba(P, n=160, devolver=False):
    """Margenes adimensionales de la familia c = 0 (positivos = la condicion se cumple):
    E1n = E[1]/D[1];  Epar, Eimpar = autovalor minimo de E frente a D (superficies libres, exterior incluido);
    Tpar, Timpar = el de T frente a RT con h(lb) = 0 (soporte compacto: exterior sin energia cinetica)."""
    la, lb = P.la, P.lb
    nodos = np.linspace(la, lb, n + 1)
    fm = formas(P, nodos)
    K = len(nodos)
    Eo = P.E_out()
    Ig, IGg, Rg = P.integrales_hueco()          # int N^3 Gcal sig^2, int G4 N, int N^3 Gcal (sig^2 + 1/r^2) en [0, la]
    E, D, T, RT = fm["E"].copy(), fm["R"].copy(), fm["T"].copy(), fm["RT"].copy()
    i0, iL = 0, 2 * (K - 1)
    E[i0, i0] += -2 * Ig
    D[i0, i0] += 2 * Rg
    T[i0, i0] += -2 * IGg
    RT[i0, i0] += 2 * IGg
    E[iL, iL] += Eo
    D[iL, iL] += abs(Eo)
    uno = np.zeros(2 * K)
    uno[0::2] = 1.0
    res = dict(E1=float(uno @ E @ uno), E1n=float((uno @ E @ uno) / (uno @ D @ uno)), E_out=Eo, gap=-2 * Ig)
    for par in ("par", "impar"):
        if la > 0:
            quita = [] if par == "par" else [i0]
        else:
            quita = [1] if par == "par" else [0]
        lib = [i for i in range(2 * K) if i not in quita]
        ev, vec = eigh(E[np.ix_(lib, lib)], D[np.ix_(lib, lib)], subset_by_index=[0, 0])
        res[f"E{par}"] = float(ev[0])
        if devolver:
            u = np.zeros(2 * K)
            u[lib] = vec[:, 0]
            res[f"modo_{par}"] = (nodos, u[0::2], u[1::2])
        libT = [i for i in lib if i != iL]
        evT = eigh(T[np.ix_(libT, libT)], RT[np.ix_(libT, libT)], eigvals_only=True, subset_by_index=[0, 0])
        res[f"T{par}"] = float(evT[0])
    res["margen"] = min(res["E1n"], res["Epar"], res["Eimpar"], res["Tpar"], res["Timpar"])
    return res


# ------------------------------------------------------------------ adaptador para las estrellas fluidas de la Fase 2h
class PerfilFluido:
    """Fondo fondo_nu.FondoNu (fluido P(Y), estrella alrededor de la garganta) visto como materia elastica con c_T = 0:
    sirve para validar la forma lagrangiana frente a la criba euleriana de la Fase 2h (con el muelle de superficie)."""

    def __init__(self, F):
        self.F = F
        self.la, self.lb = 0.0, F.ls

    def perfil(self, l):
        F = self.F
        l = np.atleast_1d(np.asarray(l, dtype=float))
        f = F.campos(l)
        Y = F.sol.sol(l)
        der = np.array([F._rhs(li, yi) for li, yi in zip(l, Y.T)])        # (N', N'', r', r'') sin normalizar
        Nu = Y[0]
        Phi1 = der[:, 0] / Nu
        Phi2 = der[:, 1] / Nu - Phi1 ** 2
        G4, H, Gc, Gcp = F.escalar(l)
        w, cs2, p = f["rp"], f["cs2"], f["p"]
        rho = w - p
        mod = modulos(rho, p, p, cs2, 0.0)
        out = dict(N=f["N"], Phi1=Phi1, Phi2=Phi2, r=f["r"], r1=der[:, 2], r2=der[:, 3], Gc=Gc, Gcp=Gcp, G4=G4,
                   sig=f["sig"], w=w, wp=-w * Phi1 * (1 + 1 / np.maximum(cs2, 1e-300)))
        out.update(mod)
        return out

    def E_out(self):
        return self.F.E_out()

    def integrales_hueco(self):
        return 0.0, 0.0, 0.0
