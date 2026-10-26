# -*- coding: utf-8 -*-
"""
BORRADOR (no validado; para la fase siguiente): solo se ha ejecutado la parte del diseno de K, que muestra que la cola
natural lejana de K tambien consume presupuesto (fase2g_resultados.md, seccion 8).

Fase 2g, parte III: analisis completo del monopolo (dos campos: fluido + escalar/metrica) y del sector impar en una
garganta de GLPV cuartica llena de fluido, con la teoria reconstruida por el metodo inverso.

Fondo (criba.FondoFluido): lapso, perfil del escalar, acoplo G4 = S(x) G4b(x) reducido en el nucleo, Gcal objetivo,
fluido P(Y) con superficie suave (eos = "exp": P = A z e^{-ks/z}, z = Y - Ys) y r de T1 con materia.
Teoria (funciones de (phi, X) a lo largo de la curva X = X_b(phi), con G3 = 0 y G4 = G4(X)):
    G4(x) analitica;  F4 = (H - Gcal)/(4 X^2);  F4_X de Gamma1 objetivo;  F4_XX = 0;
    G2, G2_X de las ligaduras con fluido;  G2_XX de k_c objetivo;  G2_XXX = 0.
Diseno de K (ausencia de fantasmas con fluido, criba.presupuesto_fantasma): Z = N xi con
    Z' = N G4 + (1 - chi)[(1 + eps) (rho~+p~) Z^2/(4 N Gcal^2) + B q] + chi N^2 K_nat/2,
chi = 1 - e^{-(l/lw)^8} (ventana), q = bump normalizado, y B tal que Z = Z_nat fuera de la ventana (ahi Gamma1 es el
natural y F4_X = 0).  Gamma1 = r Gcal H/(r' xi).
Diseno de k_c: fuera de la estrella, de P = cm2 N^2 K (como en la Fase 2f); dentro, k_c -> k_c - Lam N r^2 chi_s
(chi_s ~ 1 en la estrella): energia grande en las direcciones c != 0.

Monopolo, variables h (d = N h) y c (= delta g_ll/2) en la region con fluido [0, lc], y solo h fuera (W = 0, c = chi1 h'):
    E_in  = Q^2/(4 kap_a) - 2 N^3 Gcal (h' - sig h)^2 - k_c c^2/2 - c (e0 N h + e1 (N h)')
    T_in  = (al1 N h' + gam c)^2/(4 kap_b) - 2 G4 N h^2 - e_tt N c h
    Q     = -(4 N^2 r Gcal h')' - 2 N^2 r (rho~+p~) h + (N^2 gam c)' - N (rho~+p~) r^2 c
    E_out = -2 N^3 Gcal (h' - sig h)^2 - k_c chi1^2 h'^2/2 - chi1 h' (e0 N h + e1 (N h)'),   chi1 = 4 r Gcal/gam
    T_out = -2 G4 N h^2 - e_tt N chi1 h' h
Elementos finitos de Hermite cubicos para h y c; problema generalizado E v = w T v (T debe ser definida positiva).
"""
import os
import sys
import math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import eigh, cholesky

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
sys.path.insert(0, AQUI)
import jets as J                      # noqa: E402
import simbolico_num as SN            # noqa: E402
import criba as C                     # noqa: E402
import pickle                         # noqa: E402

NJ = 8
IMAX, JMAX = 3, 4
GX, GW = np.polynomial.legendre.leggauss(8)
_cf = pickle.load(open(os.path.join(AQUI, "resultados", "coeficientes_fluido.pkl"), "rb"))
COEF = {k: SN.compilar(SN.con_Xb(_cf[k])) for k in ("k_c", "e0", "e1", "e_tt", "G2_00", "G2_01")}
COEF["dkc_dG202"] = SN.compilar(sp.diff(SN.con_Xb(_cf["k_c"]), sp.Symbol("G2_02", real=True)))
for k in ("gam", "al1", "al0"):
    COEF[k] = SN.MONO[k]


def _G4_x(n, lam, b, lT, kmax):
    """Derivadas en x de G4(x) = S(x) (1 + x^n/(2n-1)), S = lam + (1 - lam)(1 - e^{-q^2}), q = b^2 (x^{-1/2} - 1)/lT^2."""
    x = sp.Symbol("x", positive=True)
    q = b ** 2 * (x ** sp.Rational(-1, 2) - 1) / lT ** 2
    Sx = lam + (1 - lam) * (1 - sp.exp(-q ** 2))
    G4 = Sx * (1 + x ** n / (2 * n - 1))
    out = []
    e = G4
    for k in range(kmax + 1):
        out.append(sp.lambdify(x, e, modules=J.MODULOS + ["numpy"]))
        e = sp.diff(e, x)
    fS = sp.lambdify(x, Sx, modules=J.MODULOS + ["numpy"])
    return out, fS


def _P_Y(eos, A, Ys, m, ks, kmax):
    Y = sp.Symbol("Y", positive=True)
    z = Y - Ys
    P = A * z * sp.exp(-ks / z) if eos == "exp" else A * z ** m
    out = []
    e = P
    for k in range(kmax + 1):
        out.append(sp.lambdify(Y, e, modules=J.MODULOS + ["numpy"]))
        e = sp.diff(e, Y)
    return out


class GargantaFluido:
    def __init__(self, fondo, phi1=1.0, eps=0.2, wq=1.0, lw_fac=2.5, cm2=0.8, wP=None, Lam=0.0, ls_fac=1.0,
                 rp_corte=1e-6):
        self.F = F = fondo
        assert F.ok
        self.phi1 = phi1
        self.X0 = -phi1 ** 2 / 2
        self.eps, self.wq, self.cm2, self.Lam = eps, wq, cm2, Lam
        self.lw = lw_fac * F.ls
        self.wP = wP if wP is not None else 1.5 * self.lw
        self.lsc = ls_fac * F.ls
        self.G4x, self.Sx = _G4_x(F.n, F.lam, F.b, F.lT, JMAX + IMAX + 2)
        m = F.m if F.eos != "exp" else 1.0
        self.PY = _P_Y(F.eos, F.A, F.Ys, m, F.ks, 4)
        # region con fluido apreciable: (rho~+p~) > rp_corte (rho~+p~)(0)
        lt = np.linspace(0, F.ls, 20001)
        rp = F.fluido(F.lapso(lt)[0])[0]
        dentro = np.nonzero(rp > rp_corte * rp[0])[0]
        self.lc = float(lt[dentro[-1]])
        self._disenar_K()

    # -------------------------------------------------------------- chorros de los perfiles en puntos l
    def perfiles(self, l):
        F = self.F
        l = np.atleast_1d(np.asarray(l, dtype=float))
        lj = J.var(l, NJ)
        R = J.sqrt(lj * lj + F.a ** 2)
        Phi = -F.M / R
        N = J.exp(Phi)
        b2 = F.b ** 2
        s = b2 / (lj * lj + b2)
        phip = self.phi1 * s
        x = s * s
        X = self.X0 * x
        g4 = [self.G4x[k](x) * (1.0 / self.X0 ** k) for k in range(JMAX + IMAX + 1)]
        H = g4[0] - 2 * X * g4[1]
        Sx = self.Sx(x)
        Gc = Sx * (1 - x ** F.n + F.d0 * J.exp(-lj * lj / F.wd ** 2))
        # fluido: P(Y) y derivadas en Y_b = 1/(2 N^2), cero fuera de la estrella
        Y = 1.0 / (2 * N * N)
        dentro = Y.c[0] > F.Ys
        Pf = []
        for k in range(4):
            v = np.zeros_like(Y.c)
            if np.any(dentro):
                Yd = J.Jet(Y.c[:, dentro])
                v[:, dentro] = self.PY[k](Yd).c
            Pf.append(J.Jet(v))
        return dict(l=lj, N=N, Phi=Phi, phip=phip, x=x, X=X, g4=g4, H=H, Gc=Gc, Pf=Pf, Y=Y)

    def _fondo_r(self, l, p):
        """r, u, I, J como chorros (iteracion de la integral desde los valores del ODE de criba.FondoFluido)."""
        F = self.F
        Iv, uv = F._y(l)
        rp = p["Pf"][1] / (p["N"] * p["N"])
        uj = J.const(0.0, NJ, l.size) + uv
        Ij = J.const(0.0, NJ, l.size) + Iv
        for _ in range(NJ + 1):
            Ij = (p["N"] * p["g4"][0] - 0.5 * p["N"] ** 3 * uj * rp).integral(Iv)
            uj = (2 * Ij / (p["N"] ** 3 * p["Gc"])).integral(uv)
        r = p["N"] * J.sqrt(uj)
        return dict(r=r, u=uj, I=Ij, rp=rp)

    # -------------------------------------------------------------- K de diseno (Z = N xi)
    def _Gamma1_nat(self, p, f4_0):
        X, g4 = p["X"], p["g4"]
        return g4[0] - 4 * X * g4[1] - 4 * X * X * g4[2] - 20 * X * X * f4_0

    def _escalares(self, l):
        """Valores (sin chorros) de N, G4, Gcal, rho~+p~ y K_nat, Z_nat en puntos l (para el ODE de Z)."""
        p = self.perfiles(l)
        fr = self._fondo_r(np.atleast_1d(l), p)
        X = p["X"]
        f0 = (p["H"] - p["Gc"]) / (4 * X * X)
        G1 = self._Gamma1_nat(p, f0)
        r = fr["r"]
        r1 = r.d()
        xi = r * p["Gc"] * p["H"] / (r1 * G1)
        Z = p["N"] * xi
        Kn = 2 * (Z.d() - p["N"] * p["g4"][0]) / (p["N"] * p["N"])
        return dict(N=p["N"].c[0], G4=p["g4"][0].c[0], Gc=p["Gc"].c[0], rp=fr["rp"].c[0], Kn=Kn.c[0], Zn=Z.c[0])

    def _q(self, l):
        return (2 / math.pi) * self.wq / (self.wq ** 2 + l * l)

    def _chi(self, l):
        return 1 - np.exp(-(l / self.lw) ** 8)

    def _disenar_K(self):
        F = self.F
        Lz = 4 * self.lw
        lt = np.concatenate([np.linspace(1e-4, Lz, 6001)])
        e = self._escalares(lt)
        from scipy.interpolate import CubicSpline
        sp_ = {k: CubicSpline(lt, v) for k, v in e.items()}

        def rhs(l, y, B):
            l = max(l, 1e-4)
            N, G4, Gc, rp, Kn = (float(sp_[k](l)) for k in ("N", "G4", "Gc", "rp", "Kn"))
            ch = float(self._chi(l))
            return [N * G4 + (1 - ch) * ((1 + self.eps) * rp * y[0] ** 2 / (4 * N * Gc * Gc) + B * self._q(l))
                    + ch * N * N * Kn / 2]

        def fin(B):
            s = solve_ivp(rhs, (1e-4, Lz), [float(sp_["N"](1e-4) * sp_["G4"](1e-4) * 1e-4)], args=(B,),
                          method="DOP853", rtol=1e-11, atol=1e-13)
            return s.y[0, -1] - float(sp_["Zn"](Lz))
        f0 = fin(0.0)
        self.Z_sin_B = f0
        if f0 > 0:
            self.B = None                                   # ni con B = 0 se alcanza Z_nat: fantasma
            return
        B1 = 1.0
        while fin(B1) < 0:
            B1 *= 2
        from scipy.optimize import brentq
        self.B = brentq(fin, 0.0, B1, xtol=1e-12)
        s = solve_ivp(rhs, (1e-4, Lz), [float(sp_["N"](1e-4) * sp_["G4"](1e-4) * 1e-4)], args=(self.B,),
                      method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True)
        self.Zsol = s
        self.Lz = Lz
        self.sp_ = sp_

    # -------------------------------------------------------------- todas las magnitudes en puntos l
    def V(self, l):
        F = self.F
        l = np.atleast_1d(np.asarray(l, dtype=float))
        p = self.perfiles(l)
        fr = self._fondo_r(l, p)
        lj, N, X, g4 = p["l"], p["N"], p["X"], p["g4"]
        r, rp = fr["r"], fr["rp"]
        f0 = (p["H"] - p["Gc"]) / (4 * X * X)
        # Z de diseno como chorro: iteracion de la integral del ODE desde los valores de Zsol
        G1n = self._Gamma1_nat(p, f0)
        r1 = r.d()
        Zn = N * r * p["Gc"] * p["H"] / (r1 * G1n)
        Kn = 2 * (Zn.d() - N * g4[0]) / (N * N)
        chi = 1 - J.exp(-(lj * (1.0 / self.lw)) ** 8)
        q = (2 / math.pi) * self.wq / (self.wq ** 2 + lj * lj)
        Zv = np.where(l <= self.Lz, self.Zsol.sol(np.minimum(l, self.Lz))[0], Zn.c[0])
        Zj = J.const(0.0, NJ, l.size) + Zv
        for _ in range(NJ + 1):
            Zj = (N * g4[0] + (1 - chi) * ((1 + self.eps) * rp * Zj * Zj / (4 * N * p["Gc"] * p["Gc"]) + self.B * q)
                  + chi * N * N * Kn / 2).integral(Zv)
        lejos = l > self.Lz
        if np.any(lejos):
            Zj.c[:, lejos] = Zn.c[:, lejos]
        xi = Zj / N
        G1 = r * p["Gc"] * p["H"] / (r1 * xi)
        f1 = (g4[0] - 4 * X * g4[1] - 4 * X * X * g4[2] - 20 * X * X * f0 - G1) / (8 * X * X * X)
        cero = J.const(0.0, NJ, l.size)
        phip = p["phip"]
        phij = phip.integral(self.phi1 * F.b * np.arctan(l / F.b))
        phipp = phip.d()

        def T(seq):
            out = []
            for k in range(len(seq)):
                sig = seq[k + 1] if k + 1 < len(seq) else cero
                out.append(seq[k].d() / phip + phipp * sig)
            return out
        Ph1 = p["Phi"].d()
        V = {"r0": r, "r1": r1, "r2": r1.d(), "r3": r1.d().d(), "r4": r1.d().d().d(),
             "P0": phij, "P1": phip, "P2": phipp, "P3": phipp.d(), "P4": phipp.d().d(),
             "Phi1": Ph1, "Phi2": Ph1.d(), "Phi3": Ph1.d().d(), "Phi4": Ph1.d().d().d(),
             "E": N * N, "eP": N, "Pf0": p["Pf"][0], "Pf1": p["Pf"][1]}
        f4 = [f0, f1] + [cero] * (JMAX + IMAX)
        g3 = [cero] * (JMAX + IMAX + 1)
        for nombre, seq in (("G4", g4), ("F4", f4), ("G3", g3)):
            s_ = list(seq)
            for i in range(IMAX + 1):
                for j in range(JMAX + 1):
                    V[f"{nombre}_{i}{j}"] = s_[j] if j < len(s_) else cero
                s_ = T(s_)
        g20 = COEF["G2_00"](V)
        V["G2_00"] = g20
        V["G2_01"] = COEF["G2_01"](V)
        V["G2_02"] = cero
        # k_c objetivo
        gam, al1 = COEF["gam"](V), COEF["al1"](V)
        kc0 = COEF["k_c"](V)
        dk = COEF["dkc_dG202"](V)
        e1 = COEF["e1"](V)
        Gc = p["Gc"]
        p_d = -2 * N * Gc
        A1 = al1 / gam
        K = 2 * (Zj.d() - N * g4[0]) / (N * N)
        chiP = J.exp(-(lj * (1.0 / self.wP)) ** 4)
        # P_t = (1 - chiP) P_nat + chiP cm2 N^2 K;  P = p_d + A1 e1 - (k_c/2) A1^2  ->  k_c
        P_nat = p_d + A1 * e1 - (kc0 / 2) * A1 * A1
        Pt = (1 - chiP) * P_nat + chiP * self.cm2 * N * N * K
        kcP = 2 * (p_d + A1 * e1 - Pt) / (A1 * A1)
        chis = J.exp(-(lj * (1.0 / self.lsc)) ** 8)
        kct = kcP - self.Lam * N * r * r * chis
        g22 = (kct - kc0) / dk
        V["G2_02"] = g22
        g2 = [g20, V["G2_01"], g22] + [cero] * (JMAX + IMAX - 2)
        s_ = g2
        for i in range(IMAX + 1):
            for j in range(JMAX + 1):
                V[f"G2_{i}{j}"] = s_[j] if j < len(s_) else cero
            s_ = T(s_)
        V.update(N=N, Gc=Gc, G4=g4[0], H=p["H"], rp=rp, K=K, xi=xi, G1=G1, f1=f1, l=lj,
                 Pf2=p["Pf"][2], kct=kct, P=Pt)
        return V

    def coeficientes(self, l):
        """Coeficientes de las formas E y T (valores en l) y comprobaciones."""
        V = self.V(l)
        N, r, r1, Gc, G4 = V["N"], V["r0"], V["r1"], V["Gc"], V["G4"]
        gam, al1 = COEF["gam"](V), COEF["al1"](V)
        kc, e0, e1, ett = COEF["k_c"](V), COEF["e0"](V), COEF["e1"](V), COEF["e_tt"](V)
        rp = V["rp"]
        Pf1, Pf2 = V["Pf1"], V["Pf2"]
        E = N * N
        a2 = 4 * N * N * r * Gc
        out = dict(N=N, Phi1=V["Phi1"], r=r, r1=r1, Gc=Gc, G4=G4, gam=gam, al1=al1, kc=kc, e0=e0, e1=e1, ett=ett,
                   a2=a2, a2p=a2.d(), Ng=N * N * gam, Ngp=(N * N * gam).d(), rp=rp, K=V["K"], kct=V["kct"])
        out = {k: (v.c[0] if isinstance(v, J.Jet) else v) for k, v in out.items()}
        cs2 = np.where(Pf1.c[0] > 0, Pf1.c[0] / (Pf1.c[0] + Pf2.c[0] / E.c[0] + 1e-300), 0.0)
        out["cs2"] = cs2
        out["sig"] = out["r1"] / out["r"] - out["Phi1"]
        out["Np"] = out["N"] * out["Phi1"]
        out["N_"] = out["N"]
        return out


# ------------------------------------------------------------------ elementos finitos
def _herm(t, H):
    return C.hermite(t, H)


def mallas(g, n_in=120, L=600.0, n_out=160):
    a = np.linspace(0.0, g.lc, n_in + 1)
    b = g.lc + (L - g.lc) * (np.geomspace(1.0, 1.0 + (L - g.lc), n_out + 1) - 1.0) / (L - g.lc)
    return a, b


def ensamblar(g, paridad, n_in=120, L=600.0, n_out=160):
    """Matrices E, T (grados de libertad: h, h' en todos los nodos; c, c' en los nodos interiores [0, lc])."""
    a, b = mallas(g, n_in, L, n_out)
    nodos = np.concatenate([a, b[1:]])
    Kn = len(nodos)
    Ki = len(a)
    nh = 2 * Kn
    nd = nh + 2 * Ki
    E = np.zeros((nd, nd))
    T = np.zeros((nd, nd))
    t = 0.5 * (GX + 1)
    for e in range(Kn - 1):
        x0, x1 = nodos[e], nodos[e + 1]
        H = x1 - x0
        w = 0.5 * GW * H
        l = x0 + t * H
        f = g.coeficientes(l)
        v, d1, d2 = _herm(t, H)
        ih = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        N, r, Gc, sig = f["N"], f["r"], f["Gc"], f["sig"]
        Dm = d1 - sig * v                                          # h' - sig h
        Nh_p = f["Np"] * v + N * d1                                 # (N h)'
        if e < Ki - 1:                                              # interior: h y c
            ic = [nh + 2 * e, nh + 2 * e + 1, nh + 2 * e + 2, nh + 2 * e + 3]
            idx = ih + ic
            rp = np.maximum(f["rp"], 1e-300)
            AQ = f["cs2"] / (2 * N * r * r * rp)
            AB = N / (2 * r * r * rp)
            # Q = Qh . h + Qc . c
            Qh = -f["a2"] * d2 - f["a2p"] * d1 - 2 * N * N * r * f["rp"] * v
            Qc = f["Ng"] * d1 + f["Ngp"] * v - N * f["rp"] * r * r * v
            Qv = np.vstack([Qh, Qc])
            Dv = np.vstack([Dm, 0 * Dm])
            cv = np.vstack([0 * v, v])
            Eloc = (Qv * AQ * w) @ Qv.T - (Dv * 2 * N ** 3 * Gc * w) @ Dv.T - (cv * f["kc"] / 2 * w) @ cv.T
            acop = np.vstack([f["e0"] * N * v + f["e1"] * Nh_p, 0 * v])          # (e0 N h + e1 (N h)') como vector en h
            Eloc -= 0.5 * ((cv * w) @ acop.T + (acop * w) @ cv.T)
            Wv = np.vstack([f["al1"] * N * d1, f["gam"] * v])                   # W = al1 N h' + gam c
            hv = np.vstack([v, 0 * v])
            Tloc = (Wv * AB * w) @ Wv.T - (hv * 2 * f["G4"] * N * w) @ hv.T
            cr = (cv * (-f["ett"] * N / 2) * w) @ hv.T
            Tloc += cr + cr.T
            E[np.ix_(idx, idx)] += Eloc
            T[np.ix_(idx, idx)] += Tloc
        else:                                                       # exterior: c = chi1 h'
            chi1 = 4 * r * Gc / f["gam"]
            cvec = chi1 * d1
            Eloc = -(Dm * 2 * N ** 3 * Gc * w) @ Dm.T - (cvec * f["kc"] / 2 * w) @ cvec.T
            acop = f["e0"] * N * v + f["e1"] * Nh_p
            Eloc -= 0.5 * ((cvec * w) @ acop.T + (acop * w) @ cvec.T)
            Tloc = -(v * 2 * f["G4"] * N * w) @ v.T
            cr = (cvec * (-f["ett"] * N / 2) * w) @ v.T
            Tloc += cr + cr.T
            E[np.ix_(ih, ih)] += Eloc
            T[np.ix_(ih, ih)] += Tloc
    # condiciones: paridad en la garganta, c(lc) = chi1(lc) h'(lc), h(L) = 0
    fl = g.coeficientes(np.array([g.lc]))
    chi_c = float(4 * fl["r"][0] * fl["Gc"][0] / fl["gam"][0])
    # transformacion: los grados libres son todos salvo c(lc) (sustituido por chi_c h'(lc))
    jc = nh + 2 * (Ki - 1)
    jh1 = 2 * (Ki - 1) + 1
    libres = [i for i in range(nd) if i != jc]
    Tm = np.zeros((nd, len(libres)))
    for k, i in enumerate(libres):
        Tm[i, k] = 1.0
    Tm[jc, libres.index(jh1)] = chi_c
    E2 = Tm.T @ E @ Tm
    T2 = Tm.T @ T @ Tm
    quitar = []
    if paridad == "par":
        quitar += [1, nh + 1]                      # h'(0) = 0, c'(0) = 0
    else:
        quitar += [0, nh]                          # h(0) = 0, c(0) = 0
    quitar.append(2 * (Kn - 1))                    # h(L) = 0
    keep = [k for k, i in enumerate(libres) if i not in quitar]
    return E2[np.ix_(keep, keep)], T2[np.ix_(keep, keep)], nodos


def espectro(g, paridad, n_in=120, L=600.0, n_out=160, k=6):
    E, T, nodos = ensamblar(g, paridad, n_in, L, n_out)
    E = 0.5 * (E + E.T)
    T = 0.5 * (T + T.T)
    tv = np.linalg.eigvalsh(T)
    out = dict(Tmin=float(tv[0]), Tmax=float(tv[-1]))
    if tv[0] > 0:
        w = eigh(E, T, eigvals_only=True, subset_by_index=[0, k - 1])
        out["w"] = [float(x) for x in w]
    return out
