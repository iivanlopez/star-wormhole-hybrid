# -*- coding: utf-8 -*-
"""
Fase 2h: reconstruccion de la teoria (GLPV cuartica con G3 = 0, G4 = G4(X)) sobre el fondo con fluido y zona lejana
canonica (fondo.Fondo), y coeficientes del monopolo con fluido (dos campos).

Chorros del fondo en puntos l:
  l <= lm: lapso analitico, r de T1 con materia (iteracion de la integral desde los valores del ODE);
  l >  lm: N y r del sistema lejano N'' = F_N, r'' = r' N'/N - nu (iteracion de las dos integrales), con nu como
           chorro: nu = (1 - chi) nu_int + chi (kappa/2) r phi'^2 + beta_nu 4 chi (1 - chi)/l^3, nu_int de la solucion
           interior prolongada (lapso analitico + T1).
Teoria a lo largo de la curva X = X_b(phi) (como completo.py de la Fase 2g):
  G4(x) analitica (S(x) G4b(x)); F4 = (H - Gcal)/(4X^2); F4_X de Gamma1 = r Gcal H/(r' xi), xi = Z/N con Z de diseno;
  G2, G2_X de las ligaduras con fluido; G2_XX de k_c objetivo; G3 = 0.
Diseno de K:  Z' = N G4 + (1 - chiK) [(1 + eps)(rho~+p~) Z^2/(4 N Gcal^2) + B q] + chiK N^2 K_nat/2,
  chiK = escalon((l - lw1)/(lw2 - lw1)) (exactamente 0 antes de lw1: la K natural tiene polos en la region interior,
  donde Gamma1_nat se anula), q = (2/pi) wq/(wq^2 + l^2), B tal que Z = Z_nat lejos (F4_X -> 0).  Dentro de la estrella
  N^2 K >= (1 + eps) N (rho~+p~) xi^2/(2 Gcal^2) (ausencia de fantasmas con fluido) si chiK es despreciable alli.
Diseno de k_c:  P_t = (1 - chiP) P_nat + chiP cm2 N^2 K  (P = p_d + A1 e1 - (k_c/2) A1^2, solo da k_c),
  chiP = 1 - escalon((l - wP1)/(wP2 - wP1)),
  y dentro de la estrella k_c -> k_c - Lam N r^2 chis (rigidez de las direcciones c != 0).
"""
import os
import sys
import math
import pickle
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
F2G = os.path.join(AQUI, "..", "fase2g")
for p_ in (F2F, F2G, AQUI):
    sys.path.insert(0, p_)
import jets as J                      # noqa: E402
import simbolico_num as SN            # noqa: E402
import completo as CO                 # noqa: E402  (Fase 2g: COEF, _G4_x, _P_Y)
import fondo as FH                    # noqa: E402

NJ = 8
IMAX, JMAX = CO.IMAX, CO.JMAX
COEF = CO.COEF


def _mezcla(partes, mascaras, npts):
    """Un chorro a partir de chorros parciales en subconjuntos de puntos."""
    c = None
    for jt, m in zip(partes, mascaras):
        if jt is None:
            continue
        if c is None:
            c = np.zeros((jt.c.shape[0], npts))
        c[:, m] = jt.c
    return J.Jet(c)


def escalon_jet(t):
    """escalon C^infinito como chorro (solo para 0 < t < 1)."""
    a = J.exp(-1.0 / t)
    b = J.exp(-1.0 / (1.0 - t))
    return a / (a + b)


class Teoria:
    def __init__(self, F, eps=0.2, wq=1.0, lw1=None, lw2=None, cm2=0.8, wP1=None, wP2=None, Lam=0.0, lsc=None,
                 Lz=None):
        assert F.ok and F._lejos
        self.F = F
        self.eps, self.wq, self.cm2, self.Lam = eps, wq, cm2, Lam
        self.lw1 = lw1 if lw1 is not None else F.lm + F.Delta
        self.lw2 = lw2 if lw2 is not None else 1.6 * self.lw1
        self.wP1 = wP1 if wP1 is not None else self.lw2
        self.wP2 = wP2 if wP2 is not None else 1.6 * self.wP1
        self.lw = self.lw2
        self.lsc = lsc if lsc is not None else F.ls
        self.Lz = Lz if Lz is not None else min(1.5 * self.wP2, 0.9 * F.Lmax)
        self.X0 = -F.phi1 ** 2 / 2
        self.G4x, self.Sx = CO._G4_x(F.n, F.lam, F.b, F.lT, JMAX + IMAX + 2)
        m = F.m if F.eos != "exp" else 1.0
        self.PY = CO._P_Y(F.eos, F.A, F.Ys, m, F.ks, 4)
        self._disenar_K()

    # -------------------------------------------------------------- chorros del fondo
    def _perfil_teoria(self, lj, N):
        """G4, H, Gcal (chorros) y el fluido a partir del lapso N (chorro)."""
        F = self.F
        b2 = F.b ** 2
        s = b2 / (lj * lj + b2)
        phip = F.phi1 * s
        x = s * s
        X = self.X0 * x
        g4 = [self.G4x[k](x) * (1.0 / self.X0 ** k) for k in range(JMAX + IMAX + 1)]
        H = g4[0] - 2 * X * g4[1]
        Sx = self.Sx(x)
        Gc = Sx * (1 - x ** F.n + F.d0 * J.exp(-lj * lj / F.wd ** 2))
        Y = 1.0 / (2 * N * N)
        dentro = Y.c[0] > F.Ys
        Pf = []
        for k in range(4):
            v = np.zeros_like(Y.c)
            if np.any(dentro):
                v[:, dentro] = self.PY[k](J.Jet(Y.c[:, dentro])).c
            Pf.append(J.Jet(v))
        return dict(phip=phip, x=x, X=X, g4=g4, H=H, Gc=Gc, Pf=Pf)

    def _interior(self, l):
        """Chorros interiores (lapso analitico + T1 con materia) en puntos l."""
        F = self.F
        lj = J.var(l, NJ)
        R = J.sqrt(lj * lj + F.a ** 2)
        Phi = -F.M / R
        N = J.exp(Phi)
        p = self._perfil_teoria(lj, N)
        Iv, uv = F._y(l)
        rp = p["Pf"][1] / (N * N)
        uj = J.const(0.0, NJ, l.size) + uv
        Ij = J.const(0.0, NJ, l.size) + Iv
        for _ in range(NJ + 1):
            Ij = (N * p["g4"][0] - 0.5 * N ** 3 * uj * rp).integral(Iv)
            uj = (2 * Ij / (N ** 3 * p["Gc"])).integral(uv)
        r = N * J.sqrt(uj)
        return dict(lj=lj, N=N, r=r, I=Ij, rp=rp, p=p)

    def _lejana(self, l):
        """Chorros de la zona lejana: iteracion del sistema N'' = F_N, r'' = r' N'/N - nu."""
        F = self.F
        lj = J.var(l, NJ)
        ins = self._interior(l)                         # para nu_int (solucion interior prolongada)
        rin = ins["r"]
        Phin = J.log(ins["N"])
        nu_int = Phin.d() * rin.d() - rin.d().d()
        t = (lj - F.lm) * (1.0 / F.Delta)
        tv = t.c[0]
        ch = J.const(0.0, NJ, l.size)
        m = (tv > 0) & (tv < 1)
        if np.any(m):
            ch.c[:, m] = escalon_jet(J.Jet(t.c[:, m])).c
        ch.c[0, tv >= 1] = 1.0
        Nv, Npv, rv, rpv = F.sol_lejos.sol(l)
        N = J.const(0.0, NJ, l.size) + Nv
        N.c[1] = Npv
        r = J.const(0.0, NJ, l.size) + rv
        r.c[1] = rpv
        p = self._perfil_teoria(lj, N)
        for _ in range(NJ + 2):
            Np, rpj = N.d(), r.d()
            nu = (1 - ch) * nu_int + ch * (0.5 * F.kappa) * r * p["phip"] * p["phip"] \
                + F.beta_nu * 4 * ch * (1 - ch) / (lj * lj * lj)
            Gc, G4 = p["Gc"], p["g4"][0]
            Npp = (Gc * N * (rpj * rpj - r * nu) - N * G4 + Gc.d() * (N * r * rpj - r * r * Np)) / (r * r * Gc)
            rpp = rpj * Np / N - nu
            N = Npp.integral(Npv).integral(Nv)
            r = rpp.integral(rpv).integral(rv)
            p = self._perfil_teoria(lj, N)
        return dict(lj=lj, N=N, r=r, p=p, rp=p["Pf"][1] / (N * N))

    def fondo(self, l):
        F = self.F
        l = np.atleast_1d(np.asarray(l, dtype=float))
        mi = l <= F.lm
        mf = ~mi
        out = {}
        a = self._interior(l[mi]) if np.any(mi) else None
        b = self._lejana(l[mf]) if np.any(mf) else None
        for k in ("N", "r", "rp"):
            out[k] = _mezcla([a[k] if a else None, b[k] if b else None], [mi, mf], l.size)
        pa, pb = (a["p"] if a else None), (b["p"] if b else None)
        for k in ("phip", "x", "X", "H", "Gc"):
            out[k] = _mezcla([pa[k] if pa else None, pb[k] if pb else None], [mi, mf], l.size)
        out["g4"] = [_mezcla([pa["g4"][k] if pa else None, pb["g4"][k] if pb else None], [mi, mf], l.size)
                     for k in range(len((pa or pb)["g4"]))]
        out["Pf"] = [_mezcla([pa["Pf"][k] if pa else None, pb["Pf"][k] if pb else None], [mi, mf], l.size)
                     for k in range(4)]
        out["lj"] = J.var(l, NJ)
        return out

    # -------------------------------------------------------------- K natural y diseno de Z
    def _Gamma1_nat(self, f, f0):
        X, g4 = f["X"], f["g4"]
        return g4[0] - 4 * X * g4[1] - 4 * X * X * g4[2] - 20 * X * X * f0

    def escalares(self, l):
        f = self.fondo(l)
        X = f["X"]
        f0 = (f["H"] - f["Gc"]) / (4 * X * X)
        G1 = self._Gamma1_nat(f, f0)
        r = f["r"]
        r1 = r.d()
        xi = r * f["Gc"] * f["H"] / (r1 * G1)
        Z = f["N"] * xi
        Kn = 2 * (Z.d() - f["N"] * f["g4"][0]) / (f["N"] * f["N"])
        return dict(N=f["N"].c[0], G4=f["g4"][0].c[0], Gc=f["Gc"].c[0], rp=f["rp"].c[0], Kn=Kn.c[0], Zn=Z.c[0])

    def _q(self, l):
        return (2 / math.pi) * self.wq / (self.wq ** 2 + l * l)

    def _chiK(self, l):
        return FH.escalon((np.asarray(l, dtype=float) - self.lw1) / (self.lw2 - self.lw1))

    @staticmethod
    def _escalon_chorro(lj, a, b):
        """escalon((l - a)/(b - a)) como chorro: 0 antes de a, 1 despues de b."""
        t = (lj - a) * (1.0 / (b - a))
        tv = t.c[0]
        ch = J.const(0.0, t.c.shape[0] - 1, tv.size)
        m = (tv > 0) & (tv < 1)
        if np.any(m):
            ch.c[:, m] = escalon_jet(J.Jet(t.c[:, m])).c
        ch.c[0, tv >= 1] = 1.0
        return ch

    def _disenar_K(self):
        Lz = self.Lz
        lt = np.concatenate([np.linspace(1e-4, 50, 2501), np.geomspace(50, Lz, 2500)[1:]])
        e = {}
        for i in range(0, lt.size, 500):
            ei = self.escalares(lt[i:i + 500])
            for k, v in ei.items():
                e.setdefault(k, []).append(v)
        e = {k: np.concatenate(v) for k, v in e.items()}
        self.tabla = dict(l=lt, **e)
        sp_ = {k: CubicSpline(lt, v) for k, v in e.items()}
        self.sp_ = sp_

        def rhs(l, y, B):
            l = max(l, 1e-4)
            N, G4, Gc, rp, Kn = (float(sp_[k](l)) for k in ("N", "G4", "Gc", "rp", "Kn"))
            ch = float(self._chiK(l))
            return [N * G4 + (1 - ch) * ((1 + self.eps) * rp * y[0] ** 2 / (4 * N * Gc * Gc) + B * self._q(l))
                    + ch * N * N * Kn / 2]
        z0 = float(sp_["N"](1e-4) * sp_["G4"](1e-4) * 1e-4)

        def fin(B):
            s = solve_ivp(rhs, (1e-4, Lz), [z0], args=(B,), method="DOP853", rtol=1e-11, atol=1e-13)
            if s.status != 0:
                return np.inf
            return s.y[0, -1] - float(sp_["Zn"](Lz))
        f0 = fin(0.0)
        self.exceso_sin_B = f0
        if not np.isfinite(f0) or f0 > 0:
            self.B = None
            return
        B1 = 1.0
        while fin(B1) < 0:
            B1 *= 2
        self.B = brentq(fin, 0.0, B1, xtol=1e-13, rtol=1e-12)
        self.Zsol = solve_ivp(rhs, (1e-4, Lz), [z0], args=(self.B,), method="DOP853", rtol=1e-11, atol=1e-13,
                              dense_output=True)

    # -------------------------------------------------------------- todas las magnitudes (chorros) en puntos l
    def V(self, l):
        F = self.F
        l = np.atleast_1d(np.asarray(l, dtype=float))
        f = self.fondo(l)
        lj, N, X, g4 = f["lj"], f["N"], f["X"], f["g4"]
        r, rp = f["r"], f["rp"]
        r1 = r.d()
        Gc, H = f["Gc"], f["H"]
        f0 = (H - Gc) / (4 * X * X)
        G1n = self._Gamma1_nat(f, f0)
        Zn = N * r * Gc * H / (r1 * G1n)
        Kn = 2 * (Zn.d() - N * g4[0]) / (N * N)
        chi = self._escalon_chorro(lj, self.lw1, self.lw2)
        q = (2 / math.pi) * self.wq / (self.wq ** 2 + lj * lj)
        dentro_Z = l <= self.Lz
        Zv = np.where(dentro_Z, self.Zsol.sol(np.minimum(l, self.Lz))[0], Zn.c[0])
        Zj = J.const(0.0, NJ, l.size) + Zv
        for _ in range(NJ + 1):
            Zj = (N * g4[0] + (1 - chi) * ((1 + self.eps) * rp * Zj * Zj / (4 * N * Gc * Gc) + self.B * q)
                  + chi * N * N * Kn / 2).integral(Zv)
        if np.any(~dentro_Z):
            Zj.c[:, ~dentro_Z] = Zn.c[:, ~dentro_Z]
        xi = Zj / N
        G1 = r * Gc * H / (r1 * xi)
        f1 = (g4[0] - 4 * X * g4[1] - 4 * X * X * g4[2] - 20 * X * X * f0 - G1) / (8 * X * X * X)
        # mas alla de la ventana Z = Z_nat (salvo redondeo): F4_X = 0 exactamente para no amplificar el redondeo
        lejos = l > 1.3 * self.lw2
        if np.any(lejos):
            f1.c[:, lejos] = 0.0
        cero = J.const(0.0, NJ, l.size)
        phip = f["phip"]
        phij = phip.integral(F.phi1 * F.b * np.arctan(l / F.b))
        phipp = phip.d()

        def T(seq):
            out = []
            for k in range(len(seq)):
                sig = seq[k + 1] if k + 1 < len(seq) else cero
                out.append(seq[k].d() / phip + phipp * sig)
            return out
        Phi = J.log(N)
        Ph1 = Phi.d()
        V = {"r0": r, "r1": r1, "r2": r1.d(), "r3": r1.d().d(), "r4": r1.d().d().d(),
             "P0": phij, "P1": phip, "P2": phipp, "P3": phipp.d(), "P4": phipp.d().d(),
             "Phi1": Ph1, "Phi2": Ph1.d(), "Phi3": Ph1.d().d(), "Phi4": Ph1.d().d().d(),
             "E": N * N, "eP": N, "Pf0": f["Pf"][0], "Pf1": f["Pf"][1]}
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
        gam, al1 = COEF["gam"](V), COEF["al1"](V)
        kc0 = COEF["k_c"](V)
        dk = COEF["dkc_dG202"](V)
        e1 = COEF["e1"](V)
        p_d = -2 * N * Gc
        A1 = al1 / gam
        K = 2 * (Zj.d() - N * g4[0]) / (N * N)
        chiP = 1 - self._escalon_chorro(lj, self.wP1, self.wP2)
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
        V.update(N=N, Gc=Gc, G4=g4[0], H=H, rp=rp, K=K, xi=xi, G1=G1, f1=f1, l=lj, Pf2=f["Pf"][2], kct=kct, Pt=Pt,
                 Kn=Kn, Zj=Zj)
        return V

    def coeficientes(self, l):
        V = self.V(l)
        N, r, r1, Gc, G4 = V["N"], V["r0"], V["r1"], V["Gc"], V["G4"]
        gam, al1 = COEF["gam"](V), COEF["al1"](V)
        kc, e0, e1, ett = COEF["k_c"](V), COEF["e0"](V), COEF["e1"](V), COEF["e_tt"](V)
        Pf1, Pf2 = V["Pf1"], V["Pf2"]
        E = N * N
        a2 = 4 * N * N * r * Gc
        out = dict(N=N, Phi1=V["Phi1"], r=r, r1=r1, Gc=Gc, G4=G4, gam=gam, al1=al1, kc=kc, e0=e0, e1=e1, ett=ett,
                   a2=a2, a2p=a2.d(), Ng=N * N * gam, Ngp=(N * N * gam).d(), rp=V["rp"], K=V["K"], kct=V["kct"])
        out = {k: (v.c[0] if isinstance(v, J.Jet) else v) for k, v in out.items()}
        out["cs2"] = np.where(Pf1.c[0] > 0, Pf1.c[0] / (Pf1.c[0] + Pf2.c[0] / E.c[0] + 1e-300), 0.0)
        out["sig"] = out["r1"] / out["r"] - out["Phi1"]
        out["Np"] = out["N"] * out["Phi1"]
        return out, V
