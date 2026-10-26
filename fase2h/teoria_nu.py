# -*- coding: utf-8 -*-
"""
Fase 2h: reconstruccion de la teoria sobre los fondos fondo_nu.FondoNu y diseno de K con exterior natural.

Chorros del fondo: N y r salen del sistema de T1 con materia (el mismo de FondoNu, con el lapso normalizado y el fluido
reescrito con fluido_normalizado), iterado sobre los valores del ODE en cada punto.

Diseno de K (ver diseno_ext.py y fase2h_resultados.md).  Con Z = N xi, D = Z - Z_nat (Z_nat: F4_X = 0) y
sig = (ln r/N)', en la region sin fluido
    U = U_nat(P) + (2/r^2) (2 sig D - D'),        K = K_nat + 2 D'/N^2,
asi que fuera de la esfera de fotones (sig > 0) una D >= 0 decreciente sube U, a costa de una parte del presupuesto
natural de K.  El diseno:
    Z' = N G4 + (1 - chi) [(1 + eps)(rho~+p~) Z^2/(4 N Gcal^2) + B q] + chi (1 - theta w) N^2 K_nat/2,
    chi = escalon((l - l1)/(l2 - l1))  (de la estrella al exterior natural),
    w = 1 - escalon((l - LDa)/(LD - LDa)),  q gaussiana normalizada centrada en lq de anchura wq;
tras l2, D' = -theta w N^2 K_nat/2 y B se fija con D(LD) = 0, es decir D(l) = theta int_l^LD w N^2 K_nat/2 >= 0;
mas alla de LD, Z = Z_nat y F4_X = 0 exactamente.  Dentro de la estrella N^2 K >= (1 + eps) N (rho~+p~) xi^2/(2 Gcal^2).
P (y k_c fuera del fluido): P_t = (1 - chiP) P_nat + chiP cm2 N^2 K, chiP = 1 - escalon((l - wP1)/(wP2 - wP1)).
"""
import os
import sys
import math
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2f"), os.path.join(AQUI, "..", "fase2g"), AQUI):
    sys.path.insert(0, p_)
import jets as J                      # noqa: E402
import completo as CO                 # noqa: E402
import teoria as TE                   # noqa: E402
import fondo as FH                    # noqa: E402

NJ = TE.NJ
IMAX, JMAX = TE.IMAX, TE.JMAX
COEF = TE.COEF


class TeoriaNu(TE.Teoria):
    def __init__(self, F, eps=0.05, theta=0.0, l1=None, l2=None, LDa=None, LD=None, lq=None, wq=None, cm2=1.0,
                 wP1=None, wP2=None, Lam=0.0, lsc=None, Lz=None, b_in=0.0, w_in=0.5):
        """b_in: semilla fija b_in q_in de K en la garganta (q_in gaussiana normalizada de anchura w_in), para que K(0) > 0
        (sin ella K ~ l^2 en la garganta y la matriz cinetica degenera alli)."""
        assert F.ok and F.ls > 0
        self.F = F
        self.eps, self.theta, self.cm2, self.Lam = eps, theta, cm2, Lam
        self.b_in, self.w_in = b_in, w_in
        ls = F.ls
        self.l1 = l1 if l1 is not None else ls
        self.l2 = l2 if l2 is not None else 1.15 * ls
        self.LDa = LDa if LDa is not None else 4.0 * ls
        self.LD = LD if LD is not None else 8.0 * ls
        self.lq = lq if lq is not None else 0.85 * ls
        self.wq = wq if wq is not None else 0.15 * ls
        self.wP1 = wP1 if wP1 is not None else self.LD
        self.wP2 = wP2 if wP2 is not None else 1.6 * self.wP1
        self.lsc = lsc if lsc is not None else ls
        self.Lz = Lz if Lz is not None else min(2.0 * self.wP2, 0.9 * F.Lmax)
        self.lw2 = self.LD                       # (Teoria.V usa lw2 para el corte "lejos" de F4_X)
        self.X0 = -F.phi1 ** 2 / 2
        self.G4x, self.Sx = CO._G4_x(F.n, F.lam, F.b, F.lT, JMAX + IMAX + 2)
        fl = F.fluido_normalizado()
        self.Ys_n = fl["Ys"]
        self.PY = CO._P_Y(F.eos, fl["A"], fl["Ys"], fl["m"], fl["ks"], 4)
        self._disenar()

    # -------------------------------------------------------------- chorros del fondo
    def _perfil_teoria(self, lj, N):
        F = self.F
        b2 = F.b ** 2
        s = b2 / (lj * lj + b2)
        phip = F.phi1 * s
        x = s * s
        X = self.X0 * x
        g4 = []
        for k in range(JMAX + IMAX + 1):
            res = self.G4x[k](x)
            if isinstance(res, J.Jet):
                g4.append(res * (1.0 / self.X0 ** k))
            else:                                        # derivada constante (lam = 1)
                g4.append(J.const(float(res) / self.X0 ** k, NJ, lj.c.shape[1]))
        H = g4[0] - 2 * X * g4[1]
        Gc = self.Sx(x) * (1 - x ** F.n + F.d0 * J.exp(-lj * lj / F.wd ** 2))
        Y = 1.0 / (2 * N * N)
        dentro = Y.c[0] > self.Ys_n
        Pf = []
        for k in range(4):
            v = np.zeros_like(Y.c)
            if np.any(dentro):
                res = self.PY[k](J.Jet(Y.c[:, dentro]))
                if isinstance(res, J.Jet):
                    v[:, dentro] = res.c
                else:                                    # derivada constante (m entero)
                    v[0, dentro] = float(res)
            Pf.append(J.Jet(v))
        return dict(phip=phip, x=x, X=X, g4=g4, H=H, Gc=Gc, Pf=Pf)

    def fondo(self, l):
        F = self.F
        l = np.atleast_1d(np.asarray(l, dtype=float))
        lj = J.var(l, NJ)
        Nu, Npu, rv, rpv = F.sol.sol(l)
        Nv, Npv = Nu / F._Ninf, Npu / F._Ninf
        N = J.const(0.0, NJ, l.size) + Nv
        N.c[1] = Npv
        r = J.const(0.0, NJ, l.size) + rv
        r.c[1] = rpv
        p = self._perfil_teoria(lj, N)
        joroba = J.exp(-(lj * (1.0 / F.wt)) * (lj * (1.0 / F.wt))) * (-F.nu0_actual)
        for _ in range(NJ + 2):
            Np, rpj = N.d(), r.d()
            rp = p["Pf"][1] / (N * N)
            Gc, G4 = p["Gc"], p["g4"][0]
            nu = joroba + F.f_fl * r * rp / (4 * G4) + (0.5 * F.kappa) * r * p["phip"] * p["phip"]
            Npp = (Gc * N * (rpj * rpj - r * nu) - N * G4 + 0.5 * N * r * r * rp + Gc.d() * (N * r * rpj - r * r * Np)) \
                / (r * r * Gc)
            rpp = rpj * Np / N - nu
            N = Npp.integral(Npv).integral(Nv)
            r = rpp.integral(rpv).integral(rv)
            p = self._perfil_teoria(lj, N)
        out = dict(N=N, r=r, rp=p["Pf"][1] / (N * N), lj=lj)
        out.update({k: p[k] for k in ("phip", "x", "X", "H", "Gc", "g4", "Pf")})
        return out

    # -------------------------------------------------------------- diseno de Z
    def _chi(self, l):
        return FH.escalon((np.asarray(l, dtype=float) - self.l1) / (self.l2 - self.l1))

    def _w(self, l):
        return 1.0 - FH.escalon((np.asarray(l, dtype=float) - self.LDa) / (self.LD - self.LDa))

    def _q(self, l):
        nrm = self.wq * math.sqrt(math.pi) / 2 * (1 + math.erf(self.lq / self.wq))
        return np.exp(-((np.asarray(l, dtype=float) - self.lq) / self.wq) ** 2) / nrm

    def _qin(self, l):
        return 2 / (math.sqrt(math.pi) * self.w_in) * np.exp(-(np.asarray(l, dtype=float) / self.w_in) ** 2)

    def _tabla(self, lt, trozo=500):
        e = {}
        for i in range(0, lt.size, trozo):
            ei = self.escalares(lt[i:i + trozo])
            for k, v in ei.items():
                e.setdefault(k, []).append(v)
        return {k: np.concatenate(v) for k, v in e.items()}

    def _disenar(self):
        F = self.F
        # interior (sin K natural, que tiene polos donde Gamma1_nat = 0): N, G4, Gcal; el fluido se evalua exacto
        li = np.linspace(1e-4, self.l2, 3001)
        fi = F.campos(li)
        self.sp_in = {k: CubicSpline(li, fi[k]) for k in ("N", "G4", "Gc")}
        # exterior: K y Z naturales
        a = 0.9 * self.l1
        le = np.concatenate([np.linspace(a, 4 * self.l2, 1500), np.geomspace(4 * self.l2, self.Lz, 1500)[1:]])
        te = self._tabla(le)
        self.tabla = dict(l=le, **te)
        self.sp_ex = {k: CubicSpline(le, te[k]) for k in ("N", "G4", "Kn", "Zn")}
        # D(l) = theta int_l^LD w N^2 K_nat/2  (para l >= l2)
        g = self.theta * self._w(le) * te["N"] ** 2 * te["Kn"] / 2
        cum = cumulative_trapezoid(g, le, initial=0.0)
        iLD = np.searchsorted(le, self.LD)
        Dt = np.where(le <= self.LD, np.interp(self.LD, le, cum) - cum, 0.0)
        self.sp_D = CubicSpline(le, Dt)
        self.Kn_min_ext = float(np.min(te["Kn"][le >= F.ls]))
        def rhs(l, y, B):
            l = max(l, 1e-4)
            N, G4, Gc = (float(self.sp_in[k](l)) for k in ("N", "G4", "Gc"))
            rp = float(F.fluido(F.sol.sol(l)[0])[0])
            ch = float(self._chi(l))
            out = N * G4 + (1 - ch) * ((1 + self.eps) * rp * y[0] ** 2 / (4 * N * Gc * Gc) + B * float(self._q(l))
                                       + self.b_in * float(self._qin(l)))
            if ch > 0:
                out += ch * (1 - self.theta * float(self._w(l))) * N * N * float(self.sp_ex["Kn"](l)) / 2
            return [out]
        z0 = float(self.sp_in["N"](1e-4) * self.sp_in["G4"](1e-4) * 1e-4)
        self.Z_obj = float(self.sp_ex["Zn"](self.l2) + self.sp_D(self.l2))

        def fin(B):
            s = solve_ivp(rhs, (1e-4, self.l2), [z0], args=(B,), method="DOP853", rtol=1e-11, atol=1e-13)
            if s.status != 0:
                return np.inf
            return s.y[0, -1] - self.Z_obj
        f0 = fin(0.0)
        self.exceso_sin_B = f0
        if not np.isfinite(f0) or f0 > 0:
            self.B = None
            return
        B1 = 1.0
        while fin(B1) < 0:
            B1 *= 2
        self.B = brentq(fin, 0.0, B1, xtol=1e-14, rtol=1e-13)
        self.Zsol = solve_ivp(rhs, (1e-4, self.l2), [z0], args=(self.B,), method="DOP853", rtol=1e-11, atol=1e-13,
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
        cero = J.const(0.0, NJ, l.size)
        # K y Z naturales solo donde se usan (l >= 0.9 l1): dentro tienen polos
        ext = l >= 0.9 * self.l1
        Zn = J.const(0.0, NJ, l.size)
        Kn = J.const(0.0, NJ, l.size)
        if np.any(ext):
            fe = {k: (J.Jet(v.c[:, ext]) if isinstance(v, J.Jet) else v) for k, v in f.items() if k not in ("g4", "Pf")}
            g4e = [J.Jet(v.c[:, ext]) for v in g4]
            Xe = fe["X"]
            f0e = (fe["H"] - fe["Gc"]) / (4 * Xe * Xe)
            G1e = g4e[0] - 4 * Xe * g4e[1] - 4 * Xe * Xe * g4e[2] - 20 * Xe * Xe * f0e
            Zne = fe["N"] * fe["r"] * fe["Gc"] * fe["H"] / (fe["r"].d() * G1e)
            Kne = 2 * (Zne.d() - fe["N"] * g4e[0]) / (fe["N"] * fe["N"])
            Zn.c[:, ext] = Zne.c
            Kn.c[:, ext] = Kne.c
        chi = self._escalon_chorro(lj, self.l1, self.l2)
        w = 1 - self._escalon_chorro(lj, self.LDa, self.LD)
        nrm = self.wq * math.sqrt(math.pi) / 2 * (1 + math.erf(self.lq / self.wq))
        q = J.exp(-((lj - self.lq) * (1.0 / self.wq)) * ((lj - self.lq) * (1.0 / self.wq))) * (1.0 / nrm)
        qin = J.exp(-(lj * (1.0 / self.w_in)) * (lj * (1.0 / self.w_in))) * (2 / (math.sqrt(math.pi) * self.w_in))
        # Z: diseno hasta l2; despues Z_nat + D con D' = -theta w N^2 K_nat/2; mas alla de LD, Z_nat
        a = l <= self.l2
        b = (l > self.l2) & (l < self.LD)
        Zj = J.const(0.0, NJ, l.size)
        if np.any(a):
            Zv = self.Zsol.sol(l[a])[0]
            sub = lambda v: J.Jet(v.c[:, a])            # noqa: E731
            Na, g4a, rpa, Gca, chia, qa, wa, Kna, qina = (sub(v) for v in (N, g4[0], rp, Gc, chi, q, w, Kn, qin))
            Za = J.const(0.0, NJ, int(a.sum())) + Zv
            for _ in range(NJ + 1):
                Za = (Na * g4a + (1 - chia) * ((1 + self.eps) * rpa * Za * Za / (4 * Na * Gca * Gca) + self.B * qa
                                               + self.b_in * qina)
                      + chia * (1 - self.theta * wa) * Na * Na * Kna / 2).integral(Zv)
            Zj.c[:, a] = Za.c
        if np.any(b):
            sub = lambda v: J.Jet(v.c[:, b])            # noqa: E731
            Nb, wb, Knb = sub(N), sub(w), sub(Kn)
            Db = (-self.theta * wb * Nb * Nb * Knb / 2).integral(self.sp_D(l[b]))
            Zj.c[:, b] = (sub(Zn) + Db).c
        c_ = ~(a | b)
        if np.any(c_):
            Zj.c[:, c_] = Zn.c[:, c_]
        xi = Zj / N
        G1 = r * Gc * H / (r1 * xi)
        f1 = (g4[0] - 4 * X * g4[1] - 4 * X * X * g4[2] - 20 * X * X * f0 - G1) / (8 * X * X * X)
        if np.any(c_):
            f1.c[:, c_] = 0.0                            # Z = Z_nat: F4_X = 0 exactamente
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
                 Kn=Kn, Zj=Zj, Zn=Zn, P_nat=P_nat)
        return V

    def lc(self, umbral=1e-6):
        """Borde de la region con fluido: la superficie si la densidad salta alli (m = 1); con superficie suave, el ultimo
        punto con (rho~+p~) > umbral (rho~+p~)(0) (los pesos 1/(rho~+p~) del fluido divergen en la superficie)."""
        F = self.F
        if F.eos != "exp" and F.m == 1:
            return F.ls
        lt = np.linspace(0, F.ls, 40001)
        rp = F.campos(lt)["rp"]
        return float(lt[np.nonzero(rp > umbral * rp[0])[0][-1]])
