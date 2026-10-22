# -*- coding: utf-8 -*-
"""
Reconstruccion (metodo inverso) de una garganta estatica y simetrica en la teoria cuartica mas alla de Horndeski.

Fondo:  ds^2 = -N(l)^2 dt^2 + dl^2 + r(l)^2 dOmega^2,   N = e^Phi,  Phi = -M/sqrt(l^2 + a^2)
        (masa ADM M en cada lado; N >= e^{-M/a} > 0: ningun horizonte; M = 0 es el caso ultraestatico)
        phi'(l) = phi1 s(l),  s = b^2/(l^2 + b^2)   (kink impar; cola de escalar sin masa phi' ~ D/l^2)
Teoria (invariante bajo phi -> -phi; G3 = 0; x = X/X0 = s^2, X0 = -phi1^2/2):
  G4(X)   = 1 + x^n/(2n - 1)     =>  H = G4 - 2X G4_X = 1 - x^n:  cero doble en la garganta (regularidad del
            termino cinetico del monopolo, ver monopolo_reduccion.py)
  F4      = f0(phi) + f1(phi)(X - X_b) + ...:
            Gcal = G4 - 2X G4_X - 4X^2 F4 = H + delta(l),  delta = d0 e^{-l^2/wd^2}         (fija f0)
            f1 = -lam(l)/(8 X_b^3): lam se elige para que K (cinetico del monopolo) sea el objetivo (ver abajo)
  G2      = g20 + g21 (X - X_b) + g22 (X - X_b)^2/2 + g23 (X - X_b)^3/6:
            g20, g21 de las ligaduras E_A = E_B = 0;  g22 para que P (gradiente del monopolo) sea el objetivo;
            g23 = c23 e^{-l^2/w23^2} (G2_XXX se cancela exactamente en U, ver espectro)
  f2 = F4_XX a lo largo de la curva: cf2 [e^{-(l-lf2)^2/wf2^2} + e^{-(l+lf2)^2/wf2^2}]/X_b^4 (U depende de f2 con
            peso proporcional a Phi'; G2_XXX y F4_XXX se cancelan en U)
Tercera ecuacion del fondo (combinacion de E_A y E_C, exacta):  [N r Gcal (r' - r Phi')]' = N G4, es decir
        r^2 = N^2 [ r0^2/N0^2 + 2 int_0^l I_N/(N^3 Gcal) ],   I_N = int_0^l N G4,
y E_C, E_p se cumplen entonces (Noether), lo que se comprueba con bH_gen.fondo.

Monopolo (G3 = 0, G4 = G4(X): Gamma2 = 0):  identidad exacta en la capa de masas
        N^2 K/2 = W',   W = N [xi - r Gcal (r' - r Phi')],   xi = r Gcal H/(r' Gamma1),
de modo que int N^2 K dl = 2 [W] = 6 (M+ + M-): en el caso ultraestatico (M = 0) K cambia de signo (fantasma).
  Objetivo de K: N^2 K_t/2 = (1 - chi_K) N^2 K0/2 + chi_K A/2, con K0 el de lam = 0, chi_K = e^{-(l/wK)^4} y A tal
      que W_t y W0 tengan el mismo limite;  W_t = W0 + q,  q = int_0^l chi_K (A - N^2 K0)/2,  xi_t = xi0 + q/N,
      lam = Gamma1_t - Gamma1_0 = -(r Gcal H/r') (q/N)/(xi_t xi0).
  Objetivo de P (lineal en G2_02): P_t = (1 - chi_P) P0 + chi_P cm2 N^2 K,  chi_P = e^{-(l/wP)^4}
      (cm2 = velocidad radial fisica al cuadrado, P/(K N^2), dentro de la ventana):
      G2_02 = chi_P (cm2 N^2 K - P0)/(dP/dG2_02).
Derivadas a lo largo de la curva X = X_b(phi):  (T g)_m = g_m'/phi' + phi'' g_{m+1},  G_ij = (T^i g)_j, en
aritmetica de chorros (jets.py): todas las derivadas son exactas salvo redondeo.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp

import jets as J
import simbolico_num as S

NJ = 12               # orden de los chorros
IMAX, JMAX = 3, 5      # G_ij con i <= IMAX, j <= JMAX
GL_X, GL_W = np.polynomial.legendre.leggauss(12)


class Garganta:
    def __init__(self, phi1=1.0, b=1.0, n=2, d0=2.0, wd=1.0, M=0.0, a=1.0, c23=0.0, w23=1.0, r0=1.0,
                 Lmax=400.0, wK=None, wP=None, cm2=1.0, Lq=None, cf2=0.0, wf2=1.0, lf2=0.0):
        self.p = dict(phi1=phi1, b=b, n=n, d0=d0, wd=wd, M=M, a=a, c23=c23, w23=w23, r0=r0, Lmax=Lmax, wK=wK,
                      wP=wP, cm2=cm2, cf2=cf2, wf2=wf2, lf2=lf2)
        self.cf2, self.wf2, self.lf2 = cf2, wf2, lf2
        self.phi1, self.b, self.n, self.d0, self.wd = phi1, b, n, d0, wd
        self.M, self.a, self.c23, self.w23, self.r0 = M, a, c23, w23, r0
        self.X0 = -phi1 ** 2 / 2
        self.Lmax = Lmax
        self.wK, self.wP, self.cm2 = wK, wP, cm2
        self.modo_K = False
        self.modo_P = False
        self._integrar()
        if wK is not None:
            self._preparar_K(Lq if Lq is not None else 12 * wK)
            self.modo_K = True
        if wP is not None:
            self.modo_P = True

    # ---------------------------------------------------------------- perfiles cerrados (valores o chorros)
    def perfiles(self, l):
        b2 = self.b ** 2
        s = b2 / (l * l + b2)
        phip = self.phi1 * s
        X = -phip * phip / 2
        x = X / self.X0
        n = self.n
        G4 = 1 + x ** n / (2 * n - 1)
        H = 1 - x ** n
        delta = self.d0 * J.exp(-l * l / self.wd ** 2)
        Gc = H + delta
        Phi = -self.M / J.sqrt(l * l + self.a ** 2)
        N = J.exp(Phi)
        return dict(s=s, phip=phip, X=X, x=x, G4=G4, H=H, delta=delta, Gc=Gc, Phi=Phi, N=N)

    def _integrar(self):
        N0 = math.exp(-self.M / self.a)
        self.N0 = N0

        def rhs(l, y):
            p = self.perfiles(l)
            return [p["N"] * p["G4"], 2 * y[0] / (p["N"] ** 3 * p["Gc"])]
        sol = solve_ivp(rhs, (0.0, self.Lmax), [0.0, self.r0 ** 2 / N0 ** 2], method="DOP853", rtol=1e-13,
                        atol=1e-14, dense_output=True)
        assert sol.success
        self.sol = sol
        L = self.Lmax
        lj = J.var(np.array([L]), 3)
        V = self._fondo_basico(np.array([L]))
        r, r1 = V["r0"].deriv(0)[0], V["r1"].deriv(0)[0]
        Phi1 = V["Phi"].deriv(1)[0]
        self.masa_MS = r * (1 - r1 ** 2) / 2               # masa de Misner-Sharp en L
        self.masa_komar = r ** 2 * Phi1 * V["N"].deriv(0)[0]
        self.D = r ** 2 * self.perfiles(L)["phip"]
        self.c_asint = r * r1 - L

    def _fondo_basico(self, l):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        lj = J.var(l, NJ)
        p = self.perfiles(lj)
        IN0, u0 = self.sol.sol(l)
        INj = (p["N"] * p["G4"]).integral(IN0)
        uj = (2 * INj / (p["N"] ** 3 * p["Gc"])).integral(u0)
        rj = p["N"] * J.sqrt(uj)
        return dict(l=lj, p=p, r0=rj, r1=rj.d(), Phi=p["Phi"], N=p["N"])

    # ---------------------------------------------------------------- objetivo de K (ligadura integral)
    def _K0(self, l):
        modo = self.modo_K
        self.modo_K = False
        V = self._V(l, con_P=False)
        K0 = self._mono(V, que=("K",))["K"]
        self.modo_K = modo
        return K0, V

    @staticmethod
    def _chi(lj, w):
        return J.exp(-(lj / w) ** 4)

    def _preparar_K(self, Lq):
        lmin = 1e-3 * self.wK
        bordes = np.concatenate([[0.0], np.geomspace(lmin, Lq, 400)])
        self.q_bordes = bordes
        a, b_ = bordes[:-1], bordes[1:]
        nodos = (a[:, None] + b_[:, None]) / 2 + (b_ - a)[:, None] / 2 * GL_X[None, :]
        pesos = (b_ - a)[:, None] / 2 * GL_W[None, :]
        ln = np.maximum(nodos.ravel(), lmin)
        K0, V0 = self._K0(ln)
        N2K0 = (V0["E"] * K0).deriv(0)
        chi = np.exp(-(nodos.ravel() / self.wK) ** 4)
        self.A_K = np.sum((chi * N2K0).reshape(nodos.shape) * pesos) / np.sum(chi.reshape(nodos.shape) * pesos)
        integ = (chi * (self.A_K - N2K0) / 2).reshape(nodos.shape)
        self.q_en_bordes = np.concatenate([[0.0], np.cumsum(np.sum(integ * pesos, axis=1))])
        self.q_inf = self.q_en_bordes[-1]

    def _q(self, l):
        l = np.atleast_1d(l)
        idx = np.clip(np.searchsorted(self.q_bordes, l) - 1, 0, len(self.q_bordes) - 2)
        a = self.q_bordes[idx]
        nodos = a[:, None] + (l - a)[:, None] / 2 * (1 + GL_X[None, :])
        pesos = (l - a)[:, None] / 2 * GL_W[None, :]
        ln = np.maximum(nodos.ravel(), 1e-3 * self.wK)
        K0, V0 = self._K0(ln)
        N2K0 = (V0["E"] * K0).deriv(0)
        chi = np.exp(-(nodos.ravel() / self.wK) ** 4)
        return self.q_en_bordes[idx] + np.sum((chi * (self.A_K - N2K0) / 2).reshape(nodos.shape) * pesos, axis=1)

    # ---------------------------------------------------------------- chorros de todas las cantidades en puntos l
    def _V(self, l, con_P=True, lam=None):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        assert np.all(l > 0) or (np.all(l >= 0) and not self.modo_K)
        B = self._fondo_basico(l)
        lj, p, rj = B["l"], B["p"], B["r0"]
        Phij, Nj = p["Phi"], p["N"]
        phip = p["phip"]
        phij = phip.integral(self.phi1 * self.b * np.arctan(l / self.b))
        phipp = phip.d()
        X = p["X"]
        n = self.n
        x = p["x"]
        cero = J.const(0.0, NJ, l.size)
        g4 = [p["G4"]]
        for j in range(1, JMAX + IMAX + 1):
            if j <= n:
                g4.append(math.factorial(n) / math.factorial(n - j) * x ** (n - j) / (self.X0 ** j * (2 * n - 1)))
            else:
                g4.append(cero)
        f1 = cero if lam is None else -lam / (8 * X * X * X)
        # f2 = F4_XX a lo largo de la curva: libre (U depende de el con peso ~ Phi'); f2 = cf2 e^{-(l-lf2)^2/wf2^2}/X^4
        # (simetrizado en l) para que su efecto sea O(cf2)
        if self.cf2:
            f2 = self.cf2 * (J.exp(-(lj - self.lf2) ** 2 / self.wf2 ** 2) + J.exp(-(lj + self.lf2) ** 2 / self.wf2 ** 2)) \
                / (X * X * X * X)
        else:
            f2 = cero
        f = [-p["delta"] / (4 * X * X), f1, f2] + [cero] * (JMAX + IMAX - 2)
        g3 = [cero] * (JMAX + IMAX + 1)

        def T(seq):
            out = []
            for m in range(len(seq)):
                sig = seq[m + 1] if m + 1 < len(seq) else cero
                out.append(seq[m].d() / phip + phipp * sig)
            return out

        Ph1 = Phij.d()
        V = {"r0": rj, "r1": rj.d(), "r2": rj.d().d(), "r3": rj.d().d().d(), "r4": rj.d().d().d().d(),
             "P0": phij, "P1": phip, "P2": phipp, "P3": phipp.d(), "P4": phipp.d().d(),
             "Phi1": Ph1, "Phi2": Ph1.d(), "Phi3": Ph1.d().d(), "Phi4": Ph1.d().d().d(),
             "E": Nj * Nj, "eP": Nj}
        for nombre, seq in (("G4", g4), ("F4", f), ("G3", g3)):
            s_ = seq
            for i in range(IMAX + 1):
                for j in range(JMAX + 1):
                    V[f"{nombre}_{i}{j}"] = s_[j]
                s_ = T(s_)
        g20 = S.G2_DE_LIGADURAS[0](V)
        V["G2_00"] = g20
        g21 = S.G2_DE_LIGADURAS[1](V)
        V["G2_01"] = g21
        V["G2_02"] = cero
        g22 = cero
        if con_P and self.modo_P:
            m0 = self._mono(V, que=("K", "P"))
            dP = S.MONO["dP_dG202"](self._W(V))
            g22 = self._chi(lj, self.wP) * (self.cm2 * Nj * Nj * m0["K"] - m0["P"]) / dP
        g23 = self.c23 * J.exp(-lj * lj / self.w23 ** 2)
        g2 = [g20, g21, g22, g23] + [cero] * (JMAX + IMAX - 3)
        s_ = g2
        for i in range(IMAX + 1):
            for j in range(JMAX + 1):
                V[f"G2_{i}{j}"] = s_[j]
            s_ = T(s_)
        V.update(l=lj, Gc=p["Gc"], H=p["H"], G4=p["G4"], delta=p["delta"], X=X, N=Nj, Phi=Phij)
        return V

    def V(self, l):
        if not self.modo_K:
            return self._V(l)
        l = np.atleast_1d(np.asarray(l, dtype=float))
        K0, V0 = self._K0(l)
        lj = V0["l"]
        r, r1, Nj = V0["r0"], V0["r1"], V0["N"]
        gam0 = S.MONO["gam"](V0)
        Gamma1_0 = Nj * gam0 / (4 * r * r1)
        pref = r * V0["Gc"] * V0["H"] / r1
        xi0 = pref / Gamma1_0
        chiK = self._chi(lj, self.wK)
        qj = (chiK * (self.A_K - V0["E"] * K0) / 2).integral(self._q(l))
        xit = xi0 + qj / Nj
        # fuera de la ventana q es constante (= 0 salvo redondeo); f1 = -lam/(8 X^3) amplificaria ese redondeo por
        # X^{-3} ~ l^12: se corta lam con un factor que vale 1 hasta mucho mas alla de la ventana
        corte = J.exp(-(lj / (3 * self.wK)) ** 8)
        lam = -pref * (qj / Nj) / (xit * xi0) * corte
        return self._V(l, lam=lam)

    # ---------------------------------------------------------------- comprobaciones y coeficientes
    @staticmethod
    def _W(V):
        gam = S.MONO["gam"](V)
        al1 = S.MONO["al1"](V)
        al0 = S.MONO["al0"](V)
        W = dict(V)
        gk, a1, a0 = gam, al1, al0
        for k in range(4):
            W[f"gam{k}"] = gk
            W[f"al1_{k}"] = a1
            W[f"al0_{k}"] = a0
            gk, a1, a0 = gk.d(), a1.d(), a0.d()
        return W

    @classmethod
    def _mono(cls, V, que=("K", "P", "U")):
        W = cls._W(V)
        out = {"gam": W["gam0"], "al1": W["al1_0"]}
        for k in que:
            out[k] = S.MONO[k](W)
        return out

    def residuos_fondo(self, l):
        """E_u con bH_gen (codigo independiente de las expresiones usadas en la reconstruccion)."""
        import bH_gen
        V = self.V(l)
        out = []
        for i in range(len(np.atleast_1d(l))):
            Gi = {k: float(V[k].c[0][i]) for k in bH_gen.SIMB_G}
            r_ = [float(V[f"r{k}"].deriv(0)[i]) for k in range(3)]
            P_ = [float(V[f"P{k}"].deriv(0)[i]) for k in range(3)]
            Ph = [float(V["Phi"].deriv(0)[i]), float(V["Phi1"].deriv(0)[i]), float(V["Phi2"].deriv(0)[i])]
            out.append(bH_gen.fondo(Gi, Ph, r_, P_))
        return np.array(out)

    def monopolo(self, l):
        V = self.V(l)
        out = self._mono(V)
        out["V"] = V
        return out

    def impar(self, l, LL=2):
        """Sector impar: alfa, beta, gamma (expresiones completas de impar_ecuaciones) y la accion de chi
        K_o = -alfa^2/gamma, P_o = alfa^2/beta, U_o = alfa + (alfa' + 2 alfa r'/r)^2/beta - [alfa (alfa' + 2 alfa r'/r)/beta]'."""
        V = self.V(l)
        al = S.IMPAR[LL]["alfa"](V)
        be = S.IMPAR[LL]["beta"](V)
        ga = S.IMPAR[LL]["gamma"](V)
        r = V["r0"]
        q = al.d() + 2 * al * V["r1"] / r
        return dict(alfa=al, beta=be, gamma=ga, K=-al * al / ga, P=al * al / be, U=al + q * q / be - (al * q / be).d(),
                    V=V)
