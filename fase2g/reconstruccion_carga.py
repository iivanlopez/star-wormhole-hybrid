# -*- coding: utf-8 -*-
"""
Fase 2g, parte IV: garganta de GLPV cuartica con un campo electrico radial que la atraviesa ("carga sin carga").

El campo de Maxwell no tiene grados de libertad en l = 0; tras la ley de Gauss su accion es L_em = -q2 sqrt(A B + b^2)/C
(q2 = 2 Q^2 en unidades G = c = 1: rho~ = -p~_r = p~_t = q2/r^4).  En el monopolo (estructura_em.py) la ecuacion de a
sigue siendo (N^2 W)' = 0: misma ligadura W = 0 y misma variable que en el vacio, con
    k_c -> k_c + q2 N/r^2,   e0 -> e0 + 2 q2 N/r^3,   u_d -> u_d + 4 q2 N/r^4
(la accion de la carga en (c, d) es q2 N (c^2/2 + 2 c d/r - 4 d^2/r^2)/r^2), y por T3
    P -> P - (q2 N/(2 r^2)) A1^2,
    U -> U + 4 q2 N/r^4 - (q2 N/(2 r^2)) A0^2 + A0 2 q2 N/r^3 + (1/2) [ (q2 N/r^2) A1 A0 - A1 2 q2 N/r^3 ]',
A1 = al1/gam, A0 = -Phi' A1.  K no cambia de forma (T2).  T1 con carga: [N r Gcal (r' - r Phi')]' = N G4 - N q2/r^2.
Teorema con carga (prueba d = N):  E[N] = int [ -2 N^3 Gcal (r'/r - Phi')^2 + 4 q2 N^3/r^4 ] dl.
Identidad de K con carga: N^2 K/2 = W' - N q2/r^2  =>  int N^2 K = 6 (M+ + M-) - 2 int N q2/r^2.

GargantaCargada hereda de fase2f/reconstruccion.Garganta (mismos perfiles de fondo y de teoria, mismas ventanas de K y P)
y cambia: el fondo (T1 con carga), las ligaduras de G2, G2_X (con la carga) y los coeficientes P, U del monopolo.
"""
import os
import sys
import math
import pickle
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
import jets as J                       # noqa: E402
import simbolico_num as S              # noqa: E402
import reconstruccion as R             # noqa: E402

_em = pickle.load(open(os.path.join(AQUI, "resultados", "em_simbolico.pkl"), "rb"))
_q2 = _em["q2"]
_EA = S.con_Xb(_em["fondo_eq"]["A"])
_EB = S.con_Xb(_em["fondo_eq"]["B"])
_g200, _g201 = sp.Symbol("G2_00", real=True), sp.Symbol("G2_01", real=True)
_sA = sp.solve(_EA, _g200)
_sB = sp.solve(_EB, _g201)
assert len(_sA) == 1 and len(_sB) == 1
LIG_Q = (S.compilar(_sA[0]), S.compilar(_sB[0]))


class GargantaCargada(R.Garganta):
    def __init__(self, q2=0.0, **kw):
        self.q2 = q2
        super().__init__(**kw)

    # -------------------------------------------------------------- fondo: T1 con carga
    def _integrar(self):
        N0 = math.exp(-self.M / self.a)
        self.N0 = N0
        q2 = self.q2

        def rhs(l, y):
            p = self.perfiles(l)
            N = p["N"]
            r2 = N * N * y[1]
            return [N * p["G4"] - N * q2 / r2, 2 * y[0] / (N ** 3 * p["Gc"])]

        def cero(l, y):
            return y[1]
        cero.terminal = True
        sol = solve_ivp(rhs, (0.0, self.Lmax), [0.0, self.r0 ** 2 / N0 ** 2], method="DOP853", rtol=1e-13,
                        atol=1e-14, dense_output=True, events=cero)
        self.valido = sol.status == 0
        assert sol.success
        self.sol = sol
        if not self.valido:
            return
        L = self.Lmax
        V = self._fondo_basico(np.array([L]))
        r, r1 = V["r0"].deriv(0)[0], V["r1"].deriv(0)[0]
        Phi1 = V["Phi"].deriv(1)[0]
        self.masa_MS = r * (1 - r1 ** 2) / 2 + q2 / (4 * r)       # Misner-Sharp con la carga (q2 = 2 Q^2)
        self.masa_komar = r ** 2 * Phi1 * V["N"].deriv(0)[0]
        self.D = r ** 2 * self.perfiles(L)["phip"]
        self.c_asint = r * r1 - L

    def _fondo_basico(self, l):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        lj = J.var(l, R.NJ)
        p = self.perfiles(lj)
        IN0, u0 = self.sol.sol(l)
        N = p["N"]
        uj = J.const(0.0, R.NJ, l.size) + u0
        INj = J.const(0.0, R.NJ, l.size) + IN0
        for _ in range(R.NJ + 1):
            INj = (N * p["G4"] - self.q2 / (N * uj)).integral(IN0)          # N q2/r^2 = q2/(N u)
            uj = (2 * INj / (N ** 3 * p["Gc"])).integral(u0)
        rj = N * J.sqrt(uj)
        return dict(l=lj, p=p, r0=rj, r1=rj.d(), Phi=p["Phi"], N=N)

    # -------------------------------------------------------------- ligaduras con carga y coeficientes del monopolo
    def _V(self, l, con_P=True, lam=None):
        # mismo procedimiento que la Fase 2f, pero las ligaduras incluyen la carga y P incluye su correccion
        V = super()._V(l, con_P=False, lam=lam)
        lj = V["l"]
        cero = J.const(0.0, R.NJ, lj.c.shape[1])
        V["q2"] = cero + self.q2
        g20 = LIG_Q[0](V)
        V["G2_00"] = g20
        g21 = LIG_Q[1](V)
        V["G2_01"] = g21
        V["G2_02"] = cero
        g22 = cero
        if con_P and self.modo_P:
            m0 = self._mono(V, que=("K", "P"))
            dP = S.MONO["dP_dG202"](self._W(V))
            g22 = self._chi(lj, self.wP) * (self.cm2 * V["N"] * V["N"] * m0["K"] - m0["P"]) / dP
        g23 = self.c23 * J.exp(-lj * lj / self.w23 ** 2)
        g2 = [g20, g21, g22, g23] + [cero] * (R.JMAX + R.IMAX - 3)

        def T(seq):
            phip, phipp = V["P1"], V["P2"]
            out = []
            for m in range(len(seq)):
                sig = seq[m + 1] if m + 1 < len(seq) else cero
                out.append(seq[m].d() / phip + phipp * sig)
            return out
        s_ = g2
        for i in range(R.IMAX + 1):
            for j in range(R.JMAX + 1):
                V[f"G2_{i}{j}"] = s_[j]
            s_ = T(s_)
        return V

    def _mono(self, V, que=("K", "P", "U")):
        W = self._W(V)
        out = {"gam": W["gam0"], "al1": W["al1_0"]}
        for k in que:
            out[k] = S.MONO[k](W)
        if self.q2 and ("P" in que or "U" in que):
            q2 = self.q2
            N, r = V["N"], V["r0"]
            A1 = W["al1_0"] / W["gam0"]
            A0 = -V["Phi1"] * A1
            dkc = q2 * N / (r * r)
            de0 = 2 * q2 * N / (r * r * r)
            if "P" in que:
                out["P"] = out["P"] - dkc / 2 * A1 * A1
            if "U" in que:
                out["U"] = (out["U"] + 4 * q2 * N / r ** 4 - dkc / 2 * A0 * A0 + A0 * de0
                            + 0.5 * (dkc * A1 * A0 - A1 * de0).d())
        return out

    def monopolo(self, l):
        V = self.V(l)
        out = self._mono(V)
        out["V"] = V
        return out
