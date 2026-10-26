# -*- coding: utf-8 -*-
"""
Fase 2h: fondo del hibrido denso (garganta de GLPV cuartica llena de fluido) con una zona lejana canonica.

Region interior (l <= lm): la de la Fase 2g (criba.FondoFluido): lapso Phi = -M/sqrt(l^2 + a^2), perfiles de G4 = S G4b
y Gcal = S (Hb + delta), fluido P(Y) y r de T1 con materia.  Con este lapso la K natural lejana es la de la Fase 2f,
K ~ 2 y/l^2 con y ~ 4 M^2 ln^2 l, y en los hibridos densos es negativa hasta l ~ 10^3 (las constantes de integracion
de la region interior son grandes): no se puede empalmar con ella.

Zona lejana (l > lm): se prescribe la combinacion de la condicion de energia nula
    nu = Phi' r' - r''     (K_nat ~ 2 r nu/(N r'^2) lejos; nu = (kappa/2) r phi'^2 para un escalar canonico en RG),
como nu = (1 - chi) nu_int + chi (kappa/2) r phi'^2, con chi un escalon C^infinito (exactamente 0 en lm, 1 en lm + Delta),
y N, r salen de T1 sin materia, [N r Gcal (r' - r Phi')]' = N G4, escrita como sistema de segundo orden:
    N'' = [Gcal N (r'^2 - r nu) - N G4 + Gcal' (N r r' - r^2 N')]/(r^2 Gcal),     r'' = r' N'/N - nu.
Para l <= lm (chi = 0) la solucion coincide con la interior.  Lejos el fondo es el de RG con un escalar canonico de
normalizacion kappa: K_nat > 0, y el lapso tiende a una constante N_inf (se informa; reescalar t no cambia nada).
La masa ADM sale de la asintotica (M = lim l^2 N'/N, comprobada con la de Misner-Sharp).

La clase Fondo expone la misma interfaz que criba.FondoFluido (campos, escalar, lapso, fluido, E_out, ls, M, N0...), de
modo que las condiciones necesarias de la Fase 2g (criba.criba_normalizada, presupuesto_fantasma, impar) se aplican
igual; el presupuesto de fantasma usa 3 N_inf M_ADM.
"""
import os
import sys
import math
import numpy as np
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "fase2g"))
import criba as C          # noqa: E402


def escalon(t):
    """Escalon C^infinito: 0 para t <= 0, 1 para t >= 1."""
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    m = (t > 0) & (t < 1)
    a = np.exp(-1.0 / np.where(m, t, 0.5))
    b = np.exp(-1.0 / np.where(m, 1 - t, 0.5))
    out = np.where(m, a / (a + b), out)
    return np.where(t >= 1, 1.0, out)


def escalon_d(t):
    """Derivada del escalon respecto de t."""
    t = np.asarray(t, dtype=float)
    m = (t > 0) & (t < 1)
    tt = np.where(m, t, 0.5)
    a = np.exp(-1.0 / tt)
    b = np.exp(-1.0 / (1 - tt))
    da = a / tt ** 2
    db = -b / (1 - tt) ** 2
    d = (da * (a + b) - a * (da + db)) / (a + b) ** 2
    return np.where(m, d, 0.0)


class Fondo(C.FondoFluido):
    def __init__(self, lm=None, Delta=None, kappa=1.0, phi1=1.0, beta_nu=None, **kw):
        self._lejos = False
        self.kappa, self.phi1 = kappa, phi1
        super().__init__(**kw)
        if not self.ok:
            return
        self.lm = lm if lm is not None else 1.6 * self.lT + 2 * self.ls
        self.Delta = Delta if Delta is not None else self.lm
        if beta_nu is None:
            self.beta_nu = self._ajustar_beta()
        else:
            self.beta_nu = beta_nu
        if self.ok:
            self._integrar_lejos()

    def _bump(self, l):
        """Joroba 4 chi (1 - chi)/l^3 en la zona de mezcla (anade energia nula positiva: masa)."""
        ch = escalon((np.asarray(l, dtype=float) - self.lm) / self.Delta)
        return 4 * ch * (1 - ch) / np.asarray(l, dtype=float) ** 3

    def _ajustar_beta(self):
        """Disparo en beta_nu para que la solucion lejana sea asintoticamente plana: el lapso tiene un modo que crece
        linealmente (N ~ A + B l); B = 0 equivale a que la masa de Komar l^2 Phi' y la de Misner-Sharp r (1 - r'^2)/2
        coincidan en L.  Devuelve beta_nu (o marca el fondo como invalido)."""
        from scipy.optimize import brentq

        def f(b):
            self.beta_nu = b
            self._integrar_lejos(comprobar=False)
            if not self._lejos:
                return np.nan
            N, Np, r, rp = self.sol_lejos.sol(self.Lmax)
            return float(self.Lmax ** 2 * Np / N - r * (1 - rp ** 2) / 2)
        bs = [0.0]
        fs = [f(0.0)]
        paso = 100.0
        for k in range(30):
            for sgn in (1, -1):
                b = sgn * paso * (1.6 ** k)
                v = f(b)
                if np.isfinite(v) and np.isfinite(fs[0]) and np.sign(v) != np.sign(fs[0]):
                    b0 = brentq(f, min(0.0, b), max(0.0, b), xtol=1e-6, rtol=1e-10)
                    self._lejos = False
                    return b0
        self.ok, self.motivo = False, "no se encuentra beta_nu (zona lejana no plana)"
        self._lejos = False
        return 0.0

    # -------------------------------------------------------------- nu de la solucion interior (lapso analitico + T1)
    def nu_int(self, l):
        k = C.FondoFluido.K_lejos(self, l)                          # usa el lapso analitico y el ODE interior
        return k["y"] / np.atleast_1d(l) ** 3

    def phip(self, l):
        b2 = self.b ** 2
        return self.phi1 * b2 / (np.asarray(l) ** 2 + b2)

    def _integrar_lejos(self, comprobar=True):
        self._lejos = False
        lm, D = self.lm, self.Delta
        if lm >= self.Lmax:
            raise ValueError("lm >= Lmax")
        # estado inicial en lm desde la solucion interior
        f = C.FondoFluido.campos(self, np.array([lm]))
        N0, Phi1 = f["N"][0], f["Phi1"][0]
        y0 = [N0, N0 * Phi1, f["r"][0], f["r1"][0]]

        def nu(l, r):
            t = (l - lm) / D
            ch = float(escalon(np.array([t]))[0])
            ni = float(self.nu_int(np.array([l]))[0]) if ch < 1 else 0.0
            return (1 - ch) * ni + ch * 0.5 * self.kappa * r * float(self.phip(l)) ** 2 + self.beta_nu * float(self._bump(l))

        def rhs(l, y):
            N, Np, r, rp = y
            G4, H, Gc, Gcp = C.FondoFluido.escalar(self, np.array([l]))
            G4, Gc, Gcp = G4[0], Gc[0], Gcp[0]
            n_ = nu(l, r)
            Npp = (Gc * N * (rp * rp - r * n_) - N * G4 + Gcp * (N * r * rp - r * r * Np)) / (r * r * Gc)
            rpp = rp * Np / N - n_
            return [Np, Npp, rp, rpp]
        def malo(l, y):
            return min(y[0], y[2]) - 1e-8
        malo.terminal = True
        s = solve_ivp(rhs, (lm, self.Lmax), y0, method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True, events=malo)
        if not s.success or s.status == 1:
            if comprobar:
                self.ok, self.motivo = False, "zona lejana: " + s.message
            return
        self.sol_lejos = s
        self._lejos = True
        L = self.Lmax
        N, Np, r, rp = s.sol(L)
        self.N_inf = float(N * math.exp(self.M_aprox(L, N, Np)[1]))
        self.M_adm, _ = self.M_aprox(L, N, Np)
        self.m_MS = float(r * (1 - rp ** 2) / 2)
        if comprobar:
            lt = np.geomspace(lm, L, 400)
            ff = self.campos(lt)
            if np.any(ff["r1"] <= 0):
                self.ok, self.motivo = False, "zona lejana: r' <= 0"

    def M_aprox(self, L, N, Np):
        """M ~ l^2 Phi' (con la correccion de orden 1/l del lapso de Schwarzschild en distancia propia) y
        ln(N_inf/N(L)) ~ M/L."""
        Phi1 = Np / N
        M = float(L * L * Phi1)
        return M, M / L

    # -------------------------------------------------------------- interfaz de criba.FondoFluido
    def lapso(self, l):
        if not self._lejos or np.ndim(l) == 0 and float(l) <= self.lm:
            return C.FondoFluido.lapso(self, l)
        escalar = np.ndim(l) == 0
        l = np.atleast_1d(np.asarray(l, dtype=float))
        N, Phi1 = C.FondoFluido.lapso(self, l)
        if self._lejos:
            m = l > self.lm
            if np.any(m):
                N_, Np_, r_, rp_ = self.sol_lejos.sol(l[m])
                N = N.copy()
                Phi1 = Phi1.copy()
                N[m], Phi1[m] = N_, Np_ / N_
        if escalar:
            return N[0], Phi1[0]
        return N, Phi1

    def campos(self, l):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        f = C.FondoFluido.campos(self, l)
        if self._lejos:
            m = l > self.lm
            if np.any(m):
                N_, Np_, r_, rp_ = self.sol_lejos.sol(l[m])
                G4, H, Gc, Gcp = C.FondoFluido.escalar(self, l[m])
                f = {k: (v.copy() if isinstance(v, np.ndarray) and v.shape == l.shape else v) for k, v in f.items()}
                f["N"][m] = N_
                f["Phi1"][m] = Np_ / N_
                f["r"][m] = r_
                f["r1"][m] = rp_
                f["sig"][m] = rp_ / r_ - Np_ / N_
                f["rp"][m] = 0.0
                f["cs2"][m] = 0.0
                f["p"][m] = 0.0
                a2 = 4 * N_ * N_ * r_ * Gc
                f["a2"][m] = a2
                f["a2p"][m] = a2 * (2 * Np_ / N_ + rp_ / r_ + Gcp / Gc)
        return f

    def K_lejos(self, l):
        """K natural aproximada 2 r nu/(N r'^2) con el nu efectivo (interior hasta lm, prescrito despues)."""
        l = np.atleast_1d(np.asarray(l, dtype=float))
        f = self.campos(l)
        t = (l - self.lm) / self.Delta
        ch = escalon(t)
        nu = np.where(l > self.lm, (1 - ch) * np.where(ch < 1, self.nu_int(l), 0.0)
                      + ch * 0.5 * self.kappa * f["r"] * self.phip(l) ** 2 + self.beta_nu * self._bump(l), self.nu_int(l))
        return dict(K=2 * f["r"] * nu / (f["N"] * f["r1"] ** 2), y=nu * l ** 3, nu=nu, N=f["N"], r=f["r"])

    @property
    def M_presupuesto(self):
        """3 N_inf M_ADM es el limite de W = N [xi - r Gcal (r' - r Phi')] (con N -> N_inf)."""
        return self.N_inf * self.M_adm
