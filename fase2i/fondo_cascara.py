# -*- coding: utf-8 -*-
"""
Fase 2i: fondos con una cascara elastica alrededor de una garganta de vacio (GLPV cuartica), sobre la construccion de
fase2h/fondo_nu.FondoNu (violacion de la condicion de energia nula confinada a la garganta, cola canonica ajustada
para la planitud).

Materia (solo en [la, lb], x = (l - la)/(lb - la)), prescrita como funcion de la distancia propia:
    w = rho~ + p~_r = W0 (1 - s_w x)                       (densidad finita en las dos superficies)
    p~_r = pi_p W0 x^ap (1 - x)^bp / max                  (se anula en las dos superficies libres)
    p~_t = p~_r + (r/(2 r')) (p~_r' + w Phi')             (equilibrio anisotropo)
y modulos de un solido isotropo de velocidades c_L (longitudinal radial) y c_T (cizalla) (criba_solido.modulos).
Geometria (como FondoNu, con la materia anisotropa):
    nu = Phi' r' - r'' = -nu0 e^{-(l/wt)^2} + f_m r w/(4 G4) + (kappa/2) r phi'^2
    N'' = [Gcal N (r'^2 - r nu) - N G4 + N r^2 (rho~ + p~_t)/2 + Gcal' (N r r' - r^2 N')]/(r^2 Gcal)   (T1)
    r'' = r' N'/N - nu.
T1 lleva rho~ + p~_t (de -N^2 E_A + r^2 E_C/2 con E_A^m = -r^2 rho~/(2N), E_C^m = N p~_t; generar_solido.py) y la parte
"de RG" de nu lleva rho~ + p~_r (en RG nu = r (rho~ + p~_r)/4).
"""
import math
import os
import sys
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2g"), os.path.join(AQUI, "..", "fase2h"), AQUI):
    if p_ not in sys.path:
        sys.path.insert(0, p_)
import criba as C                 # noqa: E402
import fondo_nu as FN             # noqa: E402
import criba_solido as CS         # noqa: E402


class FondoCascara(FN.FondoNu):
    def __init__(self, la=3.0, lb=20.0, W0=0.01, s_w=0.0, pi_p=0.1, ap=1.0, bp=1.0, cL2=1.0, cT2=0.3, f_m=1.0, **kw):
        self.la, self.lb = la, lb
        self.W0, self.s_w, self.pi_p, self.ap, self.bp = W0, s_w, pi_p, ap, bp
        self.cL2, self.cT2 = cL2, cT2
        xs = ap / (ap + bp)
        self._norm = xs ** ap * (1 - xs) ** bp
        kw = dict(kw)
        kw["f_fl"] = f_m
        kw["A"] = 0.0
        super().__init__(**kw)

    # -------------------------------------------------------------- materia prescrita
    def _materia(self, l):
        """w, w', p~_r, p~_r' (escalares) en l."""
        if not (self.la <= l <= self.lb):
            return 0.0, 0.0, 0.0, 0.0
        L = self.lb - self.la
        x = (l - self.la) / L
        w = self.W0 * (1 - self.s_w * x)
        wp = -self.W0 * self.s_w / L
        if 0.0 < x < 1.0:
            B = x ** self.ap * (1 - x) ** self.bp / self._norm
            Bp = B * (self.ap / x - self.bp / (1 - x))
        else:
            B = 0.0
            # pendiente en los bordes (ap = 1 o bp = 1 la dejan finita)
            if x <= 0.0:
                Bp = (1.0 / self._norm) if self.ap == 1 else 0.0
            else:
                Bp = (-1.0 / self._norm) if self.bp == 1 else 0.0
        P = self.pi_p * self.W0
        return w, wp, P * B, P * Bp / L

    def _materia_v(self, l):
        out = np.array([self._materia(li) for li in np.atleast_1d(l)])
        return out[:, 0], out[:, 1], out[:, 2], out[:, 3]

    # -------------------------------------------------------------- sistema del fondo
    def _rhs(self, l, y):
        N, Np, r, rp_ = y
        G4, Gc, Gcp = self._perfil_escalar(l)
        w, _, pr, prp = self._materia(l)
        ph = self.phi1 * self.b * self.b / (l * l + self.b * self.b)
        nu = (-self.nu0_actual * math.exp(-(l / self.wt) ** 2) + self.f_fl * r * w / (4 * G4)
              + 0.5 * self.kappa * r * ph * ph)
        if w > 0.0:
            rpt = w + r / (2 * max(rp_, 1e-7)) * (prp + w * Np / N)
        else:
            rpt = 0.0
        Npp = (Gc * N * (rp_ * rp_ - r * nu) - N * G4 + 0.5 * N * r * r * rpt + Gcp * (N * r * rp_ - r * r * Np)) \
            / (r * r * Gc)
        return [Np, Npp, rp_, rp_ * Np / N - nu]

    def _resolver(self, nu0, hasta):
        """Como FondoNu._resolver, con un evento mas: r' vuelve a anularse (ecuador; la cascara no se sostiene)."""
        from scipy.integrate import solve_ivp
        self.nu0_actual = nu0

        def malo(l, y):
            return min(y[0], y[2]) - 1e-9
        malo.terminal = True

        def explota(l, y):
            return 1e8 - abs(y[0])
        explota.terminal = True

        def ecuador(l, y):
            return y[3] - 1e-7 if l > 0.05 else 1.0
        ecuador.terminal = True
        ecuador.direction = -1
        return solve_ivp(self._rhs, (0.0, hasta), [1.0, 0.0, self.r0, 0.0], method="DOP853", rtol=1e-11, atol=1e-13,
                         dense_output=True, events=(malo, explota, ecuador))

    def _clasificar(self, s):
        if s.status == 1 and len(s.t_events[2]):
            return -1
        return super()._clasificar(s)

    def _post(self):
        super()._post()
        self.ls = self.lb
        if not self.ok:
            return
        # la cascara debe estar donde r' > 0 (fuera de la garganta)
        f = self.campos(np.array([self.la]))
        if self.la <= 0 or f["r1"][0] <= 1e-6:
            self.ok, self.motivo = False, "cascara sin r' > 0 en su borde interior"

    # -------------------------------------------------------------- perfiles para la criba
    def perfil(self, l):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        Y = self.sol.sol(l)
        der = np.array([self._rhs(li, yi) for li, yi in zip(l, Y.T)])
        Nu, Np, r, r1 = Y
        Phi1 = Np / Nu
        Phi2 = der[:, 1] / Nu - Phi1 ** 2
        G4, H, Gc, Gcp = self.escalar(l)
        w, wp, pr, prp = self._materia_v(l)
        pt = pr + r / (2 * r1) * (prp + w * Phi1)
        rho = w - pr
        out = dict(N=Nu / self._Ninf, Phi1=Phi1, Phi2=Phi2, r=r, r1=r1, r2=der[:, 3], Gc=Gc, Gcp=Gcp, G4=G4,
                   sig=r1 / r - Phi1, w=w, wp=wp, rho=rho, pr=pr, pt=pt)
        out.update(CS.modulos(rho, pr, pt, self.cL2, self.cT2))
        return out

    def integrales_hueco(self, n=4001):
        l = np.linspace(0.0, self.la, n)
        f = self.campos(l)
        G4 = f["G4"]
        g1 = f["N"] ** 3 * f["Gc"] * f["sig"] ** 2
        return (float(np.trapezoid(g1, l)), float(np.trapezoid(G4 * f["N"], l)),
                float(np.trapezoid(g1 + f["N"] ** 3 * f["Gc"] / f["r"] ** 2, l)))

    def fisica(self, n=400):
        """Condiciones sobre la materia: min rho~, min (rho~ - |p~_t|) y (rho~ - |p~_r|) (dominante), min e_v'', y
        masa propia de la cascara (int 4 pi r^2 rho dl en unidades G = c = 1: rho = rho~/16 pi)."""
        l = np.linspace(self.la, self.lb, n)
        p = self.perfil(l)
        return dict(rho_min=float(p["rho"].min()), dec_t=float(np.min(p["rho"] - np.abs(p["pt"]))),
                    dec_r=float(np.min(p["rho"] - np.abs(p["pr"]))), ev_min=float(np.min(p["Ev"] / p["w"])),
                    masa=float(np.trapezoid(p["r"] ** 2 * p["rho"], l) / 4), pt_max=float(np.max(p["pt"] / p["rho"])),
                    pt_min=float(np.min(p["pt"] / p["rho"])))
