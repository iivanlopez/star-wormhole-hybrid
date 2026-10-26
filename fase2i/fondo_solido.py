# -*- coding: utf-8 -*-
"""
Fase 2i: estrellas elasticas con un nucleo hueco alrededor de la garganta (GLPV cuartica), sobre fase2h/fondo_nu.FondoNu.

El volumen sigue el perfil de Tolman de un fluido rigido (m = 1, P = A (Y - Ys), Y = 1/(2 N^2) con el lapso sin
normalizar; superficie exterior donde N = Ns, con salto de densidad), y la materia empieza en l = la (nucleo de vacio
alrededor de la garganta, superficie interior libre):
    w = rho~ + p~_r = 2 A Y,        p~_r = A (Y - Ys) chi(l),     chi = 1 - e^{-(l - la)/lchi}     (l >= la)
    p~_t = p~_r + (r/(2 r')) [w Phi' (1 - chi) + A (Y - Ys) chi']     (equilibrio anisotropo; >= p~_r)
Con chi = 1 es el fluido de la Fase 2h (p~_t = p~_r).  Las capas interiores se sostienen con tension tangencial de
compresion (cupula).  Modulos de un solido isotropo de velocidades c_L (longitudinal radial) y c_T (cizalla)
(criba_solido.modulos).  Geometria como fondo_cascara.FondoCascara (T1 con rho~ + p~_t; parte de RG de nu con w).
"""
import math
import os
import sys
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2g"), os.path.join(AQUI, "..", "fase2h"), AQUI):
    if p_ not in sys.path:
        sys.path.insert(0, p_)
import fondo_cascara as FC         # noqa: E402
import criba_solido as CS          # noqa: E402


class FondoSolido(FC.FondoCascara):
    def __init__(self, la=1.0, lchi=1.0, A=0.03, Ns=7.0, cL2=1.0, cT2=0.2, f_m=1.0, bphi=None, **kw):
        """bphi: anchura del perfil del escalar phi' = phi1 bphi^2/(l^2 + bphi^2) en la cola canonica de nu (por defecto
        la misma b de los perfiles de G4 y Gcal)."""
        self.la, self.lchi = la, lchi
        self.bphi = bphi
        self.As, self.Nss = A, Ns
        self.Yss = 1.0 / (2 * Ns ** 2)
        self.cL2, self.cT2 = cL2, cT2
        kw = dict(kw)
        kw["f_fl"] = f_m
        kw["A"] = 0.0
        FC.FN.FondoNu.__init__(self, **kw)

    def _materia_N(self, l, N, r, rp_, Np):
        """w, p~_r, p~_t, w' (escalares) con el lapso sin normalizar N."""
        if l < self.la:
            return 0.0, 0.0, 0.0, 0.0
        Y = 1.0 / (2 * N * N)
        dY = Y - self.Yss
        if dY <= 0:
            return 0.0, 0.0, 0.0, 0.0
        w = 2 * self.As * Y
        e = math.exp(-(l - self.la) / self.lchi)
        chi, chip = 1 - e, e / self.lchi
        Pi = self.As * dY
        pr = Pi * chi
        Phi1 = Np / N
        pt = pr + r / (2 * max(rp_, 1e-7)) * (w * Phi1 * (1 - chi) + Pi * chip)
        return w, pr, pt, -2 * w * Phi1

    def _rhs(self, l, y):
        N, Np, r, rp_ = y
        G4, Gc, Gcp = self._perfil_escalar(l)
        w, pr, pt, _ = self._materia_N(l, N, r, rp_, Np)
        bp = self.b if self.bphi is None else self.bphi
        ph = self.phi1 * bp * bp / (l * l + bp * bp)
        nu = (-self.nu0_actual * math.exp(-(l / self.wt) ** 2) + self.f_fl * r * w / (4 * G4)
              + 0.5 * self.kappa * r * ph * ph)
        rpt = w + pt - pr
        Npp = (Gc * N * (rp_ * rp_ - r * nu) - N * G4 + 0.5 * N * r * r * rpt + Gcp * (N * r * rp_ - r * r * Np)) \
            / (r * r * Gc)
        return [Np, Npp, rp_, rp_ * Np / N - nu]

    def _post(self):
        FC.FN.FondoNu._post(self)
        if not self.ok:
            return
        # superficie exterior: N (sin normalizar) = Ns, con N > Ns en todo el exterior
        s = self.sol
        L = s.t[-1]
        lt = np.linspace(self.la, min(L, 400.0), 40001)
        Nu = s.sol(lt)[0]
        self.reentra = False
        if Nu[0] >= self.Nss:
            self.ok, self.motivo = False, "sin estrella: N(la) >= Ns"
            return
        k = np.nonzero(Nu >= self.Nss)[0]
        if not k.size:
            self.ok, self.motivo = False, "la estrella no termina antes de l = 400"
            return
        from scipy.optimize import brentq
        self.ls = self.lb = brentq(lambda x: s.sol(x)[0] - self.Nss, lt[k[0] - 1], lt[k[0]], xtol=1e-13)
        self.reentra = bool(np.any(Nu[k[0]:] < self.Nss * (1 - 1e-9)))
        f = self.campos(np.array([self.la]))
        if self.la <= 0 or f["r1"][0] <= 1e-6:
            self.ok, self.motivo = False, "nucleo sin r' > 0 en la superficie interior"

    def perfil(self, l):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        Y = self.sol.sol(l)
        der = np.array([self._rhs(li, yi) for li, yi in zip(l, Y.T)])
        Nu, Np, r, r1 = Y
        Phi1 = Np / Nu
        Phi2 = der[:, 1] / Nu - Phi1 ** 2
        G4, H, Gc, Gcp = self.escalar(l)
        mat = np.array([self._materia_N(li, Ni, ri, r1i, Npi) for li, Ni, ri, r1i, Npi in zip(l, Nu, r, r1, Np)])
        w, pr, pt, wp = mat.T
        rho = w - pr
        out = dict(N=Nu / self._Ninf, Phi1=Phi1, Phi2=Phi2, r=r, r1=r1, r2=der[:, 3], Gc=Gc, Gcp=Gcp, G4=G4,
                   sig=r1 / r - Phi1, w=w, wp=wp, rho=rho, pr=pr, pt=pt)
        out.update(CS.modulos(rho, pr, pt, self.cL2, self.cT2))
        return out

    def fisica(self, n=400):
        l = np.linspace(self.la, self.lb * (1 - 1e-9), n)
        p = self.perfil(l)
        return dict(dec_t=float(np.min(1 - np.abs(p["pt"]) / p["rho"])), ev_min=float(np.min(p["Ev"] / p["w"])),
                    masa=float(np.trapezoid(p["r"] ** 2 * p["rho"], l) / 4), pt_max=float(np.max(p["pt"] / p["rho"])))
