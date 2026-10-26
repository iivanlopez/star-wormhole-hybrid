# -*- coding: utf-8 -*-
"""
Fase 2h: fondos del hibrido definidos por la combinacion de la condicion de energia nula nu = Phi' r' - r''.

Motivo (fase2h/ESTADO.md y diseno_U.py): con el lapso analitico de la Fase 2g la superficie de la estrella tiene una
presion radial efectiva negativa (G^l_l ~ -2e-3) y la planitud asintotica exige una violacion de la condicion de energia
nula extendida fuera de la estrella (y ~ -10^3), que hace U < 0 en una region ancha donde el presupuesto de K no llega.
Aqui la violacion se confina a la garganta por construccion:
    nu(l) = -nu0 e^{-(l/wt)^2}  +  f_fl r (rho~+p~)/(4 G4)  +  (kappa/2) r phi'^2,
(joroba que abre la garganta, parte "de RG" de la materia y cola canonica del escalar), y N, r salen de T1 con materia:
    N'' = [Gcal N (r'^2 - r nu) - N G4 + N r^2 (rho~+p~)/2 + Gcal' (N r r' - r^2 N')]/(r^2 Gcal),
    r'' = r' N'/N - nu,
desde la garganta (r = r0 = 1, r' = N' = 0, N = 1).  El modo del lapso que crece linealmente lejos (N ~ A + B l) se
anula ajustando por biseccion la amplitud kappa de la cola canonica (que cumple la condicion de energia nula: mas kappa,
N'' mas negativa, B menor; dependencia suave, no exponencial como la de nu0).  Si hasta kappa = 0 el modo es negativo,
la planitud exigiria violar la condicion de energia nula fuera de la garganta: fondo descartado.  El fluido lee
Y = 1/(2 N^2) con este lapso sin normalizar (la normalizacion de chi se absorbe en P): la estrella es N < Ns con Ns > 1.
Las magnitudes de salida usan el lapso normalizado N/N_inf (el de la masa ADM).
Perfiles de la teoria (G4 = S G4b, Gcal = S (Hb + delta)) y fluido P(Y) como en criba.FondoFluido, salvo que la
transicion del acoplo se fija en l_T absoluto.  La clase expone la interfaz de criba.FondoFluido.
"""
import os
import sys
import math
import numpy as np
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "fase2g"))
import criba as C          # noqa: E402


class FondoNu(C.FondoFluido):
    def __init__(self, wt=1.0, nu0=1.0, f_fl=1.0, kappa=None, phi1=1.0, b=3.0, n=2, d0=0.5, wd=1.0, lam=1.0, lT=50.0,
                 A=0.0, Ns=0.8, m=1.0, eos="pot", ks=0.0, r0=1.0, Lmax=4000.0, kappa_max=1e4, tol=1e-15):
        # atributos que usan los metodos heredados (escalar, fluido, criba)
        self.M, self.a = 0.0, 1.0
        self.b, self.n, self.d0, self.wd = b, n, d0, wd
        self.lam, self.tauT, self.lT = lam, None, lT
        self.eos, self.ks, self.cola = eos, ks, "nu"
        self.A, self.Ns, self.m, self.r0 = A, Ns, m, r0
        self.Ys = 1.0 / (2 * Ns ** 2)
        self.wt, self.f_fl, self.kappa, self.phi1 = wt, f_fl, kappa, phi1
        self.Lmax = Lmax
        self.ok, self.motivo = True, ""
        self._Ninf = 1.0
        self.nu0 = self.nu0_actual = nu0
        if kappa is None:
            self.kappa = self._ajustar_kappa(kappa_max, tol=tol)
            if not self.ok:
                return
        s = self._resolver(self.nu0, hasta=Lmax)
        if s.status != 0:
            self.ok, self.motivo = False, "la solucion no llega a Lmax"
            return
        self.sol = s
        # planitud asintotica de la solucion ajustada: el cambio de clase puede darse entre colapso y explosion sin
        # pasar por una solucion plana (p. ej. N decreciente y r' < 1 lejos)
        Lf = s.t[-1]
        N, Np, r, rp = s.sol(Lf)
        MK, mMS = Lf * Lf * Np / N, r * (1 - rp * rp) / 2
        if not (0.9 < rp < 1.0 + 1e-6 and Np >= 0 and MK > 0 and abs(MK - mMS) < 0.05 * max(1.0, MK)):
            self.ok, self.motivo = False, f"no asintoticamente plano (r'(L) = {rp:.4f}, M_K = {MK:.4g}, m_MS = {mMS:.4g})"
            return
        self._Ninf = self._N_inf(s)
        if not self._Ninf > 0:
            self.ok, self.motivo = False, "N_inf no positivo"
            return
        self._post()

    # -------------------------------------------------------------- nu y sistema
    def phip(self, l):
        b2 = self.b ** 2
        return self.phi1 * b2 / (np.asarray(l) ** 2 + b2)

    def _rp(self, N):
        """rho~ + p~ con el lapso sin normalizar (N(0) = 1)."""
        return C.FondoFluido.fluido(self, np.atleast_1d(N))[0]

    def nu(self, l, N, r):
        G4, _, _, _ = C.FondoFluido.escalar(self, np.atleast_1d(l))
        rp = self._rp(N)
        return (-self.nu0_actual * np.exp(-(np.asarray(l) / self.wt) ** 2) + self.f_fl * r * rp / (4 * G4)
                + 0.5 * self.kappa * r * self.phip(l) ** 2)

    def _perfil_escalar(self, l):
        """G4, Gcal, Gcal' en un punto (aritmetica escalar; igual que criba.FondoFluido.escalar)."""
        b2 = self.b * self.b
        q = l * l + b2
        s = b2 / q
        x = s * s
        xp = -4 * l * b2 * b2 / q ** 3
        n = self.n
        xn1 = x ** (n - 1)
        G4b = 1 + x * xn1 / (2 * n - 1)
        Hb = 1 - x * xn1
        Hbp = -n * xn1 * xp
        dl = self.d0 * math.exp(-l * l / self.wd ** 2)
        dlp = -2 * l / self.wd ** 2 * dl
        lam = self.lam
        e = math.exp(-(l / self.lT) ** 4)
        S = lam + (1 - lam) * (1 - e)
        Sp = (1 - lam) * e * 4 * l ** 3 / self.lT ** 4
        return S * G4b, S * (Hb + dl), Sp * (Hb + dl) + S * (Hbp + dlp)

    def _rp_escalar(self, N):
        Y = 1.0 / (2 * N * N)
        dY = Y - self.Ys
        if dY <= 0 or self.A == 0:
            return 0.0
        if self.eos == "exp":
            ex = math.exp(-self.ks / dY)
            return 2 * Y * self.A * ex * (1 + self.ks / dY)
        if self.m == 1:
            return 2 * self.A * Y
        return 2 * self.m * self.A * Y * dY ** (self.m - 1)

    def _rhs(self, l, y):
        N, Np, r, rp_ = y
        G4, Gc, Gcp = self._perfil_escalar(l)
        rho_p = self._rp_escalar(N)
        ph = self.phi1 * self.b * self.b / (l * l + self.b * self.b)
        nu = (-self.nu0_actual * math.exp(-(l / self.wt) ** 2) + self.f_fl * r * rho_p / (4 * G4)
              + 0.5 * self.kappa * r * ph * ph)
        Npp = (Gc * N * (rp_ * rp_ - r * nu) - N * G4 + 0.5 * N * r * r * rho_p + Gcp * (N * r * rp_ - r * r * Np)) \
            / (r * r * Gc)
        return [Np, Npp, rp_, rp_ * Np / N - nu]

    def _resolver(self, nu0, hasta):
        self.nu0_actual = nu0

        def malo(l, y):
            return min(y[0], y[2]) - 1e-9
        malo.terminal = True

        def explota(l, y):
            return 1e8 - abs(y[0])
        explota.terminal = True
        return solve_ivp(self._rhs, (0.0, hasta), [1.0, 0.0, self.r0, 0.0], method="DOP853", rtol=1e-11, atol=1e-13,
                         dense_output=True, events=(malo, explota))

    def _clase(self, nu0):
        """-1 si N o r se anulan (modo B < 0 o sin garganta), +1 si N crece sin limite o el modo lineal es positivo."""
        s = self._resolver(nu0, self.Lmax)
        return self._clasificar(s)

    def _clase_kappa(self, k):
        self.kappa = k
        return self._clasificar(self._resolver(self.nu0, self.Lmax))

    def _ajustar_kappa(self, kmax, tol=1e-12):
        a, b = 0.0, kmax
        ca, cb = self._clase_kappa(a), self._clase_kappa(b)
        if ca == cb:
            self.ok, self.motivo = False, f"kappa sin cambio de clase en [0, {kmax}] (clase {ca})"
            return 0.0
        for _ in range(200):
            m = 0.5 * (a + b)
            if not a < m < b:                                 # resolucion de la maquina
                break
            cm = self._clase_kappa(m)
            if cm == ca:
                a = m
            else:
                b = m
            if b - a < tol * max(1e-3, m):
                break
        k = a if ca == 1 else b
        self.kappa = k
        return k

    def _clasificar(self, s):
        """-1: colapso (N o r -> 0, o r' < 0 al final: modo B < 0); +1: N crece sin limite con r' > 0 (modo B > 0).
        Si la solucion llega a Lmax, signo de M_Komar - m_MS (lo domina 2 (B/A) L^2 si B != 0)."""
        N, Np, r, rp = s.y[:, -1]
        L = s.t[-1]
        if s.status == 1 and len(s.t_events[0]):
            return -1
        if rp < 0:                                           # r decreciente lejos de la garganta: colapso
            return -1
        if s.status == 1:                                    # N explota con r' >= 0
            return +1
        if s.status != 0:                                    # fallo del integrador
            return -1 if (r < 1e-3 or N < 1e-6) else +1
        return 1 if L * L * Np / N - r * (1 - rp * rp) / 2 > 0 else -1

    def _ajustar_nu0(self, a, b, tol=1e-12):
        ca, cb = self._clase(a), self._clase(b)
        if ca == cb:
            self.ok, self.motivo = False, f"nu0 sin cambio de clase en [{a}, {b}] ({ca})"
            return None
        for _ in range(200):
            m = 0.5 * (a + b)
            cm = self._clase(m)
            if cm == ca:
                a = m
            else:
                b = m
            if b - a < tol * max(1.0, abs(m)):
                break
        # se devuelve el extremo del lado +1 (modo lineal B >= 0 diminuto): la solucion llega a Lmax
        return a if ca == 1 else b

    @staticmethod
    def _masa(s):
        """Masa en el extremo: con un resto del modo lineal N ~ A + B l, la masa de Komar L^2 N'/N sube en (B/A) L^2 y la
        de Misner-Sharp r (1 - r'^2)/2 baja en lo mismo; la media cancela ese resto a primer orden."""
        L = s.t[-1]
        N, Np, r, rp = s.sol(L)
        return 0.5 * (L * L * Np / N + r * (1 - rp * rp) / 2), N, L

    def _N_inf(self, s):
        M, N, L = self._masa(s)
        if not (abs(M / L) < 50):
            return float("nan")
        return float(N * math.exp(M / L))

    def _post(self):
        s = self.sol
        L = s.t[-1]
        N, Np, r, rp = s.sol(L)
        self.M_adm = float(self._masa(s)[0])
        self.M_komar = float(L * L * Np / N)
        self.m_MS = float(r * (1 - rp * rp) / 2)
        self.N0 = 1.0 / self._Ninf
        # superficie de la estrella: N (sin normalizar) = Ns
        lt = np.linspace(0, min(L, 400.0), 40001)
        Nu = s.sol(lt)[0]
        self.reentra = False
        if self.A > 0 and Nu[0] < self.Ns:
            k = np.nonzero(Nu >= self.Ns)[0]
            if k.size:
                from scipy.optimize import brentq
                self.ls = brentq(lambda x: s.sol(x)[0] - self.Ns, lt[k[0] - 1], lt[k[0]], xtol=1e-13)
                # una sola estrella: el lapso no vuelve a bajar de Ns
                self.reentra = bool(np.any(Nu[k[0]:] < self.Ns * (1 - 1e-9)))
            else:
                self.ls = float(lt[-1])
                self.ok, self.motivo = False, "la estrella no termina antes de l = 400"
        else:
            self.ls = 0.0
        lt2 = np.concatenate([np.linspace(1e-4, 20, 2000), np.geomspace(20, L, 2000)[1:]])
        f = self.campos(lt2)
        if np.any(f["r1"] <= 0):
            self.ok, self.motivo = False, "r' <= 0 fuera de la garganta"

    def fluido_normalizado(self):
        """(A, Ys, ks, m) del mismo fluido escrito con el lapso normalizado N = N_u/N_inf (Y = 1/(2 N^2) = N_inf^2 Y_u):
        P = A (Y_u - Ys)^m = A N_inf^{-2m} (Y - Ys N_inf^2)^m;  eos = "exp": A -> A/N_inf^2, ks -> ks N_inf^2."""
        c2 = self._Ninf ** 2
        if self.eos == "exp":
            return dict(A=self.A / c2, Ys=self.Ys * c2, ks=self.ks * c2, m=1.0)
        return dict(A=self.A * self._Ninf ** (-2 * self.m), Ys=self.Ys * c2, ks=self.ks, m=self.m)

    def margen_exterior(self):
        """r sig = r' - r Phi' en la superficie: > 0 si la superficie queda fuera de la esfera de fotones (minimo de r/N).
        Fuera de la estrella I = N r Gcal (r' - r Phi') crece (I' = N G4), asi que sig > 0 en todo el exterior."""
        if self.ls <= 0:
            return 1.0
        f = self.campos(np.array([self.ls * (1 + 1e-9)]))
        return float(f["r1"][0] - f["r"][0] * f["Phi1"][0])

    @property
    def M_presupuesto(self):
        return self.M_adm

    # -------------------------------------------------------------- interfaz de criba.FondoFluido
    def lapso(self, l):
        escalar = np.ndim(l) == 0
        l = np.atleast_1d(np.asarray(l, dtype=float))
        N, Np, r, rp = self.sol.sol(l)
        out = (N / self._Ninf, Np / N)
        return (out[0][0], out[1][0]) if escalar else out

    def lapso3(self, l):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        N, Np, r, rp = self.sol.sol(l)
        Npp = np.array([self._rhs(li, yi)[1] for li, yi in zip(l, self.sol.sol(l).T)])
        return np.log(N / self._Ninf), Np / N, Npp / N - (Np / N) ** 2

    def escalar(self, l):
        return C.FondoFluido.escalar(self, l)

    def fluido(self, N):
        return C.FondoFluido.fluido(self, N)

    def campos(self, l):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        Nu, Np, r, r1 = self.sol.sol(l)
        N = Nu / self._Ninf
        Phi1 = Np / Nu
        G4, H, Gc, Gcp = self.escalar(l)
        rp, cs2, p = self.fluido(Nu)
        sig = r1 / r - Phi1
        a2 = 4 * N * N * r * Gc
        a2p = a2 * (2 * Phi1 + r1 / r + Gcp / Gc)
        I = N * r * Gc * (r1 - r * Phi1)
        return dict(l=l, N=N, Phi1=Phi1, G4=G4, H=H, Gc=Gc, rp=rp, cs2=cs2, p=p, r=r, r1=r1, sig=sig, a2=a2, a2p=a2p,
                    I=I, u=(r / N) ** 2)

    def E_out(self):
        if self.ls <= 0:
            return 0.0
        return C.FondoFluido.E_out(self)

    def presion_radial(self, l):
        """G^l_l = (r'^2 - 1)/r^2 + 2 r' Phi'/r (presion radial efectiva por 1/2 en RG)."""
        f = self.campos(l)
        return (f["r1"] ** 2 - 1) / f["r"] ** 2 + 2 * f["r1"] * f["Phi1"] / f["r"]
