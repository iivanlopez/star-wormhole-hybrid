# -*- coding: utf-8 -*-
"""
Fase 2g, parte II: criba de condiciones necesarias de estabilidad para gargantas de GLPV cuartica llenas de fluido
(hibrido estrella-agujero), sin depender de las funciones libres de la teoria (familia c = 0, familia_c0.py).

Fondo (simetrico, unidades r0 = 1):
    Phi = -M/sqrt(l^2 + a^2)          (lapso N = e^Phi, minimo N0 = e^{-M/a} en la garganta, sin horizontes)
    x = s^2, s = b^2/(l^2 + b^2)      (perfil del escalar; G4 = 1 + x^n/(2n-1), H = 1 - x^n)
    Gcal = H + d0 e^{-l^2/wd^2}       (Gcal > 0: sector impar sin fantasmas)
    fluido P(Y) = A (Y - Ys)^m, Y = 1/(2 N^2):  estrella donde N < Ns; rho~ + p~ = 2 m A Y (Y - Ys)^{m-1},
                 c_s^2 = (Y - Ys)/((Y - Ys) + 2 Y (m - 1))      (m = 1: fluido rigido, c_s = 1)
    r de T1 con materia:  [N r Gcal (r' - r Phi')]' = N G4 - (1/2) N r^2 (rho~ + p~)   (I' = N G4 - N r^2 (rho~+p~)/2,
                 u' = 2 I/(N^3 Gcal), r^2 = N^2 u).
Las demas ecuaciones del fondo fijan G2, G2_X (reconstruccion); E_p se cumple por Noether.

Condiciones necesarias (si falla alguna, la garganta es inestable o tiene un fantasma, con cualesquiera G2, G3, F4_X,...):
  (E) la energia estatica de la familia h (d = N h, c = 0, fluido moviendose para cumplir la ligadura de momento),
          E[h] = int_0^ls [ c_s^2 Q^2/(2 N r^2 (rho~+p~)) - 2 N^3 Gcal (h' - sig h)^2 ] + E_out h(ls)^2,
          Q = (4 N^2 r Gcal h')' + 2 N^2 r (rho~+p~) h,  sig = r'/r - Phi',  E_out = -2 int_ls^inf N^3 Gcal sig^2,
      con h' (ls) = 0 (fuera h es constante: d = N h(ls)), debe ser >= 0 (paridades par e impar en la garganta).
      h = 1 es la prueba d = N, xi = 0:  E[1] = -2 int N^3 [Gcal sig^2 - (rho~+p~) c_s^2].
  (T) la energia cinetica de la familia con soporte compacto en la estrella,
          T[h] = int [ 8 N Gcal^2 h'^2/(rho~+p~) - 2 G4 N h^2 ],
      debe ser >= 0 (si no, el monopolo tiene un fantasma).
Discretizacion: elementos finitos de Hermite cubicos (C^1), autovalores generalizados frente a int h^2.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

GX, GW = np.polynomial.legendre.leggauss(8)


class FondoFluido:
    def __init__(self, M=2.0, a=2.0, b=4.0, n=2, d0=0.9, wd=1.4, A=1.0, Ns=0.6, m=1.5, r0=1.0, Lmax=4000.0, lam=1.0,
                 tauT=1.5, eos="pot", ks=0.0, cola="fase2f", sfac=1.0, Ax=0.0, wA=50.0):
        self.M, self.a, self.b, self.n, self.d0, self.wd = M, a, b, n, d0, wd
        self.lam, self.tauT = lam, tauT
        self.eos, self.ks = eos, ks
        self.cola = cola
        self.s_cola = sfac * max(M, a)
        self.Ax, self.wA = Ax, wA
        self.A, self.Ns, self.m, self.r0 = A, Ns, m, r0
        self.N0 = math.exp(-M / a)
        self.Ys = 1.0 / (2 * Ns ** 2)
        self.ok = True
        self.motivo = ""
        if not (self.N0 < Ns < 1.0):
            self.ok, self.motivo = False, "Ns fuera de (N0, 1)"
            return
        if cola == "fase2f":
            self.ls = math.sqrt((M / math.log(Ns)) ** 2 - a ** 2)
        else:
            from scipy.optimize import brentq
            Lb = 1.0
            while self.lapso(np.array([Lb]))[0][0] < Ns:
                Lb *= 2
            self.ls = brentq(lambda x: self.lapso(np.array([x]))[0][0] - Ns, 0.0, Lb, xtol=1e-13)
        self.lT = tauT * self.ls
        self.Lmax = Lmax
        self._integrar()

    # -------------------------------------------------------------- perfiles cerrados
    def lapso(self, l):
        """N = e^Phi y Phi'.  cola = "fase2f": Phi = -M/R (R = sqrt(l^2 + a^2));  cola = "schw": se anade
        -(M^2/(2 R^2)) ln[(1 + R^2/s^2)/(1 + a^2/s^2)], que reproduce el termino -M^2 ln l/l^2 del lapso de Schwarzschild
        en distancia propia (s >= M/2 garantiza Phi' >= 0; mismo N0 = e^{-M/a})."""
        Phi, Phi1, _ = self.lapso3(l)
        return np.exp(Phi), Phi1

    def lapso3(self, l):
        """Phi, Phi', Phi''."""
        l = np.asarray(l, dtype=float)
        M, a = self.M, self.a
        R = np.sqrt(l * l + a * a)
        Phi = -M / R
        dR = M / R ** 2                      # dPhi/dR
        d2R = -2 * M / R ** 3
        if self.cola == "schw":
            s2 = self.s_cola ** 2
            g = np.log((1 + R * R / s2) / (1 + a * a / s2))
            gp = 2 * R / (s2 + R * R)
            gpp = 2 * (s2 - R * R) / (s2 + R * R) ** 2
            Phi = Phi - M * M * g / (2 * R * R)
            dR = dR + M * M * g / R ** 3 - M * M * gp / (2 * R * R)
            d2R = d2R - 3 * M * M * g / R ** 4 + 2 * M * M * gp / R ** 3 - M * M * gpp / (2 * R * R)
            if self.Ax:
                # termino A_x/(R^2 + wA^2): ajusta el coeficiente de 1/l^2 del lapso (constante A) sin tocar el interior
                Q = R * R + self.wA ** 2
                Phi = Phi + self.Ax / Q
                dR = dR - 2 * self.Ax * R / Q ** 2
                d2R = d2R - 2 * self.Ax / Q ** 2 + 8 * self.Ax * R * R / Q ** 3
        # Phi' = (l/R) dPhi/dR, Phi'' = (a^2/R^3) dPhi/dR + (l/R)^2 d2Phi/dR2
        return Phi, l / R * dR, a * a / R ** 3 * dR + (l / R) ** 2 * d2R

    def escalar(self, l):
        """G4 = S G4b, Gcal = S (Hb + delta), H = S Hb - 2 x G4b S'/x' (H = G4 - 2X G4_X, cero doble en la garganta);
        S = lam + (1 - lam)(1 - e^{-(l/lT)^4}) reduce el acoplo G4 (y Gcal) en el nucleo (S'' (0) = 0)."""
        b2 = self.b ** 2
        s = b2 / (l * l + b2)
        x = s * s
        xp = -4 * l * b2 * b2 / (l * l + b2) ** 3
        n = self.n
        G4b = 1 + x ** n / (2 * n - 1)
        Hb = 1 - x ** n
        Hbp = -n * x ** (n - 1) * xp
        dl = self.d0 * np.exp(-l * l / self.wd ** 2)
        dlp = -2 * l / self.wd ** 2 * dl
        lam = self.lam
        e = np.exp(-(l / self.lT) ** 4)
        S = lam + (1 - lam) * (1 - e)
        Sp = (1 - lam) * e * 4 * l ** 3 / self.lT ** 4
        G4 = S * G4b
        Gc = S * (Hb + dl)
        Gcp = Sp * (Hb + dl) + S * (Hbp + dlp)
        with np.errstate(divide="ignore", invalid="ignore"):
            H = S * Hb - np.where(np.abs(xp) > 0, 2 * x * G4b * Sp / np.where(xp != 0, xp, 1.0), 0.0)
        return G4, H, Gc, Gcp

    def fluido(self, N):
        """rho~ + p~, c_s^2, p~ (cero fuera de la estrella).  eos = "pot": P = A (Y - Ys)^m;
        eos = "exp": P = A z e^{-ks/z}, z = Y - Ys (superficie C^infinito, c_s -> 1 en el interior si z >> ks)."""
        N = np.asarray(N, dtype=float)
        Y = 1.0 / (2 * N * N)
        dY = np.maximum(Y - self.Ys, 0.0)
        if self.eos == "exp":
            A, k = self.A, self.ks
            z = np.where(dY > 0, dY, 1.0)
            ex = np.where(dY > 0, np.exp(-k / z), 0.0)
            PY = A * ex * (1 + k / z)
            PYY = A * ex * k * k / z ** 3
            rp = 2 * Y * PY
            cs2 = np.where(dY > 0, PY / (PY + 2 * Y * PYY + 1e-300), 0.0)
            return np.where(dY > 0, rp, 0.0), cs2, np.where(dY > 0, A * z * ex, 0.0)
        m, A = self.m, self.A
        rp = 2 * m * A * Y * dY ** (m - 1) if m != 1 else 2 * A * Y * (dY > 0)
        cs2 = np.where(dY > 0, dY / (dY + 2 * Y * (m - 1) + 1e-300), 0.0) if m != 1 else np.where(dY > 0, 1.0, 0.0)
        p = A * dY ** m
        return np.where(dY > 0, rp, 0.0), cs2, np.where(dY > 0, p, 0.0)

    # -------------------------------------------------------------- T1 con materia
    def _integrar(self):
        def rhs(l, y):
            I, u = y
            N, _ = self.lapso(l)
            G4, _, Gc, _ = self.escalar(l)
            rp, _, _ = self.fluido(N)
            r2 = N * N * u
            return [N * G4 - 0.5 * N * r2 * rp, 2 * I / (N ** 3 * Gc)]

        def cero_u(l, y):
            return y[1]
        cero_u.terminal = True
        u0 = self.r0 ** 2 / self.N0 ** 2
        bordes = sorted({0.0, self.ls, self.Lmax})
        sols = []
        y = [0.0, u0]
        for x0, x1 in zip(bordes[:-1], bordes[1:]):
            s = solve_ivp(rhs, (x0, x1), y, method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True,
                          events=cero_u)
            if s.status == 1 or not s.success:
                self.ok, self.motivo = False, "r se anula (sin garganta)"
                return
            sols.append((x0, x1, s))
            y = s.y[:, -1]
        self.sols = sols
        # r' > 0 en l > 0 (una sola garganta) en una malla
        lt = np.concatenate([np.linspace(1e-4, 20, 4000), np.geomspace(20, self.Lmax, 2000)[1:]])
        f = self.campos(lt)
        if np.any(f["r1"] <= 0):
            self.ok, self.motivo = False, "r' <= 0 fuera de la garganta"
            return
        self.r_inf_lin = f["r"][-1] / lt[-1]

    def _y(self, l):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        out = np.empty((2, l.size))
        for x0, x1, s in self.sols:
            msk = (l >= x0) & (l <= x1)
            if np.any(msk):
                out[:, msk] = s.sol(l[msk])
        return out

    def campos(self, l):
        l = np.atleast_1d(np.asarray(l, dtype=float))
        I, u = self._y(l)
        N, Phi1 = self.lapso(l)
        G4, H, Gc, Gcp = self.escalar(l)
        rp, cs2, p = self.fluido(N)
        up = 2 * I / (N ** 3 * Gc)
        r = N * np.sqrt(u)
        r1 = r * (Phi1 + up / (2 * u))
        sig = r1 / r - Phi1
        a2 = 4 * N * N * r * Gc
        a2p = a2 * (2 * Phi1 + r1 / r + Gcp / Gc)
        return dict(l=l, N=N, Phi1=Phi1, G4=G4, H=H, Gc=Gc, rp=rp, cs2=cs2, p=p, r=r, r1=r1, sig=sig, a2=a2, a2p=a2p,
                    I=I, u=u)

    def A_lapso(self):
        """Coeficiente A del termino 1/l^2 del lapso lejos: Phi = -M/l - M^2 ln l/l^2 + A/l^2 + ..."""
        if self.cola != "schw":
            return None
        s2 = self.s_cola ** 2
        return self.M ** 2 * (0.5 * math.log(s2) + 0.5 * math.log(1 + self.a ** 2 / s2)) + self.Ax

    def C1(self, L=None):
        """C1 = lim (I - l + M ln l) (I = N r Gcal (r' - r Phi'), primer integral de T1) y el valor de Schwarzschild
        C1_S = A/M - 2M compatible con el termino A/l^2 del lapso (dC1 = C1 - C1_S controla y ~ -6 M dC1 ln l)."""
        L = self.Lmax if L is None else L
        I, u = self._y(np.array([L]))
        Phi, _, _ = self.lapso3(np.array([L]))
        # I = l - M ln l + C1 + (M^2 ln l + M^2/2 - A)/l + ...
        A = self.A_lapso()
        corr = (self.M ** 2 * math.log(L) + self.M ** 2 / 2 - (A if A is not None else 0.0)) / L
        C1 = float(I[0] - L + self.M * math.log(L) - corr)
        return C1, (A / self.M - 2 * self.M if A is not None else None)

    def K_lejos(self, l):
        """K natural lejos de la estrella (Gcal, H, G4, Gamma1 -> 1, sin materia): K = 2 r (Phi' r' - r'')/(N r'^2) y
        y = (Phi' r' - r'') l^3.  r'' se obtiene del ODE de T1 (u = r^2/N^2, I' = N G4)."""
        l = np.atleast_1d(np.asarray(l, dtype=float))
        I, u = self._y(l)
        Phi, Phi1, Phi2 = self.lapso3(l)
        N = np.exp(Phi)
        G4, H, Gc, Gcp = self.escalar(l)
        up = 2 * I / (N ** 3 * Gc)
        Ip = N * G4
        upp = 2 * (Ip * N ** 3 * Gc - I * (3 * N ** 3 * Phi1 * Gc + N ** 3 * Gcp)) / (N ** 3 * Gc) ** 2
        sq = np.sqrt(u)
        Np, Npp = N * Phi1, N * (Phi2 + Phi1 ** 2)
        r = N * sq
        r1 = Np * sq + N * up / (2 * sq)
        r2 = Npp * sq + Np * up / sq + N * upp / (2 * sq) - N * up ** 2 / (4 * u * sq)
        w = Phi1 * r1 - r2
        return dict(K=2 * r * w / (N * r1 ** 2), y=w * l ** 3, r=r, r1=r1, N=N)

    def E_out(self):
        """-2 int_ls^inf N^3 Gcal sig^2 (cola analitica ~ 2/l mas alla de Lmax)."""
        l = np.concatenate([np.linspace(self.ls, self.ls + 50, 8001), np.geomspace(self.ls + 50, self.Lmax, 4001)[1:]])
        f = self.campos(l)
        g = f["N"] ** 3 * f["Gc"] * f["sig"] ** 2
        return -2 * (np.trapezoid(g, l) + g[-1] * l[-1])


# ------------------------------------------------------------------ elementos finitos de Hermite
def hermite(t, H):
    """Valores, primeras y segundas derivadas (en l) de las 4 funciones de Hermite en t in [0, 1]."""
    v = np.array([2 * t ** 3 - 3 * t ** 2 + 1, H * (t ** 3 - 2 * t ** 2 + t), -2 * t ** 3 + 3 * t ** 2, H * (t ** 3 - t ** 2)])
    d1 = np.array([6 * t ** 2 - 6 * t, H * (3 * t ** 2 - 4 * t + 1), -6 * t ** 2 + 6 * t, H * (3 * t ** 2 - 2 * t)]) / H
    d2 = np.array([12 * t - 6, H * (6 * t - 4), -12 * t + 6, H * (6 * t - 2)]) / H ** 2
    return v, d1, d2


def malla_estrella(ls, n=240):
    """Nodos uniformes en [0, ls] (elementos muy pequenos estropean el condicionamiento: h'' ~ 1/H^2)."""
    return np.linspace(0.0, ls, n + 1)


def formas(F, nodos, que=("E", "T", "M")):
    """Matrices globales (2 grados de libertad por nodo: h, h') de las formas E (sin E_out), T y la masa M; con "R" y
    "RT" las formas positivas de referencia (mismos terminos con signo +) para normalizar los autovalores a [-1, 1]."""
    K = len(nodos)
    nd = 2 * K
    out = {k: np.zeros((nd, nd)) for k in que}
    for k in ("R", "RT"):
        if k in que:
            out[k] = np.zeros((nd, nd))
    for e in range(K - 1):
        x0, x1 = nodos[e], nodos[e + 1]
        H = x1 - x0
        t = 0.5 * (GX + 1)
        w = 0.5 * GW * H
        l = x0 + t * H
        f = F.campos(l)
        v, d1, d2 = hermite(t, H)
        idx = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        if "E" in que:
            Q = f["a2"] * d2 + f["a2p"] * d1 + 2 * f["N"] ** 2 * f["r"] * f["rp"] * v
            AQ = np.where(f["rp"] > 0, f["cs2"] / (2 * f["N"] * f["r"] ** 2 * np.maximum(f["rp"], 1e-300)), 0.0)
            B = 2 * f["N"] ** 3 * f["Gc"]
            D = d1 - f["sig"] * v
            out["E"][np.ix_(idx, idx)] += (Q * AQ * w) @ Q.T - (D * B * w) @ D.T
            if "R" in que:
                vr = v / f["r"]
                out["R"][np.ix_(idx, idx)] += (D * B * w) @ D.T + (vr * B * w) @ vr.T
        if "T" in que:
            wt = np.where(f["rp"] > 0, 8 * f["N"] * f["Gc"] ** 2 / np.maximum(f["rp"], 1e-300), 1e300)
            out["T"][np.ix_(idx, idx)] += (d1 * wt * w) @ d1.T - (v * 2 * f["G4"] * f["N"] * w) @ v.T
            if "RT" in que:
                out["RT"][np.ix_(idx, idx)] += (d1 * wt * w) @ d1.T + (v * 2 * f["G4"] * f["N"] * w) @ v.T
        if "M" in que:
            out["M"][np.ix_(idx, idx)] += (v * w) @ v.T
    return out


def presupuesto_fantasma(F):
    """Condicion necesaria de ausencia de fantasmas con cualesquiera funciones libres.
    Energia cinetica del monopolo con fluido (velocidades h_t, y_t; d = N h, c = y/gam, W = al1 N h' + y):
        T = int [ (al1 N h_t' + y_t)^2/(4 kap_b) - 2 G4 N h_t^2 - (e_tt N/gam) y_t h_t ];
    minimizando en y_t punto a punto, T >= 0 exige N^2 K >= N (rho~+p~) xi^2/(2 Gcal^2) dentro de la estrella
    (xi = 4 r Gcal H/Gamma) y K >= 0 fuera.  Con Z = N xi y N^2 K = 2 [Z' - N G4]:  Z' >= N G4 + (rho~+p~) Z^2/(4 N Gcal^2),
    Z(0) = 0.  Con T1 con materia y W(inf) = N[xi - r Gcal (r' - r Phi')] -> 3M, la condicion es
        3 M > int_0^ls [ (rho~+p~) Z*^2/(4 N Gcal^2) + N r^2 (rho~+p~)/2 ] dl
    con Z* la solucion de la igualdad (si Z* diverge dentro de la estrella hay fantasma con cualquier eleccion)."""
    def rhs(l, y):
        f = F.campos(np.array([l]))
        return [f["N"][0] * f["G4"][0] + f["rp"][0] * y[0] ** 2 / (4 * f["N"][0] * f["Gc"][0] ** 2),
                f["rp"][0] * y[0] ** 2 / (4 * f["N"][0] * f["Gc"][0] ** 2) + 0.5 * f["N"][0] * f["r"][0] ** 2 * f["rp"][0]]

    def explota(l, y):
        return y[0] - 1e8
    explota.terminal = True
    s = solve_ivp(rhs, (0.0, F.ls), [0.0, 0.0], method="DOP853", rtol=1e-10, atol=1e-12, events=explota)
    Mp = getattr(F, "M_presupuesto", F.M)          # 3 N_inf M_ADM si el lapso no tiende a 1 (fase2h/fondo.py)
    if s.status == 1:
        return dict(Zs=math.inf, gasto=math.inf, presupuesto=3 * Mp, margen_fantasma=-math.inf)
    Zs, gasto = s.y[0, -1], s.y[1, -1]
    return dict(Zs=float(Zs), gasto=float(gasto), presupuesto=3 * Mp, margen_fantasma=float(3 * Mp - gasto))


def criba(F, n=240, devolver_modos=False):
    """Autovalor minimo de E (con E_out y h'(ls) = 0) y de T (soporte compacto) por paridad, y E[1]."""
    nodos = malla_estrella(F.ls, n)
    fm = formas(F, nodos)
    K = len(nodos)
    Eo = F.E_out()
    Em = fm["E"].copy()
    Em[2 * (K - 1), 2 * (K - 1)] += Eo
    res = dict(E_out=Eo)
    # h = 1: valores 1, derivadas 0
    uno = np.zeros(2 * K)
    uno[0::2] = 1.0
    res["E1"] = float(uno @ Em @ uno)
    for par, quita0 in (("par", [1]), ("impar", [0])):
        libres = [i for i in range(2 * K) if i not in quita0 and i != 2 * (K - 1) + 1]
        ev = eigh(Em[np.ix_(libres, libres)], fm["M"][np.ix_(libres, libres)], eigvals_only=True,
                  subset_by_index=[0, 1])
        res[f"Emin_{par}"] = float(ev[0])
        libT = [i for i in range(2 * K) if i not in quita0 and i < 2 * (K - 1)]
        evT = eigh(fm["T"][np.ix_(libT, libT)], fm["M"][np.ix_(libT, libT)], eigvals_only=True, subset_by_index=[0, 0])
        res[f"Tmin_{par}"] = float(evT[0])
    return res


def energia_superficie(F):
    """Coeficiente de h'(ls)^2 en la energia de superficie de un fluido con densidad finita en ls (m = 1; 0 si la
    densidad se anula en la superficie): N^3 Phi' (4 r Gcal)^2/(2 r^2 (rho~+p~)) evaluado justo dentro."""
    if F.eos == "exp" or F.m != 1:
        return 0.0
    f = F.campos(np.array([F.ls * (1 - 1e-9)]))
    N, Phi1, r, Gc, rp = (f[k][0] for k in ("N", "Phi1", "r", "Gc", "rp"))
    return float(N ** 3 * Phi1 * 16 * r * r * Gc * Gc / (2 * r * r * rp))


def criba_normalizada(F, n=120, libre=False):
    """Margenes adimensionales (positivos = la condicion se cumple):
    E1n = E[1]/(int 2 N^3 [Gcal sig^2 + (rho~+p~) c_s^2] + |E_out|);
    Epar, Eimpar = autovalor minimo de E frente a D[h] = int 2 N^3 Gcal [(h' - sig h)^2 + h^2/r^2] + |E_out| h(ls)^2;
    Tpar, Timpar = el de T frente a RT = int [8 N Gcal^2 h'^2/(rho~+p~) + 2 G4 N h^2];  fant = (3M - gasto)/(3M + gasto).
    libre = True (Fase 2h): h'(ls) queda libre en la prueba E (superficie del fluido movil, xi(ls) = -al1 N h'(ls)/(2 kap_b);
    con h'(ls) = 0 la superficie esta fija y la prueba es mas debil).  Si la densidad salta en la superficie (m = 1),
    P_YY = A delta(Y - Ys) da a kap_a un termino C_s delta(l - ls), C_s = N r^2 (rho~+p~)/(2 Phi'), y al eliminar el
    lapso una energia de superficie N^4 W(ls-)^2/(4 C_s) = N^3 Phi' W^2/(2 r^2 (rho~+p~)), W = al1 N h' = -4 r Gcal h'."""
    nodos = malla_estrella(F.ls, n)
    fm = formas(F, nodos, que=("E", "T", "M", "R", "RT"))
    K = len(nodos)
    Eo = F.E_out()
    Em, Rm = fm["E"].copy(), fm["R"].copy()
    Em[2 * (K - 1), 2 * (K - 1)] += Eo
    Rm[2 * (K - 1), 2 * (K - 1)] += abs(Eo)
    if libre:
        Em[2 * K - 1, 2 * K - 1] += energia_superficie(F)
    uno = np.zeros(2 * K)
    uno[0::2] = 1.0
    l = np.linspace(0, F.ls, 4001)
    f = F.campos(l)
    den1 = np.trapezoid(2 * f["N"] ** 3 * (f["Gc"] * f["sig"] ** 2 + f["rp"] * f["cs2"]), l) + abs(Eo)
    out = dict(E1n=float((uno @ Em @ uno) / den1))
    for par, quita0 in (("par", [1]), ("impar", [0])):
        lib = [i for i in range(2 * K) if i not in quita0 and (libre or i != 2 * (K - 1) + 1)]
        out[f"E{par}"] = float(eigh(Em[np.ix_(lib, lib)], Rm[np.ix_(lib, lib)], eigvals_only=True,
                                    subset_by_index=[0, 0])[0])
        libT = [i for i in range(2 * K) if i not in quita0 and i < 2 * (K - 1)]
        RT = fm["RT"][np.ix_(libT, libT)]
        RT = RT + 1e-12 * np.trace(RT) / len(libT) * np.eye(len(libT))
        out[f"T{par}"] = float(eigh(fm["T"][np.ix_(libT, libT)], RT, eigvals_only=True, subset_by_index=[0, 0])[0])
    pf = presupuesto_fantasma(F)
    Mp = getattr(F, "M_presupuesto", F.M)
    out["fant"] = float((3 * Mp - pf["gasto"]) / (3 * Mp + pf["gasto"])) if np.isfinite(pf["gasto"]) else -1.0
    out["margen"] = min(out.values())
    return out


def impar(F, LL=2, L=400.0, npts=6000):
    """Sector impar l = LL con fluido (impar_fluido.py: mismas formas que en el vacio, sin funciones libres):
    alfa = N_l Gcal/N, beta = N_l lam G4/(N r^2), gamma = -N_l lam N Gcal/r^2;  K_o = -alfa^2/gamma, P_o = alfa^2/beta,
    U_o = alfa + q^2/beta - (alfa q/beta)',  q = alfa' + 2 alfa r'/r.  Numero de taquiones (par, impar) por oscilacion
    de Sturm (espectro.SL de la Fase 2f) y velocidad radial maxima c^2 = Gcal/G4."""
    import sys as _s, os as _o
    _s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), "..", "fase2f"))
    import espectro as ES
    from scipy.interpolate import CubicSpline
    Nl = LL * (LL + 1) / (2 * LL + 1)
    lam = (LL - 1) * (LL + 2)
    ls = ES.malla(1e-3, L, npts)
    f = F.campos(ls)
    G4, H, Gc, Gcp = F.escalar(ls)
    N, Phi1, r, r1 = f["N"], f["Phi1"], f["r"], f["r1"]
    al = Nl * Gc / N
    alp = Nl * (Gcp - Gc * Phi1) / N
    be = Nl * lam * G4 / (N * r * r)
    ga = -Nl * lam * N * Gc / (r * r)
    q = alp + 2 * al * r1 / r
    g = al * q / be
    gp = CubicSpline(np.log(ls), g)(np.log(ls), 1) / ls
    K = -al * al / ga
    P = al * al / be
    U = al + q * q / be - gp
    out = dict(minK=float(K.min()), minP=float(P.min()), maxc2=float(np.max(Gc / G4)))
    if K.min() > 0 and P.min() > 0:
        sl = ES.SL(ls, K, P, U)
        out["n_par"] = sl.ceros("par")[0]
        out["n_impar"] = sl.ceros("impar")[0]
    return out


def prueba_sin_compresion(F, n=4000):
    """Prueba explicita de la familia c = 0 con la superficie libre (Fase 2h): y con Q = (4 N^2 r Gcal y')' +
    2 N^2 r (rho~+p~) y = 0 en la estrella (sin compresion del fluido, Q^2/(4 kap_a) = 0), y(0) = 1, y'(0) = 0, y
    constante fuera.  Su energia es
        E_Q0 = -2 int_0^ls N^3 Gcal (y' - sig y)^2 + E_out y(ls)^2 + E_s y'(ls)^2,
    con E_s el muelle de superficie (energia_superficie; 0 si la densidad se anula en la superficie).  Sin muelle es
    negativa: el fluido solo estabiliza a traves de la compresion, y con la superficie libre puede seguir el
    desplazamiento sin comprimirse.  Devuelve las partes (y la version normalizada por la energia de la prueba h = 1)."""
    from scipy.integrate import solve_ivp

    def rhs(l, z):
        f = F.campos(np.array([l]))
        a2 = 4 * f["N"][0] ** 2 * f["r"][0] * f["Gc"][0]
        # z = (y, a2 y')
        return [z[1] / a2, -2 * f["N"][0] ** 2 * f["r"][0] * f["rp"][0] * z[0]]
    s = solve_ivp(rhs, (0.0, F.ls * (1 - 1e-9)), [1.0, 0.0], method="DOP853", rtol=1e-10, atol=1e-13, dense_output=True)
    l = np.linspace(0.0, F.ls * (1 - 1e-9), n)
    y, q = s.sol(l)
    f = F.campos(l)
    yp = q / (4 * f["N"] ** 2 * f["r"] * f["Gc"])
    vac = -2 * np.trapezoid(f["N"] ** 3 * f["Gc"] * (yp - f["sig"] * y) ** 2, l)
    Eo = F.E_out()
    Es = energia_superficie(F)
    E = vac + Eo * y[-1] ** 2 + Es * yp[-1] ** 2
    E1 = 2 * np.trapezoid(f["N"] ** 3 * (f["rp"] * f["cs2"] - f["Gc"] * f["sig"] ** 2), l) + Eo
    return dict(E_Q0=float(E), vacio=float(vac), exterior=float(Eo * y[-1] ** 2), muelle=float(Es * yp[-1] ** 2),
                y_ls=float(y[-1]), yp_ls=float(yp[-1]), E1=float(E1))
