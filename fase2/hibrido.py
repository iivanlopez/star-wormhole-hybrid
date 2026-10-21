# -*- coding: utf-8 -*-
"""
Fase 2 - Hibrido estrella-agujero de gusano: integrador en longitud propia.

Unidades geometricas G = c = 1 con longitudes en km.  Metrica (fase1, ec. 1.2):
    ds^2 = -e^{2 Phi(l)} dt^2 + dl^2 + r(l)^2 dOmega^2 ,  garganta en l = 0, simetria Z2.
Sector exotico en la normalizacion de DFKK (phi~ = sqrt(8 pi) phi, V~ = 8 pi V):
    8 pi T(fantasma) = {rho: -psi^2/2 + V~,  p_r: -psi^2/2 - V~,  p_t: +psi^2/2 - V~},  psi = dphi~/dl
Modelos de perturbacion:
    'fantasma' : escalar fantasma dinamico (DFKK)                         -> dof extra
    'ligado'   : Modelo L = sector de NOF (arXiv:2401.15868) con A = B = 0 -> sin dof propios;
                 mismo fondo, 16 pi lambda = -psi^2, unica perturbacion exotica: d(16 pi lambda).
"""
import math
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp

G_CGS, C_CGS, MSUN_G, KM = 6.67e-8, 3.0e10, 1.99e33, 1.0e5   # constantes de las tablas DFKK
MSUN_KM = G_CGS * MSUN_G / C_CGS**2 / KM                      # 1.474811 km
RHO_GEOM = G_CGS / C_CGS**2 * KM**2                           # (g cm^-3) -> km^-2
PI4, PI8 = 4.0 * math.pi, 8.0 * math.pi
RTOL, ATOL = 1e-11, 1e-14


class Politropo:
    """EoS de DFKK: P = K rho_b^(1+1/n), eps = rho_b + n P;  Gamma_1 = 1 + 1/n."""

    def __init__(self, k=0.1, gamma=2.0, nb_ch=1.0e38, mb=1.66e-24):
        self.Gam, self.n = gamma, 1.0 / (gamma - 1.0)
        K_cgs = k * C_CGS**2 * (nb_ch * mb) ** (1.0 - gamma)
        self.K = K_cgs * G_CGS ** (1 - gamma) * C_CGS ** (2 * gamma - 4) * KM ** (2 - 2 * gamma)

    def sigma(self, rho_bc):
        return self.K * rho_bc ** (1.0 / self.n)

    def inv_A(self, rho_bc):
        """Longitud 1/A de DFKK (R2, ec. 20) en km."""
        return math.sqrt(self.sigma(rho_bc) * (self.n + 1) / (PI4 * rho_bc))

    def termo(self, theta, rho_bc, sig):
        n = self.n
        thn = theta ** n
        eps = rho_bc * (1.0 + n * sig * theta) * thn
        P = sig * rho_bc * theta * thn
        deps = n * rho_bc * theta ** (n - 1.0) * (1.0 + (n + 1.0) * sig * theta)
        dP = (n + 1.0) * sig * rho_bc * thn
        return eps, P, deps, dP


@dataclass
class Fondo:
    """Configuracion estatica Z2. eos=None -> agujero de gusano puro."""
    r0: float
    eos: Politropo = None
    rho_bc: float = 0.0          # densidad barionica en la garganta [km^-2]
    V0: float = 0.0              # V~ = -V0 (1 - Lam^2 phi~^2)^2 ; V0 = 0 -> fantasma sin masa
    Lam: float = 1.0
    theta_stop: float = 1e-7     # corte de la superficie del fluido

    def __post_init__(self):
        self.fluido = self.eos is not None and self.rho_bc > 0
        if self.fluido:
            self.sig = self.eos.sigma(self.rho_bc)
            self.c1 = self.sig * (self.eos.n + 1.0)
            self.P0 = self.sig * self.rho_bc
        else:
            self.sig = self.c1 = self.P0 = 0.0
        V, _, _ = self.potencial(0.0)
        self.psi0 = math.sqrt(2.0 * (PI8 * self.P0 - V + 1.0 / self.r0**2))   # ligadura en l = 0

    # ---------------------------------------------------------------- fondo
    def potencial(self, f):
        L2 = self.Lam**2
        u = 1.0 - L2 * f * f
        return -self.V0 * u * u, 4.0 * self.V0 * L2 * f * u, 4.0 * self.V0 * L2 * (1.0 - 3.0 * L2 * f * f)

    def theta(self, Phi):
        return ((1.0 + self.c1) * math.exp(-Phi) - 1.0) / self.c1

    def bg(self, y, fluido):
        r, s, Phi, q, f, psi = y[0], y[1], y[2], y[3], y[4], y[5]
        V, dV, d2V = self.potencial(f)
        if fluido:
            th = self.theta(Phi)
            eps, P, de, dp = self.eos.termo(max(th, 0.0), self.rho_bc, self.sig)
            dth = -(1.0 + self.c1 * th) / self.c1 * q
            eps_l, P_l = de * dth, dp * dth
        else:
            eps = P = eps_l = P_l = 0.0
        X = 0.5 * psi * psi
        sl = (1.0 - s * s) / (2.0 * r) - 0.5 * r * (PI8 * eps - X + V)
        ql = PI8 * P + X - V - sl / r - q * q - s * q / r
        psil = -dV - (q + 2.0 * s / r) * psi
        return r, s, Phi, q, f, psi, sl, ql, psil, V, dV, d2V, eps, P, eps_l, P_l

    def rhs_bg(self, y, fluido):
        r, s, Phi, q, f, psi, sl, ql, psil, V, dV, d2V, eps, P, eps_l, P_l = self.bg(y, fluido)
        w = 0.5 * r * r * s
        return [s, sl, q, ql, psi, psil, w * PI8 * eps, w * (V - 0.5 * psi * psi)]

    def y0(self):
        return [self.r0, 0.0, 0.0, 0.0, 0.0, self.psi0, 0.0, 0.0]

    def ev_superficie(self):
        def ev(l, y, *a):
            return self.theta(y[2]) - self.theta_stop
        ev.terminal, ev.direction = True, -1
        return ev

    def ev_campo(self):
        """Fantasma cuartico: termina al sobrepasar 1/Lam o al retroceder (psi < 0)."""
        def sobre(l, y, *a):
            return y[4] * self.Lam - 1.0
        def bajo(l, y, *a):
            return y[5]
        sobre.terminal = bajo.terminal = True
        sobre.direction, bajo.direction = 1, -1
        return [sobre, bajo]

    def resolver(self, l_max=2.0e4, conv_tol=1e-7):
        """Integra el fondo por tramos (fluido | exterior con campo | exterior congelado)."""
        op = dict(method="DOP853", rtol=RTOL, atol=ATOL, dense_output=True)
        f_in = lambda l, y: self.rhs_bg(y, True)
        f_out = lambda l, y: self.rhs_bg(y, False)
        self.tramos, y, l = [], self.y0(), 0.0
        if self.fluido:
            sol = solve_ivp(f_in, (0.0, 1.0e4), y, events=[self.ev_superficie()], **op)
            if sol.t_events[0].size == 0:
                raise RuntimeError("no se alcanzo la superficie del fluido")
            self.ls = sol.t[-1]
            self.tramos.append((0.0, self.ls, sol.sol, True))
            y, l = sol.y[:, -1].copy(), self.ls
            self.Rs, self.Mfl, self.Msf_int = y[0], y[6], y[7]
        else:
            self.ls, self.Rs, self.Mfl, self.Msf_int = 0.0, self.r0, 0.0, 0.0
        self.l_conv = l
        if self.V0 > 0.0:                       # campo cuartico: hasta llegar al maximo 1/Lam
            def conv(l, y):
                return (1.0 / self.Lam - y[4]) - conv_tol
            conv.terminal, conv.direction = True, -1
            sol = solve_ivp(f_out, (l, l_max), y, events=[conv] + self.ev_campo(), **op)
            if sol.t_events[0].size == 0:
                raise RuntimeError("el campo no converge a 1/Lam (fondo mal ajustado)")
            self.tramos.append((l, sol.t[-1], sol.sol, False))
            y, l = sol.y[:, -1].copy(), sol.t[-1]
            y[4], y[5] = 1.0 / self.Lam, 0.0    # congelado en el maximo (V = V' = 0)
            self.l_conv = l
        sol = solve_ivp(f_out, (l, l_max), y, **op)
        self.tramos.append((l, l_max, sol.sol, False))
        self.l_max = l_max
        yf = sol.y[:, -1]
        r, s, Phi = yf[0], yf[1], yf[2]
        m = 0.5 * r * (1.0 - s * s)
        cola = -0.25 * yf[5] ** 2 * r**3          # masa del fantasma mas alla de l_max (~ -D^2/4r)
        self.M = m + cola
        self.Phi_inf = Phi - 0.5 * math.log(1.0 - 2.0 * m / r)
        self.Msf_ext = yf[7] - self.Msf_int + cola
        self.M_suma = 0.5 * self.r0 + self.Mfl + self.Msf_int + self.Msf_ext
        self.m2_inf = 8.0 * self.V0 * self.Lam**2  # masa^2 de las perturbaciones del fantasma en el infinito
        return self

    def y_bg(self, l):
        for a, b, sol, _ in self.tramos:
            if l <= b:
                return sol(l)
        return self.tramos[-1][2](l)

    def resumen(self):
        """Magnitudes en unidades de las tablas DFKK (km y masas solares)."""
        return dict(R_th=self.r0, R=self.Rs, R_prop=self.ls, M=self.M / MSUN_KM,
                    M_th=0.5 * self.r0 / MSUN_KM, M_fl=self.Mfl / MSUN_KM,
                    M_sfint=self.Msf_int / MSUN_KM, M_sfext=self.Msf_ext / MSUN_KM,
                    r_g=2.0 * self.M, M_suma=self.M_suma / MSUN_KM, Phi_inf=self.Phi_inf)

    # ---------------------------------------------------------------- perturbaciones
    def bases(self, modelo, w):
        """Datos iniciales en la garganta. Par: cumplen E3(0). Impar (solo 'ligado'): xi(0) y Psi(0)
        libres, dr = DP = dLambda = 0 en l = 0 (E3(0) se cumple identicamente)."""
        r0, P0, p0 = self.r0, self.P0, self.psi0
        if getattr(self, "paridad", "par") == "impar":
            if modelo != "ligado":
                raise ValueError("modos impares solo implementados para el Modelo L (gauge unitario)")
            return ([1.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0 / r0, 0.0, 0.0])
        g = 2.0 / r0**2 * (1.0 + w * r0**2)       # (2 dr0/r0^3)(1 + w r0^2) con dr0 = r0
        if modelo == "fantasma":
            a = [0.0, 0.0, r0, 0.0, 0.0, -g / p0]
            b = [0.0, P0, 0.0, 0.0, 0.0, PI8 * P0 / p0]
        else:
            a = [0.0, 0.0, r0, 0.0, g, 0.0]
            b = [0.0, P0, 0.0, 0.0, -PI8 * P0, 0.0]
        return (a, b) if self.fluido else (a,)

    def rhs_p(self, b, p, w, modelo, fluido):
        r, s, Phi, q, f, psi, sl, ql, psil, V, dV, d2V, eps, P, eps_l, P_l = b
        xi, DP, dr, Ps, f1, f2 = p
        e2 = math.exp(-2.0 * Phi)
        if fluido:
            G1P, hP = self.eos.Gam * P, eps + P
            xil = -DP / G1P - 2.0 * (s * xi + dr) / r
            dPe = DP - P_l * xi                                   # Euleriana
            de = hP * DP / G1P - eps_l * xi
            DPl = w * hP * e2 * xi - q * DP - hP * (Ps + ql * xi - 2.0 * q / r * (s * xi + dr))
            F = PI4 * hP * xi
            Fl = PI4 * ((eps_l + P_l) * xi + hP * xil)
            Sf = PI8 * (de - dPe)                                 # 8 pi delta(rho + p_r - 2 p_t)
        else:
            xil = DPl = F = Fl = Sf = 0.0
        if modelo == "fantasma":
            F += 0.5 * psi * f1
            Fl += 0.5 * (psil * f1 + psi * f2)
        drl = q * dr + r * F                                      # E1 (ligadura de momento)
        drll = ql * dr + q * drl + s * F + r * Fl
        drr_l = drl / r - s * dr / (r * r)
        if modelo == "fantasma":                                  # escalar fantasma perturbado
            f1l = f2
            f2l = -w * e2 * f1 - (q + 2.0 * s / r) * f2 - psi * (Ps + 2.0 * drr_l) - d2V * f1
            Sx = -4.0 * psi * f2 + 2.0 * dV * f1
        else:                                                     # E5: conservacion del sector ligado
            f1l = -(q + 2.0 * s / r) * f1 + psi * psi * (Ps + 2.0 * drr_l)
            f2l, Sx = 0.0, f1
        Psl = -2.0 * q * Ps - 2.0 * drll / r + 2.0 * sl * dr / (r * r) - 0.5 * (Sf + Sx)   # lapso
        return [xil, DPl, drl, Psl, f1l, f2l]

    def condiciones(self, modelo, w, l_end=None, rtol=1e-10, devolver=False):
        """Integra las bases desde la garganta y devuelve (F1, F2):
        F1 = Delta P en la superficie;  F2 = coef. creciente del fantasma (o u = e^Phi r^2 d(16 pi lambda))."""
        bs = self.bases(modelo, w)
        nb = len(bs)
        Y = [x for v in bs for x in v]
        sc = [1.0, max(self.P0, 1e-30), self.r0, 1.0 / self.r0, self.psi0, self.psi0 / self.r0]
        atol = [1e-13 * c for c in sc] * nb

        def rhs(l, Y, fluido):
            b = self.bg(self.y_bg(l), fluido)
            out = []
            for k in range(nb):
                out += self.rhs_p(b, Y[6 * k:6 * k + 6], w, modelo, fluido)
            return out

        op = dict(method="DOP853", rtol=rtol, atol=atol, dense_output=devolver)
        sols, l, F1 = [], 0.0, [0.0] * nb
        if self.fluido:
            sol = solve_ivp(rhs, (0.0, self.ls), Y, args=(True,), **op)
            sols.append(sol)
            Y, l = sol.y[:, -1], self.ls
            F1 = [Y[6 * k + 1] for k in range(nb)]
        if l_end is None:
            if modelo == "fantasma":
                kap = math.sqrt(max(self.m2_inf - w * math.exp(-2.0 * self.Phi_inf), 1e-30))
                l_end = min(max(self.l_conv, self.ls) + 30.0 / kap, self.l_max)
            else:
                l_end = self.l_max
        sol = solve_ivp(rhs, (l, l_end), Y, args=(False,), **op)
        sols.append(sol)
        Yf = sol.y[:, -1]
        r, s, Phi = self.y_bg(l_end)[:3]
        F2 = []
        for k in range(nb):
            p = Yf[6 * k:6 * k + 6]
            if modelo == "fantasma":
                kap = math.sqrt(max(self.m2_inf - w * math.exp(-2.0 * Phi), 1e-30))
                F2.append((p[5] + (kap + s / r) * p[4]) * r * math.exp(-kap * l_end))
            else:
                F2.append(math.exp(Phi) * r * r * p[4])
        return (F1, F2, sols) if devolver else (F1, F2)

    def det_simple(self, modelo, w, **kw):
        """Determinante por superposicion directa (valido si no hay crecimiento evanescente fuerte)."""
        F1, F2 = self.condiciones(modelo, w, **kw)
        if len(F2) == 1:
            return F2[0]
        return F1[0] * F2[1] - F1[1] * F2[0]

    def matriz_A(self, b, w, modelo, fluido):
        """Matriz del sistema lineal p' = A p, escrita a partir de rhs_p (p = xi, DP, dr, Psi, f1, f2)."""
        r, s, Phi, q, f, psi, sl, ql, psil, V, dV, d2V, eps, P, eps_l, P_l = b
        e2 = math.exp(-2.0 * Phi)
        A = np.zeros((6, 6), dtype=complex if isinstance(w, complex) else float)
        Fc, Flc, Sc = np.zeros(6), np.zeros(6), np.zeros(6)
        if fluido:
            G1P, hP = self.eos.Gam * P, eps + P
            A[0] = [-2.0 * s / r, -1.0 / G1P, -2.0 / r, 0.0, 0.0, 0.0]
            A[1] = [w * hP * e2 - hP * ql + 2.0 * hP * q * s / r, -q, 2.0 * hP * q / r, -hP, 0.0, 0.0]
            Fc[0] = PI4 * hP
            Flc[:] = PI4 * hP * A[0].real
            Flc[0] += PI4 * (eps_l + P_l)
            Sc[0], Sc[1] = PI8 * (P_l - eps_l), PI8 * (hP / G1P - 1.0)
        if modelo == "fantasma":
            Fc[4] += 0.5 * psi
            Flc[4] += 0.5 * psil
            Flc[5] += 0.5 * psi
        A[2] = r * Fc
        A[2, 2] += q
        drll = q * A[2] + s * Fc + r * Flc
        drll[2] += ql
        drr = A[2] / r
        drr[2] -= s / (r * r)
        if modelo == "fantasma":
            A[4, 5] = 1.0
            A[5] = -2.0 * psi * drr
            A[5, 3] -= psi
            A[5, 4] += -w * e2 - d2V
            A[5, 5] += -(q + 2.0 * s / r)
            Sc[4] += 2.0 * dV
            Sc[5] += -4.0 * psi
        else:
            A[4] = 2.0 * psi * psi * drr
            A[4, 3] += psi * psi
            A[4, 4] += -(q + 2.0 * s / r)
            Sc[4] += 1.0
        A[3] = -2.0 / r * drll - 0.5 * Sc
        A[3, 3] += -2.0 * q
        A[3, 2] += 2.0 * sl / (r * r)
        return A

    def det(self, modelo, w, l_end=None, rtol=1e-10, infinito="libre", escala=False):
        """Determinante por el metodo de la matriz compuesta: se evoluciona la 2-forma de Plucker
        Y = a b^T - b a^T con Y' = A Y + Y A^T, renormalizando (factor > 0) para evitar desbordes.
        D = sum_j c_j Y[1, j](l_end), con c el funcional F2 (Delta P queda congelado fuera del fluido).
        Modelo L: infinito='libre' impone u_inf = 0 (C libre); 'sujeto' impone C = e^-Phi dr = 0.
        El valor devuelto es D_true(w) / s(w), con s > 0 el producto de las renormalizaciones: conserva
        ceros y fase, pero no es analitico.  Con escala=True devuelve (D, log s); D e^{log s} es entera en w."""
        if not self.fluido:
            return self.det_simple(modelo, w, l_end=l_end, rtol=rtol)
        cplx = isinstance(w, complex)
        a, b = (np.array(v, dtype=complex if cplx else float) for v in self.bases(modelo, w))
        lig = modelo == "ligado"
        if lig:                                          # variable reescalada u = e^Phi r^2 dLambda (O(1) en el infinito)
            a[4] *= self.r0**2
            b[4] *= self.r0**2
        P0 = self.P0                                     # Delta P en unidades de P0 (escala O(1))
        a[1] /= P0
        b[1] /= P0
        iu = np.triu_indices(6, 1)                      # 15 coordenadas de Plucker independientes
        M0 = np.outer(a, b) - np.outer(b, a)
        log_s = [math.log(np.max(np.abs(M0[iu])))]      # log del producto de renormalizaciones
        Y = M0[iu] / np.max(np.abs(M0[iu]))

        def rhs(l, y, fluido):
            bg = self.bg(self.y_bg(l), fluido)
            A = self.matriz_A(bg, w, modelo, fluido)
            if lig:
                r, s, Phi, q = bg[0], bg[1], bg[2], bg[3]
                S4 = math.exp(Phi) * r * r
                A[:, 4] /= S4
                A[4, :] *= S4
                A[4, 4] += q + 2.0 * s / r
            A[:, 1] *= P0
            A[1, :] /= P0
            M = np.zeros((6, 6), dtype=y.dtype)
            M[iu] = y
            M = M - M.T
            return (A @ M + M @ A.T)[iu]

        def grande(l, y, fluido):
            return np.max(np.abs(y)) - 1e60
        grande.terminal = True

        def tramo(l, l1, Y, fluido):
            while l < l1 * (1 - 1e-14):
                sol = solve_ivp(rhs, (l, l1), Y, args=(fluido,), events=[grande], method="DOP853",
                                rtol=rtol, atol=1e-13 * np.max(np.abs(Y)))
                Y, l = sol.y[:, -1], sol.t[-1]
                m = np.max(np.abs(Y))
                log_s[0] += math.log(m)
                Y = Y / m
            return Y

        Y = tramo(0.0, self.ls, Y, True)
        if l_end is None:
            if modelo == "fantasma":
                kap = math.sqrt(max(self.m2_inf - (w.real if isinstance(w, complex) else w)
                                    * math.exp(-2.0 * self.Phi_inf), 1e-30))
                l_end = min(max(self.l_conv, self.ls) + 30.0 / kap, self.l_max)
            else:
                l_end = self.l_max
        y = tramo(self.ls, l_end, Y, False)
        Y = np.zeros((6, 6), dtype=y.dtype)
        Y[iu] = y
        Y = Y - Y.T
        r, s, Phi = self.y_bg(l_end)[:3]
        if modelo == "fantasma":
            kap = np.sqrt(self.m2_inf - w * math.exp(-2.0 * Phi) + 0j)
            kap = kap.real if not isinstance(w, complex) else kap
            D = Y[1, 5] + (kap + s / r) * Y[1, 4]
        elif infinito == "sujeto":
            D = Y[1, 2]
        else:
            # u(l) = u_inf + w C D~^2 e^{-2Phi}/r + O(r^-2), con C = dr e^{-Phi}: se resta la cola 1/r
            psi = self.y_bg(l_end)[5]
            D = Y[1, 4] - w * psi * psi * r**3 * math.exp(-Phi) * Y[1, 2]
        return (D, log_s[0]) if escala else D

    def det_gs(self, modelo, w, l_end=None, rtol=1e-10):
        """Comprobacion independiente: ortonormalizacion continua de Godunov (QR con diagonal > 0).
        sign(det) se conserva porque det(R) > 0 en cada paso. Variables sin reescalar."""
        a, b = (np.array(v, dtype=float) for v in self.bases(modelo, w))
        Q = np.stack([a, b], axis=1)
        Q = Q / np.linalg.norm(Q, axis=0)

        def rhs(l, y, fluido):
            A = self.matriz_A(self.bg(self.y_bg(l), fluido), w, modelo, fluido)
            return (A @ y.reshape(6, 2)).ravel()

        def grande(l, y, fluido):
            return np.max(np.abs(y)) - 1e8
        grande.terminal = True

        def qr_pos(M):
            q, r = np.linalg.qr(M)
            s = np.sign(np.diag(r))
            return q * s

        def tramo(l, l1, Q, fluido):
            while l < l1 * (1 - 1e-14):
                sol = solve_ivp(rhs, (l, l1), Q.ravel(), args=(fluido,), events=[grande],
                                method="DOP853", rtol=rtol, atol=1e-12)
                Q, l = qr_pos(sol.y[:, -1].reshape(6, 2)), sol.t[-1]
            return Q

        Q = tramo(0.0, self.ls, Q, True)
        if l_end is None:
            l_end = self.l_max
        Q = tramo(self.ls, l_end, Q, False)       # Delta P (fila 1) queda congelado fuera del fluido
        F1 = Q[1]
        r, s, Phi, q, f, psi = self.y_bg(l_end)[:6]
        F2 = math.exp(Phi) * r * r * Q[4] - w * psi * psi * r**3 * math.exp(-Phi) * Q[2]
        return F1[0] * F2[1] - F1[1] * F2[0]

    def w_fisico(self, w):
        """omega^2 medido en el infinito (las ecuaciones usan Phi con Phi(0) = 0)."""
        return w * math.exp(-2.0 * self.Phi_inf)


def autovalores(fondo, modelo, w_grid, **kw):
    """Barrido de det(w) y refinamiento con brentq. Devuelve (raices en w interno, valores det)."""
    from scipy.optimize import brentq
    D = [fondo.det(modelo, w, **kw) for w in w_grid]
    raices = []
    for i in range(len(w_grid) - 1):
        if D[i] == 0.0 or D[i] * D[i + 1] < 0.0:
            raices.append(brentq(lambda x: fondo.det(modelo, x, **kw), w_grid[i], w_grid[i + 1],
                                 xtol=1e-14, rtol=1e-10, maxiter=200))
    return raices, D
