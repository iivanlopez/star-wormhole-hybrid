# -*- coding: utf-8 -*-
"""
Einstein-escalar-Gauss-Bonnet con potencial escalar, en las convenciones de fase2b (longitud propia l, N = e^Phi):

  S = (1/16 pi) int sqrt(-g) [ R - (dphi)^2/2 - V(phi) + F(phi) G ]
  ds^2 = -e^{2 Phi} dt^2 + dl^2 + r^2 dOmega^2

El potencial es parte del sector escalar: T_mn -> T_mn - g_mn V/2 y box phi + F' G -> box phi - V' + F' G.  En
componentes mixtas las ecuaciones de fase2b (egb_F_gen) reciben exactamente
  E^t_t, E^th_th, E^l_l:  + V/2          E_phi:  - V'(phi)          E^l_t: nada
y las perturbadas (phi -> phi + q):  + V' q/2 en tt, thth y ll;  - V'' q en la del escalar.  Es la regla que
comprueba verificar_potencial.py (Bianchi con V en SymPy, Schwarzschild-de Sitter, propagacion de ligaduras).

Garganta lisa (r' = 0 en l = 0, r0 = 1): la ligadura ll da
    V(phi0)/2 = 1 + phi1^2/4 - 4 F'(phi0) a phi1,        a = Phi'(0),  phi1 = phi'(0),
de modo que el potencial permite a = 0 (sin gradiente del lapso) con V(phi0) = 2 + phi1^2/2.

Contenido:
  potencial(tipo, amp, **kw)          phi -> (V, V', V'', V''')
  segundas_V, est_ligadura_V          fondo con potencial
  fondo_en_V, derivadas_V, ligaduras_V, velocidades_V, ECUACIONES_V
                                      perturbaciones radiales (gancho f.ecuaciones de fase2c/estabilidad_lisa)
  LadoV                               una mitad lisa desde la garganta, con clasificacion del final
  DobleV                              las dos mitades pegadas (l de cualquier signo)
"""
import os
import sys
import math
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(AQUI, "..", "fase2b"))
sys.path.append(os.path.join(AQUI, "..", "fase2c"))
import egb_F_gen as H                                          # noqa: E402
from fondo_doble import acoplamiento                           # noqa: E402,F401  (incluye el lineal)

PHIP_MAX = 1e7      # |Phi'| o |phi'| por encima de esto: final singular
PHI_MAX = 50.0      # |phi - phi_v| por encima de esto: el escalar se escapa del vacio
DET_MIN = 1e-4


# ------------------------------------------------------------------ potenciales
def potencial(tipo, amp=1.0, **kw):
    """Devuelve phi -> (V, V', V'', V''') con V = amp U(phi).
    constante   U = 1
    masivo      U = phi^2/2                       (vacio phi = 0, masa^2 = amp)
    doble_pozo  U = (phi^2 - v^2)^2               (vacios +-v, maximo local en 0; masa^2 = 8 amp v^2)
    cuartico    U = phi^2/2 + g phi^4/4           (g > 0)
    """
    if tipo == "constante":
        return lambda x: (amp, 0.0, 0.0, 0.0)
    if tipo == "masivo":
        return lambda x: (amp * x * x / 2, amp * x, amp, 0.0)
    if tipo == "doble_pozo":
        v2 = kw["v"] ** 2
        return lambda x: (amp * (x * x - v2) ** 2, amp * 4 * x * (x * x - v2), amp * (12 * x * x - 4 * v2),
                          amp * 24 * x)
    if tipo == "cuartico":
        g = kw["g"]
        return lambda x: (amp * (x * x / 2 + g * x**4 / 4), amp * (x + g * x**3), amp * (1 + 3 * g * x * x),
                          amp * 6 * g * x)
    raise ValueError(tipo)


def forma(tipo, **kw):
    """U(phi) normalizado (amp = 1)."""
    return potencial(tipo, 1.0, **kw)


def amp_garganta(Ud, Fd, phi0, phi1, a):
    """Amplitud del potencial que cumple la ligadura ll en la garganta (r0 = 1):
    amp U(phi0)/2 = 1 + phi1^2/4 - 4 F'(phi0) a phi1."""
    U0 = Ud(phi0)[0]
    return 2.0 * (1.0 + phi1**2 / 4 - 4.0 * Fd(phi0)[1] * a * phi1) / U0


# ------------------------------------------------------------------ fondo
def segundas_V(y, Fd, Vd):
    """(r'', Phi'', phi'') y det M del fondo con potencial (y = (r, r', Phi, Phi', phi, phi', ...))."""
    r, s, Phi, q, p, u = y[:6]
    F0, F1, F2, F3 = Fd(p)
    V0, V1 = Vd(p)[:2]
    M = np.array(H.est_principal(F1, r, s, q, u))
    b = np.array(H.est_resto(F1, F2, r, s, q, u))
    b[0] += V0 / 2                     # E^t_t + V/2
    b[1] += V0 / 2                     # E^th_th + V/2
    b[2] -= V1                         # box phi - V' + F' G
    return np.linalg.solve(M, -b), np.linalg.det(M)


def est_ligadura_V(y, Fd, Vd):
    r, s, Phi, q, p, u = y[:6]
    return H.est_ligadura(Fd(p)[1], r, s, q, u) + Vd(p)[0] / 2


# ------------------------------------------------------------------ perturbaciones radiales
def fondo_en_V(f, l):
    y = f.sol.sol(l)
    x, _ = segundas_V(y, f.Fd, f.Vd)
    F0, F1, F2, F3 = f.Fd(y[4])
    return (F1, F2, F3), (y[0], y[1], x[0], y[3], x[1], y[5], x[2], math.exp(-2.0 * y[2]))


def derivadas_V(f, l, y, w):
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = fondo_en_V(f, l)
    V0, V1, V2 = f.Vd(f.sol.sol(l)[4])[:3]
    a1, d0, d1, q0, q1 = y
    M = np.array(H.lin_principal(F1, r, r1, P1, p1))
    b = np.array(H.lin_resto(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1))
    b = b + np.array([V1 * q0 / 2, V1 * q0 / 2, -V2 * q0])
    d2, a2, q2 = np.linalg.solve(M, -b)
    return [a2, d1, d2, q1, q2]


def ligaduras_V(f, l, y, w):
    (F1, F2, F3), bg = fondo_en_V(f, l)
    V1 = f.Vd(f.sol.sol(l)[4])[1]
    lt, ll = H.lin_ligaduras(F1, F2, F3, w, *bg, *y)
    return [lt, ll + V1 * y[3] / 2]


def velocidades_V(f, l):
    """Como estabilidad_F.velocidades con el fondo de las ecuaciones con potencial (V no aporta terminos en w ni
    derivadas segundas: solo cambia el fondo)."""
    import scipy.linalg as sl
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = fondo_en_V(f, l)
    M, N = matrices_MN(F1, F2, F3, r, r1, r2, P1, P2, p1, p2, e2)
    lam = sl.eig(M, N, right=False)
    return sorted(x.real for x in lam if np.isfinite(x) and abs(x) < 1e12 and abs(x.imag) < 1e-9 * max(1, abs(x)))


def matrices_MN(F1, F2, F3, r, r1, r2, P1, P2, p1, p2, e2):
    """Parte principal M (columnas d, a, q) y coeficientes de w, N, del sistema lineal en un punto del fondo."""
    M = np.array(H.lin_principal(F1, r, r1, P1, p1))
    N = np.zeros((3, 3))
    for j, (d0, q0) in ((0, (1.0, 0.0)), (2, (0.0, 1.0))):
        N[:, j] = np.array(H.lin_resto(F1, F2, F3, 1.0, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0)) \
            - np.array(H.lin_resto(F1, F2, F3, 0.0, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0))
    return M, N


def dispersion(f, l):
    """Coeficientes (c0, c1) del polinomio de dispersion det(M - lambda N) = c0 + c1 lambda (fase2d, §3):
    c0 = det M (gradiente), c1 (cinetico); lambda = -c0/c1 es el cuadrado de la velocidad radial (reloj de la garganta)."""
    (F1, F2, F3), bg = fondo_en_V(f, l)
    M, N = matrices_MN(F1, F2, F3, *bg)
    p0 = np.linalg.det(M)
    p1 = np.linalg.det(M - N)
    pm = np.linalg.det(M + N)
    c2 = (p1 + pm) / 2 - p0
    c1 = (p1 - pm) / 2
    return p0, c1, c2


ECUACIONES_V = SimpleNamespace(derivadas=derivadas_V, ligaduras=ligaduras_V, fondo_en=fondo_en_V,
                               velocidades=velocidades_V)


# ------------------------------------------------------------------ mitades lisas
@dataclass
class LadoV:
    """Una mitad (l >= 0) desde una garganta lisa (r0 = 1) con datos (phi0, phi1, a) y el potencial que cumple la
    ligadura (si Vd no se da, amp de amp_garganta con la forma Ud).

    Finales: 'AF' (r' -> 1, masa de Misner-Sharp convergida, phi en el vacio phi_v), 'sobre' (phi cruza el vacio y
    se escapa), 'bajo' (phi vuelve hacia atras sin llegar al vacio y se escapa o queda en un valle de -V),
    'horizonte' (Phi' -> -inf), 'phi' (|phi'| -> inf), 'r0' (r -> 0), 'det' (parte principal singular), 'Phi+',
    'noAF' (llega a l_max sin cumplir el criterio), 'fallo'.
    Con disparo=True la integracion se detiene en el primer cruce del vacio.  tiro() da la clasificacion binaria del
    primer suceso del escalar (cruce del vacio o vuelta), que es la que usa el disparo."""
    Fd: object
    Vd: object
    phi0: float
    phi1: float
    a: float
    phi_v: float = 0.0
    l_max: float = 400.0
    rtol: float = 1e-10
    disparo: bool = False
    l_min_vuelta: float = 1e-6

    def y0(self):
        return np.array([1.0, 0.0, 0.0, self.a, self.phi0, self.phi1, 0.0])

    def rhs(self, l, y):
        x, _ = segundas_V(y, self.Fd, self.Vd)
        return [y[1], x[0], y[3], x[1], y[5], x[2], 1.0 / y[0]]

    def resolver(self):
        y0 = self.y0()
        self.lig0 = est_ligadura_V(y0, self.Fd, self.Vd)
        x0, det0 = segundas_V(y0, self.Fd, self.Vd)
        self.rpp0, self.Ppp0, self.ppp0 = x0
        self.det0 = det0

        def ev_det(l, y):
            return segundas_V(y, self.Fd, self.Vd)[1] / abs(det0)
        ev_det.terminal = True

        def ev_r(l, y):
            return y[0] - 1e-3
        ev_r.terminal = True

        def ev_P(l, y):
            return PHIP_MAX - max(abs(y[3]), abs(y[5]))
        ev_P.terminal = True

        def ev_esc(l, y):
            return PHI_MAX - abs(y[4] - self.phi_v)
        ev_esc.terminal = True

        def ev_rp(l, y):
            return y[1] if l > 1e-9 else (1.0 if self.rpp0 > 0 else -1.0)

        def ev_vac(l, y):                      # cruce del vacio
            return y[4] - self.phi_v
        ev_vac.terminal = self.disparo

        def ev_vuelta(l, y):                   # phi' = 0 (fuera de la garganta)
            return y[5] if l > self.l_min_vuelta else (1.0 if self.phi1 != 0 or self.ppp0 >= 0 else -1.0)
        ev_vuelta.terminal = False

        sol = solve_ivp(self.rhs, (0.0, self.l_max), y0, method="DOP853", rtol=self.rtol, atol=1e-14,
                        events=[ev_det, ev_r, ev_P, ev_esc, ev_rp, ev_vac, ev_vuelta], dense_output=True)
        self.sol = sol
        self.l_end = sol.t[-1]
        self.cruces_rp = list(sol.t_events[4])
        self.cruces_vac = list(sol.t_events[5])
        self.vueltas = list(sol.t_events[6])
        r, s, Phi, q, p, u, I = sol.y[:, -1]
        self.y_fin = sol.y[:, -1]
        self.det_fin = segundas_V(sol.y[:, -1], self.Fd, self.Vd)[1] / det0
        V_fin = self.Vd(p)[0]
        if self.disparo and len(sol.t_events[5]):
            self.fin = "sobre"
        elif sol.status == 0:
            rh, sh = sol.sol(0.5 * self.l_end)[:2]
            m, mh = r * (1 - s * s) / 2, rh * (1 - sh * sh) / 2
            ok = (s > 0 and abs(m) < 0.05 * r and abs(m - mh) <= 0.05 * max(1.0, abs(m))
                  and abs(V_fin) * r * r < 1e-6 and abs(p - self.phi_v) < 1e-3)
            self.fin = "AF" if ok else "noAF"
        elif len(sol.t_events[0]):
            self.fin = "det"
        elif len(sol.t_events[1]):
            self.fin = "r0"
        elif len(sol.t_events[3]):
            self.fin = self._sobre_o_bajo()
        elif len(sol.t_events[2]):
            self.fin = "horizonte" if abs(q) >= abs(u) and q < 0 else ("phi" if abs(u) > abs(q) else "Phi+")
        elif q < -1e3:
            self.fin = "horizonte"
        elif abs(self.det_fin) < DET_MIN:
            self.fin = "det"
        else:
            self.fin = "fallo"
        self.regular = self.fin == "AF"
        self.M = r * (1 - s * s) / 2
        self.D = r * r * u
        self.Phi_inf = Phi + q * r
        self.r_fin, self.rp_fin, self.V_fin = r, s, V_fin
        return self

    def _sobre_o_bajo(self):
        """El escalar se ha escapado: 'sobre' si cruzo el vacio, 'bajo' si no."""
        return "sobre" if len(self.cruces_vac) else "bajo"

    def tiro(self):
        """Clasificacion binaria para el disparo: +1 si phi cruza el vacio antes de volverse, -1 si se vuelve antes;
        0 si ninguna de las dos cosas (final singular u otro)."""
        tv = self.cruces_vac[0] if self.cruces_vac else math.inf
        tu = min((t for t in self.vueltas if t > self.l_min_vuelta), default=math.inf)
        if tv < tu:
            return 1
        if tu < tv:
            return -1
        return 0


class _Densa:
    SIGNO = np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0, -1.0])     # r, r', Phi, Phi', phi, phi', int dl/r

    def __init__(self, der, izq):
        self.der, self.izq = der, izq

    def sol(self, l):
        if l >= 0.0:
            return self.der.sol.sol(l)
        return self.izq.sol.sol(-l) * self.SIGNO


class DobleV:
    """Agujero liso completo: l > 0 con (phi0, phi1, a), l < 0 con (phi0, -phi1, -a) reflejado.
    Compatible con fase2c/estabilidad_lisa a traves de f.ecuaciones = ECUACIONES_V."""

    def __init__(self, der, izq):
        self.der, self.izq = der, izq
        self.Fd, self.Vd = der.Fd, der.Vd
        self.sol = _Densa(der, izq)
        self.ecuaciones = ECUACIONES_V
        self.r0 = 1.0
        self.a = der.a
        self.phi0, self.phi1 = der.phi0, der.phi1
        self.rpp0 = der.rpp0
        self.regular = der.regular and izq.regular
        self.l_end = min(der.l_end, izq.l_end)
