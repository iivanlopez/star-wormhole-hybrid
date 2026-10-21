# -*- coding: utf-8 -*-
"""
Modos radiales de un agujero EsGB LISO (fondo_doble.Doble), sin cascara ni empalmes.
Perturbaciones de fase2b/estabilidad_F: estado y = (a', d, d', q, q') con A = Phi + a e^{i w^{1/2} t}, R = r + d e^{...},
phi + q e^{...}; gauge de longitud propia (g_ll = 1, shift nulo); w = omega^2 con el reloj de la garganta (Phi(0) = 0).

Gauge residual.  xi^l = f(t), xi^t = f'(t) int_0^l e^{-2 Phi} preserva g_ll = 1 y g_tl = 0; para cada w genera la
solucion exacta  g(l) = (Phi'' - w e^{-2 Phi}, r', r'', phi', phi'') (modo gauge), que no tiene parte exponencial
(q = f phi' ~ 1/r^2): cumpliria las condiciones de decaimiento en los dos infinitos para TODO w.  En la garganta
g(0) = (Phi''(0) - w, 0, r''(0), phi1, phi''(0)) tiene d'(0) = r''(0) != 0, asi que la condicion d'(0) = 0 (la garganta
perturbada, minimo de R = r + d, sigue en l = 0) elige un representante unico de cada orbita gauge.  Las magnitudes
invariantes en la garganta son
    delta r_t   = d(0)                                  (radio de la garganta; d = f r' se anula en l = 0)
    delta phi_t = q(0) - phi1 d'(0)/r''(0)              (escalar en la garganta desplazada l_t = -d'(0)/r''(0))
y el espacio de datos admisibles en l = 0 (3 dimensiones: 5 componentes menos las dos ligaduras) modulo el gauge es
el plano (delta r_t, delta phi_t).

Condicion de autovalor.  U_+ (U_-): solucion que decae en l -> +inf (-inf), q ~ e^{-k|l|}/r con k = sqrt(-w) e^{-Phi}
(rama principal, Re k > 0), integrada hacia la garganta.  El modo existe si y solo si U_+, U_- y g son linealmente
dependientes en l = 0:
    Delta(w) = delta r_t^+ delta phi_t^- - delta r_t^- delta phi_t^+  ( = det[U_+, U_-, g] / r''(0) en coordenadas
               (d, delta phi_t, d') ),
funcion analitica de w fuera del corte w > 0, sin eleccion de signo de vectores unitarios.  Cada cero con w fuera de
[0, inf) es un modo que crece en el tiempo (Re k > 0  <=>  Im omega < 0).
"""
import os
import sys
import math
import cmath
import numpy as np
from scipy.integrate import solve_ivp, quad

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fase2b"))
from types import SimpleNamespace                                           # noqa: E402
from estabilidad_F import derivadas, ligaduras, fondo_en, velocidades     # noqa: E402,F401

_POR_DEFECTO = SimpleNamespace(derivadas=derivadas, ligaduras=ligaduras, fondo_en=fondo_en, velocidades=velocidades)


def _ec(f):
    """Ecuaciones del fondo: las de fase2b (EsGB puro) salvo que el fondo traiga las suyas en f.ecuaciones
    (p. ej. EsGB + Maxwell en fase2d/canate.py).  Sin ese atributo el comportamiento es el de la Fase 2c."""
    return getattr(f, "ecuaciones", None) or _POR_DEFECTO


def matriz_ligaduras(f, l, w):
    """C (2x5): ligaduras E^l_t, E^l_l = C y (lineales y homogeneas en y)."""
    cols = []
    for k in range(5):
        e = np.zeros(5, dtype=complex if isinstance(w, complex) else float)
        e[k] = 1.0
        cols.append(_ec(f).ligaduras(f, l, e, w))
    return np.array(cols, dtype=complex if isinstance(w, complex) else float).T


def kloc(f, l, w):
    """numero de onda local k = sqrt(-w) e^{-Phi(l)}, rama principal (Re k >= 0)."""
    Phi = f.sol.sol(l)[2]
    if isinstance(w, complex):
        k = cmath.sqrt(-w) * math.exp(-Phi)
        return k if k.real >= 0 else -k
    return math.sqrt(-w) * math.exp(-Phi)


def dato_asintotico(f, L, w):
    """Solucion decreciente en l = L (L > 0: lado derecho; L < 0: izquierdo): q = 1/r, q' = (-+k - r'/r) q, d = 0;
    a' y d' de las dos ligaduras."""
    r, rp = f.sol.sol(L)[:2]
    k = kloc(f, L, w)
    s = 1.0 if L > 0 else -1.0
    q0 = 1.0 / r
    q1 = (-s * k - rp / r) * q0
    C = matriz_ligaduras(f, L, w)
    A = C[:, [0, 2]]
    rhs = -(C[:, 3] * q0 + C[:, 4] * q1)
    a1, d1 = np.linalg.solve(A, rhs)
    tipo = complex if isinstance(w, complex) else float
    return np.array([a1, 0.0, d1, q0, q1], dtype=tipo)


def fase_wkb(f, L, w):
    """int_0^|L| k dl en el lado de L (para normalizar la solucion decreciente: q(L) = e^{-int k}/r)."""
    s = 1.0 if L > 0 else -1.0
    I, _ = quad(lambda x: math.exp(-f.sol.sol(s * x)[2]), 0.0, abs(L), limit=400, epsabs=0, epsrel=1e-11)
    if isinstance(w, complex):
        k = cmath.sqrt(-w)
        k = k if k.real >= 0 else -k
        return k * I
    return math.sqrt(-w) * I


def integrar_hacia_garganta(f, L, w, rtol=1e-10, renorm=True):
    """Integra el dato decreciente desde l = L hasta 0.  Devuelve (y(0), log de la escala acumulada).
    Con renorm la solucion se divide por factores positivos cuando crece (no cambia la fase ni los ceros)."""
    y = dato_asintotico(f, L, w)
    logs = 0.0
    l = L
    tipo = complex if isinstance(w, complex) else float

    def grande(x, yy):
        return np.max(np.abs(yy)) - 1e60
    grande.terminal = True
    ev = [grande] if renorm else []
    while abs(l) > 1e-14:
        s = solve_ivp(lambda x, yy: _ec(f).derivadas(f, x, yy, w), (l, 0.0), y.astype(tipo), method="DOP853", rtol=rtol,
                      atol=1e-14 * max(np.max(np.abs(y)), 1e-300), events=ev)
        y, l = s.y[:, -1], s.t[-1]
        if s.status == -1:
            raise RuntimeError(s.message)
        if renorm:
            m = np.max(np.abs(y))
            y = y / m
            logs += math.log(m)
    return y, logs


def invariantes_garganta(f, y):
    """(delta r_t, delta phi_t) a partir de y(0) = (a', d, d', q, q')."""
    rpp0 = f.rpp0
    return y[1], y[3] - f.phi1 * y[2] / rpp0


def modo_gauge(f, w, l=0.0):
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = _ec(f).fondo_en(f, l)
    return np.array([P2 - w * e2, r1, r2, p1, p2], dtype=complex if isinstance(w, complex) else float)


def tail(f, w, l_extra=30.0):
    """longitudes de las colas: L_pm = 10 r0 + l_extra/|k_pm(inf)|, sin pasar del 90% del fondo integrado."""
    out = []
    for lado, s in ((f.der, 1.0), (f.izq, -1.0)):
        kap = abs(cmath.sqrt(-complex(w))) * math.exp(-lado.Phi_inf)
        kr = (cmath.sqrt(-complex(w)) * math.exp(-lado.Phi_inf)).real
        kk = max(min(kap, kr if kr > 0 else kap), 1e-3)
        out.append(s * min(10.0 + l_extra / kk, 0.9 * lado.l_end))
    return out


def Delta(f, w, l_extra=30.0, rtol=1e-10, renorm=True, normalizar=True, devolver_todo=False):
    """Determinante de emparejamiento Delta(w) (analitico; con renorm su modulo no tiene significado, su fase si).
    normalizar: multiplica por e^{-int_0^{L+} k - int_0^{|L-|} k} (factor analitico sin ceros que quita la rotacion de
    fase de las colas)."""
    Lp, Lm = tail(f, w, l_extra)
    yp, sp_ = integrar_hacia_garganta(f, Lp, w, rtol, renorm)
    ym, sm_ = integrar_hacia_garganta(f, Lm, w, rtol, renorm)
    rp, pp = invariantes_garganta(f, yp)
    rm, pm = invariantes_garganta(f, ym)
    D = rp * pm - rm * pp
    logesc = sp_ + sm_
    if normalizar:
        fz = fase_wkb(f, Lp, w) + fase_wkb(f, Lm, w)
        if isinstance(w, complex):
            D = D * cmath.exp(-1j * fz.imag)
            logesc -= fz.real
        else:
            logesc -= fz
    if devolver_todo:
        return D, logesc, yp, ym, (Lp, Lm)
    return D


def Delta_sin_renormalizar(f, w, l_extra=30.0, rtol=1e-10):
    """Delta con una sola integracion por lado, sin dividir por ningun factor (modulo con significado), normalizado
    con la fase WKB: Delta_true = Delta_ren * e^{log escala}.  Devuelve (Delta, log|Delta|)."""
    D, logesc, *_ = Delta(f, w, l_extra, rtol, renorm=False, normalizar=True, devolver_todo=True)
    return D, (math.log(abs(D)) + logesc if D != 0 else -math.inf)


def raices_reales(f, ws, l_extra=30.0, rtol=1e-10, refinar=True):
    """Ceros de Delta(w) en la malla real ws (w < 0): cambios de signo y brentq."""
    from scipy.optimize import brentq
    vals = [Delta(f, float(w), l_extra, rtol) for w in ws]
    out = []
    for i in range(len(ws) - 1):
        a, b = vals[i], vals[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a * b < 0:
            if not refinar:
                out.append(0.5 * (ws[i] + ws[i + 1]))
                continue
            try:
                out.append(brentq(lambda x: Delta(f, x, l_extra, rtol), ws[i], ws[i + 1], xtol=1e-12, rtol=1e-11))
            except ValueError:
                out.append(float("nan"))
    return out, vals


def hiperbolicidad(f, ls):
    """lambda_min(l) = minimo de los cuadrados de las velocidades caracteristicas (reloj de la garganta) en cada l."""
    out = []
    for l in ls:
        lam = _ec(f).velocidades(f, float(l))
        out.append(min(lam) if lam else float("nan"))
    return np.array(out)


def Delta3(f, w, coords, l_extra=30.0, rtol=1e-10, renorm=True, devolver_todo=False):
    """Forma general del emparejamiento: det de [U_+(0), U_-(0), g(0)] restringido a tres componentes coords de
    y = (a', d, d', q, q') que sean coordenadas del espacio admisible en la garganta (las otras dos despejables de las
    ligaduras).  Vale tambien cuando E^l_l no contiene a'(0) (F'(phi0) phi'(0) = 0 en la garganta), donde
    (delta r_t, delta phi_t) dejan de ser coordenadas del cociente por el gauge y Delta() no sirve.
    Para gargantas EsGB puras (fase2c) coords = (1, 3, 4) da Delta() salvo un factor."""
    Lp, Lm = tail(f, w, l_extra)
    yp, sp_ = integrar_hacia_garganta(f, Lp, w, rtol, renorm)
    ym, sm_ = integrar_hacia_garganta(f, Lm, w, rtol, renorm)
    g = modo_gauge(f, w)
    A = np.array([yp, ym, g]).T[list(coords), :]
    D = np.linalg.det(A)
    if devolver_todo:
        return D, sp_ + sm_, yp, ym, (Lp, Lm)
    return D
