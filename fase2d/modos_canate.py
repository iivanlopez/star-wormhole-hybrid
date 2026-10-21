# -*- coding: utf-8 -*-
"""
Modos radiales del agujero liso de Canate et al. con la maquinaria sin empalmes de fase2c (estabilidad_lisa, con las
ecuaciones de Maxwell de canate.py a traves del gancho f.ecuaciones).

Garganta.  F'(phi0) = 0 y Phi' = 0 en l = 0, asi que la ligadura E^l_l no contiene a'(0) y (delta r_t, delta phi_t) no
son coordenadas del cociente por el gauge (se usa Delta3).  El fondo es simetrico bajo S: l -> -l, phi -> -phi; los
sectores de S en l = 0 son
    par:   a'(0) = d'(0) = q(0) = 0       (libres d(0), q'(0), ligados por E^l_l)
    impar: d(0) = q'(0) = 0               (libres a'(0), q(0); d'(0) de E^l_t)
y el modo gauge g = (Phi'' - w, r', r'', phi', phi'') es impar.  Con U_+ la solucion decreciente en +inf (definida
modulo gauge) las condiciones de autovalor son, sin eleccion de signo y analiticas en w,
    Delta_par(w)   = phi1 a'_+(0) + w q_+(0)      (U_+ + c g par para algun c)
    Delta_impar(w) = d_+(0)                       (U_+ + c g impar; E^l_l da entonces q'_+(0) = 0)
y el determinante general de fase2c, Delta3 en coordenadas (a', q, q'), vale -2 Delta_impar Delta_par si U_- = S U_+
(se comprueba integrando U_- por separado).
Prediccion WKB: en la capa mal planteada (lambda < 0 entre los polos l1 < l < l2) la solucion oscila con
k^2 = w/lambda > 0 y fuera decae; con puntos de retorno simples en los polos,
    sqrt(|w_n|) I = (n + 1/2) pi,   I = int_{l1}^{l2} dl/sqrt(|lambda|)   (una torre por capa, pares casi degenerados).
Salida: resultados/modos_canate.json.   Uso: python modos_canate.py
"""
import os
import sys
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import canate as C                      # noqa: E402
import estabilidad_lisa as E            # noqa: E402  (fase2c, via canate.py)
from hiperbolicidad_canate import coeficientes   # noqa: E402

QS = (0.5, 1.0, 2.0, 5.0)


def capa(f):
    """polos l1 < l2 (ceros del coeficiente cinetico c1) y la integral WKB I."""
    ls = np.linspace(1e-4, 60.0, 6001)
    c1 = np.array([coeficientes(f, l)[1] for l in ls])
    pol = [brentq(lambda l: coeficientes(f, l)[1], ls[i], ls[i + 1], xtol=1e-13)
           for i in range(len(ls) - 1) if c1[i] * c1[i + 1] < 0]
    l1, l2 = pol[0], pol[1]

    def integrando(l):
        c0, c1_, _ = coeficientes(f, l)
        return math.sqrt(abs(c1_ / c0))          # 1/sqrt|lambda|, lambda = -c0/c1
    I, err = quad(integrando, l1, l2, limit=500, epsrel=1e-10)
    c0min = min(abs(coeficientes(f, l)[0]) / max(abs(np.array(coeficientes(f, l)))) for l in ls[::20])
    return l1, l2, I, c0min


def U_mas(f, w, L, rtol=1e-10):
    y, logs = E.integrar_hacia_garganta(f, L, w, rtol=rtol, renorm=True)
    return y, logs


def deltas(f, w, L, rtol=1e-10):
    """(Delta_par, Delta_impar, log de la escala, y(0)) con U_+ renormalizada (factores positivos)."""
    y, logs = U_mas(f, w, L, rtol)
    a1, d0, d1, q0, q1 = y
    return f.phi1 * a1 + w * q0, d0, logs, y


def comprobacion_general(f, w, L):
    """Delta3 en coordenadas (a', q, q') (integrando U_- desde -L por separado) frente a la factorizacion
    -2 q'_+(0) Delta_par, con q'_+(0) = (2w - 2/Q^2 - 4 Qs^2/Q^4) Q d_+(0)/Qs por E^l_l en la garganta (Q = 1:
    (2w - 2 - 4 Qs^2) d_+(0)/Qs), es decir Delta3 = -2 (2w - 2 - 4 Qs^2)/Qs Delta_impar Delta_par."""
    yp, sp_ = E.integrar_hacia_garganta(f, L, w, renorm=True)
    ym, sm_ = E.integrar_hacia_garganta(f, -L, w, renorm=True)
    g = E.modo_gauge(f, w)
    D3 = np.linalg.det(np.array([yp, ym, g]).T[[0, 3, 4], :])
    esc = ym[1] / yp[1]                                   # U_- = S U_+ salvo un factor
    Dp = f.phi1 * yp[0] + w * yp[3]
    Di = yp[1]
    fac = (2 * w - 2 - 4 * f.Qs**2) / f.Qs                # Q = 1
    espejo = np.max(np.abs(ym - esc * np.array([-yp[0], yp[1], -yp[2], -yp[3], yp[4]]))) / np.max(np.abs(ym))
    lig_q1 = abs(yp[4] - fac * yp[1]) / np.max(np.abs(yp))
    return D3, -2 * esc * fac * Di * Dp, espejo, lig_q1


def tarea(Qs):
    t0 = time.time()
    f = C.FondoCanate(Qs)
    l1, l2, I, c0min = capa(f)
    L = l2 + 15.0
    s_max = 12 * math.pi / I + 6.0                     # ~ los 12 primeros modos de la torre
    ss = np.linspace(0.02, s_max, int(s_max * I / math.pi * 40) + 200)   # ~40 puntos por semiperiodo
    ws = -ss**2
    vals = [deltas(f, float(w), max(L, 10.0 + 30.0 / math.sqrt(-w)))[:2] for w in ws]
    dp = np.array([v[0] for v in vals])
    di = np.array([v[1] for v in vals])
    raices = {}
    for nombre, arr, k in (("par", dp, 0), ("impar", di, 1)):
        rs = []
        for i in range(len(ws) - 1):
            if arr[i] * arr[i + 1] < 0:
                Lw = lambda w: max(L, 10.0 + 30.0 / math.sqrt(-w))           # noqa: E731
                w0 = brentq(lambda w: deltas(f, w, Lw(w))[k], ws[i + 1], ws[i], xtol=1e-12, rtol=1e-12)
                # comprobacion con la escala completa (log |Delta| sin renormalizar): cero simple
                vv = []
                for x in (w0 * (1 - 1e-4), w0, w0 * (1 + 1e-4)):
                    Dp_, Di_, lg, _ = deltas(f, x, Lw(x))
                    D = (Dp_, Di_)[k]
                    vv.append(math.log(abs(D)) + lg if D != 0 else -math.inf)
                rs.append(dict(w=float(w0), simple=bool(vv[1] < min(vv[0], vv[2]) - math.log(1e3)),
                               log_ratio=float(min(vv[0], vv[2]) - vv[1])))
        raices[nombre] = rs
    comp = [comprobacion_general(f, w, L) for w in (-0.5, -7.0, -40.0)]
    wkb = [-((n + 0.5) * math.pi / I) ** 2 for n in range(14)]
    return dict(Qs=Qs, l1=l1, l2=l2, I=I, c0min_rel=float(c0min), raices=raices, wkb=wkb,
                comprobacion=[dict(D3=float(a), producto=float(b), espejo=float(c), lig_q1=float(d)) for a, b, c, d in comp],
                tiempo=time.time() - t0)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(tarea, QS))
    for o in res:
        print(f"\nQs/Q = {o['Qs']}: capa mal planteada l in ({o['l1']:.6f}, {o['l2']:.6f}) Q;  I = {o['I']:.6f};  "
              f"min |det M|/escala = {o['c0min_rel']:.2e};  ({o['tiempo']:.0f} s)")
        for c in o["comprobacion"]:
            print(f"   Delta3 = {c['D3']:+.10e}  vs  factorizacion = {c['producto']:+.10e};  U_- = S U_+ a {c['espejo']:.1e};  "
                  f"q'(0) - fac d(0) = {c['lig_q1']:.1e}")
        for nombre in ("par", "impar"):
            rs = o["raices"][nombre]
            print(f"   {nombre}: {len(rs)} raices reales; todas simples: {all(r['simple'] for r in rs)}")
            for n, r in enumerate(rs[:14]):
                wk = o["wkb"][n] if n < len(o["wkb"]) else float("nan")
                print(f"      n={n:2d}: w = {r['w']:+.8e}  sqrt|w| I/pi - 1/2 = {math.sqrt(-r['w']) * o['I'] / math.pi - 0.5:+.4f}"
                      f"   (WKB: {wk:+.4e})")
    json.dump(res, open("resultados/modos_canate.json", "w", encoding="utf-8"), indent=1)
    print(f"tiempo: {time.time() - t0:.0f} s")
