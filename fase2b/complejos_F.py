# -*- coding: utf-8 -*-
"""
Modos inestables complejos de los fondos bien planteados SIN raiz real (barridos robustos), con la funcion entera Dc
(estabilidad_F.Dc, rama principal de kappa: un cero con w fuera del corte w > 0 es un modo que crece).
  1. Recuento de ceros por el principio del argumento en el semiplano superior (D(w-barra) = conj D(w)):
     rectangulo [-X a^2, X a^2] x [y0 a^2, Y a^2], muestreo adaptativo del borde (salto de fase < 0.4 rad).
  2. Localizacion: secante desde una malla de arranques (en unidades de a^2 = Phi'(0+)^2).
Casos: gaussianos 0.2821 y 0.3179 con zeta = 0; doble garganta exponencial (alpha = 0.04, phi' = -4.534) con zeta = 1.
Salida: resultados/complejos_F.json.   Uso: python complejos_F.py
"""
import sys
import json
import cmath
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import CapaF, Dc, secante

X, Y0, Y1 = 12.0, 0.01, 12.0
ARRANQUES = [complex(x, y) for x in (-3.0, -1.0, -0.3, 0.3, 1.0, 3.0, 6.0) for y in (0.1, 0.5, 2.0)]


def fondo(caso):
    tipo, al, p0, u, kw, z = caso
    f = FondoF(tipo, al, p0, u, kw=kw).resolver()
    return f, CapaF(f, zeta=z), f.y0()[3] ** 2


def _dc(args):
    caso, w = args
    f, c, _ = fondo(caso)
    return complex(Dc(f, c, w))


def borde(t, a2):
    """Parametrizacion del rectangulo (sentido antihorario), t en [0, 4)."""
    k, s = int(t) % 4, t - int(t)
    x0, x1, y0, y1 = -X * a2, X * a2, Y0 * a2, Y1 * a2
    return [complex(x0 + (x1 - x0) * s, y0), complex(x1, y0 + (y1 - y0) * s),
            complex(x1 + (x0 - x1) * s, y1), complex(x0, y1 + (y0 - y1) * s)][k]


def contar(caso, a2, ex, n=64, salto=0.4, rondas=14):
    D = {}
    pend = [float(t) for t in np.linspace(0, 4, 4 * n, endpoint=False)]
    for _ in range(rondas):
        if not pend:
            break
        for t, v in zip(pend, ex.map(_dc, [(caso, borde(t, a2)) for t in pend])):
            D[t] = v
        ts = sorted(D)
        pend = []
        for i in range(len(ts)):
            t0, t1 = ts[i], ts[(i + 1) % len(ts)]
            dph = cmath.phase(D[t1] / D[t0])
            if abs(dph) > salto:
                pend.append((t0 + ((t1 - t0) % 4) / 2) % 4)
    ts = sorted(D)
    total = sum(cmath.phase(D[ts[(i + 1) % len(ts)]] / D[ts[i]]) for i in range(len(ts)))
    return total / (2 * math.pi), len(ts), len(pend)


def _secante(args):
    caso, x0 = args
    f, c, a2 = fondo(caso)
    w, paso, dD = secante(f, c, x0 * a2)
    vec = max(abs(Dc(f, c, w * (1 + 1e-3 * e))) for e in (1, -1, 1j, -1j))
    return dict(arranque=[x0.real, x0.imag], w=[w.real, w.imag], w_a2=[(w / a2).real, (w / a2).imag], paso=paso,
                Dc=dD, Dc_vecino=vec, convergido=bool(paso < 1e-8 * abs(w) and dD < 1e-6 * vec and w.imag > 0),
                winf=[(w * math.exp(-2 * f.Phi_inf)).real, (w * math.exp(-2 * f.Phi_inf)).imag])


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    gauss = {f"{o['alpha_eta']:.4f}": o for o in json.load(open("resultados/familia_gauss.json", encoding="utf-8"))}
    casos = {"gauss 0.2821 zeta=0": ("gauss", gauss["0.2821"]["alpha"], gauss["0.2821"]["phi0"], gauss["0.2821"]["phi1"],
                                      {"gamma": 1.0}, 0.0),
             "gauss 0.3179 zeta=0": ("gauss", gauss["0.3179"]["alpha"], gauss["0.3179"]["phi0"], gauss["0.3179"]["phi1"],
                                      {"gamma": 1.0}, 0.0),
             "doble exp a=0.04 zeta=1": ("exp", 0.04, 0.0, -4.534, {}, 1.0)}
    res = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for nom, caso in casos.items():
            f, c, a2 = fondo(caso)
            n, npts, quedan = contar(caso, a2, ex)
            raices = list(ex.map(_secante, [(caso, x0) for x0 in ARRANQUES]))
            buenas = []
            for r in raices:
                if r["convergido"] and not any(abs(complex(*r["w"]) - complex(*b["w"])) < 1e-6 * abs(complex(*r["w"]))
                                               for b in buenas):
                    buenas.append(r)
            o = dict(caso=nom, a2=a2, e2Phi=math.exp(2 * f.Phi_inf), ceros_rectangulo=n, puntos=npts, sin_resolver=quedan,
                     rectangulo=[-X, X, Y0, Y1], raices=buenas)
            res.append(o)
            print(f"{nom}: a^2 = {a2:.4g}; ceros en el rectangulo (unidades de a^2) [-{X}, {X}] x [{Y0}, {Y1}]: {n:+.3f} "
                  f"({npts} puntos, {quedan} tramos sin resolver)", flush=True)
            for r in buenas:
                print(f"    w = {r['w'][0]:.6g} {r['w'][1]:+.6g} i   (w/a^2 = {r['w_a2'][0]:.4f} {r['w_a2'][1]:+.4f} i;  "
                      f"omega_inf^2 r0^2 = {r['winf'][0]:.4f} {r['winf'][1]:+.4f} i)  |Dc| = {r['Dc']:.1e} frente a "
                      f"{r['Dc_vecino']:.1e}", flush=True)
            json.dump(res, open("resultados/complejos_F.json", "w", encoding="utf-8"), indent=1)
