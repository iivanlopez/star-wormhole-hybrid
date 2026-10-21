# -*- coding: utf-8 -*-
"""
Controles de los modos complejos (complejos_F.py, modos_complejos_F.py):
  1. gaussiano 0.2821 (zeta = 0): el recuento del principio del argumento da 1 cero en [-12, 12] x [0.01, 12] (unidades
     de a^2), pero la raiz localizada (Im w/a^2 = 0.0097) cae justo por debajo del borde.  Se repite el recuento con el
     borde inferior en 0.001 y por mitades (Re w < 0 y Re w > 0) para localizar el otro cero.
  2. cuadratico 0.3087 con zeta = 2.5 y 3: resonancias cerca del corte (Im w/Re w ~ 5e-3).  Convergencia frente a la cola
     de integracion l_extra (30, 60, 120): si la raiz se mueve mas que su parte imaginaria, su signo no es fiable.
Salida: resultados/complejos_control.json.   Uso: python complejos_control.py
"""
import sys
import json
import cmath
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import CapaF, Dc
import complejos_F as CF


def contar_rect(caso, a2, ex, x0, x1, y0, y1, n=64, salto=0.4, rondas=16):
    def borde(t):
        k, s = int(t) % 4, t - int(t)
        X0, X1, Y0, Y1 = x0 * a2, x1 * a2, y0 * a2, y1 * a2
        return [complex(X0 + (X1 - X0) * s, Y0), complex(X1, Y0 + (Y1 - Y0) * s),
                complex(X1 + (X0 - X1) * s, Y1), complex(X0, Y1 + (Y0 - Y1) * s)][k]
    D, pend = {}, [float(t) for t in np.linspace(0, 4, 4 * n, endpoint=False)]
    for _ in range(rondas):
        if not pend:
            break
        for t, v in zip(pend, ex.map(CF._dc, [(caso, borde(t)) for t in pend])):
            D[t] = v
        ts = sorted(D)
        pend = [(ts[i] + ((ts[(i + 1) % len(ts)] - ts[i]) % 4) / 2) % 4 for i in range(len(ts))
                if abs(cmath.phase(D[ts[(i + 1) % len(ts)]] / D[ts[i]])) > salto]
    ts = sorted(D)
    return sum(cmath.phase(D[ts[(i + 1) % len(ts)]] / D[ts[i]]) for i in range(len(ts))) / (2 * math.pi), len(ts), len(pend)


def secante_l(args):
    caso, w0, l_extra = args
    f, c, a2 = CF.fondo(caso)
    w0 = complex(w0)
    w1 = w0 * (1 + 1e-4) + 1e-3j
    D0, D1 = Dc(f, c, w0, l_extra=l_extra), Dc(f, c, w1, l_extra=l_extra)
    for _ in range(60):
        if D1 == D0:
            break
        paso = -D1 * (w1 - w0) / (D1 - D0)
        if abs(paso) > 0.05 * abs(w1):
            paso *= 0.05 * abs(w1) / abs(paso)
        w0, D0 = w1, D1
        w1 = w1 + paso
        D1 = Dc(f, c, w1, l_extra=l_extra)
        if abs(paso) < 1e-11 * abs(w1):
            break
    return dict(l_extra=l_extra, w=[w1.real, w1.imag], paso=abs(paso), Dc=abs(D1),
                winf=[(w1 * math.exp(-2 * f.Phi_inf)).real, (w1 * math.exp(-2 * f.Phi_inf)).imag])


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    gauss = {f"{o['alpha_eta']:.4f}": o for o in json.load(open("resultados/familia_gauss.json", encoding="utf-8"))}
    g = gauss["0.2821"]
    caso_g = ("gauss", g["alpha"], g["phi0"], g["phi1"], {"gamma": 1.0}, 0.0)
    fam = {f"{o['alpha_eta']:.4f}": o for o in json.load(open("resultados/familia_cuadratica.json", encoding="utf-8"))}
    q = fam["0.3087"]
    res = {}
    with ProcessPoolExecutor(max_workers=8) as ex:
        a2 = CF.fondo(caso_g)[2]
        res["gauss_0.2821"] = {}
        for nom, rect in (("completo_y0=0.001", (-12, 12, 0.001, 12)), ("Re<0", (-12, 0, 0.001, 12)),
                          ("Re>0", (0, 12, 0.001, 12)), ("franja_Re>0_baja", (0, 12, 0.001, 0.05))):
            n, npts, quedan = contar_rect(caso_g, a2, ex, *rect)
            res["gauss_0.2821"][nom] = dict(rect=rect, ceros=n, puntos=npts, sin_resolver=quedan)
            print(f"gauss 0.2821 zeta=0, rectangulo {nom} {rect}: {n:+.3f} ceros ({npts} puntos, {quedan} sin resolver)",
                  flush=True)
        res["cuadratico_0.3087"] = []
        for z, w0 in ((2.5, 716.403 + 3.482j), (3.0, 379.981 + 1.22038j)):
            caso_q = ("cuadratico", q["alpha"], q["phi0"], q["phi1"], {}, z)
            for o in ex.map(secante_l, [(caso_q, w0, le) for le in (30.0, 60.0, 120.0)]):
                o["zeta"] = z
                res["cuadratico_0.3087"].append(o)
                print(f"cuadratico 0.3087 zeta={z:g} l_extra={o['l_extra']:g}: w = {o['w'][0]:.6f} {o['w'][1]:+.6f} i  "
                      f"(omega_inf^2 r0^2 = {o['winf'][0]:.6f} {o['winf'][1]:+.6f} i)  |Dc| = {o['Dc']:.1e}", flush=True)
    json.dump(res, open("resultados/complejos_control.json", "w", encoding="utf-8"), indent=1)
