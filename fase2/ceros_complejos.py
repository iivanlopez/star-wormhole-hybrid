# -*- coding: utf-8 -*-
"""
Localiza los ceros que el principio del argumento (contorno.py) cuenta de mas frente a las raices reales.
  (1) Barrido real fino de D(w) en [-A, B] (descarta un par de raices reales muy proximas).
  (2) Secante en el plano complejo sobre la funcion entera G(w) = D(w) e^{log s(w)} (det con escala=True),
      desde una rejilla de puntos del semiplano superior del rectangulo de contorno.py.
w en el reloj de la garganta (interno); en el infinito, omega^2 = e^{-2 Phi_inf} w.
Salida: resultados/ceros_complejos.json.   Uso: python ceros_complejos.py   (tras contorno.py)
"""
import sys
import json
import math
import cmath
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from hibrido import Politropo, Fondo, RHO_GEOM

eos = Politropo(k=0.1, gamma=2.0)
_fondos = {}


def fondo(xi0, rho):
    if (xi0, rho) not in _fondos:
        rb = rho * RHO_GEOM
        _fondos[(xi0, rho)] = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    return _fondos[(xi0, rho)]


def signo(args):
    xi0, rho, w = args
    return float(np.sign(fondo(xi0, rho).det("ligado", w)))


def secante(args):
    xi0, rho, w0, escala_max, it = args
    f = fondo(xi0, rho)
    w1 = w0 * (1.0 + 1e-3) + 1e-3j * abs(w0)
    D0, l0 = f.det("ligado", w0, escala=True)
    D1, l1 = f.det("ligado", w1, escala=True)
    for _ in range(it):
        g0 = D0 * math.exp(min(max(l0 - l1, -700.0), 700.0))
        if D1 == g0:
            break
        paso = -D1 * (w1 - w0) / (D1 - g0)
        if abs(paso) > escala_max:                     # amortiguacion de pasos enormes
            paso *= escala_max / abs(paso)
        w0, D0, l0 = w1, D1, l1
        w1 = w1 + paso
        D1, l1 = f.det("ligado", w1, escala=True)
        if abs(paso) < 1e-11 * max(abs(w1), 1e-12):
            break
    return dict(w=[w1.real, w1.imag], paso=abs(paso), D=abs(D1))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    cont = json.load(open("resultados/contorno.json", encoding="utf-8"))
    sal = []
    with ProcessPoolExecutor(max_workers=4) as ex:
        for c in cont:
            if abs(c["N"] - c["reales"]) < 0.5:
                continue
            xi0, rho, A, B, C, k = c["xi0"], c["rho"], c["A"], c["B"], c["C"], c["k"]
            # (1) barrido real fino
            g = np.concatenate([-np.logspace(math.log10(A), -6.0, 400), np.linspace(B / 200, B, 100)])
            v = np.array(list(ex.map(signo, [(xi0, rho, float(w)) for w in g], chunksize=8)))
            reales = [float(g[j]) for j in np.nonzero(v[:-1] * v[1:] < 0)[0]]
            # (2) secante desde el semiplano superior
            res = [-A * x for x in (0.9, 0.6, 0.3, 0.1, 0.03, 0.01, 0.003)] + [0.5 * B]
            ims = [C * y for y in (0.05, 0.2, 0.6)]
            arranques = [(xi0, rho, complex(x, y), 0.25 * max(A, C), 60) for x in res for y in ims]
            ceros = []
            for o in ex.map(secante, arranques):
                w = complex(*o["w"])
                dentro = -A < w.real < B and abs(w.imag) < C
                if o["paso"] < 1e-8 * max(abs(w), 1e-12) and dentro:
                    if not any(abs(w - z) < 1e-6 * max(abs(w), 1e-12) for z in ceros):
                        ceros.append(w)
            cpx = [z for z in ceros if abs(z.imag) > 1e-6 * C]
            fis = []
            for z in cpx:
                wf = z / k
                om = cmath.sqrt(wf)
                fis.append(dict(w_int=[z.real, z.imag], w_fis=[wf.real, wf.imag],
                                omega_fis=[om.real, om.imag], tau_us=1e6 / (3.0e5 * abs(om.imag))))
            sal.append(dict(xi0=xi0, rho=rho, N=c["N"], reales_finos=reales,
                            ceros_secante=[[z.real, z.imag] for z in ceros], complejos=fis))
            print(f"xi0={xi0} rho={rho:.0e}: N = {c['N']:.3f};  raices reales (barrido fino) = "
                  f"{[f'{x:.4e}' for x in reales]}", flush=True)
            print(f"   ceros por secante (w interno): {[f'{z.real:.5e}{z.imag:+.5e}j' for z in ceros]}", flush=True)
            for o in fis:
                print(f"   complejo: omega^2 (infinito) = {o['w_fis'][0]:.5e}{o['w_fis'][1]:+.5e}j km^-2;  "
                      f"omega = {o['omega_fis'][0]:.4e}{o['omega_fis'][1]:+.4e}j km^-1;  tau = {o['tau_us']:.1f} us",
                      flush=True)
            json.dump(sal, open("resultados/ceros_complejos.json", "w", encoding="utf-8"), indent=1)
