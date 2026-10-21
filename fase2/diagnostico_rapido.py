# -*- coding: utf-8 -*-
"""
Diagnostico del modo rapido del Modelo L (xi0 = 0.6, rho = 3e13, omega^2 ~ -0.474 km^-2):
  (a) convergencia con theta_stop, rtol y l_max;
  (b) mismo determinante por superposicion directa (si no desborda);
  (c) barrido moderado w in [-30, -1] (fisico) para xi0 = 0.1.
Cada resultado se imprime al terminar (as_completed).
"""
import sys
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.optimize import brentq
from hibrido import Politropo, Fondo, RHO_GEOM

eos = Politropo(k=0.1, gamma=2.0)


def fondo(xi0, rho, th=1e-7, lmax=2e4):
    rb = rho * RHO_GEOM
    return Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb, theta_stop=th).resolver(l_max=lmax)


def conv(th, rtol, lmax):
    f = fondo(0.6, 3e13, th, lmax)
    k = np.exp(2 * f.Phi_inf)
    w = brentq(lambda w: f.det("ligado", w, rtol=rtol), -0.6 * k, -0.35 * k, rtol=1e-12)
    return f"(a) theta_stop={th:g} rtol={rtol:g} l_max={lmax:g}: omega^2 = {f.w_fisico(w):.10f}"


def simple():
    f = fondo(0.6, 3e13)
    k = np.exp(2 * f.Phi_inf)
    ws = np.linspace(-0.55, -0.40, 7) * k
    dc = [f.det("ligado", w) for w in ws]
    ds = [f.det_simple("ligado", w) for w in ws]
    filas = [f"    w_fis={w / k:+.4f}: compuesto {a:+.3e}   directo {b:+.3e}" for w, a, b in zip(ws, dc, ds)]
    return "(b) signo de D (compuesto vs directo):\n" + "\n".join(filas)


def barrido(xi0, rho):
    f = fondo(xi0, rho)
    k = np.exp(2 * f.Phi_inf)
    ws = -np.logspace(np.log10(30.0), 0.0, 16)[::-1][::-1] * k
    ws = np.sort(ws)
    D = [f.det("ligado", w) for w in ws]
    r = [f.w_fisico(brentq(lambda w: f.det("ligado", w), ws[i], ws[i + 1], rtol=1e-10))
         for i in range(len(ws) - 1) if D[i] * D[i + 1] < 0]
    return f"(c) xi0={xi0} rho={rho:.0e}: raices del Modelo L con omega^2 en [-30, -1]: {r}"


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    with ProcessPoolExecutor(max_workers=8) as ex:
        fut = [ex.submit(conv, *a) for a in ((1e-7, 1e-10, 2e4), (1e-9, 1e-10, 2e4), (1e-5, 1e-10, 2e4),
                                             (1e-7, 1e-12, 2e4), (1e-7, 1e-10, 2e3))]
        fut += [ex.submit(simple), ex.submit(barrido, 0.1, 1e14), ex.submit(barrido, 0.3, 1e14)]
        for fu in as_completed(fut):
            print(fu.result(), flush=True)
