# -*- coding: utf-8 -*-
"""
Modelo L con perfil exotico compacto (Kodama: fantasma cuartico, Lambda = 1, L = 1 km):
  (1) agujero puro: u_inf(w) / dr(0) para varios w (anclaje: u_inf != 0 para todo w?)
  (2) comparacion con Ellis-Bronnikov (u_inf = pi/2 por unidad de dr(0)/r0)
  (3) con fluido (EOS2, rho_bc = 1e13 y 1e14): raices de D(w) en [-300, -1e-3] (modo rapido?)
"""
import sys
import math
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.optimize import brentq
from hibrido import Politropo, Fondo, RHO_GEOM
from validacion_estabilidad import ajustar_cuartico

eos = Politropo(k=0.1, gamma=2.0)


def puro():
    f = ajustar_cuartico()
    filas = []
    for w in (-3.0, -0.5, 0.0, 0.5, 3.0):
        F1, F2 = f.condiciones("ligado", w)
        filas.append(f"    w={w:+.2f}: u_inf/(dr0/r0) = {F2[0]:+.6f}")
    e = Fondo(r0=1.0).resolver()
    F1, F2 = e.condiciones("ligado", 0.0)
    return (f"(1) Kodama puro (r0 = {f.r0:.4f} km), Modelo L:\n" + "\n".join(filas)
            + f"\n(2) Ellis puro: u_inf/(dr0/r0) = {F2[0]:+.6f}  (pi/2 = {math.pi / 2:.6f})")


def con_fluido(rho):
    rb = rho * RHO_GEOM
    f = ajustar_cuartico(eos=eos, rho_bc=rb)
    k = math.exp(2 * f.Phi_inf)
    ws = np.sort(-np.logspace(np.log10(300.0), -3.0, 41) * k)
    D = [f.det("ligado", w) for w in ws]
    r = [f.w_fisico(brentq(lambda w: f.det("ligado", w), ws[i], ws[i + 1], rtol=1e-10))
         for i in range(len(ws) - 1) if D[i] * D[i + 1] < 0]
    rf = []
    ws2 = np.sort(-np.logspace(np.log10(300.0), -3.0, 41) * k)
    D2 = [f.det("fantasma", w) for w in ws2]
    rf = [f.w_fisico(brentq(lambda w: f.det("fantasma", w), ws2[i], ws2[i + 1], rtol=1e-10))
          for i in range(len(ws2) - 1) if D2[i] * D2[i + 1] < 0]
    return f"(3) Kodama + fluido rho_bc={rho:.0e} (r0 = {f.r0:.4f}): Modelo L omega^2<0: {r}   fantasma: {rf}"


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    with ProcessPoolExecutor(max_workers=3) as ex:
        fut = [ex.submit(puro), ex.submit(con_fluido, 1e13), ex.submit(con_fluido, 1e14)]
        for fu in as_completed(fut):
            print(fu.result(), flush=True)
