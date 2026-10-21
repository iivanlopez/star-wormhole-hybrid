# -*- coding: utf-8 -*-
"""
Test de la hipotesis H1 (fase1, sec. 5.3):
  mismo fondo estatico (DFKK: fluido EOS2 + fantasma sin masa, simetria Z2), dos dinamicas:
    'fantasma' -> escalar fantasma dinamico (DFKK)
    'ligado'   -> sector exotico de NOF sin grados de libertad (Modelo L)
  y la estrella de neutrones ordinaria con la misma rho_bc como referencia.
Salida: resultados/h1.json  (omega^2 fisicos en km^-2, medidos en el infinito)
Uso: python test_h1.py
"""
import sys
import os
import json
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from hibrido import Politropo, Fondo, autovalores, RHO_GEOM, MSUN_KM
from estrella_ns import EstrellaNS

eos = Politropo(k=0.1, gamma=2.0)
G_FANT = -np.logspace(1.0, -3.0, 25)                                   # w < 0 (continuo para w > 0)
G_LIG = np.concatenate([-np.logspace(0.0, -3.5, 15), np.linspace(-2e-4, 2e-3, 45)])


def analizar(args):
    xi0, rho = args
    rb = rho * RHO_GEOM
    out = dict(xi0=xi0, rho=rho)
    try:
        f = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    except Exception as e:  # noqa: BLE001
        out["error"] = str(e)
        return out
    s = f.resumen()
    out.update(r0=f.r0, R=s["R"], R_prop=s["R_prop"], M=s["M"], Phi_inf=f.Phi_inf,
               sl0=f.bg(f.y0(), True)[6])                                  # r_ll(0) > 0 <=> garganta
    rg, _ = autovalores(f, "fantasma", G_FANT)
    rl, _ = autovalores(f, "ligado", G_LIG)
    out["fantasma"] = sorted(f.w_fisico(w) for w in rg)
    out["ligado"] = sorted(f.w_fisico(w) for w in rl)
    ns = EstrellaNS(eos, rb)
    out["M_ns"] = ns.M / MSUN_KM
    out["ns"] = ns.modo_fundamental(-3e-3, 4e-3, 50)
    return out


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    rhos = [1e13, 3e13, 1e14, 3e14, 4e14, 5e14, 6e14, 8e14, 1e15, 2e15, 3e15]
    tareas = [(0.1, r) for r in rhos] + [(x, r) for x in (0.3, 0.6) for r in (3e13, 1e14, 4e14, 1e15)]
    t0 = time.time()
    print(f"nucleos: {os.cpu_count()}", flush=True)
    with ProcessPoolExecutor(max_workers=min(8, os.cpu_count() or 2)) as ex:
        res = list(ex.map(analizar, tareas))
    os.makedirs("resultados", exist_ok=True)
    with open("resultados/h1.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print(f"{'xi0':>4} {'rho_bc':>8} {'r0':>6} {'M':>7} {'M_ns':>7} | {'w0 fantasma':>12} | "
          f"{'w0 ligado':>12} {'w1 ligado':>12} | {'w0 NS':>12}")
    for o in res:
        if "error" in o:
            print(f"{o['xi0']:>4} {o['rho']:.1e}  ERROR: {o['error']}")
            continue
        fa = o["fantasma"][0] if o["fantasma"] else float("nan")
        li = o["ligado"] + [float("nan")] * 2
        print(f"{o['xi0']:>4} {o['rho']:.1e} {o['r0']:6.3f} {o['M']:7.4f} {o['M_ns']:7.4f} | {fa:12.5e} | "
              f"{li[0]:12.5e} {li[1]:12.5e} | {o['ns'] if o['ns'] is not None else float('nan'):12.5e}")
    print(f"tiempo: {time.time() - t0:.0f} s")
