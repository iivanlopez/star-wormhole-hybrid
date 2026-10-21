# -*- coding: utf-8 -*-
"""
Familia de Antoniou et al. (2019) con acoplamiento cuadratico F = alpha phi^2 y phi_inf = 0.
Con r0 = 1: para cada (alpha, phi0) se buscan los phi'(0+) con phi_inf = 0 (cambios de signo sobre una malla y brentq).
Validacion: Antoniou et al. encuentran agujeros sin nodos solo para 0.205 < alpha/eta0^2 < 0.480.
Salida: resultados/familia_cuadratica.json.   Uso: python familia_cuadratica.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq
from fondo_F import FondoF

US = np.logspace(-3, 1, 60)


def phi_inf(al, p0, u):
    f = FondoF("cuadratico", al, p0, u, l_max=2e3, rtol=1e-10).resolver()
    return f.phi_inf if f.regular else float("nan")


def tarea(args):
    al, p0 = args
    vals = [phi_inf(al, p0, float(u)) for u in US]
    out = []
    for i in range(len(US) - 1):
        a, b = vals[i], vals[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a * b < 0:
            try:
                u = brentq(lambda x: phi_inf(al, p0, x), US[i], US[i + 1], xtol=1e-12, rtol=1e-11)
            except ValueError:
                continue
            f = FondoF("cuadratico", al, p0, u).resolver()
            ceros = int(np.sum(np.diff(np.sign(f.sol.y[4])) != 0))
            out.append(dict(alpha=al, phi0=p0, phi1=u, regular=bool(f.regular), phi_inf=f.phi_inf,
                            alpha_eta=f.alpha_eta, eta0=f.eta0, f0=f.f0, ecuador=bool(f.ecuador), nodos=ceros,
                            M=f.M, D=f.D, Phi_inf=f.Phi_inf, rho_zeta1=f.capa(1.0), rho_zeta0=f.capa(0.0),
                            r_garganta=getattr(f, "r_garganta", 1.0)))
    return out


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    alphas = np.round(np.arange(0.05, 0.161, 0.01), 3)
    phis = (0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5)
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = [o for lst in ex.map(tarea, [(float(a), p) for a in alphas for p in phis]) for o in lst]
    json.dump(res, open("resultados/familia_cuadratica.json", "w", encoding="utf-8"), indent=1)
    for o in sorted(res, key=lambda o: (o["nodos"], o["alpha_eta"])):
        print(f"nodos={o['nodos']} alpha/r0^2={o['alpha']:.3f} phi0={o['phi0']:.2f} phi1={o['phi1']:.5f}: "
              f"alpha/eta0^2={o['alpha_eta']:.4f} f0={o['f0']:+.4g} ecuador={o['ecuador']} e^(2Phi_inf)={math.exp(2 * o['Phi_inf']):.3g} "
              f"rho(zeta=1)={o['rho_zeta1']:+.3f} D/M={o['D'] / o['M']:+.3f}")
    for n in sorted({o["nodos"] for o in res}):
        sub = [o["alpha_eta"] for o in res if o["nodos"] == n]
        print(f"nodos={n}: alpha/eta0^2 en [{min(sub):.4f}, {max(sub):.4f}]  ({len(sub)} soluciones)")
