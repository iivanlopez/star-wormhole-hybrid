# -*- coding: utf-8 -*-
"""
Muestra de la familia de Antoniou et al. (2019) con acoplamiento gaussiano F = alpha e^{-phi^2} y phi_inf = 0, para
comprobar que las conclusiones no dependen del acoplamiento.  Con r0 = 1: para cada (alpha, phi0 > 0) se buscan los
phi'(0+) < 0 con phi_inf = 0 (F'(phi0) < 0, asi que Phi'_0 > 0 exige phi'(0+) < 0).
Salida: resultados/familia_gauss.json.   Uso: python familia_gauss.py [n_procesos]
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq
from fondo_F import FondoF

US = -np.logspace(-3, 1, 60)
KW = {"gamma": 1.0}


def phi_inf(al, p0, u):
    f = FondoF("gauss", al, p0, u, l_max=2e3, rtol=1e-10, kw=KW).resolver()
    return f.phi_inf if f.regular else float("nan")


def tarea(args):
    al, p0 = args
    vals = [phi_inf(al, p0, float(u)) for u in US]
    out = []
    for i in range(len(US) - 1):
        a, b = vals[i], vals[i + 1]
        if np.isfinite(a) and np.isfinite(b) and a * b < 0:
            try:
                u = brentq(lambda x: phi_inf(al, p0, x), US[i + 1], US[i], xtol=1e-12, rtol=1e-11)
            except ValueError:
                continue
            f = FondoF("gauss", al, p0, u, kw=KW).resolver()
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
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    alphas = (0.06, 0.08, 0.10, 0.12, 0.14, 0.16)
    phis = (0.05, 0.1, 0.2, 0.3)
    with ProcessPoolExecutor(max_workers=n) as ex:
        res = [o for lst in ex.map(tarea, [(a, p) for a in alphas for p in phis]) for o in lst]
    json.dump(res, open("resultados/familia_gauss.json", "w", encoding="utf-8"), indent=1)
    for o in sorted(res, key=lambda o: (o["nodos"], o["ecuador"], o["alpha_eta"])):
        print(f"nodos={o['nodos']} alpha/r0^2={o['alpha']:.3f} phi0={o['phi0']:.2f} phi1={o['phi1']:+.5f}: "
              f"alpha/eta0^2={o['alpha_eta']:.4f} f0={o['f0']:+.4g} ecuador={o['ecuador']} "
              f"e^(2Phi_inf)={math.exp(2 * o['Phi_inf']):.3g} rho(zeta=1)={o['rho_zeta1']:+.3f} rho(zeta=0)={o['rho_zeta0']:+.3f}")
