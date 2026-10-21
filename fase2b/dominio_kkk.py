# -*- coding: utf-8 -*-
"""
Dominio de existencia de los agujeros de KKK en su plano (alpha/r0^2, f0), con phi_inf = 0 y gamma = 1.
Se barre el plano de integracion (alpha, phi'(0+)) con r0 = 1, phi0 = 0 y se transforma con la simetria de escala.
Referencias (KKK, PRD 85 044007): existen para alpha/r0^2 <~ 0.13; para alpha/r0^2 < 0.0507 cubren 1 < f0 < inf.
Salida: resultados/dominio_kkk.json.   Uso: python dominio_kkk.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from kkk import FondoKKK


def tarea(args):
    al, u0 = args
    try:
        f = FondoKKK(alpha=al, phi1=u0, l_max=2.0e3, rtol=1e-10).resolver()
        rho, lam = f.capa()
        return dict(alpha=al, phi1=u0, f0=f.f0, alpha_kkk=f.alpha_kkk, regular=bool(f.regular), l_end=f.l_end,
                    M=f.M, D=f.D, rho=rho, ligadura=f.ligadura_max)
    except Exception as e:
        return dict(alpha=al, phi1=u0, error=str(e))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    alphas = np.logspace(-3, math.log10(0.6), 34)
    u0s = -np.logspace(-2, 1.7, 40)
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(tarea, [(float(a), float(u)) for a in alphas for u in u0s], chunksize=8))
    json.dump(res, open("resultados/dominio_kkk.json", "w", encoding="utf-8"), indent=1)
    reg = [o for o in res if o.get("regular")]
    sing = [o for o in res if "f0" in o and not o.get("regular")]
    print(f"{len(reg)} regulares, {len(sing)} no regulares, {len(res) - len(reg) - len(sing)} errores")
    print(f"alpha_kkk maximo entre las regulares: {max(o['alpha_kkk'] for o in reg):.4f}")
    for lo, hi in ((0.0, 0.0507), (0.0507, 0.08), (0.08, 0.11), (0.11, 0.14)):
        sub = [o for o in reg if lo <= o["alpha_kkk"] < hi]
        if sub:
            print(f"  alpha_kkk en [{lo}, {hi}): f0 de {min(o['f0'] for o in sub):.4f} a {max(o['f0'] for o in sub):.4g}"
                  f"  ({len(sub)} soluciones; rho min = {min(o['rho'] for o in sub):.3g})")
