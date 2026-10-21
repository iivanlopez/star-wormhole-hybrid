# -*- coding: utf-8 -*-
"""
Hiperbolicidad de las perturbaciones radiales a lo largo de la CONTINUACION LISA de gargantas con mitad atractiva AF:
el dominio es (-l_s, +inf), desde el horizonte singular de la mitad repulsiva hasta el infinito de la atractiva.
Con estabilidad_F.velocidades (cuadrados de las velocidades caracteristicas, reloj de la garganta) sobre fondo_doble.Doble;
la luz tiene lambda = N^2 = e^{2 Phi}, asi que lambda/N^2 es la velocidad local relativa a la luz.
Muestra: 40 gargantas al azar del mapa (cuadratico y exponencial) con mitad atractiva AF, mas la de Antoniou de §5.
Salida: resultados/hiperbolicidad_lisa.json.   Uso: python hiperbolicidad_lisa.py
"""
import os
import sys
import json
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fondo_doble import Doble                       # noqa: E402
from estabilidad_lisa import velocidades            # noqa: E402


def tarea(args):
    tipo, al, p0, u = args
    f = Doble(tipo, al, p0, u, l_max=400.0, rtol=1e-10)
    # lado atractivo: el de a > 0 en su propia direccion; Doble pone (phi0, phi1) en l > 0 con Phi'(0) = a
    if f.a > 0:
        lat, lrep, s_rep = f.der, f.izq, -1.0
    else:
        lat, lrep, s_rep = f.izq, f.der, 1.0
    out = dict(tipo=tipo, alpha=al, phi0=p0, phi1=u, a=abs(f.a), fin_atr=lat.fin, fin_rep=lrep.fin, l_s=float(lrep.l_end))
    ls_rep = s_rep * lrep.l_end * (1 - np.logspace(-6, 0, 120)[::-1])        # densa cerca del final singular
    ls_atr = -s_rep * np.r_[np.linspace(0, 5, 100), np.logspace(np.log10(5.5), np.log10(0.95 * lat.l_end), 80)]
    vals = {}
    for nombre, ls in (("rep", ls_rep), ("atr", ls_atr)):
        lmin, cmin, neg = np.inf, np.inf, 0
        for l in ls:
            lam = velocidades(f, float(l))
            N2 = np.exp(2 * f.sol.sol(float(l))[2])
            if not lam:
                continue
            m = min(lam)
            lmin = min(lmin, m)
            cmin = min(cmin, m / N2)
            neg += m < 0
        vals[nombre] = dict(lambda_min=float(lmin), c_min=float(cmin), puntos_negativos=int(neg), n=len(ls))
    out.update(vals)
    return out


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    t0 = time.time()
    res = json.load(open("resultados/mapa_lados.json", encoding="utf-8"))
    idx = {(o["tipo"], o["alpha"], o["phi0"], o["phi1"]): o for o in res}
    cand = [k for k, o in idx.items() if o["a"] > 0 and o["fin"] == "AF"]
    rng = np.random.default_rng(7)
    sel = [cand[i] for i in rng.choice(len(cand), 40, replace=False)]
    tareas = [("cuadratico", 0.1, 0.1, 0.32989)] + sel
    with ProcessPoolExecutor(max_workers=8) as ex:
        out = list(ex.map(tarea, tareas))
    for o in out:
        print(f"{o['tipo']:10s} alpha={o['alpha']:+.3g} phi0={o['phi0']:.3g} phi1={o['phi1']:+.3g} a r0={o['a']:.3g}: "
              f"atractiva {o['fin_atr']}, c_min = {o['atr']['c_min']:.3f} ({o['atr']['puntos_negativos']} puntos lambda<0); "
              f"repulsiva {o['fin_rep']} en l_s={o['l_s']:.3g}, c_min = {o['rep']['c_min']:.3f} "
              f"({o['rep']['puntos_negativos']} puntos lambda<0)")
    hip = [o for o in out if o["atr"]["puntos_negativos"] == 0 and o["rep"]["puntos_negativos"] == 0]
    print(f"continuaciones lisas hiperbolicas en todo (-l_s, inf): {len(hip)} de {len(out)}")
    json.dump(out, open("resultados/hiperbolicidad_lisa.json", "w", encoding="utf-8"), indent=1)
    print(f"tiempo: {time.time() - t0:.0f} s")
