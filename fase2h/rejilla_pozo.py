# -*- coding: utf-8 -*-
"""
Fase 2h: rejilla de pozos profundos en la familia fondo_nu.FondoNu (criba de la Fase 2g con la superficie fija y con la
superficie libre, y margen exterior r sig(ls)).  El mejor punto con la superficie fija es el candidato del informe
(resultados/rejilla_pozo.json).
Variables: w_t en {1, 2}, nu0 en {1.5, 2.5, 4}, Ns - 1 en {1, 3, 6}, (rho~+p~)(0) en {0.03, 0.1}, lam en {1, 0.5};
el resto como el fondo de prueba de fondo_nu (b = 3, n = 2, d0 = 0.5, wd = 1, lT = 30, phi1 = 0.5, f_fl = 1, m = 1).
Uso: python rejilla_pozo.py   -> resultados/rejilla_pozo_tabla.json y resultados/rejilla_pozo.json
"""
import os
import sys
import json
import math
import time
import itertools
import warnings
from multiprocessing import Pool

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)


def evaluar(args):
    warnings.filterwarnings("ignore")
    sys.path.insert(0, os.path.join(AQUI, "..", "fase2g"))
    import optimiza_nu as O
    import criba as C
    wt, nu0, Nsm1, rp0, lam = args
    L = math.log
    v = [L(wt), L(nu0), 1.0, 0.5, L(3.0), 2, 0.5, 0.0, L(lam), L(30.0), 1.0, L(rp0), L(Nsm1)]
    x = O.parametros(v)
    try:
        import fondo_nu as FN
        F = FN.FondoNu(**x)
        if not F.ok or F.ls <= 0 or F.reentra:
            return dict(args=args, x=x, motivo=F.motivo or "sin estrella")
        fija = C.criba_normalizada(F, n=80)
        libre = C.criba_normalizada(F, n=80, libre=True)
        return dict(args=args, x=x, fija=fija, libre=libre, ext=F.margen_exterior(), M=F.M_adm, ls=F.ls, N0=F.N0)
    except Exception as e:                      # noqa: BLE001
        return dict(args=args, x=x, motivo=repr(e)[:80])


if __name__ == "__main__":
    t0 = time.time()
    rejilla = list(itertools.product((1.0, 2.0), (1.5, 2.5, 4.0), (1.0, 3.0, 6.0), (0.03, 0.1), (1.0, 0.5)))
    with Pool(8) as p:
        out = p.map(evaluar, rejilla, chunksize=1)
    ok = [o for o in out if "fija" in o]
    for o in ok:
        o["margen_fija"] = min(o["fija"]["margen"], o["ext"])
        o["margen_libre"] = min(o["libre"]["margen"], o["ext"])
    ok.sort(key=lambda o: -o["margen_fija"])
    for o in ok[:12]:
        a = o["args"]
        print(f"wt={a[0]} nu0={a[1]} Ns-1={a[2]} rp0={a[3]} lam={a[4]}: margen fija {o['margen_fija']:+.4f}, "
              f"libre {o['margen_libre']:+.4f} (Epar libre {o['libre']['Epar']:+.3f}), ext {o['ext']:+.3f} | "
              f"M={o['M']:.2f} ls={o['ls']:.2f} N0={o['N0']:.3f}")
    print(f"validos {len(ok)}/{len(rejilla)}; con margen fija > 0: {sum(o['margen_fija'] > 0 for o in ok)}; "
          f"con margen libre > 0: {sum(o['margen_libre'] > 0 for o in ok)}  [{time.time() - t0:.0f} s]")
    json.dump(out, open(os.path.join(AQUI, "resultados", "rejilla_pozo_tabla.json"), "w"), indent=1, default=float)
    json.dump(dict(x=ok[0]["x"]), open(os.path.join(AQUI, "resultados", "rejilla_pozo.json"), "w"), indent=1)
