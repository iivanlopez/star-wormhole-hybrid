# -*- coding: utf-8 -*-
"""
Busqueda de gargantas masivas estables en la familia de reconstruccion.Garganta (Fase 2f).

evalua(x) devuelve un dict con las condiciones del monopolo (K > 0, P > 0 en la malla, en la zona lejana hasta
10 veces su extremo y en la cola asintotica mas alla, velocidad fisica c^2 = P/(K N^2) <= 1)
y los angulos de Prufer (en l = L) de las soluciones de w = 0 regulares en la garganta, par e impar; el numero de
taquiones es floor(th/pi) (par: th parte de pi/2; impar: de 0).  objetivo(x) es una medida continua que se anula
solo si no hay taquiones, fantasmas, inestabilidades de gradiente ni superluminalidad en el monopolo.
Parametros x: dict con las claves de Garganta.
"""
import math
import numpy as np

import reconstruccion as R
import espectro as E

LMIN, LMAX_MALLA, NMALLA = 1e-3, 400.0, 900


def construir(x):
    n = int(x.get("n", 2))
    G40 = 2 * n / (2 * n - 1)
    b = x["b"]
    kw = dict(phi1=x.get("phi1", 1.0), b=b, n=n, d0=x["frac"] * G40, wd=x["wdb"] * b, M=x["M"], a=x["a"],
              Lmax=4000.0, wK=x["wK"], wP=x["wPK"] * x["wK"], cm2=x["cm2"], cf2=x.get("cf2", 0.0),
              wf2=x.get("wf2", 1.0), lf2=x.get("lf2", 0.0))
    return R.Garganta(**kw)


def evalua(x, malla=None, impar=False):
    out = dict(x=dict(x))
    try:
        g = construir(x)
    except Exception as e:                                  # noqa: BLE001
        out["error"] = repr(e)[:120]
        return out
    ls = E.malla(LMIN, LMAX_MALLA, NMALLA) if malla is None else malla
    K, P, U, N = E.coeficientes(g, "monopolo", ls)
    # zona lejana (l hasta 10 L): K y P deben seguir positivos (una W que rebasa 3M fuera de la malla da K < 0)
    lf = np.geomspace(ls[-1], 10 * ls[-1], 200)
    Kf, Pf, _, Nf = E.coeficientes(g, "monopolo", lf)
    c2 = np.concatenate([P / (K * N ** 2), Pf / (Kf * Nf ** 2)])
    out.update(minK=float(min(K.min(), Kf.min())), minP=float(min(P.min(), Pf.min())), maxc2=float(np.max(c2)),
               N0=float(g.N0), A_K=float(g.A_K))
    # cola (l > Lmax): signo de y(ln l) (reconstruccion.cola_K), con la dispersion de beta como margen; el desarrollo
    # exige l >> r0/N0
    cola = g.cola_K()
    out.update(cola_ymin=cola["ymin"], cola_disp=cola["dispersion_beta"], cola_Lstar=cola["Lstar"],
               cola_ok=bool(g.Lmax >= 100 * g.r0 / g.N0 and cola["ymin"] > cola["dispersion_beta"]))
    if out["minK"] > 0 and out["minP"] > 0 and out["cola_ok"]:
        sl = E.SL(ls, K, P, U)
        for par in ("par", "impar"):
            th0, thL, _ = sl.prufer(0.0, par)
            nz = int(math.floor(thL / math.pi)) - int(math.floor(th0 / math.pi))
            out[f"n_{par}"] = nz
            out[f"w_{par}"] = [float(w) for w in sl.autovalores_negativos(par, nz, tol=1e-6)] if nz else []
    if impar and "n_par" in out:
        for LL in (2, 3):
            Ko, Po, Uo, No = E.coeficientes(g, "impar", ls, LL=LL)
            out[f"impar{LL}"] = dict(minK=float(Ko.min()), minP=float(Po.min()), maxc2=float(np.max(Po / (Ko * No ** 2))))
            if Ko.min() > 0 and Po.min() > 0:
                slo = E.SL(ls, Ko, Po, Uo)
                out[f"impar{LL}"]["nodos"] = (slo.ceros("par")[0], slo.ceros("impar")[0])
    return out


def objetivo(r):
    if "error" in r:
        return 1e3
    pen = 0.0
    if r["minK"] <= 0 or r["minP"] <= 0 or not r["cola_ok"]:
        return 100.0 + 10 * (max(0, -r["minK"]) + max(0, -r["minP"])) / (1e-3 + abs(r["minK"]) + abs(r["minP"]))
    pen += 5 * max(0.0, r["maxc2"] - 1.0)
    # taquiones: suma de sqrt|w| (tasas de crecimiento) en ambas paridades
    pen += sum(math.sqrt(-w) for w in r["w_par"] + r["w_impar"])
    return pen
