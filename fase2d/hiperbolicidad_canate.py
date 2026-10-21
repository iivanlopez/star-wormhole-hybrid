# -*- coding: utf-8 -*-
"""
Hiperbolicidad de las perturbaciones radiales del agujero de Canate et al. a lo largo de toda la linea.
En el limite de alta frecuencia (w, k -> inf con w/k^2 = lambda fijo) el sistema lineal da k^2 M x = w N x, con M la
parte principal (coeficientes de d'', a'', q'') y N los coeficientes de w.  Se escribe
    P(lambda) = det(M - lambda N) = c0 + c1 lambda + c2 lambda^2       (N tiene la columna del lapso nula),
y lambda son sus raices finitas (en este fondo c2 = 0 y hay una sola: lambda = -c0/c1).  Con c0 y c1 funciones suaves
de l se localizan sin ambiguedad:
  ceros de lambda  = ceros de c0 = det M (se anula el coeficiente de gradiente; es tambien el det del sistema estatico),
  polos de lambda  = ceros de c1 (se anula el coeficiente cinetico, el de w).
lambda > 0 en todo punto <=> hiperbolico (la luz tiene lambda = N^2 = 1 en este fondo ultraestatico);
lambda < 0 en un intervalo => crecimiento sin cota e^{sqrt(|lambda|) k t}: problema mal planteado (Hadamard).
El fondo es simetrico bajo l -> -l, phi -> -phi (se comprueba), asi que basta l >= 0.
Salida: resultados/hiperbolicidad_canate.json.   Uso: python hiperbolicidad_canate.py
"""
import os
import sys
import json
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import canate as C          # noqa: E402

LS = np.unique(np.r_[np.linspace(0.0, 2.0, 801), np.linspace(2.0, 50.0, 961), np.logspace(np.log10(50.0), 6, 400)])


def coeficientes(f, l):
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = C.fondo_en_Q(f, float(l))
    M = np.array(C.H.lin_principal(F1, r, r1, P1, p1))
    N = np.zeros((3, 3))
    for j, (d0, q0) in ((0, (1.0, 0.0)), (2, (0.0, 1.0))):
        N[:, j] = np.array(C.H.lin_resto(F1, F2, F3, 1.0, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0)) \
            - np.array(C.H.lin_resto(F1, F2, F3, 0.0, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0))
    p0, pp, pm = (np.linalg.det(M - x * N) for x in (0.0, 1.0, -1.0))
    return p0, (pp - pm) / 2, (pp + pm) / 2 - p0          # c0, c1, c2


def tarea(Qs):
    f = C.FondoCanate(Qs)
    cs = np.array([coeficientes(f, l) for l in LS])
    c0, c1, c2 = cs.T
    escala = np.abs(cs).max(axis=1)
    lam = -c0 / c1
    sim = max(abs(coeficientes(f, -l)[0] / coeficientes(f, -l)[1] - c0i / c1i) / max(1.0, abs(c0i / c1i))
              for l, c0i, c1i in ((0.05, *coeficientes(f, 0.05)[:2]), (0.7, *coeficientes(f, 0.7)[:2]),
                                  (3.0, *coeficientes(f, 3.0)[:2]), (40.0, *coeficientes(f, 40.0)[:2])))
    cambios = []
    for nombre, c in (("cero", c0), ("polo", c1)):
        for i in range(len(LS) - 1):
            if c[i] * c[i + 1] < 0:
                k = 0 if nombre == "cero" else 1
                x = brentq(lambda l: coeficientes(f, l)[k], LS[i], LS[i + 1], xtol=1e-13)
                cambios.append((nombre, float(x)))
    cambios.sort(key=lambda t: t[1])
    neg = LS[lam < 0]
    return dict(Qs=Qs, c2_max_rel=float(np.max(np.abs(c2) / escala)), lambda_inf=float(lam[-1]),
                lambda_min=float(np.min(lam)), l_lambda_min=float(LS[np.argmin(lam)]), cambios=cambios,
                l_neg=[float(neg.min()), float(neg.max())] if len(neg) else None, simetria=float(sim),
                muestra=[[float(l), float(x)] for l, x in zip(LS[::20], lam[::20])])


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    t0 = time.time()
    QS = [0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0]
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(tarea, QS))
    for o in res:
        print(f"Qs/Q = {o['Qs']:6g}: lambda(l=1e6) = {o['lambda_inf']:.6f}; min lambda = {o['lambda_min']:+.4g} en "
              f"l = {o['l_lambda_min']:.4g}; lambda < 0 en l/Q de {o['l_neg']}; cambios: "
              f"{[(t, round(x, 6)) for t, x in o['cambios']]}; |c2|/escala <= {o['c2_max_rel']:.1e}; "
              f"simetria {o['simetria']:.1e}")
    mal = sum(o["l_neg"] is not None for o in res)
    print(f"fondos con lambda < 0 en algun intervalo: {mal} de {len(res)}")
    json.dump(res, open("resultados/hiperbolicidad_canate.json", "w", encoding="utf-8"), indent=1)
    print(f"tiempo: {time.time() - t0:.0f} s")
