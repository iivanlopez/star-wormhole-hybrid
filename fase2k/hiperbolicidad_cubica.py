# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): hiperbolicidad local del monopolo a lo largo de un fondo con galileon cubico y flujo
(monopolo_cubico.py).  En cada rho, con los coeficientes congelados, M(om, k) = sum C[q,f,m,n] (-i om)^m (i k)^n
(onda e^{i(k rho - om t)}); las raices om(k) de det M = 0 (problema de autovalores cuadratico 3x3 linealizado) con
k grande dan las velocidades caracteristicas radiales c = om/k en el corte phi = cte.  Velocidades complejas =
problema mal planteado (inestabilidad de gradiente).  Se compara con el cono gravitatorio radial c = -v +- Nc y se
evalua el signo de Krein de cada rama (energia en el corte tau): signo de  i conj(x)^T dM/d(om) x  sobre el vector
nulo x (dos ramas del mismo signo y slices espaciales = sin fantasma en este corte).
Uso: python hiperbolicidad_cubica.py [w a wq An wn]
"""
import os
import sys
import json
import numpy as np
import scipy.linalg as sla

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import monopolo_cubico as MC     # noqa: E402

CAMPOS = ("nu", "w", "d")


def matrices(C, val, i, k):
    """M(om) = M0 + om M1 + om^2 M2 en el punto i con numero de onda k."""
    Ms = [np.zeros((3, 3), complex) for _ in range(3)]
    for (q, f, m, n), co in C.items():
        c = val[(q, f, m, n)][i]
        Ms[m][CAMPOS.index(q), CAMPOS.index(f)] += c * (-1j) ** m * (1j * k) ** n
    return Ms


def raices(Ms):
    M0, M1, M2 = Ms
    Z, I = np.zeros((3, 3)), np.eye(3)
    A = np.block([[Z, I], [-M0, -M1]])
    B = np.block([[I, Z], [Z, M2]])
    ev, vr = sla.eig(A, B)
    ok = np.isfinite(ev) & (np.abs(ev) < 1e12)
    return ev[ok], vr[:3, ok]


def analiza(params, xs, ks=(2e3, 4e3)):
    op = MC.operador()
    C = op["C"]
    fd = MC.construir(**params)
    nombres = sorted({str(s) for s in op["usados"]})
    bgv = MC.evaluar(fd, nombres, xs)
    import sympy as sp
    val = {}
    for key, co in C.items():
        f = sp.lambdify([sp.Symbol(n) for n in nombres], co, "numpy")
        val[key] = np.broadcast_to(np.asarray(f(*[bgv[n] for n in nombres]), dtype=complex), xs.shape)
    filas = []
    for i, x in enumerate(xs):
        vel = []
        for k in ks:
            ev, _ = raices(matrices(C, val, i, k))
            vel.append(np.sort_complex(ev / k))
        # ramas que escalan con k (velocidades): las que coinciden entre los dos k
        c1 = vel[0]
        c2 = vel[1]
        prin = [c for c in c1 if np.min(np.abs(c2 - c)) < 1e-3 * (1 + abs(c))]
        filas.append(dict(rho=float(x), velocidades=[complex(c) for c in prin],
                          cono=[float(-bgv["V_0"][i] - bgv["N_0"][i]), float(-bgv["V_0"][i] + bgv["N_0"][i])],
                          G4=float(bgv["b4_0"][i]), A4=float(bgv["a4_0"][i])))
    return fd, filas


if __name__ == "__main__":
    args = [float(a) for a in sys.argv[1:]]
    claves = ("w", "a", "wq", "An", "wn", "w2", "lam")
    params = dict(zip(claves, args)) if args else dict(w=2.0, a=20.0, wq=1.0, An=1.0, wn=2.0)
    if len(args) > 5:
        params["lapso"] = "alg"
    nombre = "_".join(f"{x:g}" for x in args) if args else "base"
    xs = np.concatenate([-np.geomspace(20, 0.01, 40), np.geomspace(0.01, 20, 40)])
    fd, filas = analiza(params, xs)
    print(f"parametros {params}, V0 = {fd['V0n']:.4f}", flush=True)
    malos = 0
    for f in filas:
        vel = f["velocidades"]
        imag = max([abs(c.imag) for c in vel] + [0])
        if imag > 1e-2:
            malos += 1
        print(f"rho = {f['rho']:+9.4f}  G4 = {f['G4']:.4f}  A4 = {f['A4']:+.4f}  cono = [{f['cono'][0]:+.3f}, "
              f"{f['cono'][1]:+.3f}]  c = " + ", ".join(f"{c.real:+.4f}{c.imag:+.4f}i" for c in vel), flush=True)
    print(f"puntos con velocidades complejas (|Im c| > 0.01): {malos} de {len(filas)}")
    json.dump(dict(params=params, V0=fd["V0n"],
                   filas=[dict(f, velocidades=[[c.real, c.imag] for c in f["velocidades"]]) for f in filas]),
              open(os.path.join(AQUI, "resultados", f"hiperbolicidad_cubica_{nombre}.json"), "w"), indent=1)
