# -*- coding: utf-8 -*-
"""
Verificacion del fondo de Canate et al. y de las ecuaciones con Maxwell (canate.py):
  1. f'(phi) y f''(phi) publicados (su ec. 34) frente a las derivadas simbolicas de f (ec. 33);
  2. la solucion exacta cumple el sistema estatico de fase2b + Maxwell: derivadas segundas del sistema frente a las
     exactas y ligadura ll, a los dos lados de la garganta, con Qe^2 = Q^2 + Qs^2; la ligadura en cada punto fija Qe^2;
  3. control de Maxwell: Reissner-Nordstrom (F = 0, phi = 0) integrado con segundas_Q reproduce r' = sqrt(f(r)) y N;
  4. el modo gauge g(l) = (Phi'' - w e^{-2Phi}, r', r'', phi', phi'') es solucion exacta de derivadas_Q y cumple
     ligaduras_Q (w real y complejo, l de los dos signos);
  5. una solucion lineal generica conserva las ligaduras (identidad de Bianchi con Maxwell).
Salida: resultados/verificar_canate.json.   Uso: python verificar_canate.py
"""
import os
import sys
import json
import math
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import canate as C          # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
out = {}

# 1
d1, d2 = C.comprobar_derivadas_publicadas()
print(f"1. ec. (34) frente a d/dphi de (33): f' -> {d1},  f'' -> {d2}")
out["derivadas_publicadas"] = [str(d1), str(d2)]

# 2
res = []
for Qs in (0.1, 0.5, 1.0, 3.0):
    f = C.FondoCanate(Qs)
    emax, lmax_, qe = 0.0, 0.0, []
    for l in np.r_[-np.logspace(-3, 3, 25), 0.0, np.logspace(-3, 3, 25)]:
        y = f.sol.sol(l)
        x, det = C.segundas_Q(y, f.Fd, f.Qe2)
        ex = f.exactas_segundas(l)
        esc = np.max(np.abs(ex)) + 1e-300
        emax = max(emax, np.max(np.abs(x - ex)) / max(esc, abs(ex[0])))
        lmax_ = max(lmax_, abs(C.est_ligadura_Q(y, f.Fd, f.Qe2)) * y[0] ** 4 / f.Qe2)
        r, s, Phi, P, p, u = y[:6]
        qe.append(-r**4 * C.H.est_ligadura(f.Fd(p)[1], r, s, P, u))
    qe = np.array(qe)
    res.append(dict(Qs=Qs, err_segundas=float(emax), lig=float(lmax_), Qe2_min=float(qe.min()), Qe2_max=float(qe.max()),
                    Qe2_teo=f.Qe2))
    print(f"2. Qs/Q = {Qs}: |segundas - exactas| rel <= {emax:.1e};  |ligadura ll| r^4/Qe^2 <= {lmax_:.1e};  "
          f"Qe^2 de la ligadura en [{qe.min():.12f}, {qe.max():.12f}] (Q^2 + Qs^2 = {f.Qe2})")
out["solucion_exacta"] = res

# 3: Reissner-Nordstrom
M, Qe2 = 1.0, 0.64
Fd0 = lambda x: (0.0, 0.0, 0.0, 0.0)          # noqa: E731
r_i = 5.0
fr = 1 - 2 * M / r_i + Qe2 / r_i**2
y0 = [r_i, math.sqrt(fr), 0.0, (M / r_i**2 - Qe2 / r_i**3) / math.sqrt(fr), 0.0, 0.0]
s = solve_ivp(lambda l, y: [y[1], *C.segundas_Q(y, Fd0, Qe2)[0][:1], y[3], C.segundas_Q(y, Fd0, Qe2)[0][1], y[5],
                            C.segundas_Q(y, Fd0, Qe2)[0][2]], (0, 40), y0, method="DOP853", rtol=1e-12, atol=1e-14,
              dense_output=True)
err_rn = 0.0
for l in np.linspace(0, 40, 50):
    r, rp, Phi, P = s.sol(l)[:4]
    fr = 1 - 2 * M / r + Qe2 / r**2
    N_rn = math.sqrt(fr / (1 - 2 * M / r_i + Qe2 / r_i**2))
    err_rn = max(err_rn, abs(rp - math.sqrt(fr)), abs(math.exp(Phi) - N_rn))
print(f"3. Reissner-Nordstrom (M = 1, Qe^2 = 0.64) con segundas_Q: max |r' - sqrt(f)|, |N - N_RN| = {err_rn:.1e}")
out["reissner_nordstrom"] = float(err_rn)

# 4, 5
f = C.FondoCanate(1.0)


def g_en(l, w):
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = C.fondo_en_Q(f, l)
    return np.array([P2 - w * e2, r1, r2, p1, p2], dtype=complex if isinstance(w, complex) else float)


gauge = []
for w in (-2.5, 7.0, complex(-3.0, 1.5)):
    for l1, l2 in ((-3.0, -0.2), (-0.5, 0.7), (0.3, 6.0)):
        tipo = complex if isinstance(w, complex) else float
        ss = solve_ivp(lambda x, yy: C.derivadas_Q(f, x, yy, w), (l1, l2), g_en(l1, w).astype(tipo), method="DOP853",
                       rtol=1e-12, atol=1e-14)
        err = np.max(np.abs(ss.y[:, -1] - g_en(l2, w))) / np.max(np.abs(g_en(l2, w)))
        lig = max(np.max(np.abs(C.ligaduras_Q(f, l, g_en(l, w), w))) for l in np.linspace(l1, l2, 7))
        gauge.append(dict(w=str(w), l1=l1, l2=l2, err=float(err), lig=float(lig)))
        print(f"4. modo gauge, w = {w}, l {l1} -> {l2}: |y - g|/|g| = {err:.1e}, ligaduras <= {lig:.1e}")
out["gauge"] = gauge

rng = np.random.default_rng(5)
prop = []
for w in (-4.0, complex(-1.0, 2.0)):
    tipo = complex if isinstance(w, complex) else float
    Cm = np.array([[C.ligaduras_Q(f, 0.4, e, w)[i] for e in np.eye(5)] for i in range(2)], dtype=tipo)
    _, _, Vh = np.linalg.svd(Cm)
    y0 = (Vh[2:].conj().T @ rng.normal(size=3)).astype(tipo)          # dato en el nucleo de las ligaduras
    ss = solve_ivp(lambda x, yy: C.derivadas_Q(f, x, yy, w), (0.4, 5.0), y0, method="DOP853", rtol=1e-12, atol=1e-14,
                   dense_output=True)
    lig = max(np.max(np.abs(C.ligaduras_Q(f, l, ss.sol(l), w))) / np.max(np.abs(ss.sol(l))) for l in np.linspace(0.4, 5, 30))
    prop.append(dict(w=str(w), lig=float(lig)))
    print(f"5. solucion generica, w = {w}: ligaduras propagadas <= {lig:.1e} (relativo)")
out["propagacion"] = prop
os.makedirs("resultados", exist_ok=True)
json.dump(out, open("resultados/verificar_canate.json", "w", encoding="utf-8"), indent=1)
