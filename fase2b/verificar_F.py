# -*- coding: utf-8 -*-
"""
Verificacion del codigo de acoplamiento general F(phi) (generar_F.py -> egb_F_gen.py, fondo_F.py, estabilidad_F.py):
  1. cruce con el codigo especifico de KKK (F = alpha e^{-phi}) en puntos aleatorios: parte principal, resto y
     ligaduras de las perturbaciones; segundas derivadas y ligadura del fondo;
  2. el mismo fondo de referencia integrado con fondo_F y con kkk;
  3. test de familia (w = 0): la derivada de la familia de fondos respecto de phi'(0+) cumple las ecuaciones
     linealizadas (acoplamientos cuadratico y gaussiano, una garganta y un ecuador);
  4. redundancia de Bianchi en l = 0: la matriz de las cinco condiciones tiene un solo valor singular nulo.
Uso: python verificar_F.py
"""
import sys
import math
import numpy as np
import egb_F_gen as H
import egb_lineal_gen as L
import kkk
from fondo_F import FondoF, segundas, acoplamiento
from estabilidad_F import CapaF, matriz_garganta

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
rng = np.random.default_rng(1)


def rel(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


# ------------------------------------------------------------------ 1. cruce con el codigo de KKK
peor = {"principal": 0.0, "resto": 0.0, "ligaduras": 0.0, "segundas": 0.0, "ligadura_fondo": 0.0}
for _ in range(200):
    al, ga = rng.uniform(0.01, 0.2), 1.0
    r, r1, r2, P1, P2, p, p1, p2 = rng.uniform(0.5, 2.0), *rng.uniform(-1, 1, 7)
    e2 = rng.uniform(0.1, 2.0)
    w = rng.uniform(-50, 50)
    a1, d0, d1, q0, q1 = rng.uniform(-1, 1, 5)
    F0, F1, F2, F3 = acoplamiento("exp", al, gamma=ga)(p)
    peor["principal"] = max(peor["principal"], rel(H.lin_principal(F1, r, r1, P1, p1),
                                                   L.principal(al, ga, r, r1, r2, P1, P2, p, p1, p2, e2)))
    peor["resto"] = max(peor["resto"], rel(H.lin_resto(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1),
                                           L.resto(al, ga, w, r, r1, r2, P1, P2, p, p1, p2, e2, 0.0, a1, d0, d1, q0, q1)))
    peor["ligaduras"] = max(peor["ligaduras"],
                            rel(H.lin_ligaduras(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1),
                                L.ligaduras(al, ga, w, r, r1, r2, P1, P2, p, p1, p2, e2, 0.0, a1, d0, d1, q0, q1)))
    y = np.array([r, r1, 0.0, P1, p, p1, 0.0])
    xg, _ = segundas(y, acoplamiento("exp", al, gamma=ga))
    xk = kkk.segundas(y[:6], al, ga)[:3]
    peor["segundas"] = max(peor["segundas"], rel(xg, xk))
    peor["ligadura_fondo"] = max(peor["ligadura_fondo"],
                                 abs(H.est_ligadura(F1, r, r1, P1, p1) - kkk.ligadura(y[:6], al, ga))
                                 / max(abs(kkk.ligadura(y[:6], al, ga)), 1e-12))
print("1. cruce con el codigo de KKK (200 puntos aleatorios), diferencia relativa maxima:")
for k, v in peor.items():
    print(f"     {k:15s} {v:.1e}")

# ------------------------------------------------------------------ 2. fondo de referencia por los dos caminos
fg = FondoF("exp", 0.03, 0.0, -0.7, rtol=1e-12).resolver()
fk = kkk.FondoKKK(alpha=0.03, phi1=-0.7, rtol=1e-12).resolver()
print("2. fondo de referencia (alpha = 0.03, phi'(0+) = -0.7), fondo_F frente a kkk:")
for nom in ("f0", "Phi_inf", "phi_inf", "M", "D"):
    a, b = getattr(fg, nom), getattr(fk, nom)
    print(f"     {nom:8s} {a:+.12f}  {b:+.12f}   dif. rel. {abs(a - b) / abs(b):.1e}")

# ------------------------------------------------------------------ 3. test de familia con w = 0
print("3. test de familia (derivada respecto de phi'(0+), w = 0, diferencias centradas con extrapolacion de Richardson):")
print("   residuo relativo de las ecuaciones linealizadas")
gauss = next((("gauss", al, 0.4, u, {"gamma": 1.0}) for al in (0.1, 0.05) for u in (-0.6, -1.0, -0.3, -2.0)
              if FondoF("gauss", al, 0.4, u, kw={"gamma": 1.0}).resolver().regular), None)
CASOS = [("cuadratico", 0.09, 0.5, 0.19566, {}), ("cuadratico", 0.12, 0.03, 5.73855, {})] + ([gauss] if gauss else [])
for tipo, al, p0, u, kw in CASOS:
    h = 1e-3 * abs(u)
    fs = [FondoF(tipo, al, p0, u + s * h, rtol=1e-13, kw=kw).resolver() for s in (-2, -1, 0, 1, 2)]
    f = fs[2]
    lmax = 0.8 * min(g.l_end for g in fs)
    peor_ec, peor_lig = 0.0, 0.0
    for l in np.concatenate([np.linspace(1e-3, 0.5, 15), np.geomspace(0.5, min(lmax, 50.0), 15)]):
        ys = [g.sol.sol(l) for g in fs]
        xs = [segundas(y, g.Fd)[0] for y, g in zip(ys, fs)]
        dy = (8 * (ys[3] - ys[1]) - (ys[4] - ys[0])) / (12 * h)
        dx = (8 * (xs[3] - xs[1]) - (xs[4] - xs[0])) / (12 * h)
        r, r1, Phi, P1, p, p1 = ys[2][:6]
        r2, P2, p2 = xs[2]
        F0, F1, F2, F3 = f.Fd(p)
        e2 = math.exp(-2 * Phi)
        d0, d1, a1, q0, q1 = dy[0], dy[1], dy[3], dy[4], dy[5]
        d2, a2, q2 = dx
        M = np.array(H.lin_principal(F1, r, r1, P1, p1))
        pr = M @ np.array([d2, a2, q2])
        re = np.array(H.lin_resto(F1, F2, F3, 0.0, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1))
        peor_ec = max(peor_ec, float(np.max(np.abs(pr + re)) / max(np.max(np.abs(pr)), np.max(np.abs(re)))))
        lig = H.lin_ligaduras(F1, F2, F3, 0.0, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1)[1]
        escala = max(abs(H.lin_ligaduras(F1, F2, F3, 0.0, r, r1, r2, P1, P2, p1, p2, e2, *v)[1])
                     for v in np.eye(5) * np.array([a1, d0, d1, q0, q1]))
        peor_lig = max(peor_lig, abs(lig) / max(escala, 1e-300))
    print(f"     {tipo:10s} alpha={al} phi0={p0} phi1={u}: {'ecuador' if f.ecuador else 'garganta'}, "
          f"ecuaciones {peor_ec:.1e}, ligadura ll {peor_lig:.1e}")

# ------------------------------------------------------------------ 4. Bianchi en l = 0
print("4. redundancia de Bianchi en l = 0: valores singulares relativos de la matriz 5x5 (dos menores)")
for tipo, al, p0, u in (("cuadratico", 0.09, 0.5, 0.19566), ("cuadratico", 0.12, 0.03, 5.73855), ("exp", 0.02, 0.0, -12.34)):
    f = FondoF(tipo, al, p0, u).resolver()
    for z in (1.0, 0.0):
        c = CapaF(f, zeta=z)
        sv = []
        for w in (-1e4, -10.0, -0.1, 0.1, 10.0):
            s = np.linalg.svd(matriz_garganta(f, c, w), compute_uv=False)
            sv.append((s[-1] / s[0], s[-2] / s[0]))
        print(f"     {tipo:10s} alpha={al} phi1={u} zeta={z:g}: s5/s1 <= {max(x[0] for x in sv):.0e}, "
              f"s4/s1 >= {min(x[1] for x in sv):.0e}")
