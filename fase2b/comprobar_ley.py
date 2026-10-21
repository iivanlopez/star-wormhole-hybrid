# -*- coding: utf-8 -*-
"""
Comprueba que la raiz con C = -w/a^2 ~ 1 (zeta = 0) es un modo genuino y no un cero trivial de los menores: en una malla
fina alrededor de la raiz se miran |u_menores|, los valores singulares de la matriz 4x5 de condiciones (rango 4 si
s4/s1 no es pequeno) y el signo de D; y en la raiz, la autofuncion (decaimiento en el infinito y ligaduras).
Uso: python comprobar_ley.py
"""
import sys
import json
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from fondo_F import FondoF
from estabilidad_F import CapaF, D, vector_menores, matriz_garganta, derivadas, ligaduras

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
rob = {(o["nombre"], o["zeta"]): o for o in json.load(open("resultados/estabilidad_robusta.json", encoding="utf-8"))}
for nom in ("KKK 0.0452 1.0980", "cuadratico 0.3087", "KKK 0.0039 2.4866"):
    o = rob[(nom, 0.0)]
    f = FondoF(o["tipo"], o["alpha"], o["phi0"], o["phi1"]).resolver()
    c = CapaF(f, zeta=0.0)
    a = f.y0()[3]
    w0 = min(o["w"], key=lambda w: abs(-w / a**2 - 1))
    print(f"{nom}: a r0 = {a:.4f}, raiz w* = {w0:.6f}, C = {-w0 / a**2:.5f}")
    for x in np.linspace(0.9, 1.1, 9):
        w = w0 * x
        K = matriz_garganta(f, c, w)[[0, 1, 3, 4]]
        s = np.linalg.svd(K, compute_uv=False)
        um = np.real(vector_menores(f, c, w))
        print(f"   w/w* = {x:.3f}: D = {D(f, c, w, u=um):+.3e}  |u_m| = {np.linalg.norm(um):.3e}  s4/s1 = {s[-1] / s[0]:.2e}")
    # autofuncion en la raiz: q debe decaer como e^{-kappa l}/r en el infinito
    um = np.real(vector_menores(f, c, w0))
    um = um / np.max(np.abs(um))
    kap = math.sqrt(-w0 * math.exp(-2 * f.Phi_inf))
    l1 = min(10.0 + 30.0 / kap, 0.9 * f.l_end)
    sol = solve_ivp(lambda l, y: derivadas(f, l, y, w0), (0.0, l1), um, method="DOP853", rtol=1e-11, atol=1e-16,
                    dense_output=True)
    print(f"   autofuncion: kappa = {kap:.4g}, l1 = {l1:.4g}")
    for l in (0.0, 0.5 / a, 1.0 / a, 3.0 / a, 1.0, 0.25 * l1, 0.5 * l1, 0.75 * l1):
        y = sol.sol(l)
        r = f.sol.sol(l)[0]
        lig = ligaduras(f, l, y, w0)
        print(f"     l = {l:9.4g}: a' = {y[0]:+.3e}  d = {y[1]:+.3e}  q = {y[3]:+.3e}  q r e^(kappa l) = "
              f"{y[3] * r * math.exp(kap * l):+.3e}  ligaduras {abs(lig[0]):.1e} {abs(lig[1]):.1e}")
