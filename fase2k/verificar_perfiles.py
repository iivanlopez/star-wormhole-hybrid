# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): comprobacion del simbolo principal en forma de perfiles (simbolo_perfiles.py, usado en
criba_masa.py) frente a la ruta con las funciones de la teoria (simbolo_general.py, usado en criba_cubica.py), en un
fondo de Ellis con flujo (r = sqrt(rho^2 + 1)): s_j = P_j/D por las dos rutas.
Uso: python verificar_perfiles.py
"""
import os
import sys
import json
import numpy as np
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
sys.argv = sys.argv[:1]
import criba_cubica as CC        # noqa: E402
import criba_masa as CM          # noqa: E402

prm = dict(w=2.0, a=50.0, wq=1.0, An=1.0, wn=2.0, w2=0.2, lam=0.05)
CC.FUN = CC.construir_funciones()
fun, fP, fb2, fb2h = CC.FUN
pv = [prm[k] for k in ("w", "a", "wq", "An", "wn", "w2", "lam")]
L = CC.L
kw = dict(dense_output=True, rtol=1e-11, atol=1e-13, max_step=0.02)
sh = solve_ivp(lambda s, y: [y[1], fb2h(s, y[0], y[1], *pv, 1.0)], (0, L), [1.0, 0.0], **kw)
spp = solve_ivp(lambda s, y: [y[1], fb2(s, y[0], y[1], *pv, 1.0)], (0, L), [0.0, 0.0], **kw)
xa = np.linspace(L / 2, L, 200)
Mx = np.vstack([xa ** 2, np.ones_like(xa), 1 / xa]).T
ah = np.linalg.lstsq(Mx, sh.sol(xa)[0], rcond=None)[0]
ap = np.linalg.lstsq(Mx, spp.sol(xa)[0], rcond=None)[0]
cc = -ap[0] / ah[0]
binf = cc * ah[1] + ap[1]
V0 = 1 / np.sqrt(binf)
xs = np.array([0.05, 0.2, 0.4, 0.7, 1.0, 2.0, 5.0, 12.0])
y = (cc * sh.sol(xs) + spp.sol(xs)) / binf
vals = {u: np.broadcast_to(np.asarray(fun[u](xs, y[0], y[1], *pv, V0), float), xs.shape) for u in CC.usados}
P = [np.asarray(f(*[vals[u] for u in CC.usados]), float) for f in fP]
N0, R0, R1, V, a3N, a4N, a4 = (vals[k] for k in ("N_0", "R_0", "R_1", "V_0", "a3N_0", "a4N_0", "a4_0"))
Dn = -R0 ** 2 * (N0 ** 2 * R0 * a3N - 4 * N0 * R1 * V * a4N + 4 * R1 * V * a4) ** 2 / N0 ** 4
s_teoria = [p / Dn for p in P]
# ruta de perfiles
r = np.sqrt(xs ** 2 + 1)
z1, z2 = xs / prm["w"], xs / prm["w2"]
e1, e2 = np.exp(-z1 ** 2 / 2), np.exp(-z2 ** 2 / 2)
w1, w2, lam = prm["w"], prm["w2"], prm["lam"]
v = V0 * (z1 * e1 + lam * z2 * e2)
v1 = V0 * ((1 - z1 ** 2) * e1 / w1 + lam * (1 - z2 ** 2) * e2 / w2)
v2 = V0 * ((z1 ** 3 - 3 * z1) * e1 / w1 ** 2 + lam * (z2 ** 3 - 3 * z2) * e2 / w2 ** 2)
g = np.exp(-xs ** 2 / prm["wq"] ** 2)
u = 1 + xs ** 2 / prm["wn"] ** 2
val = dict(r0=r, r1=xs / r, r2=1 / r ** 3, r3=-3 * xs / r ** 5, v0=v, v1=v1, v2=v2,
           Q0=-2 * (1 + prm["a"] * g), Q1=-2 * prm["a"] * g * (-2 * xs / prm["wq"] ** 2),
           Nc1=-prm["An"] * 2 * xs / prm["wn"] ** 2 / u ** 2,
           Nc2=-prm["An"] * (2 / prm["wn"] ** 2 / u ** 2 - 8 * xs ** 2 / prm["wn"] ** 4 / u ** 3),
           beta0=y[0], beta1=y[1])
s_perf = [np.asarray(f(*[val[str(a_)] for a_ in CM.ARG]), float) for f in CM.fS]
RES = {}
for j, nom in enumerate(("s2", "s1", "s0")):
    rel = np.max(np.abs(s_teoria[j] - s_perf[j]) / (np.abs(s_teoria[j]) + 1e-12))
    RES[nom] = float(rel)
    print(f"{nom}: teoria {np.array2string(s_teoria[j], precision=4)}\n    perfiles {np.array2string(s_perf[j], precision=4)}"
          f"\n    diferencia relativa maxima {rel:.1e}", flush=True)
json.dump(RES, open(os.path.join(AQUI, "resultados", "verificar_perfiles.json"), "w"), indent=1)
