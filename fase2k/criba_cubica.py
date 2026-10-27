# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): criba de fondos con galileon cubico y flujo por la hiperbolicidad local del monopolo.

Fondo (como monopolo_cubico.py, parametros simbolicos y lambdify una sola vez):
    r = sqrt(rho^2 + 1),  v = V0 [z1 e^{-z1^2/2} + lam z2 e^{-z2^2/2}] (z_i = rho/w_i),  Q = -2 [1 + a e^{-rho^2/wq^2}],
    Nc = 1 + An/(1 + rho^2/wn^2),  beta de su EDO (par, acotada, beta_inf = 1 fija V0),  p de H = 0,  A3' de la
    ligadura de momento;  funciones de la teoria como funciones de Nc (f_N = f'/Nc').
Simbolo principal (simbolo_general.py): det M ~ k^4 P(c), P = P2 c^2 + P1 c + P0;  velocidades caracteristicas
radiales reales <=> P1^2 - 4 P2 P0 >= 0.  Con det M_aux < 0, el simbolo efectivo de d es -P/|D|: energia en el corte
tau positiva <=> P2 < 0 y P0 > 0.
Condiciones: (B) beta > 0 en todo el eje; (G) garganta sana (garganta_cubica.py); (H) hiperbolico en todo rho.
Uso: python criba_cubica.py [n]  |  python criba_cubica.py fijo w a wq An wn w2 lam
"""
import os
import sys
import json
import pickle
import time
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
X, B0, B1 = sp.symbols("x B0 B1", real=True)
PAR = sp.symbols("w a wq An wn w2 lam V0", real=True)
w, a, wq, An, wn, w2, lam, V0 = PAR
L = 300.0
sg = pickle.load(open(os.path.join(AQUI, "resultados", "simbolo_general.pkl"), "rb"))
P = sp.Poly(sg["P"], sp.Symbol("c"))
P2, P1, P0 = (P.coeff_monomial(sp.Symbol("c") ** j) for j in (2, 1, 0))
usados = sg["usados"]


def construir_funciones():
    r = sp.sqrt(X ** 2 + 1)
    z1, z2 = X / w, X / w2
    v_ = V0 * (z1 * sp.exp(-z1 ** 2 / 2) + lam * z2 * sp.exp(-z2 ** 2 / 2))
    Q = -2 * (1 + a * sp.exp(-X ** 2 / wq ** 2))
    Nc = 1 + An / (1 + X ** 2 / wn ** 2)
    D = lambda f: sp.diff(f, X)                     # noqa: E731
    K = 1 + r * D(D(r)) - D(r) ** 2
    S = D(Q * r ** 4 * D(v_ ** 2 / r ** 2)) / (4 * r ** 3)
    b2 = D(r) / r * B1 - K * B0 / r ** 2 - r * S

    def Dt(ex):
        return sp.diff(ex, X) + sp.diff(ex, B0) * B1 + sp.diff(ex, B1) * b2

    b4 = B0 / Nc
    a4 = Q * Nc / 2
    p = (2 * B0 * (D(r) ** 2 - 1) + 4 * B1 * r * D(r) + Q * D(r) * D(v_ ** 2 * r)) / r ** 2
    a2 = p / Nc
    a3p = 2 * v_ * D(Q * D(r)) / r
    N1 = D(Nc)
    ex = {"N_0": Nc, "N_1": N1, "R_0": r, "R_1": D(r), "R_2": D(D(r)), "V_0": v_, "V_1": D(v_)}
    for nom, f in (("a2", a2), ("b4", b4), ("a4", a4)):
        fN = Dt(f) / N1
        fNN = Dt(fN) / N1
        ex.update({f"{nom}_0": f, f"{nom}_1": Dt(f), f"{nom}N_0": fN, f"{nom}N_1": Dt(fN), f"{nom}NN_0": fNN})
    a3N = a3p / N1
    ex.update({"a3N_0": a3N, "a3N_1": Dt(a3N), "a3NN_0": Dt(a3N) / N1})
    falta = [u for u in usados if u not in ex]
    assert not falta, falta
    args = (X, B0, B1) + PAR
    t0 = time.time()
    fun = {u: sp.lambdify(args, ex[u], "numpy") for u in usados}
    sim = [sp.Symbol(u) for u in usados]
    fP = [sp.lambdify(sim, co, "numpy") for co in (P2, P1, P0)]
    fb2 = sp.lambdify(args, b2, "numpy")
    fb2h = sp.lambdify(args, D(r) / r * B1 - K * B0 / r ** 2, "numpy")
    print(f"funciones compiladas en {time.time() - t0:.1f} s", flush=True)
    return fun, fP, fb2, fb2h


FUN = None


def evalua(prm, xs=None):
    global FUN
    if FUN is None:
        FUN = construir_funciones()
    fun, fP, fb2, fb2h = FUN
    pv = [prm["w"], prm["a"], prm["wq"], prm["An"], prm["wn"], prm["w2"], prm["lam"]]
    kw = dict(dense_output=True, rtol=1e-10, atol=1e-13, max_step=0.05)
    sh = solve_ivp(lambda x, y: [y[1], fb2h(x, y[0], y[1], *pv, 1.0)], (0, L), [1.0, 0.0], **kw)
    spp = solve_ivp(lambda x, y: [y[1], fb2(x, y[0], y[1], *pv, 1.0)], (0, L), [0.0, 0.0], **kw)
    xa = np.linspace(L / 2, L, 200)
    Mx = np.vstack([xa ** 2, np.ones_like(xa), 1 / xa]).T
    ah = np.linalg.lstsq(Mx, sh.sol(xa)[0], rcond=None)[0]
    ap = np.linalg.lstsq(Mx, spp.sol(xa)[0], rcond=None)[0]
    cc = -ap[0] / ah[0]
    binf = cc * ah[1] + ap[1]
    xx = np.linspace(0, L, 30001)
    bet = cc * sh.sol(xx)[0] + spp.sol(xx)[0]
    out = dict(prm, B=bool(binf > 0 and bet.min() > 0))
    if not out["B"]:
        return out
    V0n = 1 / np.sqrt(binf)
    if xs is None:
        xs = np.concatenate([np.geomspace(0.005, 1, 120), np.linspace(1.01, 30, 300)])
    y = (cc * sh.sol(xs) + spp.sol(xs)) / binf
    vals = [fun[u](xs, y[0], y[1], *pv, V0n) for u in usados]
    vals = [np.broadcast_to(np.asarray(v_, float), xs.shape) for v_ in vals]
    p2, p1, p0 = (np.broadcast_to(np.asarray(f(*vals), float), xs.shape) for f in fP)
    esc = np.abs(p2) + np.abs(p1) + np.abs(p0) + 1e-300
    disc = (p1 ** 2 - 4 * p2 * p0) / esc ** 2
    malos = xs[disc < -1e-8]
    vmax = float(np.max(np.abs(fun["V_0"](xs, y[0], y[1], *pv, V0n))))
    out.update(V0=float(V0n), beta0=float(y[0][0]), vmax=vmax, H=bool(len(malos) == 0),
               rho_no_hip=[float(malos.min()), float(malos.max())] if len(malos) else None,
               frac_no_hip=float(np.mean(disc < -1e-8)),
               energia_tau_pos=float(np.mean((p2 < 0) & (p0 > 0))))
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "fijo":
        v = [float(s) for s in sys.argv[2:9]]
        prm = dict(zip(("w", "a", "wq", "An", "wn", "w2", "lam"), v))
        print(evalua(prm), flush=True)
        sys.exit()
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    rng = np.random.default_rng(11)
    RES = []
    t0 = time.time()
    for i in range(n):
        prm = dict(w=float(np.exp(rng.uniform(np.log(0.8), np.log(5)))), a=float(np.exp(rng.uniform(np.log(5), np.log(300)))),
                   wq=float(np.exp(rng.uniform(np.log(0.3), np.log(3)))), An=float(rng.uniform(0, 3)),
                   wn=float(np.exp(rng.uniform(np.log(0.5), np.log(5)))),
                   w2=float(np.exp(rng.uniform(np.log(0.08), np.log(0.6)))), lam=float(rng.uniform(0, 0.3)))
        try:
            c = evalua(prm)
        except Exception as ex_:                       # noqa: BLE001
            c = dict(prm, B=False, error=str(ex_))
        RES.append(c)
        if c.get("B"):
            print(f"[{time.time() - t0:6.0f} s] " + str({k_: (float(f"{x:.4g}") if isinstance(x, float) else x)
                                                       for k_, x in c.items()}), flush=True)
    nB = sum(1 for c in RES if c.get("B"))
    nH = sum(1 for c in RES if c.get("H"))
    print(f"{n} fondos: beta > 0 en {nB}; hiperbolicos en todo rho: {nH}", flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", "criba_cubica.json"), "w"), indent=1)
