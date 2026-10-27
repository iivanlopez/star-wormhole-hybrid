# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): el monopolo de una garganta con galileon cubico y flujo (G3 != 0, psi' != 0).

1. Simbolico: lagrangiana ADM dependiente del tiempo (verificada en adm_tiempo_cov.py), gauge h = 1, perturbaciones
   Nc = N + nu, v = V + w, r = R + d; ecuaciones lineales de Euler-Lagrange como operadores:
       E_q = sum_{f, m, n} C[q, f, m, n](rho) d_t^m d_rho^n f,     q, f in (nu, w, d).
2. Fondo (reconstruccion, reconstruccion_cubica.py): r = sqrt(rho^2 + 1), v = V0 x e^{-x^2/2} (x = rho/w),
   Q = 2 A4/Nc = -2 [1 + a e^{-rho^2/wq^2}], Nc = 1 + An e^{-rho^2/wn^2}, beta = Nc G4 de su EDO (par, acotada),
   p = Nc A2 de la ligadura H = 0, A3' = 2 v (Q r')'/r.  Las funciones de la teoria como funciones de Nc:
   f_N = f'/Nc', f_NN = (f_N)'/Nc' (Nc monotona en cada mitad).  Todas las derivadas analiticas salvo beta, cuyas
   derivadas altas salen de su EDO.
Uso: importar (fondo(...), coeficientes(...)).
"""
import os
import pickle
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
t, rho, e = sp.symbols("t rho e", real=True)
NOMS = ("a2", "a3", "a4", "b4")
CACHE = os.path.join(AQUI, "resultados", "monopolo_cubico_op.pkl")


def operador():
    """Coeficientes simbolicos del operador lineal; las funciones de fondo como Function(rho)."""
    if os.path.exists(CACHE):
        return pickle.load(open(CACHE, "rb"))
    G = lambda n: sp.Function(n)(rho)                   # noqa: E731
    H = lambda n: sp.Function(n)(t, rho)                # noqa: E731
    N, V, R = G("N"), G("V"), G("R")
    fun = {n: [G(n), G(n + "N"), G(n + "NN")] for n in NOMS}
    nu, w, d = H("nu"), H("w"), H("d")
    Dr = lambda x, k_=1: sp.diff(x, rho, k_)            # noqa: E731
    Dt = lambda x: sp.diff(x, t)                        # noqa: E731
    Nc, v, r = N + e * nu, V + e * w, R + e * d
    A = {n: fun[n][0] + fun[n][1] * e * nu + fun[n][2] * (e * nu) ** 2 / 2 for n in NOMS}
    Kr = -Dr(v) / Nc
    Kth = (Dt(r) - v * Dr(r)) / (r * Nc)
    R3 = (2 - 4 * r * Dr(r, 2) - 2 * Dr(r) ** 2) / r ** 2
    L = Nc * r ** 2 * (A["a2"] + A["a3"] * (Kr + 2 * Kth) + A["a4"] * (4 * Kr * Kth + 2 * Kth ** 2) + A["b4"] * R3)
    L2 = sp.expand(sp.diff(L, e, 2).subs(e, 0) / 2)
    campos = (nu, w, d)
    C = {}
    for q in campos:
        EL = sp.expand(sp.euler_equations(L2, [q], [t, rho])[0].lhs)
        ders = [a for a in EL.atoms(sp.Derivative) if a.expr in campos] + list(campos)
        for der in ders:
            if isinstance(der, sp.Derivative):
                f = der.expr
                m = sum(c_ for v_, c_ in der.variable_count if v_ == t)
                n = sum(c_ for v_, c_ in der.variable_count if v_ == rho)
            else:
                f, m, n = der, 0, 0
            co = EL.coeff(der)
            if co != 0:
                key = (str(q.func), str(f.func), m, n)
                C[key] = sp.expand(C.get(key, 0) + co)
            EL = sp.expand(EL - co * der)
        assert EL == 0, EL
    # simbolos de fondo
    bg = [N, V, R] + [x for n in NOMS for x in fun[n]]
    sims = {}
    for f in bg:
        for n in range(5):
            key = f if n == 0 else sp.Derivative(f, (rho, n))
            sims[key] = sp.Symbol(f"{f.func}_{n}")
    Cs = {k: v_.xreplace({kk: vv for kk, vv in sims.items() if isinstance(kk, sp.Derivative)}).xreplace(
        {kk: vv for kk, vv in sims.items() if not isinstance(kk, sp.Derivative)}) for k, v_ in C.items()}
    usados = sorted({s for v_ in Cs.values() for s in v_.free_symbols}, key=str)
    out = dict(C=Cs, usados=usados)
    pickle.dump(out, open(CACHE, "wb"))
    return out


# ------------------------------------------------------------------ fondo
X = sp.Symbol("x", real=True)                       # rho
B0, B1 = sp.symbols("B0 B1", real=True)             # beta, beta'
P = dict(w=2.0, a=20.0, wq=1.0, An=1.0, wn=2.0, L=400.0)


def fondo_simbolico(w, a, wq, An, wn, w2=1.0, lam=0.0, lapso="gauss"):
    r = sp.sqrt(X ** 2 + 1)
    z, z2 = X / w, X / w2
    v_ = z * sp.exp(-z ** 2 / 2) + lam * z2 * sp.exp(-z2 ** 2 / 2)       # V0 = 1; se escala despues
    Q = -2 * (1 + a * sp.exp(-X ** 2 / wq ** 2))
    Nc = 1 + An * (sp.exp(-X ** 2 / wn ** 2) if lapso == "gauss" else 1 / (1 + X ** 2 / wn ** 2))
    return r, v_, Q, Nc


def construir(w=2.0, a=20.0, wq=1.0, An=1.0, wn=2.0, L=400.0, w2=1.0, lam=0.0, lapso="gauss"):
    """Resuelve beta (V0 = 1), fija V0 con beta_inf = 1 y devuelve funciones numericas de rho para todos los
    simbolos de fondo usados por el operador."""
    r, v1, Q, Nc = fondo_simbolico(w, a, wq, An, wn, w2, lam, lapso)
    V0 = sp.Symbol("V0", positive=True)
    v_ = V0 * v1
    D = lambda f: sp.diff(f, X)                     # noqa: E731
    K = 1 + r * D(D(r)) - D(r) ** 2
    u = v_ ** 2 / r ** 2
    S = D(Q * r ** 4 * D(u)) / (4 * r ** 3)
    # EDO: (beta'/r)' + K beta/r^3 = -S  ->  beta'' = r'/r beta' - K beta/r^2 - r S
    b2 = D(r) / r * B1 - K * B0 / r ** 2 - r * S

    def Dtot(ex):
        """derivada total en rho con beta, beta' (beta'' de la EDO)."""
        return sp.diff(ex, X) + sp.diff(ex, B0) * B1 + sp.diff(ex, B1) * b2

    # integrar beta numericamente con V0 = 1 (beta proporcional a V0^2)
    fb2 = sp.lambdify((X, B0, B1), b2.subs(V0, 1), "numpy")
    fb2h = sp.lambdify((X, B0, B1), (D(r) / r * B1 - K * B0 / r ** 2), "numpy")
    kw = dict(dense_output=True, rtol=1e-12, atol=1e-14, max_step=0.05)
    sh = solve_ivp(lambda x, y: [y[1], fb2h(x, y[0], y[1])], (0, L), [1.0, 0.0], **kw)
    sp_ = solve_ivp(lambda x, y: [y[1], fb2(x, y[0], y[1])], (0, L), [0.0, 0.0], **kw)
    xa = np.linspace(L / 2, L, 300)
    Mx = np.vstack([xa ** 2, np.ones_like(xa), 1 / xa]).T
    ah = np.linalg.lstsq(Mx, sh.sol(xa)[0], rcond=None)[0]
    ap = np.linalg.lstsq(Mx, sp_.sol(xa)[0], rcond=None)[0]
    cc = -ap[0] / ah[0]
    binf = cc * ah[1] + ap[1]
    V0n = 1 / np.sqrt(binf)

    def beta_num(x):
        x = np.abs(np.asarray(x, dtype=float))
        y = cc * sh.sol(x) + sp_.sol(x)
        return y[0] / binf, y[1] / binf               # escala V0^2: beta -> beta/binf

    # signo de beta' para rho < 0 (beta par): beta'(-x) = -beta'(x)
    # funciones de la teoria a lo largo del fondo
    b4 = B0 / Nc
    a4 = Q * Nc / 2
    p = (2 * B0 * (D(r) ** 2 - 1) + 4 * B1 * r * D(r) + Q * D(r) * D(v_ ** 2 * r)) / r ** 2
    a2 = p / Nc
    a3p = 2 * v_ * D(Q * D(r)) / r                    # A3'
    N1 = D(Nc)
    base = {"N": Nc, "V": v_, "R": r}
    teor = {}
    for nom, f in (("a2", a2), ("b4", b4), ("a4", a4)):
        fN = Dtot(f) / N1
        fNN = Dtot(fN) / N1
        teor[nom] = [f, fN, fNN]
    a3N = a3p / N1
    a3NN = Dtot(a3N) / N1
    teor["a3"] = [None, a3N, a3NN]
    exprs = {}
    for nom, f in base.items():
        g = f
        for n in range(5):
            exprs[f"{nom}_{n}"] = g
            g = Dtot(g)
    for nom in NOMS:
        for j, suf in enumerate(("", "N", "NN")):
            g = teor[nom][j]
            if g is None:                              # a3 (sin derivar): solo su derivada importa (A3' = a3p)
                exprs[f"{nom}_0"] = sp.Integer(0)
                g = a3p
                for n in range(1, 5):
                    exprs[f"{nom}_{n}"] = g
                    g = Dtot(g)
                continue
            for n in range(5):
                exprs[f"{nom}{suf}_{n}"] = g
                g = Dtot(g)
    return dict(exprs=exprs, V0=V0, V0n=V0n, beta=beta_num, binf=binf,
                params=dict(w=w, a=a, wq=wq, An=An, wn=wn, w2=w2, lam=lam, lapso=lapso))


def evaluar(fd, nombres, x):
    """valores numericos en los puntos x de los simbolos de fondo pedidos."""
    b, b1 = fd["beta"](x)
    b1 = np.sign(x) * b1
    out = {}
    for n in nombres:
        ex = fd["exprs"][n].subs(fd["V0"], fd["V0n"])
        f = sp.lambdify((X, B0, B1), ex, "numpy")
        val = f(x, b, b1)
        out[n] = np.broadcast_to(np.asarray(val, dtype=float), np.shape(x)).copy()
    return out
