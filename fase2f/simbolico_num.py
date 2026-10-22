# -*- coding: utf-8 -*-
"""
Expresiones simbolicas de la Fase 2f para fondos estaticos generales (lapso N = e^Phi, E = N^2, eP = N), con
X_b = -phi'^2/2 sustituido, compiladas para valores o chorros (jets.py).  Las funciones devueltas por compilar()
reciben un dict V con los nombres de los simbolos ('r0', 'r1', 'P1', 'Phi1', 'E', 'G4_01', 'gam0', ...).

  FONDO[n]            E_A, E_B, E_C, E_p del fondo (generar_bH.py)
  G2_DE_LIGADURAS     (G2_00 de E_A = 0, G2_01 de E_B = 0)
  MONO[k]             gam, al1, al0, k_c, K, P, U del monopolo (monopolo_reduccion.py; K, P, U usan los chorros
                      abstractos gam0, gam1, ..., al1_0, ..., al0_0, ...); dP_dG202 = dP/dG2_02 (P es lineal en G2_02)
  IMPAR[l][k]         alfa, beta, gamma del sector impar (L = alfa Q^2 + beta h0^2 + gamma h1^2,
                      Q = h1_t - h0' + 2 r' h0 / r) para l = 2, 3, sin usar las ecuaciones del fondo
"""
import os
import sys
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import jets as J                                                        # noqa: E402
from chorros import rr, PP, Xs, D                                      # noqa: E402

def con_Xb(e):
    """Sustituye X_b = -phi'^2/2 (B = 1)."""
    return sp.expand(sp.expand(e).xreplace({Xs: -PP[1] ** 2 / 2}))


def compilar(expr):
    syms = sorted(expr.free_symbols, key=lambda s: s.name)
    f = sp.lambdify(syms, expr, modules=J.MODULOS + ["numpy"], cse=True)
    nombres = [s.name for s in syms]

    def g(V):
        return f(*[V[n] for n in nombres])
    g.nombres = nombres
    g.expr = expr
    return g


_d = pickle.load(open(os.path.join(AQUI, "resultados", "bH_simbolico.pkl"), "rb"))
FONDO_EXPR = {n: con_Xb(_d["fondo_eq"][n]) for n in ("A", "B", "C", "p")}
FONDO = {n: compilar(e) for n, e in FONDO_EXPR.items()}
_g200, _g201 = sp.Symbol("G2_00", real=True), sp.Symbol("G2_01", real=True)
_sA = sp.solve(FONDO_EXPR["A"], _g200)
_sB = sp.solve(FONDO_EXPR["B"], _g201)
assert len(_sA) == 1 and len(_sB) == 1
G2_DE_LIGADURAS = (compilar(_sA[0]), compilar(_sB[0]))

_m = pickle.load(open(os.path.join(AQUI, "resultados", "monopolo_reduccion.pkl"), "rb"))
MONO_EXPR = {k: con_Xb(_m[k]) for k in ("gam", "al1", "al0", "k_c")}
for k in ("K", "P", "U"):
    MONO_EXPR[k] = sp.together(con_Xb(_m[k]))
MONO_EXPR["dP_dG202"] = sp.together(sp.diff(MONO_EXPR["P"], sp.Symbol("G2_02", real=True)))
MONO = {k: compilar(e) for k, e in MONO_EXPR.items()}


def _impar(LL):
    di = pickle.load(open(os.path.join(AQUI, "resultados", f"impar_ec_L{LL}.pkl"), "rb"))
    ecs, amp, w = di["ecs"], di["amp"], di["w"]
    H0 = [amp[("h0", k)] for k in range(7)]
    H1 = [amp[("h1", k)] for k in range(7)]
    e1, e0 = ecs["h1"], ecs["h0"]
    alfa = sp.expand(e1.coeff(H0[1]) / 2)
    gam = sp.expand(e1.coeff(H1[0]).subs(w, 0) / 2)
    Q = -w * H1[0] - H0[1] + 2 * rr[1] * H0[0] / rr[0]
    bet = sp.Symbol("beta_aux")
    dif = sp.expand(e0 - (2 * alfa * 2 * rr[1] / rr[0] * Q + 2 * bet * H0[0] + 2 * D(alfa * Q, [H0, H1])))
    beta = sp.solve(dif.coeff(H0[0]), bet)[0]
    return {"alfa": con_Xb(alfa), "beta": con_Xb(beta), "gamma": con_Xb(gam)}


IMPAR_EXPR = {LL: _impar(LL) for LL in (2, 3)}
IMPAR = {LL: {k: compilar(e) for k, e in d_.items()} for LL, d_ in IMPAR_EXPR.items()}
