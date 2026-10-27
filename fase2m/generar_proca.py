# -*- coding: utf-8 -*-
"""
Fase 2m: ecuaciones del fondo y accion cuadratica (sector l = 0) de la Proca generalizada cuartica
(proca_simbolico.py) sobre un fondo estatico general con A0 y A1.

Metodo (el de fase2f/generar_bH.py, con dos invariantes):
  1. Lred(chorros, X, F) con u = (A, b, B, C, u = A0, w = A1) y sus derivadas hasta segundo orden en (t, l);
     G2(X, F), G3(X), G4(X) son funciones generales; X y F se tratan como simbolos y las derivadas usan la regla de
     la cadena  dL/du_J = L_J + L_X X_J + L_F F_J  (y la segunda variacion analoga).
  2. Fondo estatico en longitud propia:  A = N(l)^2, b = 0, B = 1, C = r(l)^2, A0 = U(l), A1 = V(l);
     X_b = U^2/(2 N^2) - V^2/2,  F_b = U'^2/(2 N^2).  Ecuaciones E_u = sum_J (-D_l)^{J_l} [dL/du_J]_fondo.
  3. Accion cuadratica  L2 = (1/2) sum_JK [d^2 L/du_J du_K]_fondo du_J du_K  con
     dA = 2 N^2 a, db = beta, dB = 2 c, dC = 2 r d, dA0 = al (u), dA1 = w;  ecuaciones lineales para modos e^{i w t}.
  D_l actua sobre los simbolos de fondo (N_k, r_k, U_k, V_k -> k+1) y G_ij -> G_(i+1)j X_b' + G_i(j+1) F_b'
  (i derivadas en X, j en F).
Comprobaciones aqui: independencia de theta, Y = 4 X F en l = 0, identidad de Noether radial y la temporal
E_b = U E_w / N^2 (fuera de la capa de masas).  Las demas (RG, Horndeski, modos gauge) en verificar_proca.py.
Uso: python generar_proca.py   (escribe resultados/proca_simbolico.pkl)
"""
import os
import sys
import json
import time
import random
import pickle
from math import comb
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import proca_simbolico as S        # noqa: E402

t, l, th = S.t, S.l, S.th
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f} s] {msg}", flush=True)


CAMPOS = ("A", "b", "B", "C", "u", "w")
ORDENES = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
funcs = {n: sp.Function(n)(t, l) for n in CAMPOS}
jet = {(n, o): sp.Symbol(f"{n}_{o[0]}{o[1]}", real=True) for n in CAMPOS for o in ORDENES}


def a_chorros(expr):
    rep = {}
    for n in CAMPOS:
        f = funcs[n]
        for (nt, nl) in sorted(ORDENES, key=lambda o: -(o[0] + o[1])):
            if nt + nl == 0:
                continue
            rep[sp.Derivative(f, *([t] * nt + [l] * nl))] = jet[(n, (nt, nl))]
    return expr.subs(rep).subs({funcs[n]: jet[(n, (0, 0))] for n in CAMPOS})


# ------------------------------------------------------------------ 1. lagrangiana reducida en chorros
Xs, Fs = sp.symbols("Xs Fs", real=True)
G2f, G3f, G4f = (sp.Function(n) for n in ("G2", "G3", "G4"))
FUNCS_G = (G2f, G3f, G4f)


def teoria(X, F):
    return G2f(Xs, Fs), G3f(Xs), G4f(Xs), sp.Derivative(G4f(Xs), Xs)


g = S.metrica_esferica(funcs["A"], funcs["b"], funcs["B"], funcs["C"])
gi, Gam, Ric, Rs = S.geometria(g)
log("geometria")
Avec = [funcs["u"], funcs["w"], 0, 0]
Lden, esc = S.lagrangiana(g, gi, Gam, Rs, Avec, teoria)
CHK = {}
dY = sp.simplify(esc["Y"] - 4 * esc["X"] * esc["F"])
assert dY == 0, dY
CHK["Y_igual_4XF"] = "identica"
log("Y = 4 X F en l = 0: identidad simbolica")
L = a_chorros(Lden)
Xj = a_chorros(esc["X"])
Fj = a_chorros(esc["F"])
random.seed(3)
prueba = {s: sp.Rational(random.randint(300, 900), 1000) for s in jet.values()}
Lnum = L.replace(lambda e: isinstance(e, sp.Derivative), lambda e: sp.Integer(1)) \
    .replace(lambda e: isinstance(e, sp.Function) and e.func in FUNCS_G, lambda e: sp.Integer(1))
d_theta = sp.N((Lnum.subs(th, sp.Rational(3, 10)) - Lnum.subs(th, sp.Rational(11, 10))).subs(prueba), 30)
assert abs(d_theta) < 1e-25, d_theta
CHK["theta"] = float(abs(d_theta))
L = L.subs(th, sp.pi / 3)
assert not any(L.has(funcs[n]) for n in CAMPOS)
log(f"lagrangiana en chorros: {sp.count_ops(L)} operaciones (theta: {float(abs(d_theta)):.1e})")

# ------------------------------------------------------------------ 2. fondo
NMAXB = 9
NN = sp.symbols(f"N0:{NMAXB}", positive=True)
rr = sp.symbols(f"r0:{NMAXB}", real=True)
UU = sp.symbols(f"U0:{NMAXB}", real=True)
VV = sp.symbols(f"V0:{NMAXB}", real=True)
fondo = {jet[("A", (0, 0))]: NN[0] ** 2, jet[("A", (0, 1))]: 2 * NN[0] * NN[1],
         jet[("A", (0, 2))]: 2 * NN[1] ** 2 + 2 * NN[0] * NN[2],
         jet[("C", (0, 0))]: rr[0] ** 2, jet[("C", (0, 1))]: 2 * rr[0] * rr[1],
         jet[("C", (0, 2))]: 2 * rr[1] ** 2 + 2 * rr[0] * rr[2], jet[("B", (0, 0))]: 1}
for k in range(3):
    fondo[jet[("u", (0, k))]] = UU[k]
    fondo[jet[("w", (0, k))]] = VV[k]
for n in CAMPOS:
    for o in ORDENES:
        fondo.setdefault(jet[(n, o)], 0)
Xb = UU[0] ** 2 / (2 * NN[0] ** 2) - VV[0] ** 2 / 2
Fb = UU[1] ** 2 / (2 * NN[0] ** 2)
GS = {}


def gsym(nombre, i, j):
    k = (nombre, i, j)
    if k not in GS:
        GS[k] = sp.Symbol(f"{nombre}_{i}{j}", real=True)
    return GS[k]


def a_simbolos_G(expr):
    rep = {}
    for d in expr.atoms(sp.Derivative):
        if isinstance(d.expr, sp.Function) and d.expr.func in FUNCS_G:
            i = sum(c for v, c in d.variable_count if v == Xs)
            j = sum(c for v, c in d.variable_count if v == Fs)
            rep[d] = gsym(d.expr.func.__name__, i, j)
    expr = expr.xreplace(rep)
    rep = {f: gsym(f.func.__name__, 0, 0) for f in expr.atoms(sp.Function) if f.func in FUNCS_G}
    out = expr.xreplace(rep)
    assert not out.atoms(sp.Derivative), out.atoms(sp.Derivative)
    assert not out.has(Xs) and not out.has(Fs)
    return out


def en_fondo(expr):
    return a_simbolos_G(expr.xreplace(fondo))


def D_bg(expr):
    """Derivada total en l de una expresion de fondo."""
    out = 0
    expr = sp.sympify(expr)
    fs = expr.free_symbols
    for seq in (NN, rr, UU, VV):
        for k in range(NMAXB - 1):
            if seq[k] in fs:
                out += sp.diff(expr, seq[k]) * seq[k + 1]
    dX = None
    for (nom, i, j), s in list(GS.items()):
        if s in fs:
            if dX is None:
                dX, dF = D_bg_simple(Xb), D_bg_simple(Fb)
            out += sp.diff(expr, s) * (gsym(nom, i + 1, j) * dX + (gsym(nom, i, j + 1) * dF if nom == "G2" else 0))
    return out


def D_bg_simple(expr):
    out = 0
    expr = sp.sympify(expr)
    for seq in (NN, rr, UU, VV):
        for k in range(NMAXB - 1):
            out += sp.diff(expr, seq[k]) * seq[k + 1]
    return out


Ld_X, Ld_F = sp.diff(L, Xs), sp.diff(L, Fs)
JETS = [jet[(n, o)] for n in CAMPOS for o in ORDENES]
dX = {j: sp.diff(Xj, j) for j in JETS}
dF = {j: sp.diff(Fj, j) for j in JETS}
Ld = {j: sp.diff(L, j) for j in JETS}
log("derivadas primeras")
prim = {}
for n in CAMPOS:
    for o in ((0, 0), (0, 1), (0, 2)):
        j = jet[(n, o)]
        prim[(n, o)] = sp.expand(en_fondo(Ld[j]) + en_fondo(Ld_X) * en_fondo(dX[j]) + en_fondo(Ld_F) * en_fondo(dF[j]))
fondo_eq = {}
for n in CAMPOS:
    fondo_eq[n] = sp.expand(prim[(n, (0, 0))] - D_bg(prim[(n, (0, 1))]) + D_bg(D_bg(prim[(n, (0, 2))])))
log("ecuaciones de fondo: " + ", ".join(f"{n}: {len(sp.Add.make_args(v))}" for n, v in fondo_eq.items()))


def rel(e, semilla=11):
    """Residuo relativo de una expresion de fondo en un punto racional aleatorio (0 si se anula identicamente)."""
    e = sp.expand(e)
    if e == 0:
        return 0.0
    rnd = random.Random(semilla)
    val = {s: sp.Rational(rnd.randint(-900, 900), 997) for s in e.free_symbols}
    for s in NN:
        val[s] = sp.Rational(rnd.randint(300, 1500), 997) if s == NN[0] else val.get(s, 0)
    num = sp.N(e.xreplace(val), 30)
    esc_ = sp.N(sum(abs(sp.N(tm.xreplace(val), 30)) for tm in sp.Add.make_args(e)), 30)
    return float(abs(num) / esc_)


# Noether radial: delta A = xi A', delta C = xi C', delta B = 2 xi' (B = 1), delta u = xi U', delta w = xi V' + V xi'
noe = (fondo_eq["A"] * 2 * NN[0] * NN[1] + fondo_eq["C"] * 2 * rr[0] * rr[1] + fondo_eq["u"] * UU[1]
       + fondo_eq["w"] * VV[1] - 2 * D_bg(fondo_eq["B"]) - D_bg(VV[0] * fondo_eq["w"]))
CHK["noether_radial"] = rel(noe)
log(f"identidad de Noether radial: {CHK['noether_radial']:.1e}")
# Noether temporal (t -> t + T(l)): delta b = -N^2 T', delta w = U T'  =>  N^2 E_b - U E_w = 0
CHK["noether_temporal"] = rel(NN[0] ** 2 * fondo_eq["b"] - UU[0] * fondo_eq["w"])
log(f"identidad temporal N^2 E_b = U E_w: {CHK['noether_temporal']:.1e}")
assert CHK["noether_radial"] < 1e-25 and CHK["noether_temporal"] < 1e-25

# ------------------------------------------------------------------ 3. accion cuadratica
PERT = ("a", "be", "c", "d", "al", "w")
pj = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in PERT for o in ORDENES}
PERT_DE = {"A": ("a", 2 * NN[0] ** 2), "b": ("be", 1), "B": ("c", 2), "C": ("d", 2 * rr[0]), "u": ("al", 1),
           "w": ("w", 1)}


def dpert(n, o):
    """Chorro (nt, nl) de f(l) p(t, l) por Leibniz en l."""
    v, f = PERT_DE[n]
    nt, nl = o
    out = 0
    fk = f
    for k in range(nl + 1):
        out += comb(nl, k) * fk * pj[(v, (nt, nl - k))]
        fk = D_bg_simple(fk)
    return out


# piezas de la regla de la cadena evaluadas una sola vez en el fondo
Ldx_b = {j: en_fondo(sp.diff(Ld[j], Xs)) for j in JETS}
Ldf_b = {j: en_fondo(sp.diff(Ld[j], Fs)) for j in JETS}
dX_b = {j: en_fondo(dX[j]) for j in JETS}
dF_b = {j: en_fondo(dF[j]) for j in JETS}
LX_b, LF_b = en_fondo(Ld_X), en_fondo(Ld_F)
LXX_b, LXF_b, LFF_b = en_fondo(sp.diff(Ld_X, Xs)), en_fondo(sp.diff(Ld_X, Fs)), en_fondo(sp.diff(Ld_F, Fs))
log("piezas de la regla de la cadena en el fondo")
INV = {jet[k]: k for k in jet}
L2 = 0
n_nz = 0
for i1, j1 in enumerate(JETS):
    for i2, j2 in enumerate(JETS):
        if i2 < i1:
            continue
        h = (en_fondo(sp.diff(Ld[j1], j2)) + Ldx_b[j1] * dX_b[j2] + Ldx_b[j2] * dX_b[j1] + Ldf_b[j1] * dF_b[j2]
             + Ldf_b[j2] * dF_b[j1] + LXX_b * dX_b[j1] * dX_b[j2] + LXF_b * (dX_b[j1] * dF_b[j2] + dX_b[j2] * dF_b[j1])
             + LFF_b * dF_b[j1] * dF_b[j2] + LX_b * en_fondo(sp.diff(dX[j1], j2)) + LF_b * en_fondo(sp.diff(dF[j1], j2)))
        hb = sp.expand(h)
        if hb == 0:
            continue
        n_nz += 1
        n1, o1 = INV[j1]
        n2, o2 = INV[j2]
        L2 += (sp.Rational(1, 2) if i1 == i2 else 1) * hb * dpert(n1, o1) * dpert(n2, o2)
    log(f"  fila {i1 + 1}/{len(JETS)}")
L2 = sp.expand(L2)
log(f"accion cuadratica: {n_nz} entradas no nulas del hessiano, {len(sp.Add.make_args(L2))} terminos")

iw = sp.Symbol("iw")
NMAX = 5
amp = {(v, k): sp.Symbol(f"{v}{k}", real=True) for v in PERT for k in range(NMAX + 1)}
modo = {pj[(v, o)]: iw ** o[0] * amp[(v, o[1])] for v in PERT for o in ORDENES}


def D_lin(expr):
    expr = sp.sympify(expr)
    out = D_bg(expr)
    fs = expr.free_symbols
    for v in PERT:
        for k in range(NMAX):
            if amp[(v, k)] in fs:
                out += sp.diff(expr, amp[(v, k)]) * amp[(v, k + 1)]
    return out


lineales = {}
for v in PERT:
    tot = 0
    for o in ORDENES:
        Q = sp.diff(L2, pj[(v, o)]).subs(modo) * (-1) ** (o[0] + o[1]) * iw ** o[0]
        for _ in range(o[1]):
            Q = D_lin(Q)
        tot += Q
    lineales[v] = sp.expand(tot)
    log(f"ecuacion lineal E_{v}: {len(sp.Add.make_args(lineales[v]))} terminos")

os.makedirs(os.path.join(AQUI, "resultados"), exist_ok=True)
with open(os.path.join(AQUI, "resultados", "proca_simbolico.pkl"), "wb") as fh:
    pickle.dump(dict(fondo_eq=fondo_eq, lineales=lineales, L2=L2, GS=GS, NN=NN, rr=rr, UU=UU, VV=VV, iw=iw,
                     amp=amp, pj=pj, CHK=CHK), fh)
CHK["tamanos"] = {"fondo": {n: len(sp.Add.make_args(v)) for n, v in fondo_eq.items()},
                  "lineales": {v: len(sp.Add.make_args(e)) for v, e in lineales.items()}, "hessiano_no_nulas": n_nz}
json.dump(CHK, open(os.path.join(AQUI, "resultados", "generar_proca.json"), "w"), indent=1)
log("guardado resultados/proca_simbolico.pkl")
