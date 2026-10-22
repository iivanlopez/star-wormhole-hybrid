# -*- coding: utf-8 -*-
"""
Ecuaciones del sector impar a partir de resultados/impar_L{L}.pkl (impar_simbolico.py): Euler-Lagrange de la
lagrangiana cuadratica para modos h0, h1 ~ e^{i omega t}, con el operador D_l sobre simbolos de fondo
(Phi_k, r_k, P_k, eP = e^{Phi0}) y G_ij, como en generar_bH.py.
Devuelve (y guarda en resultados/impar_ec_L{L}.pkl) las ecuaciones E_h0, E_h1 como combinaciones lineales de
H0_k, H1_k (k-esima derivada en l de las amplitudes) con coeficientes en w = omega^2 (E_h1 se divide por i omega si
es impar).  Uso: python impar_ecuaciones.py 2
"""
import os
import sys
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from bH_simbolico import t, l      # noqa: E402

LL = int(sys.argv[1]) if len(sys.argv) > 1 else 2
d = pickle.load(open(os.path.join(AQUI, "resultados", f"impar_L{LL}.pkl"), "rb"))
L2 = d["L2"]
Xs = sp.Symbol("Xs", real=True)
Phif, rf, Pf = (sp.Function(n)(l) for n in ("Phi", "r", "P"))
h0, h1 = sp.Function("h0")(t, l), sp.Function("h1")(t, l)
Phi = sp.symbols("Phi0:7", real=True)
rr = sp.symbols("r0:7", real=True)
PP = sp.symbols("P0:8", real=True)
eP = sp.Symbol("eP", positive=True)
G2f, G3f, G4f, F4f = (sp.Function(nm) for nm in ("G2", "G3", "G4", "F4"))
GS = {}


def gsym(nombre, i, j):
    k = (nombre, i, j)
    if k not in GS:
        GS[k] = sp.Symbol(f"{nombre}_{i}{j}", real=True)
    return GS[k]


ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
pj = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in ("h0", "h1") for o in ORD}


def a_simbolos(e):
    rep = {}
    for (fn, nm) in ((h0, "h0"), (h1, "h1")):
        for o in sorted(ORD, key=lambda o: -(o[0] + o[1])):
            if o == (0, 0):
                continue
            rep[sp.Derivative(fn, *([t] * o[0] + [l] * o[1]))] = pj[(nm, o)]
    e = e.subs(rep).subs({h0: pj[("h0", (0, 0))], h1: pj[("h1", (0, 0))]})
    # derivadas de G (respecto de Xs) antes de sustituir P(l)
    rG = {}
    for dv in e.atoms(sp.Derivative):
        if isinstance(dv.expr, sp.Function) and dv.expr.func in (G2f, G3f, G4f, F4f):
            j = sum(c for v, c in dv.variable_count if v == Xs)
            assert j == sum(c for _, c in dv.variable_count)
            rG[dv] = gsym(dv.expr.func.__name__, 0, j)
    e = e.xreplace(rG)
    e = e.xreplace({f: gsym(f.func.__name__, 0, 0) for f in e.atoms(sp.Function) if f.func in (G2f, G3f, G4f, F4f)})
    rb = {}
    for k in range(6, 0, -1):
        rb[sp.Derivative(Phif, (l, k))] = Phi[k]
        rb[sp.Derivative(rf, (l, k))] = rr[k]
        rb[sp.Derivative(Pf, (l, k))] = PP[k]
    e = e.subs(rb).subs({sp.exp(Phif): eP, sp.exp(2 * Phif): eP**2, sp.exp(-Phif): 1 / eP, sp.exp(-2 * Phif): eP**-2})
    e = e.subs({rf: rr[0], Pf: PP[0], Phif: Phi[0]})
    assert not e.atoms(sp.Derivative), e.atoms(sp.Derivative)
    assert not e.has(Phi[0]), "Phi0 sin derivar"
    return sp.expand(e)


L2s = a_simbolos(L2)
iw = sp.Symbol("iw")
NMAX = 6
amp = {(v, k): sp.Symbol(f"{v.upper()}{k}", real=True) for v in ("h0", "h1") for k in range(NMAX + 1)}
modo = {pj[(v, o)]: iw ** o[0] * amp[(v, o[1])] for v in ("h0", "h1") for o in ORD}


def D_l(e):
    out = 0
    for k in range(6):
        out += sp.diff(e, Phi[k]) * Phi[k + 1] + sp.diff(e, rr[k]) * rr[k + 1]
    for k in range(7):
        out += sp.diff(e, PP[k]) * PP[k + 1]
    out += sp.diff(e, eP) * Phi[1] * eP
    out += sp.diff(e, Xs) * (-PP[1] * PP[2])
    for (nom, i, j), s in list(GS.items()):
        if e.has(s):
            out += sp.diff(e, s) * (gsym(nom, i + 1, j) * PP[1] - gsym(nom, i, j + 1) * PP[1] * PP[2])
    for v in ("h0", "h1"):
        for k in range(NMAX):
            if e.has(amp[(v, k)]):
                out += sp.diff(e, amp[(v, k)]) * amp[(v, k + 1)]
    return out


ecs = {}
w = sp.Symbol("w", real=True)
for v in ("h0", "h1"):
    tot = 0
    for o in ORD:
        Q = sp.diff(L2s, pj[(v, o)]).subs(modo) * (-1) ** (o[0] + o[1]) * iw ** o[0]
        for _ in range(o[1]):
            Q = D_l(Q)
        tot += Q
    e = sp.expand(tot)
    pot = sorted(set(m[0] for m in sp.Poly(e, iw).monoms()))
    print(f"E_{v}: potencias de i omega {pot}")
    ecs[v] = e
# h0 y h1 tienen paridad temporal opuesta: con h0 -> H0, h1 -> iw H1 todo queda en potencias pares
par = {amp[("h1", k)]: iw * amp[("h1", k)] for k in range(NMAX + 1)}
for v in ("h0", "h1"):
    e = sp.expand(ecs[v].subs(par))
    pot = sorted(set(m[0] for m in sp.Poly(e, iw).monoms()))
    if all(p % 2 == 1 for p in pot):
        e = sp.expand(e / iw)
    e = sp.expand(e.subs(iw**2, -w))
    assert not e.has(iw), (v, pot)
    ecs[v] = e
    print(f"E_{v} con h1 -> i omega H1: amplitudes {sorted(str(s) for s in e.free_symbols if s in set(amp.values()))}")
with open(os.path.join(AQUI, "resultados", f"impar_ec_L{LL}.pkl"), "wb") as fh:
    pickle.dump(dict(ecs=ecs, GS=GS, Phi=Phi, rr=rr, PP=PP, eP=eP, Xs=Xs, amp=amp, w=w, LL=LL), fh)
print(f"guardado resultados/impar_ec_L{LL}.pkl; simbolos G: {sorted(str(s) for s in GS.values())}")
