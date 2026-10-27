# -*- coding: utf-8 -*-
"""
Fase 2j: ecuaciones radiales (sector l = 0) de la teoria quintica (q_simbolico.py: G2, G3, G4, F4, G5 y, con la opcion
F5, el termino quintico mas alla de Horndeski) con el mismo metodo que fase2f/generar_bH.py:
  1. lagrangiana reducida en chorros u = (A, b, B, C, phi) hasta segundo orden en (t, l), X como simbolo (regla de la
     cadena); 2. fondo estatico en longitud propia (A = e^{2 Phi}, b = 0, B = 1, C = r^2, phi = P(l)):
     E_u = sum_J (-D_l)^{J_l} [dL/du_J]; 3. accion cuadratica con dA = 2 e^{2 Phi} a, db = beta, dB = 2 c, dC = 2 r d,
     dphi = q y ecuaciones lineales para modos e^{i omega t}.
Comprobaciones:
  (a) identidad de Noether del fondo (reparametrizaciones de l), identica;
  (b) limite cuartico (G5 = F5 = 0): ecuaciones del fondo y accion cuadratica identicas a las de la Fase 2f;
  (c) la teoria de Lu-Pang (G2 = 8 al X^2, G3 = -8 al X, G4 = 1 + 4 al X, G5 = -4 al ln|X|) con su agujero negro
      h = 1 + (r^2/2al)(1 - sqrt(1 + 8 al M/r^3)), phi' = (sqrt(h) - 1)/r (en longitud propia) cumple las ecuaciones
      (la semilla de las gargantas de Bakopoulos-Charmousis-Kanti, arXiv:2111.09857).
Uso: python generar_q.py [F5]      (escribe resultados/q_simbolico.pkl o q5_simbolico.pkl y el .json)
"""
import os
import sys
import json
import time
import random
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import q_simbolico as S        # noqa: E402

CON_F5 = len(sys.argv) > 1 and sys.argv[1] == "F5"
NOMBRE = "q5" if CON_F5 else "q"
t, l, th = S.t, S.l, S.th
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f} s] {msg}", flush=True)


CAMPOS = ("A", "b", "B", "C", "p")
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
    e = expr.subs(rep)
    return e.subs({funcs[n]: jet[(n, (0, 0))] for n in CAMPOS})


# ------------------------------------------------------------------ 1. lagrangiana reducida en chorros
Xs = sp.Symbol("Xs", real=True)
NOMBRES = ("G2", "G3", "G4", "F4", "G5", "F5")
FN = {n: sp.Function(n) for n in NOMBRES}


def teoria(phi, X):
    G = {n: FN[n](phi, Xs) for n in NOMBRES}
    return (G["G2"], G["G3"], G["G4"], sp.Derivative(G["G4"], Xs), G["F4"], G["G5"], sp.Derivative(G["G5"], Xs),
            G["F5"])


g = S.metrica_esferica(funcs["A"], funcs["b"], funcs["B"], funcs["C"])
gi, Gam, Ric, Rs = S.geometria(g)
log("geometria")
Lden, esc = S.lagrangiana(g, gi, Gam, Ric, Rs, funcs["p"], teoria, con_F5=CON_F5)
L = a_chorros(Lden)
Xj = a_chorros(esc["X"])
assert not any(L.has(funcs[n]) for n in CAMPOS)
random.seed(3)
prueba = {s: sp.Rational(random.randint(300, 900), 1000) for s in jet.values()}
prueba[Xs] = sp.Rational(-1, 7)
Lnum = L.replace(lambda e: isinstance(e, sp.Function) and e.func in tuple(FN.values()), lambda e: sp.Integer(1)) \
    .replace(lambda e: isinstance(e, sp.Derivative), lambda e: sp.Integer(1))
d_theta = sp.N((Lnum.subs(th, sp.Rational(3, 10)) - Lnum.subs(th, sp.Rational(11, 10))).subs(prueba), 30)
assert abs(d_theta) < 1e-25, d_theta
L = L.subs(th, sp.pi / 3)
log(f"lagrangiana en chorros: {sp.count_ops(L)} operaciones (independiente de theta: {float(abs(d_theta)):.1e})")

# ------------------------------------------------------------------ 2. derivadas en el fondo estatico
Phi = sp.symbols("Phi0:9", real=True)
rr = sp.symbols("r0:9", real=True)
PP = sp.symbols("P0:10", real=True)
Ee = sp.Symbol("E", positive=True)
fondo = {jet[("A", (0, 0))]: Ee, jet[("A", (0, 1))]: 2 * Phi[1] * Ee, jet[("A", (0, 2))]: (2 * Phi[2] + 4 * Phi[1] ** 2) * Ee,
         jet[("C", (0, 0))]: rr[0] ** 2, jet[("C", (0, 1))]: 2 * rr[0] * rr[1],
         jet[("C", (0, 2))]: 2 * rr[1] ** 2 + 2 * rr[0] * rr[2],
         jet[("B", (0, 0))]: 1, jet[("p", (0, 0))]: PP[0], jet[("p", (0, 1))]: PP[1], jet[("p", (0, 2))]: PP[2]}
for n in CAMPOS:
    for o in ORDENES:
        fondo.setdefault(jet[(n, o)], 0)
GS = {}


def gsym(nombre, i, j):
    k = (nombre, i, j)
    if k not in GS:
        GS[k] = sp.Symbol(f"{nombre}_{i}{j}", real=True)
    return GS[k]


def a_simbolos_G(expr):
    rep = {}
    for d in expr.atoms(sp.Derivative):
        if isinstance(d.expr, sp.Function) and d.expr.func in tuple(FN.values()):
            i = sum(c for v, c in d.variable_count if v == PP[0])
            j = sum(c for v, c in d.variable_count if v == Xs)
            assert i + j == sum(c for _, c in d.variable_count), d
            rep[d] = gsym(d.expr.func.__name__, i, j)
    expr = expr.xreplace(rep)
    rep = {f: gsym(f.func.__name__, 0, 0) for f in expr.atoms(sp.Function) if f.func in tuple(FN.values())}
    out = expr.xreplace(rep)
    assert not out.atoms(sp.Derivative), out.atoms(sp.Derivative)
    return out


def en_fondo(expr):
    return a_simbolos_G(expr.subs(fondo))


Ldx = sp.diff(L, Xs)
Ldxx = sp.diff(Ldx, Xs)
JETS = [jet[(n, o)] for n in CAMPOS for o in ORDENES]
dX = {j: sp.diff(Xj, j) for j in JETS}
Ld = {j: sp.diff(L, j) for j in JETS}
log("derivadas primeras")
prim = {}
for n in CAMPOS:
    for o in ((0, 0), (0, 1), (0, 2)):
        j = jet[(n, o)]
        prim[(n, o)] = sp.expand(en_fondo(Ld[j] + Ldx * dX[j]))
log("primeras derivadas en el fondo")


def D_l(expr):
    out = 0
    for k in range(8):
        out += sp.diff(expr, Phi[k]) * Phi[k + 1] + sp.diff(expr, rr[k]) * rr[k + 1]
    for k in range(9):
        out += sp.diff(expr, PP[k]) * PP[k + 1]
    out += sp.diff(expr, Ee) * 2 * Phi[1] * Ee
    out += sp.diff(expr, Xs) * (-PP[1] * PP[2])
    for (nom, i, j), s in list(GS.items()):
        if expr.has(s):
            out += sp.diff(expr, s) * (gsym(nom, i + 1, j) * PP[1] - gsym(nom, i, j + 1) * PP[1] * PP[2])
    return out


fondo_eq = {}
for n in CAMPOS:
    fondo_eq[n] = sp.expand(prim[(n, (0, 0))] - D_l(prim[(n, (0, 1))]) + D_l(D_l(prim[(n, (0, 2))])))
log("ecuaciones de fondo: " + ", ".join(f"{n}: {sp.count_ops(v)}" for n, v in fondo_eq.items()))
CHK = {}
noether = sp.expand(fondo_eq["A"] * 2 * Phi[1] * Ee + fondo_eq["C"] * 2 * rr[0] * rr[1] + fondo_eq["p"] * PP[1]
                    - 2 * D_l(fondo_eq["B"]))
CHK["noether_fondo"] = "0 (identica)" if noether == 0 else str(sp.count_ops(noether))
assert noether == 0, sp.count_ops(noether)
log("(a) identidad de Noether del fondo: se anula identicamente")

# (b) limite cuartico frente a la Fase 2f
d2f = pickle.load(open(os.path.join(AQUI, "..", "fase2f", "resultados", "bH_simbolico.pkl"), "rb"))
cero5 = {s: 0 for (nom, i, j), s in GS.items() if nom in ("G5", "F5")}
for n in CAMPOS:
    dif = sp.expand(fondo_eq[n].xreplace(cero5) - d2f["fondo_eq"][n])
    assert dif == 0, (n, sp.count_ops(dif))
CHK["limite_cuartico_fondo"] = "identico a la Fase 2f"
log("(b) limite cuartico (G5 = F5 = 0): ecuaciones del fondo identicas a las de la Fase 2f")

# (c) Lu-Pang: agujero negro de la teoria 4D EGB
al_, M_, rv = sp.symbols("alpha M rv", positive=True)
h = 1 + rv ** 2 / (2 * al_) * (1 - sp.sqrt(1 + 8 * al_ * M_ / rv ** 3))
sq = sp.sqrt(h)
dl = lambda e: sp.diff(e, rv) * sq                      # noqa: E731   d/dl = sqrt(h) d/dr
ePhi = [sp.log(h) / 2]
er = [rv]
ep = [None, (sq - 1) / rv]
for k in range(1, 5):
    ePhi.append(dl(ePhi[-1]))
    er.append(dl(er[-1]))
for k in range(2, 6):
    ep.append(dl(ep[-1]))
GLP = {"G2": 8 * al_ * Xs ** 2, "G3": -8 * al_ * Xs, "G4": 1 + 4 * al_ * Xs, "F4": sp.Integer(0),
       "G5": -4 * al_ * sp.log(-Xs), "F5": sp.Integer(0)}
res_LP = {n: 0.0 for n in ("A", "B", "C", "p")}
for a_, m_, r_ in ((sp.Rational(3, 10), 1, sp.Rational(27, 10)), (sp.Rational(1, 10), 2, sp.Rational(73, 10)),
                   (sp.Rational(8, 10), sp.Rational(3, 2), sp.Rational(41, 10))):
    pt = {al_: a_, M_: m_, rv: r_}
    num = {Ee: sp.N(h.subs(pt), 40), rr[0]: sp.N(r_, 40)}
    num.update({Phi[k]: sp.N(ePhi[k].subs(pt), 40) for k in range(1, 5)})
    num.update({rr[k]: sp.N(er[k].subs(pt), 40) for k in range(1, 5)})
    num.update({PP[k]: sp.N(ep[k].subs(pt), 40) for k in range(1, 6)})
    Xn = -num[PP[1]] ** 2 / 2
    num[Xs] = Xn
    for (nom, i, j), s in GS.items():
        num[s] = 0 if i > 0 else sp.N(sp.diff(GLP[nom], Xs, j).subs({Xs: Xn, al_: a_}), 40)
    for n in res_LP:
        res_LP[n] = max(res_LP[n], abs(float(sp.N(fondo_eq[n].xreplace(num), 30))))
CHK["lu_pang"] = res_LP
assert max(res_LP.values()) < 1e-25, res_LP
log(f"(c) agujero negro de Lu-Pang: residuos de E_A, E_B, E_C, E_p {max(res_LP.values()):.1e} (30 cifras)")

# ------------------------------------------------------------------ 3. accion cuadratica y ecuaciones lineales l = 0
PERT = ("a", "be", "c", "d", "q")
pj = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in PERT for o in ORDENES}
E2 = 2 * Ee


def dpert(n, o):
    nt, nl = o
    if n == "A":
        if nl == 0:
            return E2 * pj[("a", o)]
        if nl == 1:
            return E2 * (2 * Phi[1] * pj[("a", (nt, 0))] + pj[("a", (nt, 1))])
        return E2 * ((2 * Phi[2] + 4 * Phi[1] ** 2) * pj[("a", (0, 0))] + 4 * Phi[1] * pj[("a", (0, 1))] + pj[("a", (0, 2))])
    if n == "b":
        return pj[("be", o)]
    if n == "B":
        return 2 * pj[("c", o)]
    if n == "C":
        if nl == 0:
            return 2 * rr[0] * pj[("d", o)]
        if nl == 1:
            return 2 * (rr[1] * pj[("d", (nt, 0))] + rr[0] * pj[("d", (nt, 1))])
        return 2 * (rr[2] * pj[("d", (0, 0))] + 2 * rr[1] * pj[("d", (0, 1))] + rr[0] * pj[("d", (0, 2))])
    return pj[("q", o)]


Ldj_x = {j: sp.diff(Ld[j], Xs) for j in JETS}
n_nz = 0
L2 = 0
for i1, j1 in enumerate(JETS):
    for i2, j2 in enumerate(JETS):
        if i2 < i1:
            continue
        dxx = sp.diff(dX[j1], j2)
        hh = sp.diff(Ld[j1], j2) + Ldj_x[j1] * dX[j2] + Ldj_x[j2] * dX[j1] + Ldxx * dX[j1] * dX[j2] + Ldx * dxx
        hb = sp.expand(en_fondo(hh))
        if hb == 0:
            continue
        n_nz += 1
        n1, o1 = [k for k, v in jet.items() if v == j1][0]
        n2, o2 = [k for k, v in jet.items() if v == j2][0]
        L2 += (sp.Rational(1, 2) if i1 == i2 else 1) * hb * dpert(n1, o1) * dpert(n2, o2)
    log(f"   hessiano: fila {i1 + 1}/{len(JETS)}")
L2 = sp.expand(L2)
log(f"accion cuadratica: {n_nz} entradas no nulas del hessiano, {len(sp.Add.make_args(L2))} terminos")
dif = sp.expand(L2.xreplace(cero5) - d2f["L2"])
assert dif == 0, sp.count_ops(dif)
CHK["limite_cuartico_L2"] = "identica a la Fase 2f"
log("(b') limite cuartico: accion cuadratica identica a la de la Fase 2f")

iw = sp.Symbol("iw")
NMAX = 5
amp = {(v, k): sp.Symbol(f"{v}{k}", real=True) for v in PERT for k in range(NMAX + 1)}
modo = {pj[(v, o)]: iw ** o[0] * amp[(v, o[1])] for v in PERT for o in ORDENES}


def D_lin(expr):
    out = D_l(expr)
    for v in PERT:
        for k in range(NMAX):
            s = amp[(v, k)]
            if expr.has(s):
                out += sp.diff(expr, s) * amp[(v, k + 1)]
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

with open(os.path.join(AQUI, "resultados", f"{NOMBRE}_simbolico.pkl"), "wb") as fh:
    pickle.dump(dict(fondo_eq=fondo_eq, lineales=lineales, L2=L2, GS=GS, Phi=Phi, rr=rr, PP=PP, Ee=Ee, Xs=Xs, iw=iw,
                     amp=amp, CHK=CHK), fh)
json.dump(CHK, open(os.path.join(AQUI, "resultados", f"generar_{NOMBRE}.json"), "w"), indent=1, default=str)
log(f"guardado resultados/{NOMBRE}_simbolico.pkl")
