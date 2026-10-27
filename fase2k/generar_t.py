# -*- coding: utf-8 -*-
"""
Fase 2k: ecuaciones radiales (sector l = 0) de la GLPV cuartica con simetria de desplazamiento (G2, G3, G4, F4
funciones de X) -- y, con la opcion G5, Horndeski quintico G5(X) -- sobre un fondo de metrica estatica con escalar
lineal en el tiempo:
    ds^2 = -e^{2 Phi(l)} dt^2 + dl^2 + r(l)^2 dOmega^2,    phi = qt t + psi(l),    X_b = qt^2/(2 e^{2 Phi}) - psi'^2/2.
Mismo metodo que fase2f/generar_bH.py y fase2j/generar_q.py (lagrangiana de fase2j/q_simbolico.py en chorros, X como
simbolo y regla de la cadena), con estos cambios:
  * las funciones de la teoria solo dependen de X (simbolos G_0j); la lagrangiana no depende de p_00;
  * en el fondo p_10 = qt; las derivadas de la lagrangiana en el fondo no dependen de t, asi que las ecuaciones del fondo
    siguen siendo E_u = sum_J (-D_l)^{J_l} [dL/du_J] y la accion cuadratica tiene coeficientes independientes de t;
  * D_l X_b = -qt^2 Phi'/e^{2 Phi} - psi' psi''.
Comprobaciones:
  (a) identidad de Noether del fondo (reparametrizaciones de l): E_A A' + E_C C' + E_p p' - 2 D_l E_B = 0;
  (b) corriente: E_p = -D_l J con J = dL/dp_01 - D_l dL/dp_02, y la identidad de las traslaciones de t:
      qt E_p + D_l(e^{2 Phi} E_b) = 0 (fuera de la capa de masas), de modo que J = 0 resuelve E_b = E_p = 0;
  (c) limite qt = 0 (con G_ij -> G_0j): ecuaciones del fondo y accion cuadratica de la Fase 2f / 2j.
Uso: python generar_t.py [G5]      (escribe resultados/t_simbolico.pkl o t5_simbolico.pkl y el .json)
"""
import os
import sys
import json
import time
import pickle
import random
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "fase2j"))
import q_simbolico as S        # noqa: E402

CON_G5 = len(sys.argv) > 1 and sys.argv[1] == "G5"
NOMBRE = "t5" if CON_G5 else "t"
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
NOMBRES = ("G2", "G3", "G4", "F4", "G5")
FN = {n: sp.Function(n) for n in NOMBRES}


def teoria(phi, X):
    G = {n: FN[n](Xs) for n in NOMBRES}
    G5 = G["G5"] if CON_G5 else sp.Integer(0)
    G5X = sp.Derivative(G["G5"], Xs) if CON_G5 else sp.Integer(0)
    return (G["G2"], G["G3"], G["G4"], sp.Derivative(G["G4"], Xs), G["F4"], G5, G5X, sp.Integer(0))


g = S.metrica_esferica(funcs["A"], funcs["b"], funcs["B"], funcs["C"])
gi, Gam, Ric, Rs = S.geometria(g)
log("geometria")
Lden, esc = S.lagrangiana(g, gi, Gam, Ric, Rs, funcs["p"], teoria)
L = a_chorros(Lden)
Xj = a_chorros(esc["X"])
assert not any(L.has(funcs[n]) for n in CAMPOS)
assert not L.has(jet[("p", (0, 0))]), "simetria de desplazamiento"
random.seed(3)
prueba = {s: sp.Rational(random.randint(300, 900), 1000) for s in jet.values()}
prueba[Xs] = sp.Rational(-1, 7)
Lnum = L.replace(lambda e: isinstance(e, sp.Function) and e.func in tuple(FN.values()), lambda e: sp.Integer(1)) \
    .replace(lambda e: isinstance(e, sp.Derivative), lambda e: sp.Integer(1))
d_theta = sp.N((Lnum.subs(th, sp.Rational(3, 10)) - Lnum.subs(th, sp.Rational(11, 10))).subs(prueba), 30)
assert abs(d_theta) < 1e-25, d_theta
L = L.subs(th, sp.pi / 3)
log(f"lagrangiana en chorros: {sp.count_ops(L)} operaciones (independiente de theta: {float(abs(d_theta)):.1e})")

# ------------------------------------------------------------------ 2. derivadas en el fondo
Phi = sp.symbols("Phi0:9", real=True)
rr = sp.symbols("r0:9", real=True)
PP = sp.symbols("P0:10", real=True)
Ee = sp.Symbol("E", positive=True)
qt = sp.Symbol("qt", real=True)
fondo = {jet[("A", (0, 0))]: Ee, jet[("A", (0, 1))]: 2 * Phi[1] * Ee, jet[("A", (0, 2))]: (2 * Phi[2] + 4 * Phi[1] ** 2) * Ee,
         jet[("C", (0, 0))]: rr[0] ** 2, jet[("C", (0, 1))]: 2 * rr[0] * rr[1],
         jet[("C", (0, 2))]: 2 * rr[1] ** 2 + 2 * rr[0] * rr[2],
         jet[("B", (0, 0))]: 1, jet[("p", (1, 0))]: qt, jet[("p", (0, 1))]: PP[1], jet[("p", (0, 2))]: PP[2]}
for n in CAMPOS:
    for o in ORDENES:
        fondo.setdefault(jet[(n, o)], 0)
XB = qt ** 2 / (2 * Ee) - PP[1] ** 2 / 2
DXB = -qt ** 2 * Phi[1] / Ee - PP[1] * PP[2]
assert sp.simplify(Xj.subs(fondo) - XB) == 0
GS = {}


def gsym(nombre, j):
    k = (nombre, 0, j)
    if k not in GS:
        GS[k] = sp.Symbol(f"{nombre}_0{j}", real=True)
    return GS[k]


def a_simbolos_G(expr):
    rep = {}
    for d in expr.atoms(sp.Derivative):
        if isinstance(d.expr, sp.Function) and d.expr.func in tuple(FN.values()):
            j = sum(c for v, c in d.variable_count if v == Xs)
            assert j == sum(c for _, c in d.variable_count), d
            rep[d] = gsym(d.expr.func.__name__, j)
    expr = expr.xreplace(rep)
    rep = {f: gsym(f.func.__name__, 0) for f in expr.atoms(sp.Function) if f.func in tuple(FN.values())}
    out = expr.xreplace(rep)
    assert not out.atoms(sp.Derivative), out.atoms(sp.Derivative)
    return out


def en_fondo(expr):
    return a_simbolos_G(expr.subs(fondo))


Ldx = sp.diff(L, Xs)
Ldxx = sp.diff(Ldx, Xs)
JETS = [jet[(n, o)] for n in CAMPOS for o in ORDENES if (n, o) != ("p", (0, 0))]
dX = {j: sp.diff(Xj, j) for j in JETS}
Ld = {j: sp.diff(L, j) for j in JETS}
log("derivadas primeras")
prim = {}
for n in CAMPOS:
    for o in ((0, 0), (0, 1), (0, 2)):
        if (n, o) == ("p", (0, 0)):
            prim[(n, o)] = sp.Integer(0)
            continue
        j = jet[(n, o)]
        prim[(n, o)] = sp.expand(en_fondo(Ld[j] + Ldx * dX[j]))
log("primeras derivadas en el fondo")


def D_l(expr):
    out = 0
    expr = sp.sympify(expr)
    fs = expr.free_symbols
    for k in range(8):
        if Phi[k] in fs:
            out += sp.diff(expr, Phi[k]) * Phi[k + 1]
        if rr[k] in fs:
            out += sp.diff(expr, rr[k]) * rr[k + 1]
    for k in range(9):
        if PP[k] in fs:
            out += sp.diff(expr, PP[k]) * PP[k + 1]
    if Ee in fs:
        out += sp.diff(expr, Ee) * 2 * Phi[1] * Ee
    if Xs in fs:
        out += sp.diff(expr, Xs) * DXB
    for (nom, i, j), s in list(GS.items()):
        if s in fs:
            out += sp.diff(expr, s) * gsym(nom, j + 1) * DXB
    return out


fondo_eq = {}
for n in CAMPOS:
    fondo_eq[n] = sp.expand(prim[(n, (0, 0))] - D_l(prim[(n, (0, 1))]) + D_l(D_l(prim[(n, (0, 2))])))
J = sp.expand(prim[("p", (0, 1))] - D_l(prim[("p", (0, 2))]))
log("ecuaciones de fondo: " + ", ".join(f"{n}: {sp.count_ops(v)}" for n, v in fondo_eq.items()) + f", J: {sp.count_ops(J)}")
CHK = {}
sXB = {Xs: XB}
noether = sp.expand((fondo_eq["A"] * 2 * Phi[1] * Ee + fondo_eq["C"] * 2 * rr[0] * rr[1] + fondo_eq["p"] * PP[1]
                     - 2 * D_l(fondo_eq["B"])))
CHK["noether_l"] = "0 (identica)" if noether == 0 else str(sp.count_ops(noether))
assert noether == 0, sp.count_ops(noether)
log("(a) identidad de Noether del fondo (l): se anula identicamente")
assert sp.expand(fondo_eq["p"] + D_l(J)) == 0
nt_ = sp.expand(qt * fondo_eq["p"] + D_l(Ee * fondo_eq["b"]))
log(f"(b) E_p = -D_l J identica;  qt E_p + D_l(e^(2Phi) E_b) = {nt_}")
assert nt_ == 0
# relacion algebraica: e^{2 Phi} E_b = qt J (la constante de integracion es cero)
rel = sp.expand(Ee * fondo_eq["b"] - qt * J)
CHK["Eb_igual_qJ"] = "e^{2Phi} E_b = qt J (identica)" if rel == 0 else str(sp.count_ops(rel))
log(f"(b') e^(2Phi) E_b - qt J = {'0' if rel == 0 else sp.count_ops(rel)}")
# (c) limite qt = 0 frente a la Fase 2f (cuartica) o 2j (quintica): G_ij con i > 0 -> 0
ref = pickle.load(open(os.path.join(AQUI, "..", "fase2j" if CON_G5 else "fase2f", "resultados",
                                    "q_simbolico.pkl" if CON_G5 else "bH_simbolico.pkl"), "rb"))
cero_ref = {s: 0 for e in list(ref["fondo_eq"].values()) + [ref["L2"]] for s in e.free_symbols
            if s.name[:3] in ("G2_", "G3_", "G4_", "F4_", "G5_", "F5_") and (s.name[3] != "0" or s.name[:2] == "F5")}
ren = {s: sp.Symbol(s.name, real=True) for e in ref["fondo_eq"].values() for s in e.free_symbols}
for n in CAMPOS:
    dif = sp.expand(fondo_eq[n].subs(qt, 0) - ref["fondo_eq"][n].xreplace(cero_ref))
    assert dif == 0, (n, sp.count_ops(dif))
CHK["limite_q0_fondo"] = "identico a la fase de referencia"
log("(c) limite qt = 0: ecuaciones del fondo identicas a las de la fase de referencia")

# ------------------------------------------------------------------ 3. accion cuadratica
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
PUNTO = os.path.join(AQUI, "resultados", f"{NOMBRE}_L2_punto.pkl")
if os.path.exists(PUNTO):
    L2, n_nz = pickle.load(open(PUNTO, "rb"))
    JETS_H = []
    log("accion cuadratica leida del punto de control")
else:
    JETS_H = JETS
for i1, j1 in enumerate(JETS_H):
    for i2, j2 in enumerate(JETS_H):
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
    log(f"   hessiano: fila {i1 + 1}/{len(JETS_H)}")
L2 = sp.expand(L2)
os.makedirs(os.path.join(AQUI, "resultados"), exist_ok=True)
pickle.dump((L2, n_nz), open(PUNTO, "wb"))
log(f"accion cuadratica: {n_nz} entradas no nulas del hessiano, {len(sp.Add.make_args(L2))} terminos")
dif = sp.expand(L2.subs(qt, 0) - ref["L2"].xreplace(cero_ref))
assert dif == 0, sp.count_ops(dif)
CHK["limite_q0_L2"] = "identica a la fase de referencia"
log("(c') limite qt = 0: accion cuadratica identica a la de la fase de referencia")

iw = sp.Symbol("iw")
NMAX = 5
amp = {(v, k): sp.Symbol(f"{v}{k}", real=True) for v in PERT for k in range(NMAX + 1)}
modo = {pj[(v, o)]: iw ** o[0] * amp[(v, o[1])] for v in PERT for o in ORDENES}


def D_lin(expr):
    expr = sp.sympify(expr)
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

os.makedirs(os.path.join(AQUI, "resultados"), exist_ok=True)
with open(os.path.join(AQUI, "resultados", f"{NOMBRE}_simbolico.pkl"), "wb") as fh:
    pickle.dump(dict(fondo_eq=fondo_eq, J=J, lineales=lineales, L2=L2, GS=GS, Phi=Phi, rr=rr, PP=PP, Ee=Ee, Xs=Xs,
                     qt=qt, XB=XB, DXB=DXB, iw=iw, amp=amp, CHK=CHK), fh)
json.dump(CHK, open(os.path.join(AQUI, "resultados", f"generar_{NOMBRE}.json"), "w"), indent=1, default=str)
log(f"guardado resultados/{NOMBRE}_simbolico.pkl")
