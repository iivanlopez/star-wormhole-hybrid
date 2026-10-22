# -*- coding: utf-8 -*-
"""
Ecuaciones radiales (sector l = 0) de la teoria cuartica mas alla de Horndeski a partir de la lagrangiana reducida
(bH_simbolico.py).  Genera bH_gen.py.

Metodo
  1. Lred(jets, X) con u = (A, b, B, C, phi) y sus derivadas hasta segundo orden en (t, l) como simbolos de chorro;
     G2, G3, G4, F4 son funciones generales de (phi, X) y X se trata como simbolo: las derivadas parciales usan la
     regla de la cadena d/du_J [L(J, X(J))] = L_J + L_X X_J (y la segunda variacion analoga).
  2. Fondo estatico en longitud propia: A = e^{2 Phi(l)}, b = 0, B = 1, C = r(l)^2, phi = P(l), X_b = -P'^2/2.
     Ecuaciones de fondo  E_u = sum_J (-D_l)^{J_l} [dL/du_J]_fondo  (solo chorros sin derivadas temporales).
  3. Accion cuadratica  L2 = (1/2) sum_JK [d^2 L/du_J du_K]_fondo du_J du_K  con perturbaciones lineales
     dA = 2 e^{2 Phi} a, db = beta, dB = 2 c, dC = 2 r d, dphi = q;  ecuaciones lineales = Euler-Lagrange de L2.
     Modos e^{i omega t}: D_t -> i omega, w = omega^2.  Gauge de longitud propia: beta = c = 0 despues de variar.
  D_l actua sobre los simbolos de fondo (Phi_k, r_k, P_k -> k+1), X_b -> -P1 P2 y G_ij -> G_(i+1)j P1 - G_i(j+1) P1 P2
  (i derivadas en phi, j en X).
Uso: python generar_bH.py      (escribe resultados/bH_simbolico.pkl y resultados/generar_bH.json; el codigo numerico
                                bH_gen.py lo genera despues exportar_bH.py)
"""
import os
import sys
import json
import time
import itertools
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import bH_simbolico as S        # noqa: E402

t, l, th = S.t, S.l, S.th
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f} s] {msg}", flush=True)


CAMPOS = ("A", "b", "B", "C", "p")
ORDENES = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]          # (nt, nl)
funcs = {n: sp.Function(n)(t, l) for n in CAMPOS}
jet = {(n, o): sp.Symbol(f"{n}_{o[0]}{o[1]}", real=True) for n in CAMPOS for o in ORDENES}


def a_chorros(expr):
    """Derivative(u, t^nt, l^nl) -> u_ntnl, u -> u_00 (de mayor a menor orden)."""
    rep = {}
    for n in CAMPOS:
        f = funcs[n]
        for (nt, nl) in sorted(ORDENES, key=lambda o: -(o[0] + o[1])):
            if nt + nl == 0:
                continue
            vars_ = [t] * nt + [l] * nl
            rep[sp.Derivative(f, *vars_)] = jet[(n, (nt, nl))]
    e = expr.subs(rep)
    return e.subs({funcs[n]: jet[(n, (0, 0))] for n in CAMPOS})


# ------------------------------------------------------------------ 1. lagrangiana reducida en chorros
Xs = sp.Symbol("Xs", real=True)
ps = jet[("p", (0, 0))]
G2f, G3f, G4f, F4f = (sp.Function(n) for n in ("G2", "G3", "G4", "F4"))


def teoria(phi, X):
    # se construye con phi como funcion; tras pasar a chorros phi -> p_00.  X se sustituye por el simbolo Xs.
    return (G2f(phi, Xs), G3f(phi, Xs), G4f(phi, Xs), sp.Derivative(G4f(phi, Xs), Xs), F4f(phi, Xs))


g = S.metrica_esferica(funcs["A"], funcs["b"], funcs["B"], funcs["C"])
gi, Gam, Rs = S.geometria(g)
log("geometria")
Lden, esc = S.lagrangiana(g, gi, Gam, Rs, funcs["p"], teoria)
L = a_chorros(Lden)
Xj = a_chorros(esc["X"])
assert not L.atoms(sp.Derivative) - {d for d in L.atoms(sp.Derivative) if d.expr.func in (G4f,)}
assert not any(L.has(funcs[n]) for n in CAMPOS)
# la densidad reducida no depende de theta: se comprueba en un chorro racional aleatorio y se fija theta
import random                                       # noqa: E402
random.seed(3)
prueba = {s: sp.Rational(random.randint(300, 900), 1000) for s in jet.values()}
prueba[Xs] = sp.Rational(-1, 7)
Lnum = L.replace(lambda e: isinstance(e, sp.Function) and e.func in (G2f, G3f, G4f, F4f), lambda e: sp.Integer(1)) \
    .replace(lambda e: isinstance(e, sp.Derivative), lambda e: sp.Integer(1))
d_theta = sp.N((Lnum.subs(th, sp.Rational(3, 10)) - Lnum.subs(th, sp.Rational(11, 10))).subs(prueba), 30)
assert abs(d_theta) < 1e-25, d_theta
L = L.subs(th, sp.pi / 3)
log(f"lagrangiana en chorros: {sp.count_ops(L)} operaciones (independiente de theta: {float(abs(d_theta)):.1e})")

# ------------------------------------------------------------------ 2. derivadas en el fondo estatico
Phi = sp.symbols("Phi0:6", real=True)
rr = sp.symbols("r0:6", real=True)
PP = sp.symbols("P0:7", real=True)
Ee = sp.Symbol("E", positive=True)                      # e^{2 Phi0}
fondo = {jet[("A", (0, 0))]: Ee, jet[("A", (0, 1))]: 2 * Phi[1] * Ee, jet[("A", (0, 2))]: (2 * Phi[2] + 4 * Phi[1] ** 2) * Ee,
         jet[("C", (0, 0))]: rr[0] ** 2, jet[("C", (0, 1))]: 2 * rr[0] * rr[1],
         jet[("C", (0, 2))]: 2 * rr[1] ** 2 + 2 * rr[0] * rr[2],
         jet[("B", (0, 0))]: 1, jet[("p", (0, 0))]: PP[0], jet[("p", (0, 1))]: PP[1], jet[("p", (0, 2))]: PP[2]}
for n in CAMPOS:
    for o in ORDENES:
        fondo.setdefault(jet[(n, o)], 0)
GS = {}                                                    # simbolos G_ij (i derivadas en phi, j en X)


def gsym(nombre, i, j):
    k = (nombre, i, j)
    if k not in GS:
        GS[k] = sp.Symbol(f"{nombre}_{i}{j}", real=True)
    return GS[k]


def a_simbolos_G(expr):
    """G(P0, Xs) y sus derivadas -> simbolos G_ij.  Primero las derivadas completas (xreplace estructural: replace()
    sustituiria antes la funcion interior y dejaria Derivative(simbolo, P0) sin evaluar), luego las funciones."""
    rep = {}
    for d in expr.atoms(sp.Derivative):
        if isinstance(d.expr, sp.Function) and d.expr.func in (G2f, G3f, G4f, F4f):
            i = sum(c for v, c in d.variable_count if v == PP[0])
            j = sum(c for v, c in d.variable_count if v == Xs)
            assert i + j == sum(c for _, c in d.variable_count), d
            rep[d] = gsym(d.expr.func.__name__, i, j)
    expr = expr.xreplace(rep)
    rep = {f: gsym(f.func.__name__, 0, 0) for f in expr.atoms(sp.Function) if f.func in (G2f, G3f, G4f, F4f)}
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
# primeras derivadas completas (regla de la cadena por X) en el fondo, solo chorros sin derivada temporal
prim = {}
for n in CAMPOS:
    for o in ((0, 0), (0, 1), (0, 2)):
        j = jet[(n, o)]
        prim[(n, o)] = sp.expand(en_fondo(Ld[j] + Ldx * dX[j]))
log("primeras derivadas en el fondo")


# ------------------------------------------------------------------ operador D_l sobre expresiones de fondo
def D_l(expr):
    out = 0
    for k in range(5):
        out += sp.diff(expr, Phi[k]) * Phi[k + 1] + sp.diff(expr, rr[k]) * rr[k + 1]
    for k in range(6):
        out += sp.diff(expr, PP[k]) * PP[k + 1]
    out += sp.diff(expr, Ee) * 2 * Phi[1] * Ee
    out += sp.diff(expr, Xs) * (-PP[1] * PP[2])
    for (nom, i, j), s in list(GS.items()):
        if expr.has(s):
            out += sp.diff(expr, s) * (gsym(nom, i + 1, j) * PP[1] - gsym(nom, i, j + 1) * PP[1] * PP[2])
    return out


# ecuaciones de fondo E_u = L_u - D_l L_{u_l} + D_l^2 L_{u_ll}
fondo_eq = {}
for n in CAMPOS:
    e = prim[(n, (0, 0))] - D_l(prim[(n, (0, 1))]) + D_l(D_l(prim[(n, (0, 2))]))
    fondo_eq[n] = sp.expand(e)
log("ecuaciones de fondo: " + ", ".join(f"{n}: {sp.count_ops(v)}" for n, v in fondo_eq.items()))

# ------------------------------------------------------------------ comprobaciones de las ecuaciones de fondo
sys.path.append(os.path.join(AQUI, "..", "fase2e"))
import numpy as np                     # noqa: E402
import potencial as PV                 # noqa: E402  (EsGB + V de la Fase 2e; con F = 0 es RG + escalar + V)
CHK = {}
# (a) identidad de Noether de las reparametrizaciones de l (fuera de la capa de masas): en B = 1,
#     E_A A' + E_C C' + E_p p' - 2 D_l E_B = 0
noether = fondo_eq["A"] * 2 * Phi[1] * Ee + fondo_eq["C"] * 2 * rr[0] * rr[1] + fondo_eq["p"] * PP[1] \
    - 2 * D_l(fondo_eq["B"])
random.seed(11)
val = {s: sp.Rational(random.randint(-900, 900), 997) for s in list(Phi) + list(rr) + list(PP) + list(GS.values())}
val[Ee] = sp.Rational(13, 10)
val[rr[0]] = sp.Rational(17, 10)
val[Xs] = -val[PP[1]] ** 2 / 2
for k, s in GS.items():
    val.setdefault(s, sp.Rational(random.randint(-900, 900), 997))
noether = sp.expand(noether)
if noether == 0:
    CHK["noether_fondo"] = 0.0
    log("identidad de Noether del fondo: se anula identicamente (simbolica)")
else:
    res_n = sp.N(noether.subs(val), 30)
    esc_n = sp.N(sum(abs(sp.N(tm.subs(val), 30)) for tm in sp.Add.make_args(noether)), 30)
    CHK["noether_fondo"] = float(abs(res_n) / esc_n)
    log(f"identidad de Noether del fondo (chorro aleatorio): residuo relativo {CHK['noether_fondo']:.1e}")
assert CHK["noether_fondo"] < 1e-25
# (b) limite RG + escalar canonico + potencial: G4 = 1, G2 = X - V(phi), G3 = F4 = 0
V0, V1, V2 = sp.symbols("V0 V1 V2", real=True)
rg = {s: 0 for s in GS.values()}
rg.update({gsym("G4", 0, 0): 1, gsym("G2", 0, 0): Xs - V0, gsym("G2", 0, 1): 1, gsym("G2", 1, 0): -V1,
           gsym("G2", 2, 0): -V2})
eqs_rg = [sp.expand(fondo_eq[n].subs(rg).subs(Xs, -PP[1] ** 2 / 2)) for n in ("A", "C", "p")]
sol_rg = sp.solve(eqs_rg, [Phi[2], rr[2], PP[2]], dict=True)
assert len(sol_rg) == 1
f_sol = sp.lambdify((Phi[1], rr[0], rr[1], PP[0], PP[1], V0, V1), [sol_rg[0][x] for x in (rr[2], Phi[2], PP[2])])
lig_rg = sp.expand(fondo_eq["B"].subs(rg).subs(Xs, -PP[1] ** 2 / 2))
f_lig = sp.lambdify((Ee, Phi[1], rr[0], rr[1], PP[0], PP[1], V0), lig_rg)
rng = np.random.default_rng(5)
err = 0.0
for _ in range(20):
    Ph1, r0_, r1_, p0_, p1_ = rng.uniform(-1, 1), rng.uniform(0.5, 2), rng.uniform(-0.9, 0.9), rng.uniform(-1, 1), \
        rng.uniform(-1, 1)
    m2 = rng.uniform(0.1, 2)
    Vd = PV.potencial("masivo", m2)
    Fd0 = PV.acoplamiento("cuadratico", 0.0)
    y = np.array([r0_, r1_, 0.0, Ph1, p0_, p1_])
    ref, _ = PV.segundas_V(y, Fd0, Vd)
    nuestra = np.array(f_sol(Ph1, r0_, r1_, p0_, p1_, Vd(p0_)[0], Vd(p0_)[1]), dtype=float)
    err = max(err, float(np.max(np.abs(nuestra - ref) / (1 + np.abs(ref)))))
    # la ligadura: proporcional a est_ligadura_V; se compara el cociente en dos puntos
CHK["limite_RG_segundas"] = err
log(f"limite RG + escalar + V frente a fase2e/potencial.segundas_V (20 estados): {err:.1e}")
lig_ref = []
for _ in range(6):
    Ph1, r0_, r1_, p0_, p1_ = rng.uniform(-1, 1), rng.uniform(0.5, 2), rng.uniform(-0.9, 0.9), rng.uniform(-1, 1), \
        rng.uniform(-1, 1)
    Vd = PV.potencial("masivo", 0.7)
    y = np.array([r0_, r1_, 0.0, Ph1, p0_, p1_])
    a_ = float(f_lig(1.0, Ph1, r0_, r1_, p0_, p1_, Vd(p0_)[0]))
    b_ = PV.est_ligadura_V(y, Fd0, Vd)
    lig_ref.append(a_ / (b_ * r0_**2))
CHK["limite_RG_ligadura_cociente"] = [float(x) for x in lig_ref]
log(f"ligadura l l: E_B / (est_ligadura_V r^2) en 6 estados = {[round(x, 12) for x in lig_ref]}")

# ------------------------------------------------------------------ 3. accion cuadratica y ecuaciones lineales l = 0
PERT = ("a", "be", "c", "d", "q")
pj = {(v, o): sp.Symbol(f"{v}_{o[0]}{o[1]}", real=True) for v in PERT for o in ORDENES}
E2 = 2 * Ee


def dpert(n, o):
    """chorro (n, o) de la perturbacion lineal del campo n en funcion de los chorros de (a, be, c, d, q)."""
    nt, nl = o
    if n == "A":          # dA = 2 e^{2 Phi} a
        if nl == 0:
            return E2 * pj[("a", o)]
        if nl == 1:
            return E2 * (2 * Phi[1] * pj[("a", (nt, 0))] + pj[("a", (nt, 1))])
        return E2 * ((2 * Phi[2] + 4 * Phi[1] ** 2) * pj[("a", (0, 0))] + 4 * Phi[1] * pj[("a", (0, 1))] + pj[("a", (0, 2))])
    if n == "b":
        return pj[("be", o)]
    if n == "B":
        return 2 * pj[("c", o)]
    if n == "C":          # dC = 2 r d
        if nl == 0:
            return 2 * rr[0] * pj[("d", o)]
        if nl == 1:
            return 2 * (rr[1] * pj[("d", (nt, 0))] + rr[0] * pj[("d", (nt, 1))])
        return 2 * (rr[2] * pj[("d", (0, 0))] + 2 * rr[1] * pj[("d", (0, 1))] + rr[0] * pj[("d", (0, 2))])
    return pj[("q", o)]


H2 = {}
dXX = {}
Ldj_x = {j: sp.diff(Ld[j], Xs) for j in JETS}
n_nz = 0
L2 = 0
for i1, j1 in enumerate(JETS):
    for i2, j2 in enumerate(JETS):
        if i2 < i1:
            continue
        dxx = sp.diff(dX[j1], j2)
        h = sp.diff(Ld[j1], j2) + Ldj_x[j1] * dX[j2] + Ldj_x[j2] * dX[j1] + Ldxx * dX[j1] * dX[j2] + Ldx * dxx
        hb = sp.expand(en_fondo(h))
        if hb == 0:
            continue
        n_nz += 1
        n1, o1 = [k for k, v in jet.items() if v == j1][0]
        n2, o2 = [k for k, v in jet.items() if v == j2][0]
        fac = sp.Rational(1, 2) if i1 == i2 else 1
        L2 += fac * hb * dpert(n1, o1) * dpert(n2, o2)
log(f"hessiano en el fondo: {n_nz} entradas no nulas")
L2 = sp.expand(L2)
log(f"accion cuadratica: {len(sp.Add.make_args(L2))} terminos")

# Euler-Lagrange para modos e^{i omega t}: chorro (nt, nl) -> (i omega)^nt V_nl
iw = sp.Symbol("iw")                                   # i omega;  iw^2 = -w
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
        Q = sp.diff(L2, pj[(v, o)]).subs(modo)
        Q = Q * (-1) ** (o[0] + o[1]) * iw ** o[0]
        for _ in range(o[1]):
            Q = D_lin(Q)
        tot += Q
    lineales[v] = sp.expand(tot)
    log(f"ecuacion lineal E_{v}: {len(sp.Add.make_args(lineales[v]))} terminos")

# ------------------------------------------------------------------ guardar
import pickle                                         # noqa: E402
os.makedirs(os.path.join(AQUI, "resultados"), exist_ok=True)
with open(os.path.join(AQUI, "resultados", "bH_simbolico.pkl"), "wb") as fh:
    pickle.dump(dict(fondo_eq=fondo_eq, lineales=lineales, L2=L2, GS=GS, Phi=Phi, rr=rr, PP=PP, Ee=Ee, Xs=Xs, iw=iw,
                     amp=amp, CHK=CHK), fh)
CHK["tamanos"] = {"fondo": {n: int(sp.count_ops(v)) for n, v in fondo_eq.items()},
                  "lineales": {v: len(sp.Add.make_args(e)) for v, e in lineales.items()}, "hessiano_no_nulas": n_nz}
json.dump(CHK, open(os.path.join(AQUI, "resultados", "generar_bH.json"), "w"), indent=1)
log("guardado resultados/bH_simbolico.pkl y resultados/generar_bH.json")
