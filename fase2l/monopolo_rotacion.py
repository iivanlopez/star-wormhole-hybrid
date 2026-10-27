# -*- coding: utf-8 -*-
"""
Fase 2l, parte II: el sector l = 0 de una garganta en rotacion lenta (GLPV cuartica), a O(Omega^2).

Por kk_rotacion.py, la parte de la accion con el arrastre es L_w = J w_l^2/2 (salvo derivadas totales), con
    J = (2/3) C^2 Gcal(phi, X)/sqrt(A B + b^2).
w es ciclica: su ecuacion es D_l(J w_l) = 0 y la ligadura de g_{l varphi} hace constante  pi = J w_l  (momento angular
conservado; en el exterior de Kerr lento, pi = -4 J_Kerr).  Eliminandola a pi fijo (Routh):
    L_eff = L_l0(A, b, B, C, phi) - pi^2/(2 J)        =: L_l0 + L_J,
exacto a O(Omega^2): las deformaciones l = 2 (Hartle-Thorne) no entran en la accion l = 0 a este orden (integral
angular de Y_20 Y_00 Y_00 = 0, y axial x polar = 0).
Este script:
  1. ecuaciones del fondo de L_J y su hessiano en el fondo estatico (mismas variables que fase2f/generar_bH.py);
  2. limite de RG: el exterior l = 0 de Hartle-Thorne (m = M - J^2/r^3, h0 = J^2/(r^3 (r - 2M))) cumple las
     ecuaciones de L_RG + L_J con pi^2 = 16 J^2 a O(J^2);
  3. estructura del monopolo con rotacion, en chorros exactos de la capa de masas de L_l0 + L_J (pi^2 racional):
     k_beta = 0, a solo via F = e^{2Phi} a' + beta_t, c auxiliar, la ligadura W = gam c + al1 d' + al0 d sin
     correccion de L_J y al0 = -Phi' al1;
  4. la energia de la familia c = 0:  -L_rest,tot(0, N f) = -2 N^3 Gcal (f' - sig f)^2 + pi^2 (...) f^2  (exacta), y
     las correcciones directas de K, P, U.
Uso: python monopolo_rotacion.py        -> resultados/monopolo_rotacion.pkl, .json
"""
import os
import sys
import time
import json
import pickle
import random
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
from chorros import Phi, rr, PP, Ee, Xs, D, Valores, gsym      # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


bh = pickle.load(open(os.path.join(F2F, "resultados", "bH_simbolico.pkl"), "rb"))
fe0, lin0, L20, iw, amp = (bh[k] for k in ("fondo_eq", "lineales", "L2", "iw", "amp"))
pi2 = sp.Symbol("pi2", positive=True)                     # pi^2 = (J w_l)^2

CAMPOS = ("A", "b", "B", "C", "p")
ORDENES = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
jet = {(n, o): sp.Symbol(f"{n}_{o[0]}{o[1]}", real=True) for n in CAMPOS for o in ORDENES}
JETS = list(jet.values())
Aj, bj, Bj, Cj = (jet[(n, (0, 0))] for n in ("A", "b", "B", "C"))
pt, pl = jet[("p", (1, 0))], jet[("p", (0, 1))]
Dd = Aj * Bj + bj ** 2
Xj = (Bj * pt ** 2 - Aj * pl ** 2 - 2 * bj * pt * pl) / (2 * Dd)       # X = -(1/2) g^{ab} phi_a phi_b
# Gcal como funcion de (phi, X): Gcal = G4 - 2 X G4_X - 4 X^2 F4; sus derivadas se escriben con los G_ij de la Fase 2f
Xsym = sp.Symbol("Xg", real=True)
G4f, F4f = sp.Function("G4"), sp.Function("F4")
p0 = jet[("p", (0, 0))]
Gcal_f = G4f(p0, Xsym) - 2 * Xsym * sp.Derivative(G4f(p0, Xsym), Xsym) - 4 * Xsym ** 2 * F4f(p0, Xsym)
LJ = -pi2 / 2 * 3 * sp.sqrt(Dd) / (2 * Cj ** 2 * Gcal_f)               # -pi^2/(2 J)

fondo = {jet[("A", (0, 0))]: Ee, jet[("A", (0, 1))]: 2 * Phi[1] * Ee, jet[("A", (0, 2))]: (2 * Phi[2] + 4 * Phi[1] ** 2) * Ee,
         jet[("C", (0, 0))]: rr[0] ** 2, jet[("C", (0, 1))]: 2 * rr[0] * rr[1],
         jet[("C", (0, 2))]: 2 * rr[1] ** 2 + 2 * rr[0] * rr[2],
         jet[("B", (0, 0))]: 1, jet[("p", (0, 0))]: PP[0], jet[("p", (0, 1))]: PP[1], jet[("p", (0, 2))]: PP[2]}
for nm in CAMPOS:
    for o in ORDENES:
        fondo.setdefault(jet[(nm, o)], 0)


def a_G(e):
    """G4(P0, X), F4 y derivadas -> simbolos G_ij de la Fase 2f (X ya evaluado en el fondo)."""
    rep = {}
    for d_ in e.atoms(sp.Derivative):
        f = d_.expr
        i = sum(c for v, c in d_.variable_count if v == PP[0])
        j = sum(c for v, c in d_.variable_count if v == Xsym)
        rep[d_] = gsym(f.func.__name__, i, j)
    e = e.xreplace(rep)
    rep = {f: gsym(f.func.__name__, 0, 0) for f in e.atoms(sp.Function) if f.func in (G4f, F4f)}
    return e.xreplace(rep)


def parcial(e, j):
    """derivada parcial respecto al chorro j con X = Xj(chorros) dentro de Gcal (regla de la cadena)."""
    return sp.diff(e, j) + sp.diff(e, Xsym) * sp.diff(Xj, j)


def en_fondo(e):
    e = a_G(e.subs(fondo)).subs(Xsym, -PP[1] ** 2 / 2)
    assert not e.atoms(sp.Derivative, sp.Subs)
    return sp.expand(sp.simplify(e))


# ------------------------------------------------------------------ 1. ecuaciones del fondo y hessiano de L_J
prim = {}
for nm in CAMPOS:
    for o in ((0, 0), (0, 1), (0, 2)):
        prim[(nm, o)] = en_fondo(parcial(LJ, jet[(nm, o)]))
feJ = {nm: sp.expand(prim[(nm, (0, 0))] - D(prim[(nm, (0, 1))]) + D(D(prim[(nm, (0, 2))]))) for nm in CAMPOS}
log("ecuaciones del fondo de L_J: " + ", ".join(f"{k}: {v}" for k, v in feJ.items()))
fe = {k: sp.expand(fe0[k] + feJ[k]) for k in CAMPOS}

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


presentes = [j for j in JETS if LJ.has(j) or Xj.has(j)]
L2J = 0
for i1, j1 in enumerate(presentes):
    for i2, j2 in enumerate(presentes):
        if i2 < i1:
            continue
        h = en_fondo(parcial(parcial(LJ, j1), j2))
        if h == 0:
            continue
        n1, o1 = [k for k, v in jet.items() if v == j1][0]
        n2, o2 = [k for k, v in jet.items() if v == j2][0]
        L2J += (sp.Rational(1, 2) if i1 == i2 else 1) * h * dpert(n1, o1) * dpert(n2, o2)
L2J = sp.expand(L2J)
log(f"L2 de L_J: {len(sp.Add.make_args(L2J))} terminos; chorros: "
    f"{sorted({str(s) for s in L2J.free_symbols if s in set(pj.values())})}")

# ------------------------------------------------------------------ 2. limite de RG: Hartle-Thorne l = 0
r, M, Jk, lam = sp.symbols("r M J lam", positive=True)
m_ = M - lam * Jk ** 2 / r ** 3
h0 = lam * Jk ** 2 / (r ** 3 * (r - 2 * M))
E_ = (1 - 2 * M / r) * (1 + 2 * h0)
Phi_r = sp.diff(sp.log(E_), r) / 2
r1_ = sp.sqrt(1 - 2 * m_ / r)


def dl(e):
    return sp.diff(e, r) * r1_


hjet = {Ee: E_, Phi[1]: Phi_r * r1_, rr[0]: r, rr[1]: r1_}
hjet[Phi[2]] = dl(hjet[Phi[1]])
hjet[rr[2]] = dl(r1_)
hjet[Phi[3]] = dl(hjet[Phi[2]])
hjet[rr[3]] = dl(hjet[rr[2]])
hjet[Phi[4]] = dl(hjet[Phi[3]])
hjet[rr[4]] = dl(hjet[rr[3]])
rg = {}
for e in list(fe.values()):
    for s_ in e.free_symbols:
        if s_.name[:3] in ("G2_", "G3_", "G4_", "F4_"):
            rg[s_] = 0
rg[gsym("G4", 0, 0)] = 1
for k in range(10):
    rg[PP[k]] = 0
res_HT = {}
for nm in ("A", "B", "C"):
    e = fe[nm].xreplace(rg).subs(Xs, 0).subs(pi2, 16 * lam * Jk ** 2)
    e = e.xreplace({k: v for k, v in hjet.items()})
    ser = sp.simplify(sp.series(e, lam, 0, 2).removeO())
    res_HT[nm] = str(ser)
    assert ser == 0, (nm, ser)
log("limite RG: el exterior l = 0 de Hartle-Thorne (m = M - J^2/r^3, h0 = J^2/(r^3(r-2M))) cumple E_A = E_B = E_C = 0 "
    "a O(J^2) con pi^2 = 16 J^2")

# ------------------------------------------------------------------ 3. ecuaciones lineales de L2J y estructura
NMAX = 5
modo = {pj[(v, o)]: iw ** o[0] * amp[(v, o[1])] for v in PERT for o in ORDENES}
AMPS = tuple([amp[(v, k)] for k in range(NMAX + 1)] for v in PERT)


def lineales(L2):
    out = {}
    for v in PERT:
        tot = 0
        for o in ORDENES:
            Qe = sp.diff(L2, pj[(v, o)]).subs(modo) * (-1) ** (o[0] + o[1]) * iw ** o[0]
            for _ in range(o[1]):
                Qe = D(Qe, AMPS)
            tot += Qe
        out[v] = sp.expand(tot)
    return out


linJ = lineales(L2J)
lin = {v: sp.expand(lin0[v] + linJ[v]) for v in PERT}
uni = {amp[("q", k)]: 0 for k in range(NMAX + 1)}
EbeJ = linJ["be"].subs(uni)
assert all(EbeJ.coeff(amp[(v, k)]) == 0 for v in ("c", "d", "a") for k in range(NMAX + 1))
log(f"E_beta de L_J (gauge unitario) = {sp.factor(EbeJ)}  -> no corrige W = gam c + al1 d' + al0 d")
Ebe = sp.expand(lin["be"].subs(uni))
k_beta = Ebe.coeff(amp[("be", 0)])
resto = sp.expand(Ebe - k_beta * amp[("be", 0)])
gam = sp.expand(resto.coeff(amp[("c", 0)]) / iw)
al1 = sp.expand(resto.coeff(amp[("d", 1)]) / iw)
al0 = sp.expand(resto.coeff(amp[("d", 0)]) / iw)
assert sp.expand(resto - iw * (gam * amp[("c", 0)] + al1 * amp[("d", 1)] + al0 * amp[("d", 0)])) == 0
idl = sp.expand(al0 + Phi[1] * al1)
log(f"al0 + Phi' al1 fuera de la capa de masas: {'0 (identidad)' if idl == 0 else str(len(sp.Add.make_args(idl))) + ' terminos'}")


def chorro_rot(semilla, orden=3, pi2v=sp.Rational(3, 7)):
    v = Valores(semilla)
    v[pi2] = pi2v
    v[Ee] = sp.Rational(121, 100)
    v[Phi[1]] = v.r()
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[rr[1]] = v.r()
    v[PP[0]] = v.r()
    v[PP[1]] = v.r()
    v[Xs] = -v[PP[1]] ** 2 / 2
    g200 = gsym("G2", 0, 0)
    eB = v.ev(fe["B"], salvo=(g200,))
    c1 = eB.coeff(g200)
    v[g200] = -sp.expand(eB - c1 * g200) / c1
    eqs = [fe["A"], fe["C"], fe["p"]]
    for k in range(2, orden + 1):
        inc = [Phi[k], rr[k], PP[k]]
        ev = [v.ev(e) for e in eqs]
        Mk = sp.Matrix([[x.coeff(u) for u in inc] for x in ev])
        bk = sp.Matrix([sp.expand(x - sum(x.coeff(u) * u for u in inc)) for x in ev])
        assert all(not x.free_symbols for x in bk) and Mk.det() != 0
        for u, s_ in zip(inc, Mk.LUsolve(-bk)):
            v[u] = sp.Rational(s_)
        if k < orden:
            eqs = [D(e) for e in eqs]
    return v


Ec = sp.expand(lin["c"].subs(uni))
Ea = sp.expand(lin["a"].subs(uni))
cA1 = Ec.coeff(amp[("a", 1)])
cB0 = sp.expand(Ec.coeff(amp[("be", 0)]) / iw)
CHK = {}
for sem in range(3):
    v = chorro_rot(700 + sem)
    assert v.ev(k_beta) == 0
    assert v.ev(Ea.coeff(amp[("a", 0)])) == 0 and v.ev(Ec.coeff(amp[("a", 0)])) == 0
    assert v.ev(cA1 - Ee * cB0) == 0
    assert v.ev(idl) == 0
    # sin rotacion la misma comprobacion falla si se omite L_J (control: el chorro no es de la capa de masas de L_l0)
    kb0 = v.ev(sp.expand(lin0["be"].subs(uni)).coeff(amp[("be", 0)]))
    CHK[f"control_kbeta_sin_LJ_{sem}"] = str(kb0)
log("capa de masas de L_l0 + L_J (3 chorros exactos, pi^2 = 3/7): k_beta = 0, a solo via F, al0 = -Phi' al1; "
    f"control (sin L_J, mismo chorro): k_beta = {[CHK[k][:12] for k in CHK]}")

# ------------------------------------------------------------------ 4. familia c = 0 y correcciones de K, P, U
Lrest = sp.expand(L20.xreplace({pj[(v, o)]: 0 for v in ("a", "be", "q") for o in ORDENES}) +
                  L2J.xreplace({pj[(v, o)]: 0 for v in ("a", "be", "q") for o in ORDENES}))
LrJ = sp.expand(L2J.xreplace({pj[(v, o)]: 0 for v in ("a", "be", "q") for o in ORDENES}))
log(f"L_rest de L_J (a = beta = q = 0): {sp.factor(LrJ.subs(Xs, -PP[1] ** 2 / 2))}")
N = sp.sqrt(Ee)
f = sp.symbols("f0:4", real=True)


def Df(e):
    out = D(e)
    for k in range(3):
        out += sp.diff(e, f[k]) * f[k + 1]
    return sp.expand(out)


d0 = N * f[0]
d1 = Df(d0)
d2 = Df(d1)
cero_c = {pj[("c", o)]: 0 for o in ORDENES}
sub_d = {pj[("d", (0, 0))]: d0, pj[("d", (0, 1))]: d1, pj[("d", (0, 2))]: d2, pj[("d", (1, 0))]: 0,
         pj[("d", (2, 0))]: 0, pj[("d", (1, 1))]: 0}
menosL = sp.expand(-Lrest.xreplace(cero_c).xreplace(sub_d).subs(Xs, -PP[1] ** 2 / 2))
Gc = gsym("G4", 0, 0) + gsym("G4", 0, 1) * PP[1] ** 2 - gsym("F4", 0, 0) * PP[1] ** 4
sig = rr[1] / rr[0] - Phi[1]
centrif = sp.simplify(menosL + 2 * N ** 3 * Gc * (f[1] - sig * f[0]) ** 2)
log(f"-L_rest,tot(c = 0, d = N f) + 2 N^3 Gcal (f' - sig f)^2 = {centrif}")
esperado = 9 * pi2 * N ** 3 * f[0] ** 2 / (Gc * rr[0] ** 6)
assert sp.simplify(centrif - esperado) == 0
log("=> E[N f] = int [ -2 N^3 Gcal (f' - sig f)^2 + 9 pi^2 N^3 f^2/(Gcal r^6) ]   (exacta, puntual)")
# kinetica de la familia y correcciones de K, P, U con c = c*(d) = -(al1 d' + al0 d)/gam = -(al1/gam) N (d/N)'
cin = sp.expand(LrJ.subs(Xs, -PP[1] ** 2 / 2))
k_t = [s_ for s_ in cin.free_symbols if s_.name in ("d_10", "c_10", "d_20", "c_20", "d_11", "c_11")]
log(f"chorros temporales en L_rest de L_J: {k_t}  (=> dK directo = 0)")
gS, a1S = sp.symbols("gam al1", real=True)
c_coef = {o: cin.coeff(pj[("c", (0, 0))], 2) for o in [0]}
kcc = sp.factor(cin.coeff(pj[("c", (0, 0))], 2))
kcd = sp.factor(sp.expand(cin.coeff(pj[("c", (0, 0))], 1)).coeff(pj[("d", (0, 0))], 1))
kdd = sp.factor(cin.coeff(pj[("d", (0, 0))], 2))
log(f"L_rest de L_J = kcc c^2 + kcd c d + kdd d^2:  kcc = {kcc},  kcd = {kcd},  kdd = {kdd}")
assert sp.simplify(cin - kcc * pj[("c", (0, 0))] ** 2 - kcd * pj[("c", (0, 0))] * pj[("d", (0, 0))]
                   - kdd * pj[("d", (0, 0))] ** 2) == 0
# forma cerrada: L_J(c, d) = -(3 pi^2 N/4) sqrt(1 + 2c) / ((r^2 + 2 r d)^2 Gcal(phi, X_b/(1 + 2c)))  a segundo orden
eps_, cS, dS, Xv = sp.symbols("eps_ cS dS Xv")
Xb = -PP[1] ** 2 / 2


Xc = Xb / (1 + 2 * eps_ * cS)
G4t = sum(gsym("G4", 0, k) * (Xv - Xb) ** k / sp.factorial(k) for k in range(4))
F4t = sum(gsym("F4", 0, k) * (Xv - Xb) ** k / sp.factorial(k) for k in range(4))
Gcx = (G4t - 2 * Xv * sp.diff(G4t, Xv) - 4 * Xv ** 2 * F4t).subs(Xv, Xc)
cerr = -sp.Rational(3, 4) * pi2 * N * sp.sqrt(1 + 2 * eps_ * cS) / ((rr[0] ** 2 + 2 * rr[0] * eps_ * dS) ** 2 * Gcx)
seg = sp.diff(cerr, eps_, 2).subs(eps_, 0) / 2
assert sp.simplify(seg - cin.xreplace({pj[("c", (0, 0))]: cS, pj[("d", (0, 0))]: dS})) == 0
log("L_rest(L_J) = segundo orden de -(3 pi^2 N/4) sqrt(1+2c)/((r^2+2rd)^2 Gcal(phi, X_b/(1+2c)))  (comprobado)")
# con c = -(al1/gam)(d' - Phi' d):  -dP d'^2 - dU d^2 (salvo derivada total)
x = a1S / gS
dP = sp.factor(-kcc * x ** 2)
log(f"dP = -kcc (al1/gam)^2 = {dP};   dU = -[kdd + kcc x^2 Phi'^2 - kcd x Phi'] - (1/2)(kcd x - 2 kcc x^2 Phi')'  "
    "(x = al1/gam)")
out = dict(feJ=feJ, L2J=L2J, gam=gam, al1=al1, al0=al0, kcc=kcc, kcd=kcd, kdd=kdd, centrif=centrif)
pickle.dump(out, open(os.path.join(AQUI, "resultados", "monopolo_rotacion.pkl"), "wb"))
CHK.update({"feJ": {k: str(v) for k, v in feJ.items()}, "Hartle_Thorne": res_HT, "al0+Phi1*al1_fuera_capa": str(idl == 0),
            "centrifugo": str(centrif), "kcc": str(kcc), "kcd": str(kcd), "kdd": str(kdd), "dP": str(dP)})
json.dump(CHK, open(os.path.join(AQUI, "resultados", "monopolo_rotacion.json"), "w"), indent=1, ensure_ascii=False)
log("guardado resultados/monopolo_rotacion.pkl, .json")
