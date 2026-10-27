# -*- coding: utf-8 -*-
"""
Fase 2l, parte I: rotacion lenta en la GLPV cuartica (y Horndeski quintico) por reduccion de Kaluza-Klein.

Metrica estacionaria-axisimetrica a primer orden en la rotacion, con el sector l = 0 general y el arrastre w(t, l):
    ds^2 = -A dt^2 + 2 b dt dl + B dl^2 + C [ dtheta^2 + sin^2 theta (dvarphi - eps w dt)^2 ],    phi = phi(t, l).
(g_{l varphi} = 0 fija el gauge de la componente axial l = 1; w es su unica funcion.)  Se calcula, de forma exacta,
    L_w = (1/2) int_0^pi dtheta [ sqrt(-g) L ]_{eps^2}            (normalizacion de la lagrangiana l = 0 de la Fase 2f:
                                                                   int_0^pi sin theta dtheta / 2 = 1)
y se comprueba que, salvo derivadas totales,
    L_w = (1/2) J(A, b, B, C, phi; derivadas) w_l^2 ,
verificando que la ecuacion de Euler-Lagrange de w es exactamente E_w = -D_l (J w_l) (identidad en el espacio de
chorros, funciones de la teoria generales).  J es el coeficiente de arrastre.  Despues:
  - limite de RG (Hartle): J = (2/3) C^2 / sqrt(A B + b^2);
  - J en el fondo estatico frente a r^4 Gcal / N (Gcal: coeficiente impar de la Fase 2g/2j).
Uso: python kk_rotacion.py [quintico]      -> resultados/kk_cuartico.pkl (o kk_quintico.pkl)
"""
import os
import sys
import time
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "fase2j"))
import q_simbolico as Q          # noqa: E402

QUINT = len(sys.argv) > 1 and sys.argv[1].startswith("q")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


t, l, th = Q.t, Q.l, Q.th
COORD = Q.COORD
eps = sp.Symbol("eps")
CAMPOS = ("A", "b", "B", "C", "p", "w")
FN = {n: sp.Function(n)(t, l) for n in CAMPOS}
A, b, B, C, p, w = (FN[n] for n in CAMPOS)
s = sp.sin(th)


def trunc(e, n=2):
    e = sp.expand(e)
    return sum(e.coeff(eps, k) * eps ** k for k in range(n + 1))


# ------------------------------------------------------------------ geometria truncada a O(eps^2)
g = sp.Matrix([[-A + eps ** 2 * C * s ** 2 * w ** 2, b, 0, -eps * C * s ** 2 * w],
               [b, B, 0, 0], [0, 0, C, 0], [-eps * C * s ** 2 * w, 0, 0, C * s ** 2]])
D2 = A * B + b ** 2
# inversa exacta (forma de Kaluza-Klein): bloque (t, l) = inversa 2D; g^{t phi} = w g^{tt}_2D, ...
gi = sp.Matrix(4, 4, lambda i, j: 0)
g2i = sp.Matrix([[-B / D2, b / D2], [b / D2, A / D2]])
for i in range(2):
    for j in range(2):
        gi[i, j] = g2i[i, j]
gi[0, 3] = gi[3, 0] = eps * w * g2i[0, 0]
gi[1, 3] = gi[3, 1] = eps * w * g2i[0, 1]
gi[2, 2] = 1 / C
gi[3, 3] = 1 / (C * s ** 2) + eps ** 2 * w ** 2 * g2i[0, 0]
assert sp.simplify(g * gi - sp.eye(4)) == sp.zeros(4, 4)
log("inversa (Kaluza-Klein) comprobada")
n = 4
dg = [[[sp.diff(g[a, c], COORD[k]) for k in range(n)] for c in range(n)] for a in range(n)]
Gam = [[[trunc(sum(gi[a, d] * (dg[d][c][k] + dg[d][k][c] - dg[c][k][d]) for d in range(n)) / 2)
         for k in range(n)] for c in range(n)] for a in range(n)]
log("Christoffel")
Ric = sp.zeros(n, n)
for c in range(n):
    for d in range(c, n):
        v = sum(sp.diff(Gam[a][c][d], COORD[a]) - sp.diff(Gam[a][c][a], COORD[d]) for a in range(n))
        v += sum(Gam[a][a][e] * Gam[e][c][d] - Gam[a][d][e] * Gam[e][c][a] for a in range(n) for e in range(n))
        Ric[c, d] = Ric[d, c] = trunc(v)
log("Ricci")
Rs = trunc(sum(gi[a, c] * Ric[a, c] for a in range(n) for c in range(n)))
log("escalar de Ricci")

d1 = [sp.diff(p, x) for x in COORD]
d2 = sp.Matrix(n, n, lambda m, k: trunc(sp.diff(p, COORD[m], COORD[k]) - sum(Gam[q][m][k] * d1[q] for q in range(n))))
up = [trunc(sum(gi[m, k] * d1[k] for k in range(n))) for m in range(n)]
Smat = sp.Matrix(n, n, lambda a, k: trunc(sum(gi[a, m] * d2[m, k] for m in range(n))))
X = trunc(-sum(up[m] * d1[m] for m in range(n)) / 2)
assert not X.has(eps)
box = trunc(sum(Smat[a, a] for a in range(n)))
S2 = (Smat * Smat).applyfunc(trunc)
trS2 = trunc(sum(S2[a, a] for a in range(n)))
pSp = trunc(sum(up[m] * d2[m, k] * up[k] for m in range(n) for k in range(n)))
pS2p = trunc(sum(up[m] * d2[m, k] * Smat[k, j] * up[j] for m in range(n) for k in range(n) for j in range(n)))
E4 = trunc(2 * X * (box ** 2 - trS2) + 2 * box * pSp - 2 * pS2p)
ESC = {"box": box, "R": Rs, "K4": trunc(box ** 2 - trS2), "E4": E4}
if QUINT:
    trS3 = trunc(sum((S2 * Smat)[a, a] for a in range(n)))
    Rup = sp.Matrix(n, n, lambda a, k: trunc(sum(gi[a, m] * Ric[m, k] for m in range(n))))
    Gphi = trunc(sum(Rup[a, k] * Smat[k, a] for a in range(n) for k in range(n)) - Rs * box / 2)
    ESC["Gphi"] = Gphi
    ESC["K5"] = trunc(box ** 3 - 3 * box * trS2 + 2 * trS3)
log("invariantes del escalar")

# ------------------------------------------------------------------ chorros
ORD = [(i, j) for i in range(7) for j in range(7) if i + j <= 6]
jet = {(nm, o): sp.Symbol(f"{nm}_{o[0]}{o[1]}", real=True) for nm in CAMPOS for o in ORD}


def a_chorros(e):
    rep = {}
    for nm in CAMPOS:
        f = FN[nm]
        for (i, j) in sorted(ORD, key=lambda o: -(o[0] + o[1])):
            if i + j == 0 or i + j > 4:
                continue
            rep[sp.Derivative(f, *([t] * i + [l] * j))] = jet[(nm, (i, j))]
    e = e.subs(rep)
    return e.subs({FN[nm]: jet[(nm, (0, 0))] for nm in CAMPOS})


Ssym, Csym = sp.symbols("S_th C_th")


def integra_theta(e):
    """(1/2) int_0^pi dtheta e, con e = polinomio en sin, cos (posibles potencias negativas de sin que deben cancelarse)."""
    e = sp.expand(e.subs({sp.sin(th): Ssym, sp.cos(th): Csym}))
    pol = sp.Poly(e, Csym)
    tot = 0
    for (k,), co in pol.terms():
        tot += co * (1 - Ssym ** 2) ** (k // 2) * Csym ** (k % 2)
    tot = sp.expand(tot)
    out = 0
    grupos = {}
    for term in sp.Add.make_args(tot):
        dep, ind = 1, 1
        for fac in sp.Mul.make_args(term):
            if fac.has(Ssym) or fac.has(Csym):
                dep *= fac
            else:
                ind *= fac
        grupos.setdefault(dep, []).append(ind)
    for dep, lst in grupos.items():
        co = sp.expand(sum(lst))
        if co == 0:
            continue
        kS = sp.degree(dep, Ssym) if dep.has(Ssym) else 0
        if dep.has(Ssym) and not sp.Poly(dep * Ssym ** 20, Ssym, Csym).is_monomial:
            raise ValueError(dep)
        pw = sp.Poly(dep * Ssym ** 20, Ssym, Csym).monoms()[0]
        kS, kC = pw[0] - 20, pw[1]
        if kC % 2 == 1:
            continue                    # int_0^pi sin^k cos = 0
        assert kS >= 0, ("potencia negativa de sin que no se cancela", kS)
        out += co * sp.integrate(sp.sin(th) ** kS, (th, 0, sp.pi)) / 2
    return sp.expand(out)


sqrtg = sp.sqrt(D2) * C * s
ESC2 = {}
for k, v in ESC.items():
    dens = sp.expand(v).coeff(eps, 2)
    ESC2[k] = sp.factor(integra_theta(a_chorros(sqrtg * dens)))
    log(f"[sqrt(-g) {k}]_eps2 integrado en theta: {len(sp.Add.make_args(ESC2[k]))} terminos")
# sin eps, para comprobar la normalizacion: (1/2) int sqrt(-g) R = lagrangiana l = 0 de la Fase 2f
R0 = integra_theta(a_chorros(sqrtg * sp.expand(Rs).coeff(eps, 0)))
g0 = Q.metrica_esferica(A, b, B, C)
ref = a_chorros(Q.lagrangiana(g0, *Q.geometria(g0), p, lambda ph, XX: (0, 0, 1, 0, 0, 0, 0, 0))[0])
ref = sp.expand(ref.subs(th, sp.pi / 3))
assert sp.simplify(R0 - ref) == 0
log("normalizacion: (1/2) int sqrt(-g) R dtheta (eps^0) = lagrangiana de la Fase 2f con G4 = 1")

# ------------------------------------------------------------------ funciones de la teoria y derivada total
Xj = a_chorros(X)
GN = ("G2", "G3", "G4", "F4", "G5")


def G(nm, i, j):
    return sp.Symbol(f"{nm}_{i}{j}", real=True)


L_w = -G("G3", 0, 0) * ESC2["box"] + G("G4", 0, 0) * ESC2["R"] + G("G4", 0, 1) * ESC2["K4"] + G("F4", 0, 0) * ESC2["E4"]
if QUINT:
    L_w += G("G5", 0, 0) * ESC2["Gphi"] - G("G5", 0, 1) / 6 * ESC2["K5"]
L_w = sp.expand(sp.simplify(L_w))
log(f"L_w: {len(sp.Add.make_args(L_w))} terminos")
import re                                   # noqa: E402
GP = re.compile(r"^(G2|G3|G4|F4|G5)_(\d)(\d)$")


def Dtot(e, k):
    """derivada total en t (k = 0) o en l (k = 1) en el espacio de chorros, G_ij(phi, X) por la regla de la cadena."""
    out = 0
    fs = e.free_symbols
    DX = None
    for (nm, (i, j)), sy in jet.items():
        if sy in fs:
            o2 = (i + 1, j) if k == 0 else (i, j + 1)
            out += sp.diff(e, sy) * jet[(nm, o2)]
    for sy in fs:
        m = GP.match(sy.name)
        if m:
            if DX is None:
                DX = Dtot(Xj, k)
            nm, i, j = m.group(1), int(m.group(2)), int(m.group(3))
            pk = jet[("p", (1, 0) if k == 0 else (0, 1))]
            out += sp.diff(e, sy) * (G(nm, i + 1, j) * pk + G(nm, i, j + 1) * DX)
    return sp.expand(out)


def euler(Lag, nm):
    tot = 0
    for (i, j) in ORD:
        sy = jet[(nm, (i, j))]
        if not Lag.has(sy):
            continue
        Qe = sp.diff(Lag, sy) * (-1) ** (i + j)
        for _ in range(i):
            Qe = Dtot(Qe, 0)
        for _ in range(j):
            Qe = Dtot(Qe, 1)
        tot += Qe
    return sp.expand(tot)


# w solo a traves de derivadas (simetria varphi -> varphi + f(t))
wj = {o: jet[("w", o)] for o in ORD}
assert not L_w.has(wj[(0, 0)]) and not L_w.has(wj[(1, 0)]) and not L_w.has(wj[(2, 0)])
mon = sorted({str(sy) for sy in L_w.free_symbols if sy.name.startswith("w_")})
log(f"chorros de w presentes: {mon}  (ni w, ni w_t, ni w_tt: invariancia gauge)")
# forma cuadratica en (w_l, w_tl, w_ll): J/2 = a - D_l(b)/2 - D_t(c)/2  (b w_l w_ll, c w_l w_tl)
wl, wtl, wll = wj[(0, 1)], wj[(1, 1)], wj[(0, 2)]
pol = sp.Poly(L_w, wl, wtl, wll)
coefs = {mo: sp.together(co) for mo, co in pol.terms()}
log("monomios: " + ", ".join(str(mo) for mo in coefs))
assert set(coefs) <= {(2, 0, 0), (1, 1, 0), (1, 0, 1)}, set(coefs)
a_ = sp.sympify(coefs.get((2, 0, 0), 0))
bb = sp.sympify(coefs.get((1, 0, 1), 0))
cc = sp.sympify(coefs.get((1, 1, 0), 0))
Jk = sp.together(2 * a_ - Dtot(bb, 1) - Dtot(cc, 0))
log(f"J: {len(sp.Add.make_args(Jk))} terminos")
if not QUINT:
    Ew = euler(L_w, "w")
    res = sp.simplify(Ew + Dtot(Jk * wl, 1))
    assert res == 0, res
    log("E_w = -D_l (J w_l)  IDENTICAMENTE: L_w = J w_l^2 / 2 + derivadas totales")
else:
    # con solo los monomios w_l^2, w_l w_ll, w_l w_tl la forma J w_l^2/2 + divergencia es una identidad algebraica:
    # b w_l w_ll = D_l(b w_l^2)/2 - (D_l b) w_l^2/2,  c w_l w_tl = D_t(c w_l^2)/2 - (D_t c) w_l^2/2
    resto = sp.expand(L_w - Jk * wl ** 2 / 2 - Dtot(bb * wl ** 2, 1) / 2 - Dtot(cc * wl ** 2, 0) / 2)
    import random
    rnd = random.Random(5)
    for k in range(4):                     # 4 chorros aleatorios, 50 cifras
        val = {sy: sp.Rational(rnd.randint(200, 900), 997) for sy in resto.free_symbols}
        r_ = sp.N(resto.xreplace(val), 50)
        esc = max(abs(sp.N(tm.xreplace(val), 50)) for tm in sp.Add.make_args(resto))
        assert abs(r_) < 1e-40 * (1 + esc), (r_, esc)
    log("L_w = J w_l^2/2 + D_l(b w_l^2)/2 + D_t(c w_l^2)/2  (comprobado en 4 chorros aleatorios, 50 cifras)")

# ------------------------------------------------------------------ simplificacion y comprobaciones
Jsim = Jk if QUINT else sp.factor(Jk)
log(f"J factorizado: {Jsim}" if len(str(Jsim)) < 1500 else f"J factorizado ({len(str(Jsim))} caracteres)")
Ar, br, Br, Cr = (jet[(nm, (0, 0))] for nm in ("A", "b", "B", "C"))
if not QUINT:
    Gcov = G("G4", 0, 0) - 2 * Xj * G("G4", 0, 1) - 4 * Xj ** 2 * G("F4", 0, 0)
    assert sp.simplify(Jk - sp.Rational(2, 3) * Cr ** 2 * Gcov / sp.sqrt(Ar * Br + br ** 2)) == 0
    log("J = (2/3) C^2 Gcal(phi, X) / sqrt(A B + b^2),  Gcal = G4 - 2 X G4_X - 4 X^2 F4  (covariante, exacto)")
gr = {sy: 0 for sy in Jk.free_symbols if GP.match(sy.name)}
gr[G("G4", 0, 0)] = 1
J_RG = sp.simplify(Jk.xreplace(gr))
assert sp.simplify(J_RG - sp.Rational(2, 3) * Cr ** 2 / sp.sqrt(Ar * Br + br ** 2)) == 0
log("limite RG: J = (2/3) C^2/sqrt(A B + b^2)  (Hartle: (r^4 j wbar')' = 0)")
# fondo estatico en longitud propia
Phi1, r1, r2, P1, P2 = sp.symbols("Phi1 r1 r2 P1 P2", real=True)
r0, Nn = sp.symbols("r0 N", positive=True)
est = {}
for sy in Jk.free_symbols:
    if sy.name.startswith(("A_", "b_", "B_", "C_", "p_", "w_")):
        est[sy] = 0
est.update({Ar: Nn ** 2, jet[("A", (0, 1))]: 2 * Phi1 * Nn ** 2, Br: 1, Cr: r0 ** 2,
            jet[("C", (0, 1))]: 2 * r0 * r1, jet[("C", (0, 2))]: 2 * r1 ** 2 + 2 * r0 * r2,
            jet[("p", (0, 1))]: P1, jet[("p", (0, 2))]: P2})
if QUINT:
    est[jet[("A", (0, 2))]] = sp.Symbol("A2", real=True)
J_est = sp.simplify(Jk.xreplace(est))
Gc4 = G("G4", 0, 0) + G("G4", 0, 1) * P1 ** 2 - G("F4", 0, 0) * P1 ** 4
dif = sp.simplify(J_est - sp.Rational(2, 3) * r0 ** 4 * Gc4 / Nn)
log(f"J estatico = (2/3) r^4 Gcal_4 / N + {dif}")
if QUINT:
    Xb_ = -P1 ** 2 / 2
    Gc5 = G("G4", 0, 0) - 2 * Xb_ * G("G4", 0, 1) - 4 * Xb_ ** 2 * G("F4", 0, 0)         + Xb_ * (Phi1 * P1 * G("G5", 0, 1) + G("G5", 1, 0))
    dif5 = sp.simplify(J_est - sp.Rational(2, 3) * r0 ** 4 * Gc5 / Nn)
    log(f"J estatico - (2/3) r^4 Gcal_5/N = {dif5}   (Gcal_5 = G4 - 2X G4X - 4X^2 F4 + X(Phi' phi' G5X + G5phi), KMS)")
    dif = dif5
out = dict(J=Jk, J_est=J_est, dif_est=dif, L_w=L_w, coefs=coefs, Xj=Xj, quintico=QUINT)
os.makedirs(os.path.join(AQUI, "resultados"), exist_ok=True)
nom = "kk_quintico" if QUINT else "kk_cuartico"
pickle.dump(out, open(os.path.join(AQUI, "resultados", nom + ".pkl"), "wb"))
json.dump({"J_estatico": str(J_est), "J_menos_(2/3)r^4Gcal4/N": str(dif), "J": str(Jsim)[:4000],
           "E_w=-D_l(J w_l)": True, "limite_RG": "J = (2/3) C^2/sqrt(AB+b^2)"},
          open(os.path.join(AQUI, "resultados", nom + ".json"), "w"), indent=1, ensure_ascii=False)
log("guardado " + nom)
