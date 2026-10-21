# -*- coding: utf-8 -*-
"""
Deduccion analitica de la ley de §12: limite de zona cercana (a r0 >> 1) de las ecuaciones de §2 (egb_F_gen) con la
capa de §3 (estabilidad_F.CapaF / matriz_garganta), en SymPy.

Escalas (r0 = 1).  X = a l, Z = 1 + X.  Fondo: r = 1, r' = S(X)/a, Phi' = a P(X), phi' = U(X), F'(phi) = G(X)/a;
F'' = F2 y F''' = F3 constantes en la zona.  Perturbaciones: a' = a^2 A(X), d = D(X)/a, q = Q(X), w = -a^2 W.

 1. Fondo de zona cercana a orden dominante (se comprueba en las cuatro ecuaciones de fondo):
      P = 1/Z (lapso de Rindler N = Z),  U = u/Z,  G = G0 + F2 u ln Z con G0 = (1 + u^2/4)/(4u)  (ligadura en l = 0),
      S = [Z - (1 + 4 F2 u^2 ln Z)/Z]/2.
 2. Ecuaciones lineales a orden dominante y su solucion general (se comprueba que las anula):
      A = A0/Z^2,  Q = u A0/((W - 1) Z) + c+ Z^s + c- Z^(-s)  (s^2 = W),
      V = (4 + u^2) Q + 8 u D = [2u^2 (Q + Z Q') - u A0 (4 Z^2 + u^2)/Z]/(1 - W) - 16 F2 u^2 ln Z Q.
    El desplazamiento de la capa respecto del fondo estatico (modo gauge de la Fase 2c) es A0 = f (W - 1), c+- = 0.
 3. Empalme con la zona lejana: la solucion que decae en el infinito llega a la zona cercana como
    alpha rho^(-s) + beta rho^s (rho = Z/a, distancia al horizonte frustrado), asi que c+/c- = O((a r0)^(-2s)) y a orden
    dominante c+ = 0.  La condicion de contorno es local.
 4. Empalmes de la capa a orden dominante: E^th_th en a^2 y E_phi en a (E^t_t y las ligaduras se anulan
    identicamente).  Sistema 2x2 en (f, c-), cuyo determinante es un polinomio P(s; zeta, u, F0, F2).
 5. Limites: zeta = 0, coeficiente principal (zeta_c) y limite de agujero negro (u, F0, F2 u^2 -> 0).

Salida: resultados/ley_analitica.json (polinomio y limites).   Uso: python ley_analitica.py   (~3 min)
"""
import sys
import json
import sympy as sp
import egb_F_gen as H

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

X = sp.symbols("X", nonnegative=True)
a, u, F2, F3, F0, zeta, mu, W, f, c = sp.symbols("a u F2 F3 F0 zeta mu W f c", real=True)
s = sp.symbols("s", positive=True)
Z = 1 + X
L = sp.log(Z)
G0 = (1 + u**2 / 4) / (4 * u)


def orden_dominante(e, nmax=4, nmin=-3):
    e = sp.expand(e)
    for n in range(nmax, nmin - 1, -1):
        cc = sp.simplify(sp.limit(e / a**n, a, sp.oo))
        if cc != 0:
            return n, cc
    return None, sp.Integer(0)


# ---------------------------------------------------------------- 1. fondo de zona cercana
P, U, G = 1 / Z, u / Z, G0 + F2 * u * L
S = sp.Rational(1, 2) * (Z - (1 + 4 * F2 * u**2 * L) / Z)
e2 = 1 / Z**2
r, r1, r2 = sp.Integer(1), S / a, sp.diff(S, X)
P1, P2 = a * P, a**2 * sp.diff(P, X)
p1, p2 = U, a * sp.diff(U, X)
F1 = G / a

print("1. Fondo de zona cercana: orden dominante de las ecuaciones de fondo (tt, thth, phi, ligadura ll)")
Mb = sp.Matrix(H.est_principal(F1, r, r1, P1, p1))
eqs_fondo = list(Mb * sp.Matrix([r2, P2, p2]) + sp.Matrix(H.est_resto(F1, F2, r, r1, P1, p1)))
eqs_fondo.append(H.est_ligadura(F1, r, r1, P1, p1))
for nombre, e in zip(("tt", "thth", "phi", "ligadura"), eqs_fondo):
    e = sp.expand(e)
    ords = [sp.simplify(sp.limit(e / a**n, a, sp.oo)) if n == 2 else None for n in (2,)]
    lead2 = ords[0]
    lead1 = sp.simplify(sp.limit((e - lead2 * a**2) / a, a, sp.oo))
    lead0 = sp.simplify(sp.limit(e - lead2 * a**2 - lead1 * a, a, sp.oo))
    print(f"   {nombre:9s}: a^2 -> {lead2},  a^1 -> {lead1},  a^0 -> {'0' if lead0 == 0 else 'no nulo (correccion O(1/a) de Phi)'}")

# ---------------------------------------------------------------- 2. ecuaciones lineales a orden dominante
Af, Df, Qf = [sp.Function(n)(X) for n in ("A", "D", "Q")]
a1, d0, d1, q0, q1 = a**2 * Af, Df / a, sp.diff(Df, X), Qf, a * sp.diff(Qf, X)
d2, a2, q2 = a * sp.diff(Df, X, 2), a**3 * sp.diff(Af, X), a**2 * sp.diff(Qf, X, 2)
w = -a**2 * W
Mp = sp.Matrix(H.lin_principal(F1, r, r1, P1, p1))
Rs = sp.Matrix(H.lin_resto(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1))
lineales = list(Mp * sp.Matrix([d2, a2, q2]) + Rs)
lineales += list(H.lin_ligaduras(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1))
nombres_lin = ("E^t_t", "E^th_th", "E_phi", "E^l_t (ligadura)", "E^l_l (ligadura)")

A0 = sp.symbols("A0")
cp, cm = sp.symbols("c_p c_m")
Qsol = u * A0 / ((W - 1) * Z) + cp * Z**s + cm * Z**(-s)
Asol = A0 / Z**2
Vsol = (2 * u**2 * (Qsol + Z * sp.diff(Qsol, X)) - u * A0 * (4 * Z**2 + u**2) / Z) / (1 - W) - 16 * F2 * u**2 * L * Qsol
Dsol = (Vsol - (4 + u**2) * Qsol) / (8 * u)
print("\n2. Ecuaciones lineales a orden dominante y solucion general de zona cercana:")
for nombre, e in zip(nombres_lin, lineales):
    n, cc = orden_dominante(e)
    res = cc.subs({Af: Asol, Df: Dsol, Qf: Qsol}).doit()
    res = sp.simplify(sp.expand(res.subs(W, s**2)))
    print(f"   {nombre:18s}: orden a^{n};  residuo de la solucion general = {res}")

# modo de desplazamiento (gauge): A0 = f (W - 1), c+- = 0  ->  coincide con (Phi'' - w e^-2Phi, r', phi')
g = {Asol: (f * (W - 1) / Z**2), Qsol: f * U}
Dg = sp.simplify(Dsol.subs({A0: f * (W - 1), cp: 0, cm: 0}))
print(f"   modo de desplazamiento: D = {sp.simplify(Dg)}  (f S = {sp.simplify(f * S)});"
      f"  Q = {sp.simplify(Qsol.subs({A0: f * (W - 1), cp: 0, cm: 0}))}")

# ---------------------------------------------------------------- 4. empalmes en l = 0+
sol = {A0: f * (W - 1), cp: 0, cm: c}
Q = Qsol.subs(sol)
Afin = Asol.subs(sol)
D = Dsol.subs(sol)
y0 = dict(a1=a**2 * Afin, d0=D / a, d1=sp.diff(D, X), q0=Q, q1=a * sp.diff(Q, X))
y0 = {k: sp.simplify(v.subs(X, 0)) for k, v in y0.items()}
bg0 = dict(r=sp.Integer(1), r1=sp.Integer(0), r2=sp.simplify(r2.subs(X, 0)), P1=a, P2=-a**2, p1=u, p2=-a * u,
           e2=sp.Integer(1))
F1_0 = G0 / a


def capa(zeta):
    """CapaF a orden dominante: kappa = 4 zeta a, Lambda_s = -2a; rho y beta de los empalmes del fondo."""
    kappa, Lam = 4 * zeta * a, -2 * a
    M0 = sp.Matrix(H.est_principal(F1_0, bg0["r"], bg0["r1"], bg0["P1"], bg0["p1"]))
    dd = M0 * sp.Matrix([0, 2 * a, 2 * u])
    rho = -dd[0] + 2 * kappa * F0 - Lam
    beta = (dd[2] + 4 * kappa * F1_0) / (2 * rho)
    assert sp.simplify(dd[1] + Lam) == 0
    return kappa, rho, beta


def fuentes(zeta, d0, q0):
    kappa, rho, beta = capa(zeta)
    Stt = rho * (2 * d0 - beta * q0) + 2 * kappa * (F1_0 * q0 - 2 * F0 * d0)
    Sthth = -2 * kappa * w * (F0 * d0 + F1_0 * q0)
    dR = -4 * d0 - 4 * w * d0
    sphi = -2 * kappa * (F2 * q0 * 2 + F1_0 * dR) + 2 * rho * beta * (-2 * d0 + mu * q0)
    return Stt, Sthth, sphi


rr, rr1, rr2, PP1, PP2, pp1, pp2, ee2 = (bg0[k] for k in ("r", "r1", "r2", "P1", "P2", "p1", "p2", "e2"))
M0 = sp.Matrix(H.lin_principal(F1_0, rr, rr1, PP1, pp1))
Pp, pp = sp.symbols("Pp pp")
R = sp.Matrix(H.lin_resto(F1_0, F2, F3, w, rr, rr1, rr2, PP1, Pp, pp1, pp, ee2, *(y0[k] for k in ("a1", "d0", "d1", "q0", "q1"))))
R0 = R.subs({Pp: 0, pp: 0})
delta = M0 * sp.Matrix([2 * y0["d1"], 2 * y0["a1"], 2 * y0["q1"]]) + (R.subs({Pp: 1, pp: 0}) - R0) * 2 * PP1 \
    + (R.subs({Pp: 0, pp: 1}) - R0) * 2 * pp1
lig0 = H.lin_ligaduras(F1_0, F2, F3, w, rr, rr1, rr2, PP1, PP2, pp1, pp2, ee2, *(y0[k] for k in ("a1", "d0", "d1", "q0", "q1")))
St = fuentes(zeta, y0["d0"], y0["q0"])
filas = dict(lig_t=lig0[0], lig_l=lig0[1], J_tt=delta[0] - St[0], J_thth=delta[1] - St[1], J_phi=delta[2] - St[2])
print("\n4. Empalmes y ligaduras en l = 0+ a orden dominante, sobre f (desplazamiento) + c (onda escalar Z^-s):")
coef = {}
for k, e in filas.items():
    n, cc = orden_dominante(e)
    cc = sp.expand(sp.simplify(cc.subs(W, s**2)))
    cf, cs = sp.factor(cc.coeff(f)), sp.factor(cc.coeff(c))
    coef[k] = (cf, cs)
    print(f"   {k:7s} (orden a^{n}):  coef f = {cf};   coef c = {cs}")
# las ligaduras y E^t_t se anulan identicamente; quedan J_thth y J_phi
tf, tc = coef["J_thth"]
pf, pc = coef["J_phi"]
det = sp.factor(sp.together(tf * pc - tc * pf))
Ps = sp.factor(sp.numer(det))
print("\n   determinante 2x2 =", det)

# ---------------------------------------------------------------- 5. limites
Ppol = sp.Poly(sp.expand(Ps), s)
lider = sp.factor(Ppol.LC())
print("\n5. Polinomio P(s) de grado", Ppol.degree())
print("   zeta = 0:", sp.factor(Ps.subs(zeta, 0)))
print("   coeficiente principal:", lider, " -> zeta_c = ", sp.solve(sp.factor(lider / zeta), zeta))
lim_bh = sp.factor(Ps.subs({F0: 0, F2: 0, u: 0}))
print("   agujero negro (F0 = F2 = 0, u -> 0; zeta != 0):", lim_bh)
print("     -> C = s^2 =", sorted({sp.simplify(x**2) for x in sp.solve(lim_bh, s)}, key=str))
rigido = {sp.simplify(x**2) for x in sp.solve(sp.Eq(tf, 0), s)}
print("   desplazamiento rigido (sin onda escalar; solo J_thth): C =", rigido)

json.dump(dict(polinomio=str(sp.expand(Ps)), variable="s", C="s^2 (rama Re s > 0)", parametros=["u", "F0", "F2", "zeta"],
               zeta_0=str(sp.factor(Ps.subs(zeta, 0))), coeficiente_principal=str(lider),
               agujero_negro=str(lim_bh),
               coef=dict((k, [str(v[0]), str(v[1])]) for k, v in coef.items())),
          open("resultados/ley_analitica.json", "w", encoding="utf-8"), indent=1)
print("\nGuardado: resultados/ley_analitica.json")
