# -*- coding: utf-8 -*-
"""
El coeficiente cinetico de las perturbaciones radiales en EsGB + potencial, con lapso variable (SymPy).
Generaliza la Parte III de la Fase 2d (fase2d/teorema_chi.py), que suponia N = 1.

1. Polinomio de dispersion de las ecuaciones generadas (fase2b/egb_F_gen.py), en cualquier punto:
       det(M - lam N) = (2/r^5) (G - lam e^{-2 Phi} K),     lam = N^2 G/K,
       K = A^2 r^2 - 16 F'^2 B (4 A r'' + B C),               (la misma forma que en la Fase 2d)
       G = A^2 r^2 - 16 F'^2 B [4 A Phi' r' + B (1 - 4 F_l Phi')],
   con A = r - 4 F_l r', B = 1 - r'^2, C = 1 - 4 F_ll, F_l = F'(phi) phi', F_ll = F'' phi'^2 + F' phi''.
2. La fila tt del fondo (con V) es exactamente  A r'' + 2 B F_ll = (B - r^2 rho)/2,  rho = phi'^2/4 + V/2.
   Eliminando F_ll:  K = A^2 r^2 - 16 F'^2 B (6 A r'' + r^2 rho).
3. chi = r r' + 2 F_l (1 - 3 r'^2) es r^2/2 por el coeficiente de Phi' en la ligadura ll:
       2 chi Phi' = r^2 num,     num = (1 - r'^2)/r^2 + phi'^2/4 - V/2.
   Identidad:  K = -a chi' + b chi + c num,   a = 48 F'^2 A B^2 / (r' (r r' - 4 F_l)),  b y c explicitos.
4. En un cero regular de chi (Phi' finito, luego num = 0):
       K = -(16 F'^2 (1 - r'^2)^2 / r'^2) chi'      (el lapso desaparece).
5. Garganta lisa (r0 = 1): chi(0) = 2 F1 phi1,  chi'(0) = 2 F1 a phi1 - phi1^2/4.
Uso: python teorema_chi_lapso.py
"""
import os
import sys
import json
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import potencial as P               # noqa: E402,F401  (rutas a fase2b)
import egb_F_gen as H               # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
out = {}


def racional(e):
    return sp.nsimplify(e, rational=True)


F1, F2, F3, w = sp.symbols("F1 F2 F3 w", real=True)
r, r1, r2, P1, P2, p1, p2, e2, lam = sp.symbols("r r1 r2 P1 P2 p1 p2 e2 lam", real=True)

# ------------------------------------------------------------------ 1. dispersion en un punto generico
M = sp.Matrix(H.lin_principal(F1, r, r1, P1, p1)).applyfunc(racional)
Nm = sp.zeros(3, 3)
for j, (d0, q0) in ((0, (1, 0)), (2, (0, 1))):
    c = [sp.diff(e, w) for e in H.lin_resto(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0)]
    Nm[:, j] = sp.Matrix(c).applyfunc(racional)
pol = sp.Poly(sp.expand((M - lam * Nm).det()), lam)
c0, c1 = pol.coeff_monomial(1), pol.coeff_monomial(lam)
Fl, Fll = F1 * p1, F1 * p2 + F2 * p1**2
A, B, C = r - 4 * Fl * r1, 1 - r1**2, 1 - 4 * Fll
K = A**2 * r**2 - 16 * F1**2 * B * (4 * A * r2 + B * C)
G = A**2 * r**2 - 16 * F1**2 * B * (4 * A * P1 * r1 + B * (1 - 4 * Fl * P1))
chk1 = [sp.simplify(pol.degree() - 1), sp.simplify(c0 - 2 * G / r**5), sp.simplify(c1 + 2 * e2 * K / r**5)]
print("1. det(M - lam N) = (2/r^5)(G - lam e2 K):  grado en lam - 1, c0 - 2G/r^5, c1 + 2 e2 K/r^5 =", chk1)
out["dispersion"] = [str(x) for x in chk1]

# ------------------------------------------------------------------ 2. fila tt con potencial
V0, V1 = sp.symbols("V0 V1", real=True)
Mest = sp.Matrix(H.est_principal(F1, r, r1, P1, p1)).applyfunc(racional)
best = sp.Matrix(H.est_resto(F1, F2, r, r1, P1, p1)).applyfunc(racional) + sp.Matrix([V0 / 2, V0 / 2, -V1])
tt = sp.expand((Mest * sp.Matrix([r2, P2, p2]) + best)[0])
rho = p1**2 / 4 + V0 / 2
chk2 = sp.simplify(tt * r**2 / 2 - (A * r2 + 2 * B * Fll - (B - r**2 * rho) / 2))
print("2. (r^2/2) E^t_t - [A r'' + 2 B F_ll - (B - r^2 rho)/2] =", chk2)
p2_tt = sp.solve(sp.Eq(A * r2 + 2 * B * Fll, (B - r**2 * rho) / 2), p2)[0]
K_tt = sp.simplify(K.subs(p2, p2_tt))
chk2b = sp.simplify(K_tt - (A**2 * r**2 - 16 * F1**2 * B * (6 * A * r2 + r**2 * rho)))
print("   K con F_ll de la fila tt - [A^2 r^2 - 16 F'^2 B (6 A r'' + r^2 rho)] =", chk2b)
out["tt"] = [str(chk2), str(chk2b)]

# ------------------------------------------------------------------ 3. identidad K = -a chi' + b chi + c num
# variables independientes del punto: r, r', r'', F', phi', V (phi'' sale de la fila tt); se usa X = F_l como
# variable en lugar de phi' (F' = X/phi')
X, u = sp.symbols("X u", real=True, nonzero=True)          # X = F_l, u = phi'
f = X / u
Ab, Bb = r - 4 * X * r1, 1 - r1**2
rhob = u**2 / 4 + V0 / 2
Fll_b = ((Bb - r**2 * rhob) / 2 - Ab * r2) / (2 * Bb)        # F_ll de la fila tt
Kb = Ab**2 * r**2 - 16 * f**2 * Bb * (6 * Ab * r2 + r**2 * rhob)
chi = r * r1 + 2 * X * (1 - 3 * r1**2)
dchi = r1**2 + r * r2 + 2 * Fll_b * (1 - 3 * r1**2) - 12 * X * r1 * r2     # chi' (F_l' = F_ll)
num = Bb / r**2 + u**2 / 4 - V0 / 2
a = 48 * f**2 * Ab * Bb**2 / (r1 * (r * r1 - 4 * X))
resto = sp.factor(sp.together(sp.expand(Kb + a * dchi)))
assert sp.simplify(sp.diff(resto, r2)) == 0, "a no elimina r''"
resto = sp.factor(sp.together(resto.subs(r2, 0)))
# resto es lineal en V0 a traves de num: resto = b chi + c num
cc = sp.simplify(sp.diff(resto, V0) / sp.diff(num, V0))
bb = sp.simplify((resto - cc * num) / chi)
chk3 = sp.simplify(Kb - (-a * dchi + bb * chi + cc * num))
print("3. K = -a chi' + b chi + c num:  residuo", chk3, ";  b contiene V0:", bb.has(V0))
print("   c =", sp.factor(cc))
print("   b =", sp.factor(bb))
out["identidad"] = dict(residuo=str(chk3), b=str(sp.factor(bb)), c=str(sp.factor(cc)))

# ------------------------------------------------------------------ 4. en un cero regular de chi
X0 = sp.solve(sp.Eq(chi, 0), X)[0]                           # F_l en el cero de chi
V00 = sp.solve(sp.Eq(num, 0), V0)[0]                         # num = 0
a0 = sp.simplify(a.subs(X, X0))
chk4 = sp.simplify(a0 - 16 * (X0 / u)**2 * Bb**2 / r1**2)
K0 = sp.simplify((Kb + 16 * f**2 * Bb**2 * dchi / r1**2).subs(V0, V00).subs(X, X0))
print("4. en chi = 0: a - 16 F'^2 B^2/r'^2 =", chk4, ";  K + 16 F'^2 B^2 chi'/r'^2 (con num = 0) =", K0)
out["cero_chi"] = [str(chk4), str(K0)]

# ------------------------------------------------------------------ 5. garganta lisa
aa, ph1 = sp.symbols("a phi1", real=True)
chi_g = (r * r1 + 2 * F1 * p1 * (1 - 3 * r1**2)).subs({r: 1, r1: 0, p1: ph1})
# chi'(0) con la fila tt y la ligadura V0/2 = 1 + phi1^2/4 - 4 F1 a phi1
V0g = 2 * (1 + ph1**2 / 4 - 4 * F1 * aa * ph1)
dchi_g = sp.simplify(dchi.subs({X: F1 * ph1, u: ph1, V0: V0g}).subs({r: 1, r1: 0}))
print("5. garganta: chi(0) =", chi_g, ";  chi'(0) =", sp.expand(dchi_g))
out["garganta"] = dict(chi0=str(chi_g), dchi0=str(sp.expand(dchi_g)))

# ------------------------------------------------------------------ 6. phi'(0) = 0 con a != 0: chi''(0)
# desarrollo de Taylor de tercer orden en l (V0 = 2 por la ligadura); chi = chi''(0) l^2/2 + O(l^3)
l = sp.Symbol("l")
F1s, F2s, F3s, F4s, V1s, V2s, V3s = sp.symbols("F1 F2 F3 F4 V1 V2 V3", real=True)
r2s, r3s, r4s, P2s, P3s, P4s, q2s, q3s, q4s = sp.symbols("r2 r3 r4 P2 P3 P4 q2 q3 q4", real=True)
rr = 1 + r2s * l**2 / 2 + r3s * l**3 / 6 + r4s * l**4 / 24
PP = aa * l + P2s * l**2 / 2 + P3s * l**3 / 6 + P4s * l**4 / 24
dd = q2s * l**2 / 2 + q3s * l**3 / 6 + q4s * l**4 / 24                  # phi - phi0
Fp_ = F1s + F2s * dd + F3s * dd**2 / 2 + F4s * dd**3 / 6
Fpp_ = F2s + F3s * dd + F4s * dd**2 / 2
V_ = 2 + V1s * dd + V2s * dd**2 / 2 + V3s * dd**3 / 6
Vp_ = V1s + V2s * dd + V3s * dd**2 / 2
Mt = sp.Matrix(H.est_principal(Fp_, rr, rr.diff(l), PP.diff(l), dd.diff(l))).applyfunc(racional)
bt = sp.Matrix(H.est_resto(Fp_, Fpp_, rr, rr.diff(l), PP.diff(l), dd.diff(l))).applyfunc(racional) \
    + sp.Matrix([V_ / 2, V_ / 2, -Vp_])
eqs = Mt * sp.Matrix([rr.diff(l, 2), PP.diff(l, 2), dd.diff(l, 2)]) + bt
ser = [sp.expand(sp.series(e, l, 0, 3).removeO()) for e in eqs]
sol6 = {}
for k, inc in enumerate(([r2s, P2s, q2s], [r3s, P3s, q3s], [r4s, P4s, q4s])):
    s = sp.solve([sp.expand(e.coeff(l, k).subs(sol6)) for e in ser], inc, dict=True)[0]
    sol6.update({x: sp.simplify(v) for x, v in s.items()})
chi6 = rr * rr.diff(l) + 2 * Fp_ * dd.diff(l) * (1 - 3 * rr.diff(l) ** 2)
c6 = sp.expand(sp.series(chi6, l, 0, 4).removeO())
chk6 = [sp.simplify(c6.coeff(l, k).subs(sol6)) for k in (0, 1, 2)]
chi3 = sp.factor(sp.simplify(c6.coeff(l, 3).subs(sol6)))
print("6. phi'(0) = 0, a libre: coeficientes de l^0, l^1, l^2 de chi =", chk6)
print("   coeficiente de l^3, c3(a) =", chi3)
print("   c3(a) - c3(0) =", sp.factor(sp.simplify(chi3 - chi3.subs(aa, 0))))
out["phi1_cero"] = dict(chi012=[str(c) for c in chk6], c3=str(chi3))

# ------------------------------------------------------------------ 7. materia adicional sin grados de libertad radiales
# (Maxwell, NLED, ...: no entra en M ni en N).  Con densidad rho_m y presion radial p_m, la fila tt y la ligadura son
#   A r'' + 2 B F_ll = (B - r^2 (rho + rho_m))/2,     2 chi Phi' = r^2 (num + p_m)
# y en un cero regular de chi:  K = -16 F'^2 B^2 chi'/r'^2 - 8 F'^2 B^2 r^2 (rho_m + p_m)/r'^2.
rm, pm = sp.symbols("rho_m p_m", real=True)
rho_t = rhob + rm
Fll_m = ((Bb - r**2 * rho_t) / 2 - Ab * r2) / (2 * Bb)
Km = Ab**2 * r**2 - 16 * f**2 * Bb * (6 * Ab * r2 + r**2 * rho_t)
dchim = r1**2 + r * r2 + 2 * Fll_m * (1 - 3 * r1**2) - 12 * X * r1 * r2
V0m = sp.solve(sp.Eq(num + pm, 0), V0)[0]
K7 = sp.simplify((Km + 16 * f**2 * Bb**2 * dchim / r1**2 + 8 * f**2 * Bb**2 * r**2 * (rm + pm) / r1**2)
                 .subs(V0, V0m).subs(X, X0))
print("7. materia adicional: K + 16 F'^2 B^2 chi'/r'^2 + 8 F'^2 B^2 r^2 (rho_m + p_m)/r'^2 en chi = 0 =", K7)
out["materia_adicional"] = str(K7)
json.dump(out, open(os.path.join(AQUI, "resultados", "teorema_chi_lapso.json"), "w", encoding="utf-8"), indent=1,
          ensure_ascii=False)
