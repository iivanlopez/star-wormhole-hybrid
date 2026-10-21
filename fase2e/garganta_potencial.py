# -*- coding: utf-8 -*-
"""
Garganta lisa de EsGB con potencial escalar (SymPy, ecuaciones de fase2b/egb_F_gen.py + potencial.py).
r0 = 1, r'(0) = 0; datos u = phi'(0), a = Phi'(0); F_k = F^(k)(phi0), V_k = V^(k)(phi0).

1. Ligadura ll:  V0/2 = 1 + u^2/4 - 4 F1 a u.   La fila thth en la garganta no contiene terminos de Gauss-Bonnet
   salvo el factor (1 - 4 F1 a u), y con la ligadura da
       N''/N(0) = Phi'' + a^2 = -(1 - 4 F1 a u)(1 + r''(0)) - u^2/2 .
   Con a = 0 (lapso sin gradiente, solo posible con V0 = 2 + u^2/2):  N''/N(0) = -(1 + r'') - u^2/2 < -1  si r'' > 0:
   el lapso tiene un maximo estricto en la garganta, para cualquier F y V.
2. Energia mecanica del escalar E = phi'^2/2 - V(phi) (el escalar estatico es una particula en el potencial -V):
       E(0) = -2 + 8 F1 a u      y      E' = -(Phi' + 2 r'/r) phi'^2 - F_l G ,  F_l = F'(phi) phi',
   con el invariante de Gauss-Bonnet como derivada total, G = -(8/(r^2 N)) [N' (1 - r'^2)]'.
3. Clase A (escalar par, u = 0, a = 0, V0 = 2):
       phi''(0) = (V1 - 8 F1)/(1 - 16 F1^2),  r''(0) = -2 F1 phi''(0),  Phi''(0) = -1 - r''(0),
       det M(0) = 2 (1 - 16 F1^2),   lambda(0) = 1/(1 - 96 F1^2 r''(0)/(1 - 16 F1^2))   (velocidad radial^2).
4. Clase B (escalar impar o "kink": F y V pares, phi0 = 0, F1 = V1 = 0, a = 0, V0 = 2 + u^2/2):
       r''(0) = -u^2 (2 F2 + 1/4),  Phi''(0) = -1 - u^2/4 + 2 F2 u^2,  phi''(0) = 0,  lambda(0) en forma cerrada.
5. Datos generales (u, a): r'', Phi'', phi'' y lambda(0) en forma cerrada (comprobados frente a potencial.segundas_V).
6. Identidades globales (X = F_l (1 - r'^2), rho = phi'^2/4 + V/2):
     fila tt:   r r'' + 2 X' = (1 - r'^2)/2 - r^2 rho/2,  integrada  r r' + 2 X = (1/2) int_0^l (1 + r'^2 - r^2 rho)
     Komar:     (r^2 N' - 2 N X - 4 r r' F_l N')' = -4 N' X - r^2 N V/2
7. Limite de Nariai (numerico): V = 2 constante y F = eps phi: N -> cos(l), horizonte en pi/2, error O(eps^2).
Uso: python garganta_potencial.py
"""
import os
import sys
import json
import numpy as np
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import potencial as P                      # noqa: E402
import egb_F_gen as H                      # noqa: E402
from egb_simbolico import tensores, gauss_bonnet   # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
out = {}

F1, F2, F3, V1, V2, u, a = sp.symbols("F1 F2 F3 V1 V2 u a", real=True)
V0 = 2 * (1 + u**2 / 4 - 4 * F1 * a * u)                       # ligadura ll


def racional(e):
    """egb_F_gen escribe (1/4) como 0.25: se pasan los coeficientes a racionales exactos."""
    return sp.nsimplify(e, rational=True)


def garganta(F1, F2, V0, V1, u, a):
    """(r'', Phi'', phi''), M, b en la garganta r0 = 1 con potencial."""
    M = sp.Matrix(H.est_principal(F1, 1, 0, a, u)).applyfunc(racional)
    b = (sp.Matrix(H.est_resto(F1, F2, 1, 0, a, u)) + sp.Matrix([V0 / 2, V0 / 2, -V1])).applyfunc(racional)
    return -M.LUsolve(b), M, b


def lambda_garganta(F1, F2, F3, x, u, a):
    """velocidad radial^2 (reloj de la garganta) en l = 0: det(M - lam N) = c0 + c1 lam."""
    r2, P2, p2 = x
    lam = sp.Symbol("lam")
    M = sp.Matrix(H.lin_principal(F1, 1, 0, a, u)).applyfunc(racional)
    N = sp.zeros(3, 3)
    for j, (d0, q0) in ((0, (1, 0)), (2, (0, 1))):
        c = [sp.diff(e, sp.Symbol("w")) for e in H.lin_resto(F1, F2, F3, sp.Symbol("w"), 1, 0, r2, a, P2, u, p2, 1,
                                                                 0, d0, 0, q0, 0)]
        N[:, j] = sp.Matrix(c).applyfunc(racional)
    pol = sp.Poly(sp.expand((M - lam * N).det()), lam)
    c = pol.all_coeffs()[::-1] + [0, 0]
    return c[0], c[1], c[2]


# ------------------------------------------------------------------ 1. thth en la garganta
x, M, b = garganta(F1, F2, V0, V1, u, a)
rpp, Ppp, ppp = [sp.simplify(v) for v in x]
NppN = sp.simplify(Ppp + a**2)
chk1 = sp.simplify(NppN + (1 - 4 * F1 * a * u) * (1 + rpp) + u**2 / 2)
print("1. Garganta lisa con potencial (r0 = 1)")
print("   ligadura ll en la garganta con V0 = 2(1 + u^2/4 - 4 F1 a u):",
      sp.simplify(H.est_ligadura(F1, 1, 0, a, u) + V0 / 2))
print("   N''/N(0) + (1 - 4 F1 a u)(1 + r'') + u^2/2 =", chk1)
out["NppN_identidad"] = str(chk1)

# ------------------------------------------------------------------ 2. energia mecanica y GB como derivada total
E0 = sp.simplify(u**2 / 2 - V0)
print("2. Energia mecanica del escalar en la garganta: E(0) =", sp.expand(E0))
t, l, th, ph = sp.symbols("t l theta varphi", real=True)
Nf, rf = sp.Function("N")(l), sp.Function("r")(l)
g = sp.diag(-Nf**2, 1, rf**2, rf**2 * sp.sin(th) ** 2)
gi, Gam, Riem, Ric, Rs = tensores(g, [t, l, th, ph])
GB = gauss_bonnet(g, gi, Riem, Ric, Rs)
total = -8 / (rf**2 * Nf) * sp.diff(Nf.diff(l) * (1 - rf.diff(l) ** 2), l)
chk2 = sp.simplify(GB - total)
print("   G + (8/(r^2 N)) [N' (1 - r'^2)]' =", chk2)
out["GB_derivada_total"] = str(chk2)

# ------------------------------------------------------------------ 3. clase A
xA, MA, bA = garganta(F1, F2, 2, V1, 0, 0)
rA, PA, pA = [sp.simplify(v) for v in xA]
detA = sp.factor(MA.det())
c0A, c1A, c2A = lambda_garganta(F1, F2, F3, xA, 0, 0)
lamA = sp.simplify(-c0A / c1A)
chk3 = [sp.simplify(pA - (V1 - 8 * F1) / (1 - 16 * F1**2)), sp.simplify(rA + 2 * F1 * pA), sp.simplify(PA + 1 + rA),
        sp.simplify(lamA - 1 / (1 - 96 * F1**2 * rA / (1 - 16 * F1**2))), sp.simplify(c2A)]
print("3. Clase A (u = 0, a = 0, V0 = 2)")
print("   r''(0) =", sp.factor(rA), "  Phi''(0) =", sp.factor(PA), "  phi''(0) =", sp.factor(pA))
print("   det M(0) =", detA, "   lambda(0) =", sp.factor(lamA))
print("   comprobaciones (deben ser 0):", chk3)
out["claseA"] = dict(rpp=str(sp.factor(rA)), Ppp=str(sp.factor(PA)), ppp=str(sp.factor(pA)), detM=str(detA),
                     lam=str(sp.factor(lamA)), chk=[str(c) for c in chk3])
# ventana hiperbolica con V1 = 0 (garganta en un maximo de V): lambda(0) > 0 y r'' > 0
lam0 = sp.factor(lamA.subs(V1, 0))
print("   con V1 = 0: r''(0) =", sp.factor(rA.subs(V1, 0)), "; lambda(0) =", lam0)
sol_borde = [s for s in sp.solve(sp.denom(sp.together(lam0)), F1) if s.is_real and s > 0]
print("   lambda(0) -> infinito (c1 = 0) en |F1| =", [float(s) for s in sol_borde])
out["claseA"]["F1_borde_V1_0"] = [float(s) for s in sol_borde]

# ------------------------------------------------------------------ 4. clase B (kink)
xB, MB, bB = garganta(0, F2, 2 + u**2 / 2, 0, u, 0)
rB, PB, pB = [sp.simplify(v) for v in xB]
c0B, c1B, c2B = lambda_garganta(0, F2, F3, xB, u, 0)
lamB = sp.factor(sp.simplify(-c0B / c1B))
print("4. Clase B (phi0 = 0, F1 = V1 = 0, a = 0, V0 = 2 + u^2/2)")
print("   r''(0) =", sp.factor(rB), "  Phi''(0) =", sp.factor(PB), "  phi''(0) =", pB)
print("   det M(0) =", sp.factor(MB.det()), "   lambda(0) =", lamB, "   c2 =", sp.simplify(c2B))
out["claseB"] = dict(rpp=str(sp.factor(rB)), Ppp=str(sp.factor(PB)), ppp=str(pB), lam=str(lamB))

# ------------------------------------------------------------------ 5. datos generales: comprobacion numerica
rng = np.random.default_rng(3)
f_x = sp.lambdify((F1, F2, V1, u, a), [rpp, Ppp, ppp])
err = 0.0
for _ in range(30):
    tipoF = rng.choice(["cuadratico", "exp", "gauss", "lineal"])
    alpha = rng.uniform(-0.5, 0.5)
    phi0, uu, aa = rng.uniform(-1, 1), rng.uniform(-2, 2), rng.uniform(-1, 1)
    Fd = P.acoplamiento(tipoF, alpha)
    Ud = P.forma("cuartico", g=rng.uniform(0, 2))
    amp = P.amp_garganta(Ud, Fd, phi0, uu, aa)
    Vd = (lambda x, Ud=Ud, amp=amp: tuple(amp * c for c in Ud(x)))
    y0 = np.array([1.0, 0.0, 0.0, aa, phi0, uu])
    xs, _ = P.segundas_V(y0, Fd, Vd)
    Fk = Fd(phi0)
    ref = np.array(f_x(Fk[1], Fk[2], Vd(phi0)[1], uu, aa), dtype=float)
    err = max(err, float(np.max(np.abs(xs - ref) / (1 + np.abs(ref)))), abs(P.est_ligadura_V(y0, Fd, Vd)))
print(f"5. Formulas generales frente a potencial.segundas_V (30 gargantas aleatorias): error max {err:.1e}")
out["comprobacion_numerica"] = err
c0g, c1g, c2g = lambda_garganta(F1, F2, F3, x, u, a)
print("   coeficiente c2 del polinomio de dispersion en la garganta (debe ser 0):", sp.simplify(c2g))

# ------------------------------------------------------------------ 6. identidades globales (fila tt y Komar)
l_ = sp.Symbol("l")
rf, Pf, pf = (sp.Function(n)(l_) for n in ("r", "Phi", "p"))
Ff, Vf = sp.Function("F"), sp.Function("V")
F1f, F2f = Ff(pf).diff(pf), Ff(pf).diff(pf, 2)
d = [f.diff(l_) for f in (rf, Pf, pf)]
dd = [f.diff(l_, 2) for f in (rf, Pf, pf)]
Mf = sp.Matrix(H.est_principal(F1f, rf, d[0], d[1], d[2])).applyfunc(racional)
bf = sp.Matrix(H.est_resto(F1f, F2f, rf, d[0], d[1], d[2])).applyfunc(racional) \
    + sp.Matrix([Vf(pf) / 2, Vf(pf) / 2, -Vf(pf).diff(pf)])
Ef = Mf * sp.Matrix(dd) + bf                                     # (E^t_t, E^th_th, E_phi)
llf = racional(H.est_ligadura(F1f, rf, d[0], d[1], d[2])) + Vf(pf) / 2
Nf_ = sp.exp(Pf)
Fl_ = Ff(pf).diff(l_)
Xf = Fl_ * (1 - d[0] ** 2)
rho_f = d[2] ** 2 / 4 + Vf(pf) / 2
# (a) fila tt integrada: r r'' + 2 X' = (1 - r'^2)/2 - r^2 rho/2  <=>  (r^2/2) E^t_t = 0
chk_tt = sp.simplify(rf**2 * Ef[0] / 2 - (rf * dd[0] + 2 * Xf.diff(l_) - (1 - d[0] ** 2) / 2 + rf**2 * rho_f / 2))
# (b) Komar: -(r^2 N/2)(E^t_t - E^l_l - 2 E^th_th) = (r^2 N')' - (2 N X + 4 r r' F_l N')' + 4 N' X + r^2 N V/2
Ckom = Ef[0] - llf - 2 * Ef[1]
kom = (rf**2 * Nf_.diff(l_)).diff(l_) - (2 * Nf_ * Xf + 4 * rf * d[0] * Fl_ * Nf_.diff(l_)).diff(l_) \
    + 4 * Nf_.diff(l_) * Xf + rf**2 * Nf_ * Vf(pf) / 2
chk_kom = sp.simplify(-rf**2 * Nf_ * Ckom / 2 - kom)
print("6. fila tt: (r^2/2) E^t_t - [r r'' + 2 X' - (1 - r'^2)/2 + r^2 rho/2] =", chk_tt)
print("   Komar: -(r^2 N/2)(E^t_t - E^l_l - 2 E^th_th) - [(r^2 N')' - (2NX + 4 r r' F_l N')' + 4 N' X + r^2 N V/2] =",
      chk_kom)
out["identidades_globales"] = dict(tt=str(chk_tt), komar=str(chk_kom))

# ------------------------------------------------------------------ 7. limite de Nariai (numerico)
# V = 2 constante (Lambda = 1) y F = eps phi: la garganta de clase A tiende a dS2 x S2, N = cos(l), horizonte en pi/2
nariai = []
for eps in (1e-2, 1e-3, 1e-4):
    Ld = P.LadoV(P.acoplamiento("lineal", eps), P.potencial("constante", 2.0), 0.0, 0.0, 0.0, phi_v=-1.0, l_max=5.0,
                 rtol=1e-12).resolver()
    ls_ = np.linspace(0.0, 0.95 * Ld.l_end, 200)
    errN = max(abs(np.exp(Ld.sol.sol(x)[2]) - np.cos(x)) for x in ls_)
    nariai.append(dict(eps=eps, rpp0=Ld.rpp0, fin=Ld.fin, l_s=Ld.l_end, err_N=float(errN)))
    print(f"7. Nariai, F = {eps:g} phi, V = 2: r''(0) = {Ld.rpp0:.2e}, final '{Ld.fin}' en l = {Ld.l_end:.6f} "
          f"(pi/2 = {np.pi / 2:.6f}); max |N - cos l| = {errN:.1e}")
out["nariai"] = nariai
json.dump(out, open(os.path.join(AQUI, "resultados", "garganta_potencial.json"), "w", encoding="utf-8"), indent=1,
          ensure_ascii=False)
