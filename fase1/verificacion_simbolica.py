# -*- coding: utf-8 -*-
"""
Verificacion simbolica y numerica de fase1_marco_teorico.md

Uso:  python fase1/verificacion_simbolica.py
Requiere: sympy (>=1.12)
"""
import sys
import sympy as sp

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

t, r, th, ph, l, a = sp.symbols("t r theta phi ell a", real=True, positive=True)
RESULTS = []


def check(nombre, expr):
    ok = sp.simplify(expr) == 0
    RESULTS.append(ok)
    print(f"[{'OK   ' if ok else 'FALLO'}] {nombre}")
    if not ok:
        print("        residuo:", sp.simplify(expr))


def einstein_mixed(g, X):
    """G^a_b para una metrica g en coordenadas X (convencion MTW)."""
    n = len(X)
    gi = g.inv()
    Gam = [[[sum(gi[i, d] * (sp.diff(g[d, j], X[k]) + sp.diff(g[d, k], X[j]) - sp.diff(g[j, k], X[d]))
                 for d in range(n)) / 2 for k in range(n)] for j in range(n)] for i in range(n)]
    Ric = sp.Matrix(n, n, lambda j, k: sum(
        sp.diff(Gam[i][j][k], X[i]) - sp.diff(Gam[i][j][i], X[k])
        + sum(Gam[i][i][d] * Gam[d][j][k] - Gam[i][k][d] * Gam[d][j][i] for d in range(n))
        for i in range(n)))
    Rmix = sp.simplify(gi * Ric)
    return sp.simplify(Rmix - sp.eye(n) * Rmix.trace() / 2), Gam


# ---------------------------------------------------------------- 1. Morris-Thorne en r
Phi, b = sp.Function("Phi")(r), sp.Function("b")(r)
g = sp.diag(-sp.exp(2 * Phi), 1 / (1 - b / r), r**2, r**2 * sp.sin(th)**2)
G, _ = einstein_mixed(g, [t, r, th, ph])
rho, pr, pt = -G[0, 0] / (8 * sp.pi), G[1, 1] / (8 * sp.pi), G[2, 2] / (8 * sp.pi)
bp, Pp = b.diff(r), Phi.diff(r)
check("(2.1) G_tt = b'/r^2", 8 * sp.pi * rho - bp / r**2)
check("(2.1) G_rr", 8 * sp.pi * pr - (-b / r**3 + 2 * (1 - b / r) * Pp / r))
check("(2.1) G_thth", 8 * sp.pi * pt - (1 - b / r) * (Phi.diff(r, 2) + Pp**2 - (bp * r - b) / (2 * r * (r - b)) * Pp
                                                     - (bp * r - b) / (2 * r**2 * (r - b)) + Pp / r))
check("(2.2) conservacion (Bianchi)", pr.diff(r) + (rho + pr) * Pp - 2 / r * (pt - pr))
m = sp.Function("m")(r)
check("(2.3) Phi' = (m+4 pi r^3 p_r)/(r(r-2m))",
      (Pp - (m + 4 * sp.pi * r**3 * pr) / (r * (r - 2 * m))).subs(b, 2 * m).doit())
bp0 = sp.symbols("bp0", real=True)  # b'(r0) libre; b(r0)=r0 anula (1-b/r)Phi' si Phi' es finito
check("(2.5) rho+p_r en garganta = (b'(r0)-1)/(8 pi r0^2)",
      (rho + pr).subs(bp, bp0).subs(b, r) - (bp0 - 1) / (8 * sp.pi * r**2))

# ---------------------------------------------------------------- 2. gauge de longitud propia
R, Pl = sp.Function("R")(l), sp.Function("Phi")(l)
gl = sp.diag(-sp.exp(2 * Pl), 1, R**2, R**2 * sp.sin(th)**2)
Gl, Gaml = einstein_mixed(gl, [t, l, th, ph])
rl, rll, pl1, pll = R.diff(l), R.diff(l, 2), Pl.diff(l), Pl.diff(l, 2)
rho_l, pr_l, pt_l = -Gl[0, 0] / (8 * sp.pi), Gl[1, 1] / (8 * sp.pi), Gl[2, 2] / (8 * sp.pi)
check("(2.9) 8 pi rho", 8 * sp.pi * rho_l - ((1 - rl**2) / R**2 - 2 * rll / R))
check("(2.9) 8 pi p_r", 8 * sp.pi * pr_l - (-(1 - rl**2) / R**2 + 2 * rl * pl1 / R))
check("(2.9) 8 pi p_t", 8 * sp.pi * pt_l - (rll / R + pll + pl1**2 + rl * pl1 / R))

# Resultado 4: en la garganta Z2 (r_l = Phi_l = 0)
Rll, Pll, r0 = sp.symbols("R_ll Phi_ll r_0", real=True)
garg = lambda e: e.subs({R.diff(l, 2): Rll, Pl.diff(l, 2): Pll}).subs({R.diff(l): 0, Pl.diff(l): 0}).subs(R, r0)
check("Ligadura en garganta: 8 pi p_r = -1/r0^2", garg(8 * sp.pi * pr_l) + 1 / r0**2)
check("Resultado 4: Phi_ll0 = 4 pi (rho + p_r + 2 p_t)_0",
      Pll - 4 * sp.pi * garg(rho_l + pr_l + 2 * pt_l))

# (2.11) fantasma sin masa + plasma
rs0, P0 = sp.symbols("rho_s0 P_0", real=True)
X0 = P0 + 1 / (8 * sp.pi * r0**2)                      # de la ligadura 8 pi (P0 - X0) = -1/r0^2
rll0 = 1 / (2 * r0) - 4 * sp.pi * r0 * (rs0 - X0)      # (2.10) con r_l = 0
pll0 = 8 * sp.pi * (P0 + X0) - rll0 / r0               # (2.10) con Phi_l = 0
check("(2.11) r_ll0 = 1/r0 - 4 pi r0 (rho_s0 - P0)", rll0 - (1 / r0 - 4 * sp.pi * r0 * (rs0 - P0)))
check("(2.11) Phi_ll0 = 4 pi (rho_s0 + 3 P0)", pll0 - 4 * sp.pi * (rs0 + 3 * P0))
D = sp.symbols("D")
check("(2.11) Box phi = 0 => phi_l = D/(r^2 e^Phi)",
      sp.diff(sp.exp(Pl) * R**2 * (D / (R**2 * sp.exp(Pl))), l))

# (4.3) k^mu k^nu nabla_mu nabla_nu F en gauge l (base ortonormal)
F = sp.Function("F")(l)
X = [t, l, th, ph]
hess = lambda i, j: sp.diff(F, X[i], X[j]) - sum(Gaml[k][i][j] * sp.diff(F, X[k]) for k in range(4))
check("(4.3) kk nabla nabla F = F_ll - Phi_l F_l",
      hess(0, 0) / sp.exp(2 * Pl) + hess(1, 1) - (F.diff(l, 2) - pl1 * F.diff(l)))

# ---------------------------------------------------------------- 3. capa delgada (3.2)-(3.4)
ba, Pa = sp.Function("b")(a), sp.Function("Phi")(a)
s = sp.sqrt(1 - ba / a)
sigma = -s / (2 * sp.pi * a)
Pst = s * (Pa.diff(a) + 1 / a) / (4 * sp.pi)
rho_pr = (ba.diff(a) / a**2 - ba / a**3 + 2 * (1 - ba / a) * Pa.diff(a) / a) / (8 * sp.pi)
Xi = 2 / s * rho_pr                                    # (3.4) con a_punto = 0
check("(3.4) d sigma/da + 2(sigma+P)/a = Xi  (limite estatico)",
      sp.simplify((sigma.diff(a) + 2 * (sigma + Pst) / a - Xi) * s))
M = sp.symbols("M", positive=True)
Psch = Pst.subs({ba: 2 * M}).subs(Pa, sp.log(1 - 2 * M / a) / 2).doit()
check("(3.2) limite Schwarzschild (Visser 1989)", Psch - (1 - M / a) / (4 * sp.pi * a * sp.sqrt(1 - 2 * M / a)))

# ---------------------------------------------------------------- 4. escalas numericas (SI)
Gn, c, GMsun = 6.67430e-11, 299792458.0, 1.32712440018e20
Msun, Rsun, Mearth, AU = GMsun / 6.67430e-11, 6.957e8, 5.9722e24, 1.495978707e11
kB, mu_u, mu, Tc, Pc = 1.380649e-23, 1.66053906660e-27, 0.85, 1.57e7, 2.5e16
Rs = 2 * GMsun / c**2
M0_1m = c**2 / (2 * Gn)
tau = lambda r0_: c**4 / (8 * sp.pi.evalf() * Gn * r0_**2)
print("\nEscalas (SI):")
print(f"  R_s(Sol) = 2GM/c^2            = {Rs:.4e} m")
print(f"  compacidad 2GM/(R c^2)        = {Rs / Rsun:.3e}")
print(f"  M0 por metro de r0            = {M0_1m:.3e} kg = {M0_1m / Msun:.3e} Msun = {M0_1m / Mearth:.1f} Mtierra")
print(f"  tau0(1 m), tau0(1 km)         = {float(tau(1.0)):.3e} Pa, {float(tau(1e3)):.3e} Pa")
print(f"  tau0(1 km)/P_c(Sol)           = {float(tau(1e3)) / Pc:.2e}")
print(f"  r0 con tau0 = P_c             = {float(sp.sqrt(c**4 / (8 * sp.pi * Gn * Pc))):.3e} m = {float(sp.sqrt(c**4 / (8 * sp.pi * Gn * Pc))) / AU:.0f} UA")
print(f"  rho_max flare-out (1 km)      = {c**2 / (8 * float(sp.pi) * Gn * 1e6):.3e} kg/m^3")
print(f"  r_B/(eps r0) = mu m_u c^2/2kT = {mu * mu_u * c**2 / (2 * kB * Tc):.3e}")
print(f"  f ~ c/(2 pi r0) para 1 km     = {c / (2 * float(sp.pi) * 1e3):.3e} Hz")
print(f"  Omega_sol r0/c para 1 km      = {2.865e-6 * 1e3 / c:.1e}")
print(f"\n{sum(RESULTS)}/{len(RESULTS)} comprobaciones simbolicas superadas.")
