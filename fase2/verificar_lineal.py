# -*- coding: utf-8 -*-
"""
Fase 2 - Verificacion simbolica de las ecuaciones perturbadas (gauge de longitud
propia con shift nulo: g_ll = 1, g_tl = 0).

Metrica:  ds^2 = -e^{2A(t,l)} dt^2 + dl^2 + R(t,l)^2 dOmega^2,
          A = Phi(l) + eps*a(t,l),  R = r(l) + eps*d(t,l)   (a = dPhi, d = dr)

Comprueba (G = c = 1):
  [1] R_ll exacto = -A_ll - A_l^2 - 2 R_ll/R            (ecuacion del lapso, sin d/dt)
  [2] G_tl exacto = -(2/R)(dR_l/dt - A_l dR/dt)         (ligadura de momento, E1)
  [3] G^l_l exacto = -(1-R_l^2)/R^2 + 2R_l A_l/R - e^{-2A}(2R_tt/R + R_t^2/R^2 - 2A_t R_t/R)   (E3)
  [4] G^t_t lineal no contiene derivadas temporales (ligadura hamiltoniana, E2)
  [5] E1 + matter  =>  delta m = -4 pi r^2 [p_r d + (rho+P) r_l xi]  (masa de Misner-Sharp)
Uso: python fase2/verificar_lineal.py
"""
import sys
import sympy as sp

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

t, l, th, ph, eps = sp.symbols("t ell theta phi epsilon", real=True)
Phi, r = sp.Function("Phi")(l), sp.Function("r")(l)
a, d = sp.Function("a")(t, l), sp.Function("d")(t, l)
OK = []


def check(nombre, expr):
    ok = sp.simplify(expr) == 0
    OK.append(ok)
    print(f"[{'OK   ' if ok else 'FALLO'}] {nombre}")
    if not ok:
        print("        residuo:", sp.simplify(expr))


def tensores(g, X):
    n = len(X)
    gi = sp.simplify(g.inv())
    Gam = [[[sp.simplify(sum(gi[i, s] * (sp.diff(g[s, j], X[k]) + sp.diff(g[s, k], X[j])
                                          - sp.diff(g[j, k], X[s])) for s in range(n)) / 2)
             for k in range(n)] for j in range(n)] for i in range(n)]
    Ric = sp.Matrix(n, n, lambda j, k: sp.simplify(sum(
        sp.diff(Gam[i][j][k], X[i]) - sp.diff(Gam[i][j][i], X[k])
        + sum(Gam[i][i][s] * Gam[s][j][k] - Gam[i][k][s] * Gam[s][j][i] for s in range(n))
        for i in range(n))))
    Rs = sp.simplify((gi * Ric).trace())
    G = sp.simplify(Ric - g * Rs / 2)          # G_{mu nu}
    Gmix = sp.simplify(gi * G)                  # G^mu_nu
    return Ric, G, Gmix


# Metrica general (A, R funciones de t y l) -> formas exactas
A, R = sp.Function("A")(t, l), sp.Function("R")(t, l)
X = [t, l, th, ph]
g = sp.diag(-sp.exp(2 * A), 1, R**2, R**2 * sp.sin(th)**2)
Ric, G, Gmix = tensores(g, X)
At, Al, Rt, Rl = A.diff(t), A.diff(l), R.diff(t), R.diff(l)

check("[1] R_ll = -A_ll - A_l^2 - 2R_ll/R (sin derivadas temporales)",
      Ric[1, 1] - (-A.diff(l, 2) - Al**2 - 2 * R.diff(l, 2) / R))
check("[2] G_tl = -(2/R)(R_tl - A_l R_t)", G[0, 1] - (-(2 / R) * (R.diff(t, l) - Al * Rt)))
check("[3] G^l_l completo",
      Gmix[1, 1] - (-(1 - Rl**2) / R**2 + 2 * Rl * Al / R
                    - sp.exp(-2 * A) * (2 * R.diff(t, 2) / R + Rt**2 / R**2 - 2 * At * Rt / R)))
Gtt_static = -((1 - Rl**2) / R**2 - 2 * R.diff(l, 2) / R)
check("[4] G^t_t = -[(1-R_l^2)/R^2 - 2R_ll/R] - e^{-2A} R_t^2/R^2  (sin terminos lineales en d/dt)",
      (Gmix[0, 0] - Gtt_static + sp.exp(-2 * A) * Rt**2 / R**2))

# Linealizacion explicita de G^l_l alrededor del fondo estatico
sub = {A: Phi + eps * a, R: r + eps * d}
Gll = Gmix[1, 1].subs(sub).doit()
Gll1 = sp.simplify(sp.diff(Gll, eps).subs(eps, 0))
rl, Pl = r.diff(l), Phi.diff(l)
dl, al = d.diff(l), a.diff(l)
esperado = (2 * rl * dl / r**2 + 2 * (1 - rl**2) * d / r**3 + 2 * Pl * dl / r
            + 2 * rl * al / r - 2 * rl * Pl * d / r**2 - 2 * sp.exp(-2 * Phi) * d.diff(t, 2) / r)
check("[3'] delta G^l_l lineal (E3)", Gll1 - esperado)

# [5] masa de Misner-Sharp perturbada usando E1 y la ecuacion rr de fondo
rho, P, pr, xi = sp.symbols("rho P p_r xi", real=True)
dr_l = Pl * d + 4 * sp.pi * r * (rho + P) * xi        # E1 (fluido) resuelto para d_l
dm = d / 2 * (1 - rl**2) - r * rl * dr_l               # delta[(R/2)(1-R_l^2)]
# fondo: 8 pi p_r = -(1-r_l^2)/r^2 + 2 r_l Phi_l / r  =>  r r_l Phi_l = 4 pi r^2 p_r + (1-r_l^2)/2
dm = dm.subs(Pl, (4 * sp.pi * r**2 * pr + (1 - rl**2) / 2) / (r * rl))
check("[5] delta m = -4 pi r^2 [p_r dr + (rho+P) r_l xi]",
      dm - (-4 * sp.pi * r**2 * (pr * d + (rho + P) * rl * xi)))

print(f"\n{sum(OK)}/{len(OK)} comprobaciones superadas.")
