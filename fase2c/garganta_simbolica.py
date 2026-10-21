# -*- coding: utf-8 -*-
"""
Resultados analiticos de la Fase 2c con las ecuaciones de fase2b/egb_F_gen.py (SymPy).

1. Garganta lisa (r0 = 1, r'(0) = 0; u = phi'(0), F_k = F^(k)(phi0)):
     Phi'(0) = P = (1 + u^2/4)/(4 F1 u)            (ligadura ll)
     r''(0)    = 1 - 2 F2 u^2/(1 + 4 F1^2 u^2)
     Phi''(0)  = -P^2 - (u^2/4)(1 - r''(0))
     N''/N(0)  = Phi'' + Phi'^2 = -F2 u^4 / (2 (1 + 4 F1^2 u^2))      (N = e^Phi: lapso)
2. Colapso del lapso: en un punto generico (r, r', Phi' = P, phi') el determinante de la parte principal es afin en P,
     det M = d0 + d1 P,   d1 = 128 F1^2 (1 - r'^2) [F1 phi' (1 + 3 r'^2) - r r'] / r^5,
   y el numerador de Cramer de Phi'' es -d1 P^3 + O(P^2); por tanto, si d1 != 0,
     Phi'' = -Phi'^2 + O(Phi')     (|Phi'| -> inf),
   como en relatividad general: la ecuacion de Riccati Psi' = -Psi^2 + ... lleva Psi = Phi' -> -inf a distancia finita
   ~ 1/|Psi| y N = e^Phi se anula linealmente.
Uso: python garganta_simbolica.py
"""
import os
import sys
import sympy as sp

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fase2b"))
import egb_F_gen as H      # noqa: E402

# ------------------------------------------------------------------ 1. garganta lisa
F1, F2, u = sp.symbols("F1 F2 u", real=True, nonzero=True)
P0 = (1 + u**2 / 4) / (4 * F1 * u)
lig = sp.simplify(H.est_ligadura(F1, 1, 0, P0, u))
M0 = sp.Matrix(H.est_principal(F1, 1, 0, P0, u))
b0 = sp.Matrix(H.est_resto(F1, F2, 1, 0, P0, u))
x = sp.simplify(-M0.LUsolve(b0))
rpp, Ppp, ppp = [sp.simplify(sp.factor(v)) for v in x]
print("1. Garganta lisa")
print("   ligadura ll en la garganta:", lig)
print("   det M(0) =", sp.factor(sp.simplify(M0.det())))
print("   r''(0)   : comprobacion", sp.simplify(rpp - (1 - 2 * F2 * u**2 / (1 + 4 * F1**2 * u**2))))
print("   Phi''(0) : comprobacion", sp.simplify(Ppp - (-P0**2 - u**2 / 4 * (1 - rpp))))
NppN = sp.simplify(Ppp + P0**2)
print("   N''/N(0) : comprobacion", sp.simplify(NppN + F2 * u**4 / (2 * (1 + 4 * F1**2 * u**2))))
print("   N''(0)/P^2 =", sp.factor(sp.nsimplify(sp.simplify(NppN / P0**2))))

# ------------------------------------------------------------------ 2. colapso del lapso
F1, F2, r, r1, p1, P = sp.symbols("F1 F2 r r1 p1 P", real=True)
M = sp.Matrix(H.est_principal(F1, r, r1, P, p1))
b = sp.Matrix(H.est_resto(F1, F2, r, r1, P, p1))
detM = sp.Poly(sp.expand(M.det()), P)
d1 = 128 * F1**2 * (1 - r1**2) * (F1 * p1 * (1 + 3 * r1**2) - r * r1) / r**5
M2 = M.copy()
M2[:, 1] = -b                                   # Cramer para Phi'' (segunda incognita)
num = sp.Poly(sp.expand(M2.det()), P)
print("2. Colapso del lapso")
print("   grado de det M en P:", detM.degree(), "; coeficiente de P - d1 =", sp.simplify(detM.coeff_monomial(P) - d1))
print("   grado del numerador de Phi'' en P:", num.degree(), "; coeficiente de P^3 + d1 =",
      sp.simplify(num.coeff_monomial(P**3) + d1))
print("   => Phi'' = -P^2 + O(P) si d1 != 0")
