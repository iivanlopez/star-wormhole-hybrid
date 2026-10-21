# -*- coding: utf-8 -*-
"""
Coeficientes del polinomio de dispersion de las perturbaciones radiales en un fondo ULTRAESTATICO (Phi = 0) de EsGB con
un campo electromagnetico (Maxwell o NLED) radial, a partir de las ecuaciones generadas de fase2b (egb_F_gen.py):

  det(M - lambda N) = c0 + c1 lambda,   lambda = -c0/c1 = G/K,
  G = (2/r^5)^{-1} c0,   K = -(2/r^5)^{-1} c1        (coeficientes de gradiente y cinetico).

Resultado (comprobado aqui en SymPy), con B = 1 - r'^2, A = r - 4 F' phi' r', C = 1 - 4 (F' phi')' :
  G = A^2 r^2 - 16 F'^2 B^2 = (A r - 4 F' B)(A r + 4 F' B)
  K = A^2 r^2 - 16 F'^2 B (4 A r'' + B C)
Reduccion con las ecuaciones de fondo ultraestaticas (T^t_t = T^r_r, 8 pi T = diag(-eps, -eps, P, P)):
  escalar:  (r^2 phi')' = 0                    -> phi' = D/r^2
  ll:       1 - r'^2 = r^2 eps - D^2/(4 r^2)
  thth:     r''/r = P - D^2/(4 r^4)
  tt:       (B h)' = -r^2 P/2,  h = dF/dl = F' phi'   (unica ecuacion con F)
y con ellas  F' = h r^2/D  y  B C = B - 4 B h' = B - 4 (B h)' + 4 B' h = B + 2 r^2 P - 8 r' r'' h.
Uso: python cinetico_simbolico.py
"""
import os
import sys
import sympy as sp

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fase2b"))
import egb_F_gen as H      # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

F1, F2, F3, r, r1, r2, p1, p2, lam = sp.symbols("F1 F2 F3 r r1 r2 p1 p2 lambda", real=True)
# fondo ultraestatico: Phi' = Phi'' = 0, e^{-2Phi} = 1
P1 = P2 = 0
e2 = 1
M = sp.Matrix(H.lin_principal(F1, r, r1, P1, p1))
N = sp.zeros(3, 3)
cero = H.lin_resto(F1, F2, F3, 0, r, r1, r2, P1, P2, p1, p2, e2, 0, 0, 0, 0, 0)
for j, (d0, q0) in ((0, (1, 0)), (2, (0, 1))):
    con = sp.Matrix(H.lin_resto(F1, F2, F3, 1, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0))
    sin = sp.Matrix(H.lin_resto(F1, F2, F3, 0, r, r1, r2, P1, P2, p1, p2, e2, 0, d0, 0, q0, 0))
    N[:, j] = sp.expand(con - sin)
print("N (coeficientes de w) =", N)
Pl = sp.expand((M - lam * N).det())
c0 = sp.factor(Pl.coeff(lam, 0))
c1 = sp.factor(Pl.coeff(lam, 1))
c2 = sp.simplify(Pl.coeff(lam, 2))
print("c2 =", c2)
A = r - 4 * F1 * p1 * r1
B = 1 - r1**2
C = 1 - 4 * (F1 * p2 + F2 * p1**2)
G = A**2 * r**2 - 16 * F1**2 * B**2
K = A**2 * r**2 - 16 * F1**2 * B * (4 * A * r2 + B * C)
print("c0 - (2/r^5) G =", sp.simplify(c0 - 2 * G / r**5))
print("c1 + (2/r^5) K =", sp.simplify(c1 + 2 * K / r**5))
print("G factoriza:", sp.factor(G - (A * r - 4 * F1 * B) * (A * r + 4 * F1 * B)))

# ------------------------------------------------------------------ ecuaciones de fondo ultraestaticas
eps, P, D, h, hp = sp.symbols("epsilon P D h hp", real=True)
Fd1 = sp.Symbol("F1")
Mest = sp.Matrix(H.est_principal(F1, r, r1, 0, p1))
best = sp.Matrix(H.est_resto(F1, F2, r, r1, 0, p1))
# filas tt, thth, phi con (r'', 0, phi'') y Maxwell/NLED: b_tt + eps, b_thth - P
x = sp.Matrix([r2, 0, p2])
ec = Mest * x + best + sp.Matrix([eps, -P, 0])
lig = H.est_ligadura(F1, r, r1, 0, p1) + eps
print("\nfila tt   :", sp.factor(sp.expand(ec[0])))
print("fila thth :", sp.factor(sp.expand(ec[1])))
print("fila phi  :", sp.factor(sp.expand(ec[2])))
print("ligadura ll:", sp.factor(sp.expand(lig)))
# comprobacion de la forma (B h)' = -r^2 P/2 de la fila tt: con h = F1 p1, h' = F2 p1^2 + F1 p2, phi' = D/r^2,
# eps de la ligadura y P de la fila thth
sus = {p1: D / r**2, p2: -2 * D * r1 / r**3}
eps_l = sp.solve(lig.subs(sus), eps)[0]
P_t = sp.solve(ec[1].subs(sus), P)[0]
tt = sp.expand(ec[0].subs(sus).subs(eps, eps_l))
hval = F1 * D / r**2
hder = (F2 * D / r**2) * (D / r**2) + F1 * (-2 * D * r1 / r**3)
Bv = 1 - r1**2
derivada_Bh = sp.expand(-2 * r1 * r2 * hval + Bv * hder)       # (B h)'
print("fila tt - (2/r^2)[(Bh)' + r^2 P/2]* :", sp.simplify(tt - sp.expand(4 * (derivada_Bh + r**2 * P_t / 2) / r**2)))
