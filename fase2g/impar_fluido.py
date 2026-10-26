# -*- coding: utf-8 -*-
"""
Fase 2g: sector impar (axial) con el fluido P(Y).  El fluido irrotacional no tiene grados de libertad impares; su
accion cuadratica en (h0, h1) (gauge de Regge-Wheeler, g_{t varphi} = h0 W, g_{l varphi} = h1 W, W = sin th dP_l/dth)
es algebraica.  Se calcula aqui (SymPy, exacta) y se comprueba que
    L2^fluido = N_l [ -rho~ h0^2/N - N p~ h1^2 ],      N_l = l(l+1)/(2l+1),
que cancela exactamente los terminos 2 E_A, -2 E_B de la parte gravitatoria (T4 de la Fase 2f, con
E_A^g = r^2 rho~/(2N) y E_B^g = -N r^2 p~/2 en la capa de masas con fluido).  El sector impar con fluido es, por tanto,
    alfa = N_l Gcal/N,   beta = N_l lam_l G4/(N r^2),   gamma = -N_l lam_l N Gcal/r^2      (lam_l = (l-1)(l+2)),
igual que en el vacio y sin funciones libres: K_o = N_l r^2 Gcal/(lam_l N^3), P_o = N_l r^2 Gcal^2/(lam_l N G4),
c^2 = Gcal/G4.  Uso: python impar_fluido.py
"""
import sympy as sp

t, l, th = sp.symbols("t l theta", real=True)
eps = sp.Symbol("eps")
N, r, h0, h1 = sp.symbols("N r h0 h1", positive=True)
P0, P1, P2 = sp.symbols("P0 P1 P2", real=True)              # P, P_Y, P_YY en Y_b = 1/(2 N^2)
for LL in (2, 3):
    W = sp.simplify(sp.sin(th) * sp.diff(sp.legendre(LL, sp.cos(th)), th))
    g = sp.diag(-N ** 2, 1, r ** 2, r ** 2 * sp.sin(th) ** 2)
    g[0, 3] = g[3, 0] = eps * h0 * W
    g[1, 3] = g[3, 1] = eps * h1 * W
    gi = g.inv()
    Y = -gi[0, 0] / 2                                   # chi = t
    Yb = 1 / (2 * N ** 2)
    P = P0 + P1 * (Y - Yb) + P2 * (Y - Yb) ** 2 / 2
    dens = sp.sqrt(-g.det()) * P
    d2 = sp.simplify(sp.diff(dens, eps, 2).subs(eps, 0) / 2)
    L2 = sp.simplify(sp.integrate(sp.expand(d2), (th, 0, sp.pi)))
    Nl = sp.Rational(LL * (LL + 1), 2 * LL + 1)
    rho = P1 / N ** 2 - P0
    esperado = 2 * Nl * (-rho * h0 ** 2 / N - N * P0 * h1 ** 2) / 2
    # la integral angular de sqrt(-g) incluye sin(theta); la normalizacion de la Fase 2f usa N_l = int sin th (P_l')^2 / 2
    assert sp.simplify(L2 - esperado) == 0, (LL, sp.simplify(L2), esperado)
    print(f"l = {LL}: L2^fluido = N_l [ -rho~ h0^2/N - N p~ h1^2 ]  (exacta)")
