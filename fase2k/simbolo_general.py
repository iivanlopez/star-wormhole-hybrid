# -*- coding: utf-8 -*-
"""
Fase 2k (parte II): simbolo principal del monopolo en un punto cualquiera del fondo (con desplazamiento v != 0).
Como simbolo_garganta.py pero con valores generales: los coeficientes congelados del operador de monopolo_cubico.py
(C[q,f,m,n]) evaluados con simbolos para el fondo y sus derivadas; M(om, k) = sum C (-i om)^m (i k)^n.
Partes principales (om = c k, k -> oo):  det M ~ k^4 P(c),  det M_aux ~ k^2 D(c)  (M_aux: bloque de nu y w).
Se guardan P y D (polinomios en c) para la criba (criba_cubica.py).
Uso: python simbolo_general.py
"""
import os
import sys
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import monopolo_cubico as MC     # noqa: E402

op = MC.operador()
C = op["C"]
om, k, c = sp.symbols("omega k c")
CAMPOS = ("nu", "w", "d")
M = sp.zeros(3, 3)
for (q, f, m, n), co in C.items():
    M[CAMPOS.index(q), CAMPOS.index(f)] += co * (-sp.I * om) ** m * (sp.I * k) ** n
M = M.applyfunc(sp.expand)


def lider(expr):
    Pk = sp.Poly(sp.expand(expr.subs(om, c * k)), k)
    return Pk.degree(), sp.expand(Pk.coeff_monomial(k ** Pk.degree()))


g, P = lider(M.det(method="berkowitz"))
gx, Dx = lider(M[:2, :2].det())
print(f"det M ~ k^{g}, det M_aux ~ k^{gx}", flush=True)
Pc = sp.Poly(P, c)
print("grado de P en c:", Pc.degree(), flush=True)
usados = sorted({str(s) for s in P.free_symbols | Dx.free_symbols} - {"c"})
print("simbolos de fondo en P y D:", usados, flush=True)
pickle.dump(dict(P=P, D=Dx, gP=g, gD=gx, usados=usados),
            open(os.path.join(AQUI, "resultados", "simbolo_general.pkl"), "wb"))
print("D =", sp.factor(Dx), flush=True)
