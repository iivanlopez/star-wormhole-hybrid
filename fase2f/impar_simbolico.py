# -*- coding: utf-8 -*-
"""
Sector impar (axial) de la teoria cuartica mas alla de Horndeski sobre un fondo estatico esferico (SymPy).

Fondo: ds^2 = -e^{2 Phi(l)} dt^2 + dl^2 + r(l)^2 dOmega^2, phi = P(l).
Perturbacion axial en el gauge de Regge-Wheeler, m = 0:
  g_{t varphi} = eps h0(t, l) W(theta),  g_{l varphi} = eps h1(t, l) W(theta),  W = sin(theta) dP_L(cos theta)/dtheta.
Se calcula la lagrangiana sqrt(-g) L hasta O(eps^2) con la inversa en serie exacta (g = g0 + eps g1:
g^-1 = g0^-1 - eps g0^-1 g1 g0^-1 + eps^2 g0^-1 g1 g0^-1 g1 g0^-1), se integra en theta y se obtienen las ecuaciones
de Euler-Lagrange para modos e^{i omega t} con el operador D_l de generar_bH.py (simbolos de fondo y G_ij).
X = -(1/2) g^{ll} P'^2 cambia solo a O(eps^2): L2 = (1/2) d^2/deps^2 L(eps, X_b) + L_X X_2.
Salida: resultados/impar_L{L}.pkl.   Uso: python impar_simbolico.py 2   (y 3, 4)
"""
import os
import sys
import time
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from bH_simbolico import COORD, t, l, th, vp      # noqa: E402

T0 = time.time()
LL = int(sys.argv[1]) if len(sys.argv) > 1 else 2


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


eps = sp.Symbol("eps")
Xs = sp.Symbol("Xs", real=True)
Phif, rf, Pf = (sp.Function(n)(l) for n in ("Phi", "r", "P"))
h0, h1 = sp.Function("h0")(t, l), sp.Function("h1")(t, l)
W = sp.simplify(sp.sin(th) * sp.diff(sp.legendre(LL, sp.cos(th)), th))
g0 = sp.diag(-sp.exp(2 * Phif), 1, rf**2, rf**2 * sp.sin(th) ** 2)
g1 = sp.zeros(4, 4)
g1[0, 3] = g1[3, 0] = h0 * W
g1[1, 3] = g1[3, 1] = h1 * W
g = g0 + eps * g1
g0i = sp.diag(*[1 / g0[i, i] for i in range(4)])
gi = g0i - eps * g0i * g1 * g0i + eps**2 * g0i * g1 * g0i * g1 * g0i


def trunc(e, n=2):
    e = sp.expand(e)
    return sum(e.coeff(eps, k) * eps**k for k in range(n + 1))


gi = gi.applyfunc(trunc)
n = 4
dg = [[[sp.diff(g[a, b], COORD[c]) for c in range(n)] for b in range(n)] for a in range(n)]
Gam = [[[trunc(sum(gi[a, d] * (dg[d][b][c] + dg[d][c][b] - dg[b][c][d]) for d in range(n)) / 2)
         for c in range(n)] for b in range(n)] for a in range(n)]
log("Christoffel")
Ric = sp.zeros(n, n)
for b_ in range(n):
    for d in range(b_, n):
        v = sum(sp.diff(Gam[a][b_][d], COORD[a]) - sp.diff(Gam[a][b_][a], COORD[d]) for a in range(n))
        v += sum(Gam[a][a][e] * Gam[e][b_][d] - Gam[a][d][e] * Gam[e][b_][a] for a in range(n) for e in range(n))
        Ric[b_, d] = Ric[d, b_] = trunc(v)
Rs = trunc(sum(gi[a, b_] * Ric[a, b_] for a in range(n) for b_ in range(n)))
log("Ricci")
# escalar: phi = P(l)
d1 = [0, sp.diff(Pf, l), 0, 0]
d2 = sp.Matrix(n, n, lambda m, k: trunc(sp.diff(Pf, COORD[m], COORD[k]) - Gam[1][m][k] * d1[1]))
up = [trunc(gi[m, 1] * d1[1]) for m in range(n)]
S = sp.Matrix(n, n, lambda a, k: trunc(sum(gi[a, m] * d2[m, k] for m in range(n))))
Xexact = trunc(-up[1] * d1[1] / 2)
box = trunc(sum(S[a, a] for a in range(n)))
trS2 = trunc(sum(S[a, k] * S[k, a] for a in range(n) for k in range(n)))
pSp = trunc(sum(up[m] * d2[m, k] * up[k] for m in range(n) for k in range(n)))
pS2p = trunc(sum(up[m] * d2[m, k] * S[k, j] * up[j] for m in range(n) for k in range(n) for j in range(n)))
log("escalares del campo phi")
G2f, G3f, G4f, F4f = (sp.Function(nm) for nm in ("G2", "G3", "G4", "F4"))
G2, G3, G4, F4 = (f(Pf, Xs) for f in (G2f, G3f, G4f, F4f))
G4X = sp.Derivative(G4f(Pf, Xs), Xs)
# el factor X explicito de E4 es el exacto (cambia a O(eps^2)); las funciones G se evaluan en Xs y su dependencia
# en X se recupera con la regla de la cadena L_X X_2
E4 = trunc(2 * Xexact * (box**2 - trS2) + 2 * box * pSp - 2 * pS2p)
Lden = G2 - G3 * box + G4 * Rs + G4X * (box**2 - trS2) + F4 * E4
# sqrt(-g) hasta O(eps^2): det(1 + eps M) = 1 + eps tr M + eps^2 [(tr M)^2 - tr(M^2)]/2, con tr M = 0
Mg = g0i * g1
sq = sp.exp(Phif) * rf**2 * sp.sin(th) * (1 - eps**2 * (Mg * Mg).trace() / 4)
Ltot = trunc(sp.expand(sq * Lden))
X2 = sp.expand(Xexact).coeff(eps, 2)
L2 = sp.expand(Ltot.coeff(eps, 2) + sp.diff(sp.expand(sq * Lden).coeff(eps, 0), Xs) * X2)
log(f"lagrangiana a O(eps^2): {len(sp.Add.make_args(L2))} terminos")
L2 = sp.integrate(L2, (th, 0, sp.pi))
L2 = sp.expand(L2)
log(f"integrada en theta: {len(sp.Add.make_args(L2))} terminos")
with open(os.path.join(AQUI, "resultados", f"impar_L{LL}.pkl"), "wb") as fh:
    pickle.dump(dict(L2=L2, LL=LL), fh)
log(f"guardado resultados/impar_L{LL}.pkl")
