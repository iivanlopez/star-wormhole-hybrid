# -*- coding: utf-8 -*-
"""
Fase 2g, parte II: estructura del sector l = 0 con fluido (GLPV cuartica + fluido P(Y)) y reduccion a dos campos.

Gauges: unitario para el escalar (q = 0) y comovil para el fluido (x = delta chi = 0); ambos fijan por completo las
reparametrizaciones zeta (de l) y tau (de t).  Quedan a, beta, c, d.  Se comprueba, con chorros exactos en la capa de
masas con fluido (aritmetica racional), que la accion cuadratica tiene la forma
    L2 = -(N^2 a' + beta_t) W + kap_a a^2 + a M - kap_b beta^2 + L_rest(c, d)        (salvo derivadas totales)
con
    W     = gam c + al1 d' + al0 d                    (los coeficientes de la Fase 2f, sin cambio: al0 = -Phi' al1)
    kap_b = r^2 (rho~ + p~)/(2N)                      (energia cinetica del fluido)
    kap_a = N r^2 (rho~ + p~)/(2 c_s^2)               (compresibilidad)
    M     = -N (rho~ + p~)(r^2 c + 2 r d)             (acoplo del lapso a la variacion de volumen)
    L_rest(c, d) = L2(a = beta = 0).
Consecuencias (ecuaciones de beta y a):  beta = W_t/(2 kap_b) y a = -[(N^2 W)' + M]/(2 kap_a), y con el desplazamiento
del fluido xi (dl/dt = -beta):  W = -2 kap_b xi = -(r^2/N)(rho~ + p~) xi  (ligadura de momento con el fluido).
Accion reducida en (xi, d):
    L_red = kap_b xi_t^2 - Q^2/(4 kap_a) + L_rest(c*, d),   Q = (N^2 W)' + M,   c* = (W - al1 d' - al0 d)/gam.
Prueba d = N, xi = 0:  W = 0, c* = 0, Q = M(c = 0, d = N) = -2 N^2 r (rho~ + p~) y
    E[N, 0] = int [ Q^2/(4 kap_a) - L_rest(0, N) ] = -2 int N^3 [ Gcal (r'/r - Phi')^2 - (rho~ + p~) c_s^2 ] dl
            = -2 int N^3 [ Gcal (r'/r - Phi')^2 - 16 pi Gamma_1 p ] dl       (rho~ + p~ = 16 pi (rho + p)).
La energia cinetica de esta prueba es la del vacio (el fluido no se mueve): int K N^2 (con K de la Fase 2f).
Escribe resultados/reduccion_fluido.json.
Uso: python reduccion_fluido.py
"""
import os
import sys
import json
import time
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
sys.path.insert(0, F2F)
import chorros as CH                                       # noqa: E402
from chorros import Phi, rr, PP, Ee, eP, Xs              # noqa: E402

T0 = time.time()
RES = {}


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


d_ = pickle.load(open(os.path.join(AQUI, "resultados", "fluido_simbolico.pkl"), "rb"))
L2g, L2m, fe, Pf = d_["L2_g"], d_["L2_m"], d_["fondo_eq"], d_["Pf"]
mr = pickle.load(open(os.path.join(F2F, "resultados", "monopolo_reduccion.pkl"), "rb"))
NPF = len(Pf)
Xb = {Xs: -PP[1] ** 2 / 2}
N = sp.sqrt(Ee)


def S(n):
    return sp.Symbol(n, real=True)


def Dbg(e, amps=()):
    out = CH.D(e, amps)
    for k in range(NPF - 1):
        if e.has(Pf[k]):
            out += sp.diff(e, Pf[k]) * (-Pf[k + 1] * Phi[1] / Ee)
    return sp.expand(out)


def capa(semilla, orden=5):
    """Chorro exacto en la capa de masas con fluido (como generar_fluido.chorro_capa_fluido)."""
    v = CH.Valores(semilla)
    v[Ee] = sp.Rational(121, 100)
    v[eP] = sp.Rational(11, 10)
    v[Phi[1]] = v.r()
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[rr[1]] = v.r()
    v[PP[0]] = v.r()
    v[PP[1]] = v.r()
    v[Xs] = -v[PP[1]] ** 2 / 2
    for k in range(NPF):
        v[Pf[k]] = v.r()
    g200 = CH.gsym("G2", 0, 0)
    eB = v.ev(fe["B"], salvo=(g200,))
    c1 = eB.coeff(g200)
    v[g200] = -sp.expand(eB - c1 * g200) / c1
    eqs = [fe["A"], fe["C"], fe["p"]]
    for k in range(2, orden + 1):
        inc = [Phi[k], rr[k], PP[k]]
        ev = [v.ev(e) for e in eqs]
        Mk = sp.Matrix([[x.coeff(u) for u in inc] for x in ev])
        bk = sp.Matrix([sp.expand(x - sum(x.coeff(u) * u for u in inc)) for x in ev])
        sol = Mk.LUsolve(-bk)
        for u, s in zip(inc, sol):
            assert s.is_Rational
            v[u] = s
        if k < orden:
            eqs = [Dbg(e) for e in eqs]
    return v


# ------------------------------------------------------------------ lagrangiana en el gauge q = x = 0, chorros (t, l)
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
cero = {S(f"{v}_{o[0]}{o[1]}"): 0 for v in ("q", "x") for o in ORD}
L = sp.expand((L2g + L2m.subs(eP, N)).xreplace(cero).xreplace(Xb))
NMAX = 6
CAMPOS = ("a", "be", "c", "d")
J = {f: {(i, j): S(f"{f}{i}{j}") for i in range(NMAX + 1) for j in range(NMAX + 1)} for f in CAMPOS}
L = L.xreplace({S(f"{f}_{o[0]}{o[1]}"): J[f][o] for f in CAMPOS for o in ORD})
FILAS = [tuple(J[f][(i, j)] for j in range(NMAX + 1)) for f in CAMPOS for i in range(NMAX + 1)]


def Dl(e, extra=()):
    return Dbg(e, tuple(FILAS) + tuple(extra))


def Dt(e):
    out = 0
    for f in CAMPOS:
        for (i, j), s in J[f].items():
            if i < NMAX and e.has(s):
                out += sp.diff(e, s) * J[f][(i + 1, j)]
    return sp.expand(out)


def euler_lagrange(Lg, f, extra=()):
    tot = 0
    for (i, j), s in J[f].items():
        if i + j > 3 or not Lg.has(s):
            continue
        Q = sp.expand(sp.diff(Lg, s)) * (-1) ** (i + j)
        for _ in range(i):
            Q = Dt(Q)
        for _ in range(j):
            Q = Dl(Q, extra)
        tot += Q
    return sp.expand(tot)


EL = {f: euler_lagrange(L, f) for f in CAMPOS}
log("ecuaciones de Euler-Lagrange (gauge q = x = 0): " + ", ".join(f"{f}: {len(sp.Add.make_args(e))}" for f, e in EL.items()))

# ------------------------------------------------------------------ hipotesis
gam = sp.expand(mr["gam"].xreplace(Xb))
al1 = sp.expand(mr["al1"].xreplace(Xb))
al0 = sp.expand(mr["al0"].xreplace(Xb))
rho = Pf[1] / Ee - Pf[0]
pre = Pf[0]
rp = Pf[1] / Ee                                        # rho~ + p~
cs2 = Pf[1] / (Pf[1] + Pf[2] / Ee)
kap_b = rr[0] ** 2 * rp / (2 * N)
kap_a = N * rr[0] ** 2 * rp / (2 * cs2)
c, d, a, b = J["c"], J["d"], J["a"], J["be"]
W = gam * c[(0, 0)] + al1 * d[(0, 1)] + al0 * d[(0, 0)]
M = -N * rp * (rr[0] ** 2 * c[(0, 0)] + 2 * rr[0] * d[(0, 0)])
F = Ee * a[(0, 1)] + b[(1, 0)]
Lrest = sp.expand(L.xreplace({s: 0 for f in ("a", "be") for s in J[f].values()}))
Lhyp = sp.expand(-F * W + kap_a * a[(0, 0)] ** 2 + a[(0, 0)] * M - kap_b * b[(0, 0)] ** 2 + Lrest)
ELh = {f: euler_lagrange(Lhyp, f) for f in CAMPOS}
for sem in range(3):
    v = capa(1300 + sem, orden=6)
    for f in CAMPOS:
        dif = sp.expand(EL[f] - ELh[f])
        val = v.ev(dif)
        # los chorros de los campos quedan libres: el resultado debe ser identicamente 0 en ellos
        assert val == 0, (f, val)
RES["forma"] = ("L2 = -(N^2 a' + beta_t) W + kap_a a^2 + a M - kap_b beta^2 + L_rest(c, d): ecuaciones de Euler-Lagrange "
                "identicas en 3 chorros exactos en la capa de masas con fluido")
log("forma comprobada: L2 = -(N^2 a' + beta_t) W + kap_a a^2 + a M - kap_b beta^2 + L_rest(c, d)  (E-L exactas, 3 chorros)")

# ------------------------------------------------------------------ prueba d = N, xi = 0: energia estatica
# L_rest(c = 0, d) es la de la Fase 2f (el fluido no aporta d^2 ni d d'): -L_rest(0, N) = -2 N^3 Gcal (r'/r - Phi')^2
Lr0 = sp.expand(Lrest.xreplace({s: 0 for s in J["c"].values()}))
Lm_d = sp.expand(L2m.subs(eP, N).xreplace(cero).xreplace({S(f"{v}_{o[0]}{o[1]}"): 0 for v in ("a", "be", "c")
                                                         for o in ORD}))
assert Lm_d == 0
G4, G401, F400 = S("G4_00"), S("G4_01"), S("F4_00")
Gc = G4 + G401 * PP[1] ** 2 - F400 * PP[1] ** 4
dN = {(0, 0): N, (0, 1): N * Phi[1], (0, 2): N * (Phi[2] + Phi[1] ** 2)}
menosL = sp.expand(-Lr0.xreplace({d[o]: dN.get(o, 0) for o in d}))
assert sp.simplify(menosL + 2 * N ** 3 * Gc * (rr[1] / rr[0] - Phi[1]) ** 2) == 0
Q0 = M.xreplace({c[(0, 0)]: 0, d[(0, 0)]: N})
EN_fluido = sp.simplify(Q0 ** 2 / (4 * kap_a) - 2 * N ** 3 * rp * cs2)
assert EN_fluido == 0
RES["prueba_N"] = "E[N, xi=0] = -2 int N^3 [Gcal (r'/r - Phi')^2 - (rho~ + p~) c_s^2] (exacta)"
log("prueba d = N, xi = 0:  -L_rest(0, N) = -2 N^3 Gcal (r'/r - Phi')^2 (el fluido no aporta d^2) y "
    "Q^2/(4 kap_a) = 2 N^3 (rho~ + p~) c_s^2  =>  E = -2 int N^3 [Gcal (r'/r - Phi')^2 - (rho~ + p~) c_s^2]")

# ------------------------------------------------------------------ guardar la lagrangiana para la reduccion numerica
with open(os.path.join(AQUI, "resultados", "reduccion_fluido.pkl"), "wb") as fh:
    pickle.dump(dict(Lrest=Lrest, J=J, gam=gam, al1=al1, al0=al0, kap_a=kap_a, kap_b=kap_b, M=M), fh)
json.dump(RES, open(os.path.join(AQUI, "resultados", "reduccion_fluido.json"), "w"), indent=1)
log("guardado resultados/reduccion_fluido.pkl y .json")
