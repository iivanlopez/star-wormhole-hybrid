# -*- coding: utf-8 -*-
"""
Fase 2l, parte V: rotacion lenta en Horndeski quintico (con F4) a partir de resultados/kk_quintico.pkl.

1. Forma estatica de J: J = (2/3) r^4 Gt/N con Gt = G4 - 2X G4_X - 4X^2 F4 + X (phi' r'/r G5_X + G5_phi)
   (el coeficiente de al1 de la Fase 2j, con r'/r; no el Gcal del teorema, con Phi').
2. Termino de Routh L_J = -pi^2/(2J) sobre la prueba d = N f, c = a = beta = 0 (estatica): J depende de C y C', asi que
   la energia centrifuga tiene ahora un termino en d'.  Se da la densidad exacta para d = N f.
3. Ligadura: terminos cruzados de L2_J entre (a, beta) y (c, d) en el gauge unitario, en 3 fondos racionales al azar
   (funciones de la teoria por Taylor en X alrededor de X_b; phi no se perturba).  Si se anulan, W no cambia.
Uso: python quintico_rotacion.py
"""
import os
import re
import json
import pickle
import random
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
d = pickle.load(open(os.path.join(AQUI, "resultados", "kk_quintico.pkl"), "rb"))
J, Xj = d["J"], d["Xj"]
S = lambda n: sp.Symbol(n, real=True)                       # noqa: E731
GP = re.compile(r"^(G2|G3|G4|F4|G5)_(\d)(\d)$")
pi2 = sp.Symbol("pi2", positive=True)
r0, N = sp.symbols("r0 N", positive=True)
Phi1, r1, r2, P1, P2, A2 = sp.symbols("Phi1 r1 r2 P1 P2 A2", real=True)
out = {}

# ------------------------------------------------------------------ 1-2. prueba estatica d = N f
est = {s: 0 for s in J.free_symbols if s.name[:2] in ("A_", "b_", "B_", "C_", "p_")}
est.update({S("A_00"): N ** 2, S("A_01"): 2 * Phi1 * N ** 2, S("B_00"): 1, S("p_01"): P1, S("p_02"): P2})
Cs, C1s = S("C_00"), S("C_01")
est.pop(Cs)
est.pop(C1s)
Js0 = J.xreplace(est)
Xb = -P1 ** 2 / 2
G = lambda n, i, j: S(f"{n}_{i}{j}")                        # noqa: E731
Gt = G("G4", 0, 0) - 2 * Xb * G("G4", 0, 1) - 4 * Xb ** 2 * G("F4", 0, 0)     + Xb * (P1 * (C1s / (2 * Cs)) * G("G5", 0, 1) + G("G5", 1, 0))
Js = sp.Rational(2, 3) * Cs ** 2 * Gt / N
rnd = random.Random(9)
for k in range(5):                              # forma cerrada frente a la J completa: 5 puntos, 50 cifras
    val = {x: sp.Rational(rnd.randint(200, 900), 997) for x in Js0.free_symbols | Js.free_symbols}
    a_, b_ = sp.N(Js0.xreplace(val), 50), sp.N(Js.xreplace(val), 50)
    assert abs(a_ - b_) < 1e-40 * abs(b_), (a_, b_)
print("J estatico(C, C') = (2/3) C^2 Gt/N con Gt = G4 - 2X G4X - 4X^2 F4 + X(phi' C'/(2C) G5X + G5phi): "
      "comprobado en 5 puntos con 50 cifras")
out["J_estatico"] = "J = (2/3) C^2 Gt/N,  Gt = G4 - 2X G4X - 4X^2 F4 + X(phi' C'/(2C) G5X + G5phi)"
eps, f0, f1 = sp.symbols("eps f0 f1", real=True)
dd, dd1 = N * f0, N * (f1 + Phi1 * f0)                     # d = N f, d' = N (f' + Phi' f)
sub = {Cs: r0 ** 2 + 2 * r0 * eps * dd, C1s: 2 * r0 * r1 + 2 * eps * (r1 * dd + r0 * dd1)}
LJ = -pi2 / (2 * Js)
seg = sp.simplify(sp.diff(LJ.xreplace(sub), eps, 2).subs(eps, 0) / 2)
dens = sp.factor(-seg)
Gt0 = Gt.xreplace({Cs: r0 ** 2, C1s: 2 * r0 * r1})
print(f"-L2_J(d = N f) = {dens}")
# comparacion con el caso cuartico (G5 = 0): 9 pi^2 N^3 f^2/(Gcal r^6)
cuart = sp.simplify(dens.xreplace({G("G5", 0, 1): 0, G("G5", 1, 0): 0}) - 9 * pi2 * N ** 3 * f0 ** 2 /
                    (Gt0.xreplace({G("G5", 0, 1): 0, G("G5", 1, 0): 0}) * r0 ** 6))
print(f"limite cuartico: diferencia = {cuart}")
assert cuart == 0
y = sp.Symbol("y")                                          # y = X phi' G5X / (2 r Gt): peso del termino en r'
print(f"densidad con f = 1 (prueba d = N): {sp.factor(dens.subs({f0: 1, f1: 0}))}")
out["densidad_dN"] = str(dens)

# ------------------------------------------------------------------ 3. terminos cruzados (a, beta) x (c, d)
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]
pert = {v: {o: S(f"{v}_{o[0]}{o[1]}") for o in ORD} for v in ("a", "be", "c", "d")}


def fondo_aleatorio(sem):
    rnd = random.Random(sem)
    q = lambda a=-900, b=900: sp.Rational(rnd.randint(a, b), 997)       # noqa: E731
    v = dict(N=sp.Rational(rnd.randint(40, 99), 100), Phi1=q(), Phi2=q(), r0=sp.Rational(rnd.randint(700, 2000), 1000),
             r1=q(), r2=q(), P1=q(300, 900), P2=q())
    gv = {}
    for s in J.free_symbols:
        m = GP.match(s.name)
        if m:
            for k in range(3):
                gv[(m.group(1), int(m.group(2)), int(m.group(3)) + k)] = q()
    return v, gv


def cruzados(sem):
    v, gv = fondo_aleatorio(sem)
    Nn, E = v["N"], v["N"] ** 2
    # fondo + eps * perturbacion (mismas convenciones que fase2f/generar_bH.dpert)
    bg = {"A": [E, 2 * v["Phi1"] * E, (2 * v["Phi2"] + 4 * v["Phi1"] ** 2) * E], "C":
          [v["r0"] ** 2, 2 * v["r0"] * v["r1"], 2 * v["r1"] ** 2 + 2 * v["r0"] * v["r2"]]}
    a, be, c, dd_ = pert["a"], pert["be"], pert["c"], pert["d"]
    rep = {}
    E2 = 2 * E
    rep[S("A_00")] = E + eps * E2 * a[(0, 0)]
    rep[S("A_10")] = eps * E2 * a[(1, 0)]
    rep[S("A_01")] = bg["A"][1] + eps * E2 * (2 * v["Phi1"] * a[(0, 0)] + a[(0, 1)])
    for o in ORD:
        rep[S(f"b_{o[0]}{o[1]}")] = eps * be[o]
        rep[S(f"B_{o[0]}{o[1]}")] = (1 if o == (0, 0) else 0) + eps * 2 * c[o]
    rep[S("C_00")] = bg["C"][0] + eps * 2 * v["r0"] * dd_[(0, 0)]
    rep[S("C_10")] = eps * 2 * v["r0"] * dd_[(1, 0)]
    rep[S("C_01")] = bg["C"][1] + eps * 2 * (v["r1"] * dd_[(0, 0)] + v["r0"] * dd_[(0, 1)])
    rep.update({S("p_01"): v["P1"], S("p_02"): v["P2"], S("p_10"): 0, S("p_11"): 0, S("p_20"): 0})
    Xnew = Xj.xreplace(rep)
    Xb_ = -v["P1"] ** 2 / 2
    grep = {}
    for s in J.free_symbols:
        m = GP.match(s.name)
        if m:
            nm, i, j = m.group(1), int(m.group(2)), int(m.group(3))
            grep[s] = sum(gv[(nm, i, j + k)] * (Xnew - Xb_) ** k / sp.factorial(k) for k in range(3))
    Je = J.xreplace(rep).xreplace(grep)
    L = -1 / (2 * Je)                                      # pi^2 = 1
    L2 = sp.expand(sp.diff(L, eps, 2).subs(eps, 0) / 2)
    nulos, no_nulos = [], []
    for x in list(a.values()) + list(be.values()):
        for y_ in list(c.values()) + list(dd_.values()):
            co = sp.expand(L2.coeff(x).coeff(y_)) if L2.has(x) and L2.has(y_) else 0
            (nulos if co == 0 else no_nulos).append((str(x), str(y_), str(co)))
    jets_a_be = sorted({str(s) for s in L2.free_symbols if s.name.startswith(("a_", "be_"))})
    # forma cerrada del termino que modifica W: L2_J ⊃ kapJ beta d_t + N^2 kapJ a d' (= -F kapJ d, salvo derivadas)
    Xv, P1v = -v["P1"] ** 2 / 2, v["P1"]
    g_ = lambda nm, i, j: gv.get((nm, i, j), 0)            # noqa: E731
    Gt_v = g_("G4", 0, 0) - 2 * Xv * g_("G4", 0, 1) - 4 * Xv ** 2 * g_("F4", 0, 0)         + Xv * (P1v * v["r1"] / v["r0"] * g_("G5", 0, 1) + g_("G5", 1, 0))
    kap = sp.Rational(3, 4) * Xv * P1v * g_("G5", 0, 1) / (Nn * v["r0"] ** 5 * Gt_v ** 2)
    L2e = sp.expand(L2)
    c_bd = sp.expand(L2e.coeff(be[(0, 0)]).coeff(dd_[(1, 0)]))
    c_ad1 = sp.expand(L2e.coeff(a[(0, 0)]).coeff(dd_[(0, 1)]))
    return no_nulos, jets_a_be, dict(beta_dt=c_bd, kapJ=kap, a_dl=c_ad1, N2kapJ=E * kap,
                                     ok=bool(sp.simplify(c_bd - kap) == 0 and sp.simplify(c_ad1 - E * kap) == 0))


for sem in range(3):
    nn, jab, kp = cruzados(40 + sem)
    print(f"[fondo {sem}] chorros de (a, beta) en L2_J: {jab};  terminos cruzados (a, beta) x (c, d) no nulos: {nn}")
    print(f"[fondo {sem}] beta d_t: {kp['beta_dt']} frente a kapJ = 3 X phi' G5X/(4 N r^5 Gt^2) = {kp['kapJ']};  "
          f"a d': {kp['a_dl']} frente a N^2 kapJ = {kp['N2kapJ']};  coinciden: {kp['ok']}")
    out[f"cruzados_{sem}"] = dict(no_nulos=nn, chorros_a_beta=jab, kapJ={k: str(x) for k, x in kp.items()})
json.dump(out, open(os.path.join(AQUI, "resultados", "quintico_rotacion.json"), "w"), indent=1)
