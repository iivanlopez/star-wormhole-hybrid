# -*- coding: utf-8 -*-
"""
Fase 2l, parte IV: la ecuacion de c en la garganta con rotacion (el analogo de la obstruccion de la carga, Fase 2g §7).

La ecuacion de c del problema l = 0 es  E_c = k_c c + e0 d + e1 d' + e_tt d_tt - gam F = 0  (fase2g/coeficientes_fluido.py,
sin fluido), y en una garganta simetrica gam(0) = 0.  Con carga, e0_tot(0) = 4 q^2 N0/r0^3 != 0 obligaba a d(0) = 0.
Con rotacion, L_J aporta k_c -> k_c + 2 kcc y e0 -> e0 + kcd (monopolo_rotacion.py).  Aqui se evaluan k_c y e0 totales
en una garganta simetrica (r' = Phi' = phi'' = 0) que cumple las ecuaciones del fondo con L_J: se resuelven E_A, E_B,
E_C, E_p y k_c(0) = 0 para (G2_00, G2_02, Phi'', r''), con el resto de funciones de la teoria y (N0, r0, phi'(0)) racionales al azar (con la simetria phi -> -phi: G2, G4, F4
pares y G3 impar en phi, phi(0) = 0) y
pi^2 simbolico, imponiendo las condiciones de regularidad de K de la Fase 2f (H(0) = 0 y k_c(0) = 0; esta ultima
fija G2_XX(0)).  Control: con L_em = -q^2 sqrt(AB+b^2)/C en lugar de L_J se recupera e0_tot(0) = 4 q^2 N0/r0^3.
Uso: python garganta_rotacion.py        -> resultados/garganta_rotacion.json
"""
import os
import sys
import json
import pickle
import random
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
F2F = os.path.join(AQUI, "..", "fase2f")
F2G = os.path.join(AQUI, "..", "fase2g")
sys.path.insert(0, F2F)
from chorros import Phi, rr, PP, Ee, Xs, gsym      # noqa: E402

bh = pickle.load(open(os.path.join(F2F, "resultados", "bH_simbolico.pkl"), "rb"))
mr = pickle.load(open(os.path.join(AQUI, "resultados", "monopolo_rotacion.pkl"), "rb"))
cf = pickle.load(open(os.path.join(F2G, "resultados", "coeficientes_fluido.pkl"), "rb"))
pi2 = sp.Symbol("pi2", positive=True)
q2 = sp.Symbol("q2", positive=True)


def por_nombre(e):
    """unifica simbolos por nombre (los pickles de distintas fases crean simbolos con supuestos distintos)."""
    canon = {s.name: s for s in list(Phi) + list(rr) + list(PP) + [Ee, Xs, pi2, q2]}
    return e.xreplace({s: canon.get(s.name, sp.Symbol(s.name, real=True)) for s in e.free_symbols})


fe0 = {k: por_nombre(v) for k, v in bh["fondo_eq"].items()}
feJ = {k: por_nombre(v) for k, v in mr["feJ"].items()}
Pf0 = sp.Symbol("Pf0", real=True)
kc_v = por_nombre(cf["k_c"].xreplace({s: 0 for s in cf["k_c"].free_symbols if s.name == "Pf0"}))
e0_v = por_nombre(cf["e0"].xreplace({s: 0 for s in cf["e0"].free_symbols if s.name == "Pf0"}))
kcc, kcd = por_nombre(mr["kcc"]), por_nombre(mr["kcd"])
N = sp.sqrt(Ee)
# control: carga (Fase 2g): L_em = -q^2 sqrt(AB + b^2)/C -> fondo: E_A = -q^2/(2 N r^2)... se obtiene igual que L_J
Aj, Bj, Cj = sp.symbols("A_00 B_00 C_00", positive=True)
Lem = -q2 * sp.sqrt(Aj * Bj) / Cj
feEM = {"A": sp.diff(Lem, Aj).subs({Aj: Ee, Bj: 1, Cj: rr[0] ** 2}),
        "B": sp.diff(Lem, Bj).subs({Aj: Ee, Bj: 1, Cj: rr[0] ** 2}),
        "C": sp.diff(Lem, Cj).subs({Aj: Ee, Bj: 1, Cj: rr[0] ** 2}), "b": 0, "p": 0}
kc_em, e0_em = q2 * N / rr[0] ** 2, 2 * q2 * N / rr[0] ** 3


def en_garganta(fuente, kc_extra, e0_extra, semilla):
    rnd = random.Random(semilla)
    v = {Phi[1]: 0, rr[1]: 0, PP[2]: 0, Ee: sp.Rational(rnd.randint(30, 99), 100) ** 2,
         rr[0]: sp.Rational(rnd.randint(600, 2000), 1000), PP[1]: sp.Rational(rnd.randint(300, 900), 997),
         PP[0]: 0}
    v[Xs] = -v[PP[1]] ** 2 / 2
    inc = [gsym("G2", 0, 0), gsym("G2", 0, 2), Phi[2], rr[2]]
    eqs = [sp.expand(fe0[k] + fuente[k]) for k in ("A", "B", "C", "p")]
    todos = set().union(*[e.free_symbols for e in eqs], kc_v.free_symbols, e0_v.free_symbols,
                        sp.sympify(kc_extra).free_symbols, sp.sympify(e0_extra).free_symbols)
    for s in sorted(todos, key=str):
        if s not in v and s not in inc and s not in (pi2, q2) and s.name[:3] in ("G2_", "G3_", "G4_", "F4_"):
            i = int(s.name[3])
            # simetria phi -> -phi de la garganta (phi impar): G2, G4, F4 pares y G3 impar en phi
            v[s] = 0 if (i % 2 == 1) != (s.name[:2] == "G3") else sp.Rational(rnd.randint(-900, 900), 997)
    v[gsym("G4", 0, 0)] = 1 + sp.Rational(rnd.randint(0, 500), 997)
    # regularidad de K en la garganta (Fase 2f §5): H(0) = G4 + G4_X phi'^2 = 0;  Gcal(0) = -F4 phi'^4 > 0
    v[gsym("G4", 0, 1)] = -v[gsym("G4", 0, 0)] / v[PP[1]] ** 2
    v[gsym("F4", 0, 0)] = -abs(v.get(gsym("F4", 0, 0), sp.Rational(1, 3))) - sp.Rational(1, 10)
    for s in sorted(todos, key=str):
        if s not in v and s not in inc and s not in (pi2, q2):
            if s.name.startswith(("Phi", "r", "P")):
                v[s] = sp.Rational(rnd.randint(-900, 900), 997)              # derivadas de orden alto: no entran
    ev = [sp.expand(e.xreplace(v)) for e in eqs] + [sp.expand((kc_v + kc_extra).xreplace(v))]   # k_c(0) = 0
    sol = sp.solve(ev, inc, dict=True)
    assert len(sol) == 1, sol
    v.update(sol[0])
    kc = sp.simplify((kc_v + kc_extra).xreplace(v))
    e0 = sp.simplify((e0_v + e0_extra).xreplace(v))
    return kc, e0, v


if __name__ == "__main__":
    out = {}
    for sem in range(3):
        kc, e0, v = en_garganta(feEM, kc_em, e0_em, 100 + sem)
        ref = 4 * q2 * sp.sqrt(v[Ee]) / v[rr[0]] ** 3
        print(f"[carga, chorro {sem}] k_c(0) = {kc},  e0(0) = {e0}  (4 q^2 N0/r0^3 = {ref});  "
              f"e0 = 4 q^2 N0/r0^3: {sp.simplify(e0 - ref) == 0}")
        out[f"carga_{sem}"] = dict(kc=str(kc), e0=str(e0), coincide=bool(sp.simplify(e0 - ref) == 0))
        kc, e0, v = en_garganta(feJ, 2 * kcc, kcd, 200 + sem)
        G = (gsym("G4", 0, 0) + gsym("G4", 0, 1) * PP[1] ** 2 - gsym("F4", 0, 0) * PP[1] ** 4).xreplace(v)
        cand = 6 * pi2 * sp.sqrt(v[Ee]) / (G * v[rr[0]] ** 5)
        print(f"[rotacion, chorro {sem}] k_c(0) = {kc},  e0(0) = {e0};  e0/(6 pi^2 N0/(Gcal r0^5)) = "
              f"{sp.simplify(e0 / cand)}")
        out[f"rotacion_{sem}"] = dict(kc=str(kc), e0=str(e0), cociente_6pi2N_Gr5=str(sp.simplify(e0 / cand)))
    json.dump(out, open(os.path.join(AQUI, "resultados", "garganta_rotacion.json"), "w"), indent=1)
    # forma cerrada en la garganta: con H(0) = 0, E_B,tot = 0 da S(0) = -2 E_B^J/(N0 r0^2) y e0_tot = -4 E_B^J/r0 + kcd
    th = {Phi[1]: 0, rr[1]: 0, PP[2]: 0, PP[0]: 0, Xs: -PP[1] ** 2 / 2}
    e0s = sp.simplify((-4 * feJ["B"] / rr[0] + kcd).xreplace(th))
    Xv, GG, GX, GXX = sp.symbols("X Gcal Gcal_X Gcal_XX", real=True)
    Gfun = sp.Function("Gf")
    Lc = -sp.Rational(3, 4) * pi2 * sp.sqrt(Aj * Bj) / (Cj ** 2 * Gfun(-PP[1] ** 2 / (2 * Bj)))
    EB = sp.diff(Lc, Bj).subs({Aj: Ee, Bj: 1, Cj: rr[0] ** 2})
    eps_, dS, cS = sp.symbols("eps_ dS cS")
    Lcd = -sp.Rational(3, 4) * pi2 * N * sp.sqrt(1 + 2 * eps_ * cS) / ((rr[0] ** 2 + 2 * rr[0] * eps_ * dS) ** 2
                                                                      * Gfun(-PP[1] ** 2 / (2 * (1 + 2 * eps_ * cS))))
    kcd_c = sp.diff(Lcd, eps_, 2).subs(eps_, 0).diff(cS).diff(dS) / 2
    e0c = (-4 * EB / rr[0] + kcd_c).doit()
    Xb_ = -PP[1] ** 2 / 2
    e0c = e0c.replace(lambda e: isinstance(e, sp.Subs), lambda e: e.doit())
    e0c = e0c.subs(sp.Derivative(Gfun(Xb_), Xb_), GX) if False else e0c
    rep = {}
    for d_ in e0c.atoms(sp.Derivative):
        rep[d_] = GX
    e0c = sp.simplify(e0c.xreplace(rep).xreplace({Gfun(Xb_): GG}))
    print("e0_tot(0) (forma cerrada) =", sp.factor(e0c), "  con X = -phi'^2/2, Gcal_X = dGcal/dX en la garganta")
    # comprobacion de la forma cerrada frente a la expresion de los G_ij (Gcal_X = G4_X - 2 G4_X - 2X G4_XX - 8X F4 - 4X^2 F4_X)
    Xg = -PP[1] ** 2 / 2
    Gc_ = gsym("G4", 0, 0) - 2 * Xg * gsym("G4", 0, 1) - 4 * Xg ** 2 * gsym("F4", 0, 0)
    GcX = -gsym("G4", 0, 1) - 2 * Xg * gsym("G4", 0, 2) - 8 * Xg * gsym("F4", 0, 0) - 4 * Xg ** 2 * gsym("F4", 0, 1)
    dif = sp.simplify(e0s - e0c.xreplace({GG: Gc_, GX: GcX}))
    print("forma cerrada = expresion completa:", dif == 0)
    out["e0_forma_cerrada"] = str(sp.factor(e0c))
    json.dump(out, open(os.path.join(AQUI, "resultados", "garganta_rotacion.json"), "w"), indent=1)
