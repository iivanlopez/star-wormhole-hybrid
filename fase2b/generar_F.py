# -*- coding: utf-8 -*-
"""
Ecuaciones de Einstein-escalar-Gauss-Bonnet con acoplamiento GENERAL F(phi) (Antoniou et al. 2019, arXiv:1904.13091):
    S = (1/16 pi) int sqrt(-g) [ R - (dphi)^2/2 + F(phi) G ]          (F incluye la constante de acoplamiento)
en el gauge de longitud propia con shift nulo.  Se usa el campo auxiliar chi = F'(phi) (egb_simbolico, argumento Fp):
las ecuaciones solo contienen chi y sus primeras derivadas, que luego se sustituyen por
    chi = F1 + eps F2 q,   chi_l = F2 p' + eps (F3 p' q + F2 q'),   chi_t = eps F2 q_t,
con F1, F2, F3 = F', F'', F''' evaluadas en el fondo.  Genera egb_F_gen.py:
    est_principal(F1, r, r1, P1, p1)            matriz de (r'', Phi'', phi'') en (E^t_t, E^th_th, E_phi)
    est_resto(F1, F2, r, r1, P1, p1)            resto de esas tres ecuaciones
    est_ligadura(F1, r, r1, P1, p1)             E^l_l (sin derivadas segundas)
    lin_principal / lin_resto / lin_ligaduras   lo mismo para las perturbaciones (a, d, q) ~ e^{i omega t}
Uso: python generar_F.py
"""
import sys
import time
import sympy as sp
from egb_simbolico import ecuaciones

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
t0 = time.time()
t, l, th, ph, eps = sp.symbols("t l theta varphi epsilon", real=True)
w = sp.Symbol("w", real=True)
X = [t, l, th, ph]
F1, F2, F3 = sp.symbols("F1 F2 F3", real=True)
r_, r1, r2, P1, P2, p_, p1, p2, e2 = sp.symbols("r r1 r2 P1 P2 p p1 p2 e2", real=True)
A0, A1, A2, D0, D1, D2, Q0, Q1, Q2 = sp.symbols("a0 a1 a2 d0 d1 d2 q0 q1 q2")
At, Dt, Qt, Alt, Dlt, Qlt, Att, Dtt, Qtt = sp.symbols("At Dt Qt Alt Dlt Qlt Att Dtt Qtt")
Phi, r, p = (sp.Function(n)(l) for n in ("Phi", "r", "p"))


def chi_sub(expr, chi, lineal):
    """Sustituye chi y sus derivadas (primeras) por F1, F2, F3 y las del escalar."""
    if lineal:
        q = sp.Function("q")(t, l)
        rep = {chi.diff(l): F2 * p.diff(l) + eps * (F3 * p.diff(l) * q + F2 * q.diff(l)),
               chi.diff(t): eps * F2 * q.diff(t), chi: F1 + eps * F2 * q}
    else:
        rep = {chi.diff(l): F2 * p.diff(l), chi: F1}
    for k in list(expr.atoms(sp.Derivative)):
        if k.expr == chi and sum(c for _, c in k.variable_count) > 1:
            raise RuntimeError("aparecen derivadas segundas de chi")
    return expr.subs(rep)


def plano(expr):
    a, d, q = (sp.Function(n)(t, l) for n in ("a", "d", "q"))
    rep = {}
    for f, s0, s1, s2, st, slt, stt in ((a, A0, A1, A2, At, Alt, Att), (d, D0, D1, D2, Dt, Dlt, Dtt),
                                         (q, Q0, Q1, Q2, Qt, Qlt, Qtt)):
        rep.update({f.diff(t, 2): stt, f.diff(t, l): slt, f.diff(l, t): slt, f.diff(l, 2): s2, f.diff(t): st,
                    f.diff(l): s1})
    rep.update({r.diff(l, 2): r2, Phi.diff(l, 2): P2, p.diff(l, 2): p2, r.diff(l): r1, Phi.diff(l): P1, p.diff(l): p1})
    e = expr.subs(rep).subs({a: A0, d: D0, q: Q0, r: r_, p: p_})
    e = sp.expand(e, power_exp=True).subs(sp.exp(-2 * Phi), e2).subs(sp.exp(2 * Phi), 1 / e2)
    e = e.subs(Phi, -sp.log(e2) / 2)
    return sp.expand(sp.powsimp(sp.expand(e, power_exp=True)).subs(sp.sin(th), 1))


# ----------------------------------------------------------------- estatico
chi_s = sp.Function("chi")(l)
E, Ephi, GB, (g, gi, Gam) = ecuaciones(Phi, sp.Integer(0), r, p, X, sp.Integer(1), None, Fp=chi_s)
est = {"tt": gi[0, 0] * E[0, 0], "thth": gi[2, 2] * E[2, 2], "phi": Ephi, "ll": gi[1, 1] * E[1, 1]}
est = {k: plano(chi_sub(v, chi_s, False)) for k, v in est.items()}
print(f"estatico: {time.time() - t0:.0f} s", flush=True)
seg_s = [r2, P2, p2]
Ms = sp.Matrix([[sp.simplify(est[k].coeff(s)) for s in seg_s] for k in ("tt", "thth", "phi")])
Rs = [sp.expand(est[k] - sum(Ms[i, j] * seg_s[j] for j in range(3))) for i, k in enumerate(("tt", "thth", "phi"))]
assert not any(x.has(s) for x in Rs for s in seg_s)
assert not any(est["ll"].has(s) for s in seg_s)
permit_s = {F1, F2, r_, r1, P1, p1}
for k, v in est.items():
    assert v.free_symbols <= permit_s | {r2, P2, p2}, (k, v.free_symbols - permit_s)

# ----------------------------------------------------------------- lineal
a, d, q = (sp.Function(n)(t, l) for n in ("a", "d", "q"))
chi = sp.Function("chi")(t, l)
E, Ephi, GB, (g, gi, Gam) = ecuaciones(Phi + eps * a, sp.Integer(0), r + eps * d, p + eps * q, X, sp.Integer(1), None,
                                       Fp=chi)
mix = {"tt": gi[0, 0] * E[0, 0], "ll": gi[1, 1] * E[1, 1], "thth": gi[2, 2] * E[2, 2], "lt": gi[1, 1] * E[1, 0],
       "phi": Ephi}
lin = {k: sp.diff(chi_sub(v, chi, True), eps).subs(eps, 0) for k, v in mix.items()}
planos = {k: plano(v) for k, v in lin.items()}
print(f"lineal: {time.time() - t0:.0f} s", flush=True)
par = {At: 0, Dt: 0, Qt: 0, Alt: 0, Dlt: 0, Qlt: 0, Att: -w * A0, Dtt: -w * D0, Qtt: -w * Q0}
imp = {At: A0, Dt: D0, Qt: Q0, Alt: A1, Dlt: D1, Qlt: Q1, Att: 0, Dtt: 0, Qtt: 0}
for k in ("tt", "thth", "phi", "ll"):
    assert sp.expand(planos[k] - planos[k].subs({At: 0, Dt: 0, Qt: 0, Alt: 0, Dlt: 0, Qlt: 0})) == 0, k
assert sp.expand(planos["lt"].subs({s: 0 for s in (At, Dt, Qt, Alt, Dlt, Qlt, Att, Dtt, Qtt)})) == 0
assert sp.expand(planos["lt"] - planos["lt"].subs({Att: 0, Dtt: 0, Qtt: 0})) == 0
ec = {k: sp.expand(planos[k].subs(par)) for k in ("tt", "thth", "phi", "ll")}
ec["lt"] = sp.expand(planos["lt"].subs(imp))
permit = {F1, F2, F3, r_, r1, r2, P1, P2, p1, p2, e2, A0, A1, D0, D1, Q0, Q1, w}
for k, v in ec.items():
    assert v.free_symbols <= permit | {A2, D2, Q2}, (k, v.free_symbols - permit)
    assert not v.has(A0), k
seg = [D2, A2, Q2]
Mp = sp.Matrix([[sp.simplify(ec[k].coeff(s)) for s in seg] for k in ("tt", "thth", "phi")])
resto = [sp.expand(ec[k] - sum(Mp[i, j] * seg[j] for j in range(3))) for i, k in enumerate(("tt", "thth", "phi"))]
assert not any(x.has(s) for x in resto for s in seg)
assert not any(ec[k].has(s) for k in ("ll", "lt") for s in seg)
assert sp.simplify(Mp - Ms) == sp.zeros(3, 3)            # misma parte principal que el fondo
print(f"comprobaciones de estructura: OK ({time.time() - t0:.0f} s)", flush=True)

bg_s = "F1, r, r1, P1, p1"
cuerpo = ['# -*- coding: utf-8 -*-', '"""Generado por generar_F.py: no editar a mano."""', "import math", "",
          f"def est_principal({bg_s}):",
          "    return [" + ", ".join("[" + ", ".join(sp.pycode(Ms[i, j]) for j in range(3)) + "]" for i in range(3)) + "]", "",
          f"def est_resto(F1, F2, r, r1, P1, p1):",
          "    return [" + ", ".join(sp.pycode(x) for x in Rs) + "]", "",
          f"def est_ligadura({bg_s}):",
          f"    return {sp.pycode(est['ll'])}", "",
          "def lin_principal(F1, r, r1, P1, p1):",
          "    return est_principal(F1, r, r1, P1, p1)", "",
          "def lin_resto(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1):",
          "    return [" + ", ".join(sp.pycode(x) for x in resto) + "]", "",
          "def lin_ligaduras(F1, F2, F3, w, r, r1, r2, P1, P2, p1, p2, e2, a1, d0, d1, q0, q1):",
          f"    return [{sp.pycode(ec['lt'])}, {sp.pycode(ec['ll'])}]", ""]
open("egb_F_gen.py", "w", encoding="utf-8").write("\n".join(cuerpo))
print(f"escrito egb_F_gen.py ({time.time() - t0:.0f} s)")
