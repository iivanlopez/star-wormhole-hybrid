# -*- coding: utf-8 -*-
"""
Linealiza las ecuaciones de EdGB (egb_simbolico.py, verificadas en verificar_egb.py) en el gauge de longitud propia
con shift nulo:
    A = Phi(l) + eps a(t,l),   B = 0,   R = r(l) + eps d(t,l),   phi = p(l) + eps q(t,l),    (a, d, q) ~ e^{i w^{1/2} t}
y genera egb_lineal_gen.py con:
    principal(bg)          -> matriz 3x3 de los coeficientes de (d'', a'', q'') en (E^t_t, E^th_th, E_phi)
    resto(bg, pert, w)     -> terminos sin derivadas segundas de la perturbacion en esas tres ecuaciones
    ligaduras(bg, pert, w) -> (E^l_t / (i omega), E^l_l) lineales
bg = (r, r1, r2, P1, P2, p, p1, p2, e2) con r1 = r', P1 = Phi', p1 = phi', e2 = e^{-2 Phi};
pert = (a, a1, d, d1, q, q1);  w = omega^2 (medido con el tiempo coordenado t).
Uso: python generar_lineal.py
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
al, ga, w = sp.symbols("alpha gamma w", real=True)
X = [t, l, th, ph]
Phi, r, p = (sp.Function(n)(l) for n in ("Phi", "r", "p"))
a, d, q = (sp.Function(n)(t, l) for n in ("a", "d", "q"))
E, Ephi, GB, (g, gi, Gam) = ecuaciones(Phi + eps * a, sp.Integer(0), r + eps * d, p + eps * q, X, al, ga)
print(f"tensores: {time.time() - t0:.0f} s", flush=True)
mix = {"tt": gi[0, 0] * E[0, 0], "ll": gi[1, 1] * E[1, 1], "thth": gi[2, 2] * E[2, 2], "lt": gi[1, 1] * E[1, 0],
       "phi": Ephi}
lin = {k: sp.diff(v, eps).subs(eps, 0) for k, v in mix.items()}
print(f"linealizado: {time.time() - t0:.0f} s", flush=True)

# simbolos planos
r_, r1, r2, P1, P2, p_, p1, p2, e2 = sp.symbols("r r1 r2 P1 P2 p p1 p2 e2", real=True)
A0, A1, A2, D0, D1, D2, Q0, Q1, Q2 = sp.symbols("a0 a1 a2 d0 d1 d2 q0 q1 q2")
At, Dt, Qt, Alt, Dlt, Qlt, Att, Dtt, Qtt = sp.symbols("At Dt Qt Alt Dlt Qlt Att Dtt Qtt")


def plano(expr):
    rep = {}
    for f, sym0, sym1, sym2, symt, symlt, symtt in ((a, A0, A1, A2, At, Alt, Att), (d, D0, D1, D2, Dt, Dlt, Dtt),
                                                     (q, Q0, Q1, Q2, Qt, Qlt, Qtt)):
        rep.update({f.diff(t, 2): symtt, f.diff(t, l): symlt, f.diff(l, t): symlt, f.diff(l, 2): sym2,
                    f.diff(t): symt, f.diff(l): sym1})
    rep.update({r.diff(l, 2): r2, Phi.diff(l, 2): P2, p.diff(l, 2): p2, r.diff(l): r1, Phi.diff(l): P1, p.diff(l): p1})
    e = expr.subs(rep)
    e = e.subs({a: A0, d: D0, q: Q0, r: r_, p: p_})
    e = sp.expand(e, power_exp=True).subs(sp.exp(-2 * Phi), e2).subs(sp.exp(2 * Phi), 1 / e2)
    e = e.subs(Phi, -sp.log(e2) / 2)          # Phi solo aparece en exponenciales: e^{-2 Phi} = e2
    return sp.expand(sp.powsimp(sp.expand(e, power_exp=True)))


planos = {k: plano(v) for k, v in lin.items()}
print(f"aplanado: {time.time() - t0:.0f} s", flush=True)
for k, v in planos.items():
    extra = v.free_symbols - {r_, r1, r2, P1, P2, p_, p1, p2, e2, A0, A1, A2, D0, D1, D2, Q0, Q1, Q2, At, Dt, Qt, Alt, Dlt,
                              Qlt, Att, Dtt, Qtt, al, ga, th}
    print(f"  {k}: simbolos inesperados = {extra}", flush=True)
# dependencia temporal e^{i omega t}: d_t -> i omega (con omega^2 = w); se trabaja con amplitudes reales
#   ecuaciones pares en d_t (tt, ll, thth, phi): d_tt -> -w;   lt: lineal en d_t -> se divide por i omega
par = {At: 0, Dt: 0, Qt: 0, Alt: 0, Dlt: 0, Qlt: 0, Att: -w * A0, Dtt: -w * D0, Qtt: -w * Q0}
imp = {At: A0, Dt: D0, Qt: Q0, Alt: A1, Dlt: D1, Qlt: Q1, Att: 0, Dtt: 0, Qtt: 0}
cero_t = {s: 0 for s in (At, Dt, Qt, Alt, Dlt, Qlt, Att, Dtt, Qtt)}
impares = {s: 0 for s in (Att, Dtt, Qtt)}
for k in ("tt", "thth", "phi", "ll"):          # ecuaciones pares en el tiempo: sin derivadas temporales impares
    assert sp.expand(planos[k] - planos[k].subs({At: 0, Dt: 0, Qt: 0, Alt: 0, Dlt: 0, Qlt: 0})) == 0, k
assert sp.expand(planos["lt"].subs(cero_t)) == 0                       # E^l_t: todo termino lleva una d/dt
assert sp.expand(planos["lt"] - planos["lt"].subs(impares)) == 0       # ... y solo una
print("  paridad temporal de las ecuaciones: OK", flush=True)
ec = {k: sp.expand(planos[k].subs(par)) for k in ("tt", "thth", "phi", "ll")}
ec["lt"] = sp.expand(planos["lt"].subs(imp))
for k in ec:
    ec[k] = sp.expand(ec[k].subs(sp.sin(th), 1))
    resid = ec[k].free_symbols - {r_, r1, r2, P1, P2, p_, p1, p2, e2, A0, A1, A2, D0, D1, D2, Q0, Q1, Q2, al, ga, w}
    assert not resid, (k, resid)
# comprobaciones de estructura
for k in ("ll", "lt"):
    print(f"  {k}: contiene derivadas segundas? {any(ec[k].has(s) for s in (A2, D2, Q2))}", flush=True)
print(f"  'a' sin derivar aparece en: {[k for k in ec if ec[k].has(A0)]}", flush=True)
seg = [D2, A2, Q2]
Mp = sp.Matrix([[sp.simplify(ec[k].coeff(s)) for s in seg] for k in ("tt", "thth", "phi")])
resto = [sp.expand(ec[k] - sum(Mp[i, j] * seg[j] for j in range(3))) for i, k in enumerate(("tt", "thth", "phi"))]
for i, k in enumerate(("tt", "thth", "phi")):
    assert not any(resto[i].has(s) for s in seg), k
args_bg = "r, r1, r2, P1, P2, p, p1, p2, e2"
args_pt = "a0, a1, d0, d1, q0, q1"
cuerpo = ['# -*- coding: utf-8 -*-', '"""Generado por generar_lineal.py: no editar a mano."""', "import math", "",
          f"def principal(alpha, gamma, {args_bg}):",
          "    return [" + ", ".join("[" + ", ".join(sp.pycode(Mp[i, j]) for j in range(3)) + "]" for i in range(3)) + "]", "",
          f"def resto(alpha, gamma, w, {args_bg}, {args_pt}):",
          "    return [" + ", ".join(sp.pycode(x) for x in resto) + "]", "",
          f"def ligaduras(alpha, gamma, w, {args_bg}, {args_pt}):",
          f"    return [{sp.pycode(ec['lt'])}, {sp.pycode(ec['ll'])}]", ""]
open("egb_lineal_gen.py", "w", encoding="utf-8").write("\n".join(cuerpo))
print(f"escrito egb_lineal_gen.py ({time.time() - t0:.0f} s)", flush=True)
print("matriz principal:", Mp, flush=True)
