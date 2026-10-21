# -*- coding: utf-8 -*-
"""
Verificacion de egb_simbolico.py:
  [1] Invariante de Gauss-Bonnet de Schwarzschild: G = 48 M^2 / r^6.
  [2] Identidad de Bianchi generalizada  nabla^m E_mn + (1/2) E_phi d_n phi = 0, fuera de la capa de masas.
      Se evalua sobre "jets" aleatorios (valores arbitrarios de A, B, R, phi y de todas sus derivadas en un
      punto): la identidad debe cumplirse para cualquier jet.  Comprueba signos y factores del tensor de GB.
  [3] Caso estatico en las coordenadas de KKK (ds^2 = -e^{2 nu} dt^2 + f dl^2 + (l^2 + r0^2) dOmega^2):
      las ecuaciones ll, tt y del dilaton coinciden con las ec. (4) de Kanti-Kleihaus-Kunz (PRD 85, 044007),
      comparadas sobre jets aleatorios.
Uso: python verificar_egb.py
"""
import sys
import random
import sympy as sp
from egb_simbolico import ecuaciones, divergencia

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
t, l, th, ph = sp.symbols("t l theta varphi", real=True)
X = [t, l, th, ph]
al, ga, M, r0 = sp.symbols("alpha gamma M r_0", positive=True)
OK = []
random.seed(12345)


def informe(nombre, ok, extra=""):
    ok = bool(ok)
    OK.append(ok)
    print(f"[{'OK   ' if ok else 'FALLO'}] {nombre} {extra}", flush=True)


def jet(expr, extra=None):
    """Sustituye todas las derivadas y funciones por valores aleatorios (de mayor a menor orden)."""
    sub = dict(extra or {})
    atoms = sorted(expr.atoms(sp.Derivative), key=lambda d: -sum(c for _, c in d.variable_count))
    for d in atoms:
        sub.setdefault(d, sp.Rational(random.randint(-900, 900), 1000))
    fun = [f for f in expr.atoms(sp.Function) if not isinstance(f, sp.Derivative) and f.func not in (sp.sin, sp.cos, sp.exp, sp.log, sp.sqrt)]
    for f in fun:
        sub.setdefault(f, sp.Rational(random.randint(200, 900), 1000))
    return sub


# [1] Schwarzschild
A = sp.log(1 - 2 * M / l) / 2
E, Ephi, GB, _ = ecuaciones(A, -A, l, sp.Integer(0), X, al, ga)
informe("Gauss-Bonnet de Schwarzschild = 48 M^2/r^6", sp.simplify(GB - 48 * M**2 / l**6) == 0)

# [2] Bianchi sobre jets aleatorios (funciones genericas de t y l)
A, B, R, phi = (sp.Function(n)(t, l) for n in ("A", "B", "R", "phi"))
E, Ephi, GB, (g, gi, Gam) = ecuaciones(A, B, R, phi, X, al, ga)
for nu, nom in ((0, "t"), (1, "l")):
    div = divergencia(E, g, gi, Gam, X, nu)
    res = div + Ephi * sp.diff(phi, X[nu]) / 2
    for prueba in range(2):
        sub = jet(res, {al: sp.Rational(7, 10) + prueba, ga: 1 + sp.Rational(prueba, 3), th: sp.Rational(11, 10)})
        val, esc = sp.N(res.subs(sub), 25), sp.N(div.subs(sub), 25)
        informe(f"Bianchi ({nom}), jet {prueba + 1}", abs(val) < 1e-18 * max(1, abs(esc)),
                f"(residuo {float(val):.1e}; |nabla^m E_m{nom}| = {float(abs(esc)):.2e})")

# [3] Estatico, coordenadas de KKK
nu, f, p = sp.Function("nu")(l), sp.Function("f")(l), sp.Function("p")(l)
r2 = l**2 + r0**2
E, Ephi, GB, _ = ecuaciones(nu, sp.log(f) / 2, sp.sqrt(r2), p, X, al, ga)
F = sp.exp(-ga * p)
nup, fp, pp, ppp, nupp = nu.diff(l), f.diff(l), p.diff(l), p.diff(l, 2), nu.diff(l, 2)
kkk = {
    "ll  (4b, nu')": (nup, (r2 * f - l**2) / (2 * l * r2) + pp**2 * r2 / (8 * l)
                      + 2 * al * ga * F / (l * r2 * f) * nup * pp * (r2 * f - 3 * l**2), E[1, 1]),
    "tt  (4a, f')": (fp, -f * (r2 * f + l**2 - 2 * r2) / (l * r2) + r2 * f * pp**2 / (4 * l)
                     + 2 * al * ga * F / (l * r2) * (2 * (r2 * f - l**2) * (ga * pp**2 - ppp)
                                                    + pp * (fp / f * (r2 * f - 3 * l**2) + 4 * l * r0**2 / r2)), E[0, 0]),
    "dil (4c, phi'')": (ppp, -nup * pp - pp * (4 * l * f - r2 * fp) / (2 * r2 * f)
                        + 4 * al * ga * F / (r2**2 * f) * (-2 * (r2 * f - l**2) * (nup**2 + nupp)
                                                          + nup * (fp / f * (r2 * f - 3 * l**2) + 4 * l * r0**2 / r2)), Ephi),
}
for nombre, (var, rhs_kkk, nuestra) in kkk.items():
    # ambas ecuaciones son lineales en 'var': se despeja var de la nuestra y de la de KKK y se comparan sobre jets
    sol = sp.solve(sp.Eq(nuestra, 0), var)
    ref = sp.solve(sp.Eq(var, rhs_kkk), var)
    ok = len(sol) == 1 and len(ref) == 1
    peor = 0.0
    for prueba in range(3):
        # jets sobre la union de atomos (si la diferencia se simplifica a 0 exacto, no tiene atomos)
        sub = jet(sol[0] + ref[0], {al: sp.Rational(1, 20) + sp.Rational(prueba, 50), ga: 1, r0: 1,
                                    l: sp.Rational(3 + prueba, 10)}) if ok else {}
        if ok:
            d = float(abs(sp.N((sol[0] - ref[0]).subs(sub), 25)))
            s = float(abs(sp.N(ref[0].subs(sub), 25)))
            peor = max(peor, d / max(s, 1e-30))
    informe(f"estatico KKK, ec. {nombre}", ok and peor < 1e-15, f"(diferencia relativa max {peor:.1e})")
print(f"\n{sum(OK)}/{len(OK)} comprobaciones superadas", flush=True)
