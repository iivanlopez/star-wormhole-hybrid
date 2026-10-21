# -*- coding: utf-8 -*-
"""
Verificacion de las ecuaciones EsGB + potencial (potencial.py):
  1. Identidad de Bianchi generalizada con V (SymPy, jets aleatorios, egb_simbolico de fase2b):
       E_mn -> E_mn + g_mn V/2,  E_phi -> E_phi - V'(phi),   nabla^m E_mn + (1/2) E_phi d_n phi = 0.
     Control: con el signo de V' cambiado la identidad falla (la prueba es sensible al signo relativo).
  2. Las ecuaciones estaticas de potencial.py frente a las covariantes con V (acoplamiento exponencial y potencial
     cuartico), evaluadas en estados aleatorios: (r'', Phi'', phi'') y ligadura ll.
  3. Schwarzschild-de Sitter (V = 2 Lambda constante, phi constante): r' = sqrt(1 - 2M/r - Lambda r^2/3) y N.
  4. Propagacion de la ligadura ll del fondo con escalar y potencial no triviales (Bianchi con V).
  5. Perturbaciones: el modo gauge g = (Phi'' - w e^{-2 Phi}, r', r'', phi', phi'') es solucion exacta de
     derivadas_V y cumple ligaduras_V; una solucion lineal generica conserva las ligaduras; la derivada de la familia
     de fondos respecto de un dato inicial cumple el sistema lineal con w = 0.
Salida: resultados/verificar_potencial.json.   Uso: python verificar_potencial.py
"""
import os
import sys
import json
import math
import random
from types import SimpleNamespace
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import potencial as P                                            # noqa: E402
from egb_simbolico import ecuaciones, divergencia                # noqa: E402  (fase2b, via potencial)

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
out = {}
random.seed(2026)
rng = np.random.default_rng(7)

# ------------------------------------------------------------------ 1. Bianchi con V
t, l, th, ph = sp.symbols("t l theta varphi", real=True)
X = [t, l, th, ph]
al, ga, m2, lam = sp.symbols("alpha gamma m2 lambda", positive=True)


def jet(expr, extra=None):
    sub = dict(extra or {})
    atoms = sorted(expr.atoms(sp.Derivative), key=lambda d: -sum(c for _, c in d.variable_count))
    for d in atoms:
        sub.setdefault(d, sp.Rational(random.randint(-900, 900), 1000))
    fun = [f for f in expr.atoms(sp.Function) if not isinstance(f, sp.Derivative)
           and f.func not in (sp.sin, sp.cos, sp.exp, sp.log, sp.sqrt)]
    for f in fun:
        sub.setdefault(f, sp.Rational(random.randint(200, 900), 1000))
    return sub


A, B, R, phi = (sp.Function(n)(t, l) for n in ("A", "B", "R", "phi"))
E, Ephi, GB, (g, gi, Gam) = ecuaciones(A, B, R, phi, X, al, ga)
Vs = m2 * phi**2 / 2 + lam * phi**4
EV = E + g * Vs / 2
bianchi = []
RAPIDO = len(sys.argv) > 1 and sys.argv[1] == "rapido"          # sin la identidad de Bianchi (~9 min)
for nu, nom in (() if RAPIDO else ((0, "t"), (1, "l"))):
    div = divergencia(EV, g, gi, Gam, X, nu)
    for signo, nombre in ((-1, "correcto (-V')"), (+1, "control (+V')")):
        res = div + (Ephi + signo * sp.diff(Vs, phi)) * sp.diff(phi, X[nu]) / 2
        for prueba in range(2):
            sub = jet(res, {al: sp.Rational(7, 10) + prueba, ga: 1, m2: sp.Rational(13, 10),
                            lam: sp.Rational(3, 10) + prueba, th: sp.Rational(11, 10)})
            val, esc = float(abs(sp.N(res.subs(sub), 30))), float(abs(sp.N(div.subs(sub), 30)))
            bianchi.append(dict(signo=nombre, componente=nom, jet=prueba + 1, residuo=val, escala=esc))
            print(f"1. Bianchi con V, {nombre}, componente {nom}, jet {prueba + 1}: residuo {val:.1e} "
                  f"(|nabla^m E_m{nom}| = {esc:.2e})", flush=True)
out["bianchi"] = bianchi

# ------------------------------------------------------------------ 2. estatico: covariante con V frente al modulo
Phi_s, r_s, p_s = (sp.Function(n)(l) for n in ("Phi", "r", "p"))
Es, Ephis, _, (gs, gis, _) = ecuaciones(Phi_s, sp.Integer(0), r_s, p_s, X, al, ga)
Vp = m2 * p_s**2 / 2 + lam * p_s**4
mix = [gis[0, 0] * Es[0, 0] + Vp / 2, gis[2, 2] * Es[2, 2] + Vp / 2, Ephis - sp.diff(Vp, p_s)]
lig = gis[1, 1] * Es[1, 1] + Vp / 2
segs = [r_s.diff(l, 2), Phi_s.diff(l, 2), p_s.diff(l, 2)]
r_, r1_, P1_, p_, p1_ = sp.symbols("r r1 P1 p p1")
rep = {r_s.diff(l): r1_, Phi_s.diff(l): P1_, p_s.diff(l): p1_, r_s: r_, p_s: p_, sp.sin(th): 1}
x1, x2, x3 = sp.symbols("x1 x2 x3")
mix_n = [e.subs(dict(zip(segs, (x1, x2, x3)))).subs(rep).subs(Phi_s, 0) for e in mix]
# la ligadura covariante trae derivadas segundas que solo se cancelan al simplificar: se sustituyen tambien y se
# comprueba abajo que el resultado no depende de ellas
lig_n = lig.subs(dict(zip(segs, (x1, x2, x3)))).subs(rep).subs(Phi_s, 0)
assert not any(e.has(sp.Derivative) or e.has(Phi_s) for e in mix_n + [lig_n])
f_mix = sp.lambdify((x1, x2, x3, r_, r1_, P1_, p_, p1_, al, ga, m2, lam), mix_n)
f_lig = sp.lambdify((x1, x2, x3, r_, r1_, P1_, p_, p1_, al, ga, m2, lam), lig_n)
err2 = 0.0
for _ in range(20):
    rr, r1, P1, pp, p1 = rng.uniform(0.5, 3.0), rng.uniform(-0.9, 0.9), rng.uniform(-2, 2), rng.uniform(-1, 1), \
        rng.uniform(-2, 2)
    alpha, mm, ll_ = rng.uniform(-0.5, 0.5), rng.uniform(0.1, 3), rng.uniform(0, 1)
    Fd = P.acoplamiento("exp", alpha)
    Vd = lambda x, mm=mm, ll_=ll_: (mm * x * x / 2 + ll_ * x**4, mm * x + 4 * ll_ * x**3,  # noqa: E731
                                   mm + 12 * ll_ * x * x, 24 * ll_ * x)
    y = np.array([rr, r1, 0.0, P1, pp, p1])
    xs, _ = P.segundas_V(y, Fd, Vd)
    res = np.array(f_mix(*xs, rr, r1, P1, pp, p1, alpha, 1.0, mm, ll_), dtype=float)
    lg = abs(f_lig(*xs, rr, r1, P1, pp, p1, alpha, 1.0, mm, ll_) - P.est_ligadura_V(y, Fd, Vd))
    lg_x = abs(f_lig(*rng.normal(size=3), rr, r1, P1, pp, p1, alpha, 1.0, mm, ll_) - P.est_ligadura_V(y, Fd, Vd))
    err2 = max(err2, float(np.max(np.abs(res))), float(lg), float(lg_x))
print(f"2. estatico covariante con V frente a potencial.py (20 estados): residuo max {err2:.1e}", flush=True)
out["estatico_covariante"] = err2

# ------------------------------------------------------------------ 3. Schwarzschild-de Sitter
Mm, Lam = 1.0, 0.01
Fd0 = P.acoplamiento("cuadratico", 0.3)          # phi = 0 constante y F'(0) = 0: el acoplamiento no actua
Vd0 = P.potencial("constante", 2 * Lam)
r_i = 4.0


def fsds(r):
    return 1 - 2 * Mm / r - Lam * r * r / 3


y0 = [r_i, math.sqrt(fsds(r_i)), 0.0, (Mm / r_i**2 - Lam * r_i / 3) / math.sqrt(fsds(r_i)), 0.0, 0.0]


def rhs0(x, y):
    s, _ = P.segundas_V(y, Fd0, Vd0)
    return [y[1], s[0], y[3], s[1], y[5], s[2]]


s0 = solve_ivp(rhs0, (0, 6), y0, method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
err3 = 0.0
lig3 = 0.0
for x in np.linspace(0, 6, 40):
    r, rp, Phi, Pp, pp, p1 = s0.sol(x)
    err3 = max(err3, abs(rp - math.sqrt(fsds(r))), abs(math.exp(Phi) - math.sqrt(fsds(r) / fsds(r_i))))
    lig3 = max(lig3, abs(P.est_ligadura_V(s0.sol(x), Fd0, Vd0)))
print(f"3. Schwarzschild-de Sitter (M = 1, Lambda = 0.01, V = 2 Lambda): max |r' - sqrt(f)|, |N - N_SdS| = "
      f"{err3:.1e}; ligadura {lig3:.1e}", flush=True)
out["schwarzschild_de_sitter"] = dict(err=err3, ligadura=lig3)

# ------------------------------------------------------------------ 4. propagacion de la ligadura del fondo


def dato_consistente(Fd, Vd, r, r1, Phi, p, p1):
    """Phi' de la ligadura ll (lineal en Phi')."""
    c0 = P.est_ligadura_V(np.array([r, r1, Phi, 0.0, p, p1]), Fd, Vd)
    c1 = P.est_ligadura_V(np.array([r, r1, Phi, 1.0, p, p1]), Fd, Vd) - c0
    return np.array([r, r1, Phi, -c0 / c1, p, p1])


def escala_ligadura(y, Fd, Vd):
    """Suma de los valores absolutos de los terminos de E^l_l (para el residuo relativo)."""
    r, s, Phi, q, p, u = y[:6]
    F1 = Fd(p)[1]
    return (abs((s * s - 1) / r**2) + abs(2 * q * s / r) + u * u / 4 + abs(4 * F1 * q * u * (1 - 3 * s * s) / r**2)
            + abs(Vd(p)[0] / 2))


prop4 = []
for tipoF, alpha, tipoV, amp, kw, dato in (
        ("cuadratico", -0.2, "masivo", 0.8, {}, (1.3, 0.35, 0.0, 0.4, -0.6)),
        ("exp", 0.15, "doble_pozo", 0.3, {"v": 1.2}, (2.0, 0.6, 0.0, 0.4, -0.3)),
        ("lineal", 0.3, "cuartico", 0.5, {"g": 2.0}, (2.0, 0.7, 0.0, 0.2, 0.2))):
    Fd = P.acoplamiento(tipoF, alpha)
    Vd = P.potencial(tipoV, amp, **kw)
    y0 = dato_consistente(Fd, Vd, *dato)

    def rhs4(x, y, Fd=Fd, Vd=Vd):
        s, _ = P.segundas_V(y, Fd, Vd)
        return [y[1], s[0], y[3], s[1], y[5], s[2]]
    s4 = solve_ivp(rhs4, (0, 3), y0, method="DOP853", rtol=1e-12, atol=1e-14, dense_output=True)
    assert s4.status == 0
    lg = max(abs(P.est_ligadura_V(s4.sol(x), Fd, Vd)) / escala_ligadura(s4.sol(x), Fd, Vd)
             for x in np.linspace(0, s4.t[-1], 50))
    prop4.append(dict(F=tipoF, alpha=alpha, V=tipoV, amp=amp, l_fin=float(s4.t[-1]), ligadura_rel=float(lg)))
    print(f"4. ligadura ll propagada, F = {tipoF} ({alpha}), V = {tipoV} ({amp}): max |E^l_l|/escala = {lg:.1e} "
          f"en [0, {s4.t[-1]:.2f}]", flush=True)
out["propagacion_fondo"] = prop4

# ------------------------------------------------------------------ 5. perturbaciones
# fondo de prueba: una mitad lisa de clase A (garganta en la ladera de un doble pozo) integrada en l > 0
Fd = P.acoplamiento("cuadratico", 0.06)
Ud = P.forma("doble_pozo", v=1.0)
amp = P.amp_garganta(Ud, Fd, 0.5, 0.0, 0.0)
Vd = P.potencial("doble_pozo", amp, v=1.0)
lado = P.LadoV(Fd, Vd, 0.5, 0.0, 0.0, phi_v=1.0, l_max=3.0, rtol=1e-13).resolver()
print(f"5. fondo de prueba: clase A, F = 0.06 phi^2, doble pozo (v = 1, amp = {amp:.4f}), phi0 = 0.5; "
      f"r''(0) = {lado.rpp0:.4f}, Phi''(0) = {lado.Ppp0:.4f}, phi''(0) = {lado.ppp0:.4f}, final '{lado.fin}' en "
      f"l = {lado.l_end:.3f}, ligadura en la garganta {lado.lig0:.1e}", flush=True)
fb = P.DobleV(lado, lado)


def g_en(x, w):
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = P.fondo_en_V(fb, x)
    return np.array([P2 - w * e2, r1, r2, p1, p2], dtype=complex if isinstance(w, complex) else float)


gauge = []
l_hi = min(2.5, 0.9 * lado.l_end)
for w in (-2.5, 7.0, complex(-3.0, 1.5)):
    for l1, l2 in ((-l_hi, -0.2), (-0.5, 0.7), (0.3, l_hi)):
        tipo = complex if isinstance(w, complex) else float
        ss = solve_ivp(lambda x, yy: P.derivadas_V(fb, x, yy, w), (l1, l2), g_en(l1, w).astype(tipo),
                       method="DOP853", rtol=1e-12, atol=1e-14)
        err = np.max(np.abs(ss.y[:, -1] - g_en(l2, w))) / np.max(np.abs(g_en(l2, w)))
        lg = max(np.max(np.abs(P.ligaduras_V(fb, x, g_en(x, w), w))) for x in np.linspace(l1, l2, 7))
        gauge.append(dict(w=str(w), l1=l1, l2=l2, err=float(err), lig=float(lg)))
        print(f"5a. modo gauge, w = {w}, l {l1:.2f} -> {l2:.2f}: |y - g|/|g| = {err:.1e}, ligaduras <= {lg:.1e}",
              flush=True)
out["gauge"] = gauge

prop = []
for w in (-4.0, complex(-1.0, 2.0)):
    tipo = complex if isinstance(w, complex) else float
    Cm = np.array([[P.ligaduras_V(fb, 0.4, e, w)[i] for e in np.eye(5)] for i in range(2)], dtype=tipo)
    _, _, Vh = np.linalg.svd(Cm)
    y0 = (Vh[2:].conj().T @ rng.normal(size=3)).astype(tipo)
    ss = solve_ivp(lambda x, yy: P.derivadas_V(fb, x, yy, w), (0.4, l_hi), y0, method="DOP853", rtol=1e-12,
                   atol=1e-14, dense_output=True)
    lg = max(np.max(np.abs(P.ligaduras_V(fb, x, ss.sol(x), w))) / np.max(np.abs(ss.sol(x)))
             for x in np.linspace(0.4, l_hi, 30))
    prop.append(dict(w=str(w), lig=float(lg)))
    print(f"5b. solucion lineal generica, w = {w}: ligaduras propagadas <= {lg:.1e} (relativo)", flush=True)
out["propagacion_lineal"] = prop

# derivada de la familia de fondos (w = 0): datos fuera de la garganta, se varia phi' y se recalcula Phi'
Fd, Vd = P.acoplamiento("gauss", 0.2), P.potencial("cuartico", 0.4, g=0.5)
eps = 1e-6
sols = []
for sgn in (+1, -1):
    y0 = dato_consistente(Fd, Vd, 2.0, 0.6, 0.0, 0.35, 0.3 + sgn * eps)

    def rhs6(x, y, Fd=Fd, Vd=Vd):
        s, _ = P.segundas_V(y, Fd, Vd)
        return [y[1], s[0], y[3], s[1], y[5], s[2]]
    sols.append(solve_ivp(rhs6, (0, 2.0), y0, method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True))
    assert sols[-1].status == 0


fam = SimpleNamespace(Fd=Fd, Vd=Vd,                              # fondo central: media de los dos (O(eps^2))
                      sol=SimpleNamespace(sol=lambda x: np.r_[(sols[0].sol(x) + sols[1].sol(x)) / 2, 0.0]))


def dlin(x):
    d = (sols[0].sol(x) - sols[1].sol(x)) / (2 * eps)            # (dr, dr', dPhi, dPhi', dphi, dphi')
    return np.array([d[3], d[0], d[1], d[4], d[5]])              # (a', d, d', q, q')


y6 = solve_ivp(lambda x, yy: P.derivadas_V(fam, x, yy, 0.0), (0.0, 2.0), dlin(0.0), method="DOP853", rtol=1e-12,
               atol=1e-14, dense_output=True)
err6 = max(np.max(np.abs(y6.sol(x) - dlin(x))) / np.max(np.abs(dlin(x))) for x in np.linspace(0, 2, 21))
# solo E^l_l: E^l_t es proporcional a d/dt y para una perturbacion estatica se anula identicamente
lg6 = max(abs(P.ligaduras_V(fam, x, dlin(x), 0.0)[1]) / np.max(np.abs(dlin(x))) for x in np.linspace(0, 2, 21))
print(f"5c. derivada de la familia (F gaussiano, V cuartico) frente al sistema lineal con w = 0: "
      f"error relativo {err6:.1e}; ligadura E^l_l {lg6:.1e}", flush=True)
out["derivada_familia"] = dict(err=float(err6), lig=float(lg6))

if not RAPIDO:
    os.makedirs(os.path.join(AQUI, "resultados"), exist_ok=True)
    json.dump(out, open(os.path.join(AQUI, "resultados", "verificar_potencial.json"), "w", encoding="utf-8"),
              indent=1)
