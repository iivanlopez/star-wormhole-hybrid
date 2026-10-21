# -*- coding: utf-8 -*-
"""
¿Hay ventana de estabilidad con tension dependiente del escalar?  Analisis de la relacion de zona cercana
P(s; zeta, u, F0, F2, l1, s2) = 0 de tension_cercana.py, con t = l1 u = Lambda_s'(phi0) phi1 / (a^2 r0) y
s2 = (Lambda_s'' + rho beta mu)/a.
  1. Identidades (SymPy):  u^2 P(1) = 8 u^2 [ (2t - zeta(u^2+4))^2 + 16 F0 zeta u^2 (1 + s2 + 8 F2 (1 - zeta)) ],
     coeficiente principal  u^2 zeta (u^2+4)^2 [zeta(u^2+4) - 8(1 - F0)]  (zeta = 0: P(1) = 32 t^2 u^2, LC = -32 u^4).
     => para 0 <= zeta < zeta_c y 1 + s2 + 8F2(1 - zeta) >= 0 hay una raiz real s >= 1: C >= 1.
  2. Minimo global de la abscisa espectral max Re s sobre (t, s2) (evolucion diferencial + Nelder-Mead) para cinco
     fondos y zeta en [0, 3]; el minimizador se comprueba con raices a 80 cifras (mpmath.polyroots).
Salida: resultados/ventanas_tension.json.   Uso: python ventanas_tension.py
"""
import sys
import json
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.optimize import differential_evolution, minimize

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

d = json.load(open("resultados/tension_cercana.json", encoding="utf-8"))
u, F0, F2, z, l1, s2, t = sp.symbols("u F0 F2 zeta l1 s2 t", real=True)
s = sp.symbols("s", positive=True)
P = sp.sympify(d["polinomio"], locals=dict(s=s, u=u, F0=F0, F2=F2, zeta=z, l1=l1, s2=s2))
Pt = sp.expand(P.subs(l1, t / u) * u**2)
P1 = sp.factor(Pt.subs(s, 1))
LC = sp.factor(sp.Poly(Pt, s).LC())
cuadrado = 8 * u**2 * ((2 * t - z * (u**2 + 4))**2 + 16 * F0 * z * u**2 * (1 + s2 + 8 * F2 * (1 - z)))
print("1. Identidades")
print("   u^2 P(1) =", P1)
print("   u^2 P(1) - 8u^2[(2t - zeta(u^2+4))^2 + 16 F0 zeta u^2 (1 + s2 + 8F2(1 - zeta))] =", sp.simplify(P1 - cuadrado))
print("   coeficiente principal (u^2 LC) =", LC)
print("   zeta = 0:  u^2 P(1) =", sp.factor(Pt.subs({z: 0, s: 1})), ";  u^2 LC =", sp.factor(sp.Poly(Pt.subs(z, 0), s).LC()))

coefs = sp.Poly(Pt, s).all_coeffs()
cf = sp.lambdify((u, F0, F2, z, t, s2), coefs)


def alfa(par, tt, ss):
    c = np.array(cf(*par, tt, ss), dtype=float)
    while len(c) > 1 and abs(c[0]) < 1e-13 * np.abs(c).max():
        c = c[1:]
    r = [x for x in np.roots(c) if abs(x) > 1e-9]
    return max(x.real for x in r) if r else -np.inf


def alfa_mp(par, tt, ss):
    mp.mp.dps = 80
    sub = dict(zip((u, F0, F2, z, t, s2), [sp.Rational(repr(float(x))) for x in (*par, tt, ss)]))
    c = [mp.mpf(sp.Rational(x.subs(sub))) for x in coefs]
    while c and c[0] == 0:
        c = c[1:]
    r = [x for x in mp.polyroots(c, maxsteps=500, extraprec=400) if abs(x) > mp.mpf("1e-30")]
    return float(max(mp.re(x) for x in r))


FONDOS = {"KKK referencia": (-0.7, 0.03, 0.03), "cuadratico 0.3087": (0.1957, 0.0225, 0.18),
          "cuadratico 0.3465": (0.0723, 0.0081, 0.18), "gaussiano 0.2326": (-0.7611, 0.05985, -0.1191),
          "KKK f0=1.06 (u=-3)": (-3.0, 0.003, 0.003)}
ZETAS = (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
salida = dict(identidad_P1=str(P1), coeficiente_principal=str(LC), minimos=[])
print("\n2. Minimo de max Re s sobre (t, s2), t en [-60, 60], s2 = 10 sinh(y) con y en [-12, 12]:")
for nom, (uu, f0, f2) in FONDOS.items():
    zc = 2 * (1 - f0) / (1 + uu * uu / 4)
    print(f"\n   {nom}: u = {uu}, F0 = {f0}, F2 = {f2:+}, zeta_c = {zc:.3f}")
    for zz in ZETAS:
        par = (uu, f0, f2, zz)
        fobj = lambda x: alfa(par, x[0], 10 * np.sinh(x[1]))
        de = differential_evolution(fobj, bounds=[(-60, 60), (-12, 12)], seed=1, tol=1e-10, maxiter=300, popsize=40,
                                    polish=False)
        nm = minimize(fobj, de.x, method="Nelder-Mead", bounds=[(-60, 60), (-12, 12)],
                      options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
        x = nm.x if nm.fun < de.fun else de.x
        tt, ss = x[0], 10 * np.sinh(x[1])
        a64, a80 = fobj(x), alfa_mp(par, tt, ss)
        # sin acotar t: tendencia asintotica (se informa aparte)
        nl = minimize(fobj, x, method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
        tl, sl = nl.x[0], 10 * np.sinh(nl.x[1])
        al = alfa_mp(par, tl, sl)
        salida["minimos"].append(dict(fondo=nom, zeta=zz, zeta_c=zc, t=tt, s2=ss, alfa_float=a64, alfa_80=a80,
                                      sin_acotar=dict(t=tl, s2=sl, alfa_80=al)))
        print(f"      zeta = {zz:4.2f}{' (> zeta_c)' if zz > zc else '           '}: min max Re s = {a80:+.5f}"
              f" (doble precision: {a64:+.5f}) en t = {tt:+.4g}, s2 = {ss:+.4g};   sin acotar t: {al:+.5f}"
              f" en t = {tl:+.3g}")
json.dump(salida, open("resultados/ventanas_tension.json", "w", encoding="utf-8"), indent=1)
print("\nGuardado: resultados/ventanas_tension.json")
