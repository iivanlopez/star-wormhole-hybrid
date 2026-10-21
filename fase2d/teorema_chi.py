# -*- coding: utf-8 -*-
"""
Fase 2d, Parte III: demostracion general del cambio de signo del coeficiente cinetico K de los agujeros de gusano
ULTRAESTATICOS de EsGB sostenidos por un campo electromagnetico (Maxwell o NLED), para cualquier carga escalar D.

Identidad en la capa de masas (con la ecuacion de fondo (B h)' = -r^2 P/2, r^2 P = r r'' + D^2/(4 r^2)):
    K = -a chi' + b chi,      chi = r r' - 2 h (3 r'^2 - 1),      h = dF/dl,  B = 1 - r'^2,  A = r - 4 h r',
    a = 48 h^2 r^4 B^2 A / (D^2 r' (r r' - 4 h)),
    b = r^2 [D^2 (r^2 r' - 2 h r (r'^2 + 3) + 8 h^2 r' (5 - 3 r'^2)) + 32 h^2 r^2 r' B^2] / (D^2 r' (r r' - 4 h)).
chi es (r^2/2) x (coeficiente de Phi' en la ligadura ll): sus ceros son los "horizontes" de la Fase 2c.
Teorema: en un lado donde r' > 0 y h < 0 (si c = 0, en los dos lados) a > 0 y b > 0 (todos los terminos del
numerador de b son positivos); chi(0) = 2c <= 0, chi'(0) = -D^2/4 < 0 y chi -> +oo en el infinito asintoticamente
plano.  Luego chi tiene un minimo interior negativo l_m, y alli K(l_m) = b chi(l_m) < 0.  Ademas, en el primer cero
ascendente de chi, K = -a chi' <= 0.
Parte A: SymPy (identidad, a, b, chi en la garganta y en el infinito, G - K).  Parte B: comprobacion numerica en las
familias de cinetico_nled.py (hiperbolas, superelipses; la hiperbola L = 1 es Maxwell/Ellis) y nled_concretas.py
(Born-Infeld, Euler-Heisenberg).  Salida: resultados/teorema_chi.json.   Uso: python teorema_chi.py
"""
import os
import sys
import json
import numpy as np
import sympy as sp
from scipy.integrate import cumulative_trapezoid, solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cinetico_nled as CN      # noqa: E402
import nled_concretas as NC     # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ------------------------------------------------------------------ Parte A: identidad simbolica
r, s, rho, h, D, c, rho0, beta, sig = sp.symbols("r s rho h D c rho0 beta sigma", real=True)
B = 1 - s**2
hp = (-(r * rho + D**2 / (4 * r**2)) / 2 + 2 * s * rho * h) / B        # h' de (B h)' = -r^2 P/2
A = r - 4 * h * s
C = 1 - 4 * hp
F1 = h * r**2 / D                                                     # F'(phi) = h / phi',  phi' = D/r^2
K = A**2 * r**2 - 16 * F1**2 * B * (4 * A * rho + B * C)              # cinetico_simbolico.py
G = A**2 * r**2 - 16 * F1**2 * B**2
chi = r * s - 2 * h * (3 * s**2 - 1)
chip = s**2 + r * rho - 2 * hp * (3 * s**2 - 1) - 12 * h * s * rho    # r' = s, s' = r'' = rho
chip_forma = s**2 + (3 * s**2 - 1) * D**2 / (4 * r**2 * B) + (2 * s * rho / B) * (r * s - 4 * h)
a = sp.factor(sp.cancel(-sp.diff(K, rho) / sp.diff(chip, rho)))
b = sp.factor(sp.cancel((K + a * chip) / chi))
print("PARTE A (SymPy)")
print("  a =", a)
print("  b =", b)
print("  a no depende de r'':", sp.simplify(sp.diff(a, rho)) == 0)
print("  K - (-a chi' + b chi) =", sp.simplify(K - (-a * chip + b * chi)))
print("  chi' en la capa de masas - forma compacta =", sp.simplify(chip - chip_forma))
num_b, den_b = sp.fraction(sp.together(b))
num_b_agrupado = r**2 * (D**2 * (r**2 * s - 2 * h * r * (s**2 + 3) + 8 * h**2 * s * (5 - 3 * s**2))
                         + 32 * h**2 * r**2 * s * B**2)
print("  denominador de b:", sp.factor(den_b))
print("  numerador de b - r^2 [D^2 (r^2 r' - 2 h r (r'^2 + 3) + 8 h^2 r' (5 - 3 r'^2)) + 32 h^2 r^2 r' B^2] =",
      sp.simplify(num_b - num_b_agrupado), " -> b > 0 si r' in (0, 1) y h < 0")
# garganta: r = 1, r' = 0, h = c
print("  garganta: chi(0) =", sp.simplify(chi.subs({s: 0, r: 1, h: c})),
      ", chi'(0) =", sp.simplify(chip.subs({s: 0, r: 1, h: c, rho: rho0})))
# infinito asintoticamente plano: B = beta/r^2 (r'' = beta/r^3), B h -> -sigma
chi_inf = chi.subs({h: -sig * r**2 / beta, s: sp.sqrt(1 - beta / r**2)})
print("  infinito: lim chi/r^2 =", sp.limit(chi_inf / r**2, r, sp.oo), " (sigma = |B h| en el infinito, beta > 0)")
# en un cero de chi (h = r r'/(2 (3 r'^2 - 1))): K = -a chi' y a/(A^2 r^2) = 4 r^2/D^2
h0 = r * s / (2 * (3 * s**2 - 1))
print("  en chi = 0: a/(A^2 r^2) =", sp.simplify((a / (A**2 * r**2)).subs(h, h0)),
      "; b*chi =", sp.simplify((K + a * chip).subs(h, h0)))
# G - K = 64 F'^2 B (A r'' - B h'),  A r'' - B h' = r''(r - 6 h r') + r^2 P/2
print("  G - K - 64 F'^2 B (A r'' - B h') =", sp.simplify(G - K - 64 * F1**2 * B * (A * rho - B * hp)),
      "; A r'' - B h' - [r''(r - 6 h r') + r^2 P/2] =",
      sp.simplify(A * rho - B * hp - (rho * (r - 6 * h * s) + (r * rho + D**2 / (4 * r**2)) / 2)))

# para evaluar con B = 1 - r'^2 << 1 sin redondeo, B como simbolo propio
Bs = sp.Symbol("B", positive=True)
a_B = 48 * h**2 * r**4 * Bs**2 * (r - 4 * h * s) / (D**2 * s * (r * s - 4 * h))
b_B = r**2 * (D**2 * (r**2 * s - 2 * h * r * (s**2 + 3) + 8 * h**2 * s * (5 - 3 * s**2))
              + 32 * h**2 * r**2 * s * Bs**2) / (D**2 * s * (r * s - 4 * h))
print("  a, b con B explicita - a, b:", sp.simplify(a_B.subs(Bs, B) - a), sp.simplify(b_B.subs(Bs, B) - b))
fa = sp.lambdify((r, s, h, Bs, D), a_B, "numpy")
fb = sp.lambdify((r, s, h, Bs, D), b_B, "numpy")

# ------------------------------------------------------------------ Parte C: garganta y limite interior (SymPy)
# (1) desarrollo en la garganta (c = 0, forma par r = 1 + rho0 l^2/2 + r4 l^4/24): chi = -(D^2/4) l + k3 l^3 + ...
#     y, en el cero de chi, q = 4 r^2 chi'/D^2 -> 2, es decir K/norma = -q/(2 + q) -> -1/2
l, r4 = sp.symbols("l r4", real=True)
rl = 1 + rho0 * l**2 / 2 + r4 * l**4 / 24
sl, rhol = sp.diff(rl, l), sp.diff(rl, l, 2)
Il = sp.integrate(sp.series(rl * rhol + D**2 / (4 * rl**2), l, 0, 5).removeO(), (l, 0, l))
hl = sp.series(-Il / (2 * (1 - sl**2)), l, 0, 5).removeO()
chil = sp.expand(sp.series(rl * sl - 2 * hl * (3 * sl**2 - 1), l, 0, 4).removeO())
k3 = sp.factor(chil.coeff(l, 3))
print("\nPARTE C (SymPy)")
print("  garganta: chi =", sp.factor(chil.coeff(l, 1)), "* l +", k3, "* l^3 + O(l^5)")
lz = D / (2 * sp.sqrt(k3))
print("  en el cero l_chi = D/(2 sqrt(k3)): chi' =", sp.simplify(sp.diff(chil, l).subs(l, lz)),
      " -> q = 4 r^2 chi'/D^2 = 2 + O(D^2), K/norma = -q/(2 + q) = -1/2 + O(D^2)")
# (2) limite interior D -> oo (x = D^2 l, r'' = D^2 kappa, r = 1): K >= 0 <=> B Z^2 - 8 U^2 - 96 kappa U^2 Z >= 0,
#     Z = 1 + 4 s U/B, dU/ds = 1/2 + 1/(8 kappa); la trayectoria rasante (kappa maxima) es exacta:
Ui = s * B / (2 * (1 - 3 * s**2))
kap = (1 - 3 * s**2)**2 / (24 * s**2 * B)
Zi = 1 + 4 * s * Ui / B
print("  interior: U = r'(1 - r'^2)/(2(1 - 3 r'^2)), kappa = (1 - 3 r'^2)^2/(24 r'^2 B):")
print("    dU/ds - (1/2 + 1/(8 kappa)) =", sp.simplify(sp.diff(Ui, s) - (sp.Rational(1, 2) + 1 / (8 * kap))),
      "; K =", sp.simplify(B * Zi**2 - 8 * Ui**2 - 96 * kap * Ui**2 * Zi),
      "; chi = r' - 2 U (1 - 3 r'^2)/B =", sp.simplify(s - 2 * Ui * (1 - 3 * s**2) / B))


# ------------------------------------------------------------------ Parte B: comprobacion numerica
def forma_nled(tipo, par, Dv, lmax=400.0):
    """forma de las NLED concretas (como nled_concretas.analiza), en la malla de cinetico_nled."""
    m = NC.modelo(tipo, par, Dv)
    if m is None:
        return None
    eps, P = m
    if P(1.0) - Dv * Dv / 4 <= 0:
        return None
    sol = solve_ivp(lambda l, y: [y[1], y[0] * (P(y[0]) - Dv * Dv / (4 * y[0] ** 4))], (0, lmax), [1.0, 0.0],
                    method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True)
    l = malla(lmax)
    rr, rp = sol.sol(l)
    rpp = rr * (np.array([P(v) for v in rr]) - Dv * Dv / (4 * rr**4))
    Bv = rr * rr * np.array([eps(v) for v in rr]) - Dv * Dv / (4 * rr * rr)
    if np.min(Bv[1:]) <= 0:
        return None
    return l, rr, rp, rpp, Bv


def malla(lmax):
    """la de cinetico_nled mas un tramo geometrico fino junto a la garganta (resuelve chi < 0 cuando D es pequena)."""
    m = CN.malla(lmax)
    return np.concatenate(([0.0], np.geomspace(1e-9, 0.05, 6000, endpoint=False), m[m >= 0.05]))


def cruces(x, y, sube=None):
    """ceros de y (interpolados); sube = True/False filtra los ascendentes/descendentes."""
    i = np.nonzero(np.sign(y[:-1]) * np.sign(y[1:]) < 0)[0]
    if sube is not None:
        i = i[(y[i + 1] > y[i]) == sube]
    return [float(x[k] - y[k] * (x[k + 1] - x[k]) / (y[k + 1] - y[k])) for k in i]


def campos(l, rr, rp, rpp, Bv, Dv, ce):
    """h, K, G, norma, chi y chi' (capa de masas) en el lado l >= 0 con asimetria efectiva ce."""
    r2P = rr * rpp + Dv**2 / (4 * rr**2)
    I = cumulative_trapezoid(r2P, l, initial=0.0)
    hv = (ce - I / 2) / Bv
    hpv = (-r2P / 2 + 2 * rp * rpp * hv) / Bv
    Av = rr - 4 * hv * rp
    Cv = 1 - 4 * hpv
    F1v = hv * rr**2 / Dv
    E = 4 * Av * rpp + Bv * Cv
    Kv = Av**2 * rr**2 - 16 * F1v**2 * Bv * E
    Gv = Av**2 * rr**2 - 16 * F1v**2 * Bv**2
    nrm = Av**2 * rr**2 + 16 * F1v**2 * Bv * np.abs(E)
    chiv = rr * rp - 2 * hv * (3 * rp**2 - 1)
    chipv = rp**2 + (3 * rp**2 - 1) * Dv**2 / (4 * rr**2 * Bv) + (2 * rp * rpp / Bv) * (rr * rp - 4 * hv)
    return dict(h=hv, K=Kv, G=Gv, norma=nrm, chi=chiv, chip=chipv)


def lado(l, rr, rp, rpp, Bv, Dv, ce):
    """lado l >= 0 con asimetria efectiva ce (el lado izquierdo de c equivale al derecho de -c)."""
    f = campos(l, rr, rp, rpp, Bv, Dv, ce)
    hv, Kv, Gv, nrm, chiv, chipv = f["h"], f["K"], f["G"], f["norma"], f["chi"], f["chip"]
    j = rp > 1e-6                                        # a, b ~ 1/r': se evaluan lejos de r' = 0
    with np.errstate(all="ignore"):
        av = fa(rr[j], rp[j], hv[j], Bv[j], Dv)
        bv = fb(rr[j], rp[j], hv[j], Bv[j], Dv)
        escala = np.abs(av * chipv[j]) + np.abs(bv * chiv[j]) + np.abs(Kv[j])
        err_id = float(np.max(np.abs(Kv[j] - (-av * chipv[j] + bv * chiv[j])) / escala))
    # chi' de la capa de masas integrado frente a chi (valida la reconstruccion del fondo)
    err_chip = float(np.max(np.abs(chiv - chiv[0] - cumulative_trapezoid(chipv, l, initial=0.0)))
                     / np.max(np.abs(chiv[l < 0.5 * l[-1]])))
    creciente = bool(np.all(hv[1:] * rp[1:] < 0))       # lado "creciente": h r' < 0  ->  a, b > 0
    kn = Kv / nrm
    # primer minimo de chi: primer cero ascendente de chi' (chi'(0) = -D^2/4 < 0)
    zm = cruces(l, chipv, sube=True)
    lm = zm[0] if zm else None
    if lm is not None:
        bm = float(fb(np.interp(lm, l, rr), np.interp(lm, l, rp), np.interp(lm, l, hv), np.interp(lm, l, Bv), Dv))
        chim = float(np.interp(lm, l, chiv))
        Km = float(np.interp(lm, l, Kv))
        nm = float(np.interp(lm, l, nrm))
        minimo = dict(l_m=lm, chi_m=chim, b_m=bm, K_norma_m=Km / nm, err_K_igual_b_chi=abs(Km - bm * chim) / nm)
    else:
        minimo = None
    zc = cruces(l, chiv, sube=True)
    lchi = zc[0] if zc else None
    k_en_chi = float(np.interp(lchi, l, kn)) if lchi is not None else None
    neg = np.nonzero(Kv < 0)[0]
    zK = cruces(l, Kv)
    G_en_ceros_K = [float(np.interp(z, l, Gv)) for z in zK]
    return dict(creciente=creciente,
                a_b_positivas=bool(np.all(av > 0) and np.all(bv > 0)) if creciente else None,
                chi_garganta=float(chiv[0]), chip_garganta=float(chipv[0]), chi_final=float(chiv[-1]),
                minimo=minimo, l_chi=lchi, K_norma_en_l_chi=k_en_chi, min_K_norma=float(np.min(kn)),
                l_Kneg=[float(l[neg[0]]), float(l[neg[-1]])] if len(neg) else None,
                ceros_K=zK, G_en_ceros_K=G_en_ceros_K, err_identidad=err_id, err_chip=err_chip)


def caso(etiqueta, l, rr, rp, rpp, Bv, Dv, cv):
    out = []
    for nombre, ce in (("der", cv), ("izq", -cv)):
        if ce > 0:                                       # lado "decreciente" de c != 0: el teorema usa el otro
            continue
        d = lado(l, rr, rp, rpp, Bv, Dv, ce)
        d.update(etiqueta=etiqueta, D=Dv, c=cv, lado=nombre)
        out.append(d)
    return out


if __name__ == "__main__":
    print("\nPARTE B (numerica): lado(s) con h r' < 0")
    res = []
    for L in (0.03, 0.3, 1.0, 3.0, 30.0, 300.0):
        for Dv in (0.01, 0.1, 0.5, 2.0, 7.0, 30.0, 100.0):
            for cv in (0.0, 0.3):
                l = malla(400.0 * max(1.0, L))
                rr, rp, rpp, Bv = CN.forma("hiperbola", L, l)
                res += caso(f"hiperbola L={L:g}" + (" (Maxwell)" if L == 1.0 else ""), l, rr, rp, rpp, Bv, Dv, cv)
    for n in (2, 5):
        for Dv in (0.01, 0.1, 0.5, 2.0, 7.0, 30.0, 100.0):
            for cv in (0.0, 0.3):
                l = malla(400.0)
                with np.errstate(divide="ignore"):
                    rr, rp, rpp, Bv = CN.forma("superelipse", n, l)
                res += caso(f"superelipse n={n}", l, rr, rp, rpp, Bv, Dv, cv)
    for tipo, pars in (("BI", (0.5, 2.0, 20.0)), ("EH", (-0.1, 0.0, 0.3))):
        for par in pars:
            for Dv in (0.1, 0.5, 1.0, 3.0):
                f = forma_nled(tipo, par, Dv)
                if f is not None:
                    res += caso(f"{tipo} {par:g}", *f, Dv, 0.0)
    ok = [o for o in res if o["creciente"]]
    mins = [o for o in ok if o["minimo"] is not None]
    print(f"  configuraciones (lados): {len(res)}; con h r' < 0 en todo el lado: {len(ok)}")
    print(f"  a > 0 y b > 0 en todos los puntos con r' > 1e-6: {all(o['a_b_positivas'] for o in ok)}")
    print(f"  max |K - (-a chi' + b chi)|/(|a chi'| + |b chi| + |K|) = {max(o['err_identidad'] for o in ok):.1e}")
    print(f"  max |chi - chi(0) - int chi'_capa de masas|/max|chi| = {max(o['err_chip'] for o in ok):.1e}")
    print(f"  chi'(0) = -D^2/4: max error relativo = "
          f"{max(abs(o['chip_garganta'] + o['D']**2 / 4) / (o['D']**2 / 4) for o in ok):.1e}")
    print(f"  chi(0) = 2c_lado: max error = {max(abs(o['chi_garganta'] + 2 * o['c']) for o in ok):.1e}")
    print(f"  chi > 0 al final de la malla: {sum(o['chi_final'] > 0 for o in ok)}/{len(ok)}")
    print(f"  con minimo interior de chi: {len(mins)}/{len(ok)}; chi < 0 alli: "
          f"{sum(o['minimo']['chi_m'] < 0 for o in mins)}; K < 0 alli: {sum(o['minimo']['K_norma_m'] < 0 for o in mins)}; "
          f"max |K - b chi|/norma alli = {max(o['minimo']['err_K_igual_b_chi'] for o in mins):.1e}")
    print(f"  con cero ascendente de chi: {sum(o['l_chi'] is not None for o in ok)}/{len(ok)}; "
          f"K <= 0 en ese cero: {sum(o['K_norma_en_l_chi'] is not None and o['K_norma_en_l_chi'] <= 1e-9 for o in ok)}")
    print(f"  K < 0 en algun punto: {sum(o['min_K_norma'] < 0 for o in ok)}/{len(ok)}")
    print(f"  G > 0 en todos los ceros de K: {sum(all(g > 0 for g in o['G_en_ceros_K']) for o in ok)}/{len(ok)}")
    print("\n  muestra (lado creciente, c = 0):")
    for o in ok:
        if o["c"] == 0.0 and o["lado"] == "der" and (o["D"] in (0.1, 2.0, 30.0) or o["etiqueta"][:2] in ("BI", "EH")):
            lk, m = o["l_Kneg"], o["minimo"]
            print(f"    {o['etiqueta']:22s} D = {o['D']:6g}: K < 0 en l = [{lk[0]:.3g}, {lk[1]:.3g}]; "
                  f"min de chi en l = {m['l_m']:.3g} (K/norma = {m['K_norma_m']:+.3f}); "
                  f"cero de chi en l = {o['l_chi']:.3g} (K/norma = {o['K_norma_en_l_chi']:+.3f}); "
                  f"min K/norma = {o['min_K_norma']:+.3f}")
    json.dump(res, open("resultados/teorema_chi.json", "w", encoding="utf-8"), indent=1)
