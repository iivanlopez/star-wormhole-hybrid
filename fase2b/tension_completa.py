# -*- coding: utf-8 -*-
"""
Cascara con tension dependiente del escalar en el calculo completo (estabilidad_F.CapaF con dLam = Lambda_s'(phi0),
d2Lam = Lambda_s''(phi0) y mu = 0), frente a la relacion de zona cercana (tension_cercana.py).

Correspondencia (r0 = 1): Lambda_s' = a^2 t / u (t = l1 u) y Lambda_s'' + rho beta mu = a s2 (con mu = 0, d2Lam = a s2).
  1. Validacion: con dLam = d2Lam = 0 se recupera la raiz de KKK (zeta = 1, w = -391.88143); con dLam, d2Lam != 0 el
     empalme tt sigue siendo redundante (Bianchi): residuo relativo del empalme tt sobre el vector nulo de los otros cuatro.
  2. Casos (fondo, zeta, t, s2): raices predichas con Re s > 0; en el calculo completo, las reales con raices_menores en
     una malla local y en un barrido global de w en [-1e6, -1e-4], y las complejas por secante sobre Dc desde la
     prediccion.
  3. (python tension_completa.py sobre_zc) cuatro casos con zeta = 2.5 > zeta_c, con y sin tension;
     (python tension_completa.py extremos) pendientes t = -60 a -1000 con zeta = 3 y el s2 optimo de la zona cercana.
Salida: resultados/tension_completa.json (y tension_completa_sobre_zc.json).   Uso: python tension_completa.py [sobre_zc]
"""
import sys
import json
import math
import cmath
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import sympy as sp
from fondo_F import FondoF
from estabilidad_F import CapaF, vector_garganta, raices_menores, secante, Dc

FONDOS = {"KKK referencia": ("exp", 0.03, 0.0, -0.7), "cuadratico 0.3465": None}
CASOS = [("KKK referencia", 0.0, 0.0, 0.0), ("KKK referencia", 0.0, 2.0, 0.0), ("KKK referencia", 0.0, 6.0, 0.0),
         ("KKK referencia", 0.0, -2.0, 0.0), ("KKK referencia", 0.0, 2.0, 100.0),
         ("KKK referencia", 1.0, 0.0, 0.0), ("KKK referencia", 1.0, "t*", 10.0), ("KKK referencia", 1.0, "t*", 100.0),
         ("KKK referencia", 1.0, "t*", 1000.0), ("KKK referencia", 1.0, 6.0, 0.0), ("KKK referencia", 1.0, -2.0, 0.0),
         ("cuadratico 0.3465", 0.0, 2.0, 0.0), ("cuadratico 0.3465", 0.0, 6.0, 0.0),
         ("cuadratico 0.3465", 1.0, "t*", 100.0), ("cuadratico 0.3465", 1.0, "t*", 1000.0)]
# por encima de zeta_c: el mejor punto de la caja |t| <= 60 de ventanas_tension.py (resonancia mas debil)
CASOS_SOBRE_ZC = [("KKK referencia", 2.5, -60.0, 82.51), ("cuadratico 0.3465", 2.5, -60.0, 8284.0),
                  ("KKK referencia", 2.5, 0.0, 0.0), ("cuadratico 0.3465", 2.5, 0.0, 0.0)]
# pendientes extremas con zeta = 3 (cuadratico 0.3465): s2 optimo de la zona cercana para cada t
CASOS_EXTREMOS = [("cuadratico 0.3465", 3.0, -60.0, "opt"), ("cuadratico 0.3465", 3.0, -300.0, "opt"),
                  ("cuadratico 0.3465", 3.0, -1000.0, "opt"), ("KKK referencia", 3.0, -300.0, "opt")]


def parametros():
    tr = json.load(open("resultados/trayectorias_robustas.json", encoding="utf-8"))
    o = next(x for x in tr if x["nombre"] == "cuadratico 0.3465")
    FONDOS["cuadratico 0.3465"] = (o["tipo"], o["alpha"], o["phi0"], o["phi1"])


def polinomio():
    d = json.load(open("resultados/tension_cercana.json", encoding="utf-8"))
    u, F0, F2, z, l1, s2 = sp.symbols("u F0 F2 zeta l1 s2", real=True)
    s = sp.symbols("s", positive=True)
    P = sp.Poly(sp.sympify(d["polinomio"], locals=dict(s=s, u=u, F0=F0, F2=F2, zeta=z, l1=l1, s2=s2)), s)
    return sp.lambdify((u, F0, F2, z, l1, s2), P.all_coeffs())


def prediccion(cf, u, F0, F2, z, t, s2):
    c = np.array(cf(u, F0, F2, z, t / u, s2), dtype=float)
    r = np.roots(c)
    return [complex(x * x) for x in r if x.real > 1e-9 and abs(x) > 1e-9]


def s2_optimo(cf, u, F0, F2, z, t):
    """s2 que minimiza max Re s de la zona cercana con t fijo (malla en asinh y refinamiento acotado)."""
    from scipy.optimize import minimize_scalar

    def g(y):
        r = [x for x in np.roots(np.array(cf(u, F0, F2, z, t / u, 10 * np.sinh(y)), dtype=float)) if abs(x) > 1e-9]
        return max(x.real for x in r)
    ys = np.linspace(-14, 14, 561)
    y0 = ys[int(np.argmin([g(y) for y in ys]))]
    return 10 * np.sinh(minimize_scalar(g, bounds=(y0 - 0.05, y0 + 0.05), method="bounded").x)


def caso(args):
    nom, z, t, s2, tipo, al, p0, p1 = args
    cf = polinomio()
    f = FondoF(tipo, al, p0, p1).resolver()
    a = f.y0()[3]
    F0, F1, F2, _ = f.Fd(p0)
    u = p1
    if t == "t*":
        t = z * (u * u + 4) / 2
    if s2 == "opt":
        s2 = s2_optimo(cf, u, F0, F2, z, t)
    capa = CapaF(f, zeta=z, mu=0.0, dLam=a * a * t / u, d2Lam=a * s2)
    Wp = prediccion(cf, u, F0, F2, z, t, s2)
    out = dict(nombre=nom, zeta=z, t=t, s2=s2, a=a, W_pred=[[w.real, w.imag] for w in Wp], reales=[], complejas=[])
    # reales: malla local alrededor de cada prediccion real y barrido global
    ws = -np.logspace(-4, 6, 241)
    for w in Wp:
        if abs(w.imag) < 1e-9 * abs(w):
            ws = np.concatenate([ws, -a * a * w.real * np.linspace(0.5, 2.0, 61)])
    ws = np.unique(ws)[::-1]
    rr = [x for x in raices_menores(f, capa, ws) if np.isfinite(x)]
    for x in rr:
        Dx = abs(Dc(f, capa, x))
        Dv = max(abs(Dc(f, capa, x * (1 + 1e-4))), abs(Dc(f, capa, x * (1 - 1e-4))))
        out["reales"].append(dict(w=x, C=-x / a**2, Dc=Dx, Dc_vecino=Dv))
    for w in Wp:
        if abs(w.imag) >= 1e-9 * abs(w):
            if w.imag > 0:
                continue
            res = secante(f, capa, -a * a * w)
            wn = complex(res[0])
            out["complejas"].append(dict(W_pred=[w.real, w.imag], W=[(-wn / a**2).real, (-wn / a**2).imag],
                                         Dc=abs(Dc(f, capa, wn)), Dc_vecino=abs(Dc(f, capa, wn * (1 + 1e-3)))))
    return out


def validacion():
    out = {}
    f = FondoF("exp", 0.03, 0.0, -0.7).resolver()
    capa = CapaF(f, zeta=1.0)
    w0 = -391.88143094
    out["compat_Dc"] = [float(np.real(Dc(f, capa, w0 * (1 - 1e-6)))), float(np.real(Dc(f, capa, w0 * (1 + 1e-6))))]
    a = f.y0()[3]
    res = []
    for dl, d2 in ((a * a * 3.0, a * 5.0), (-a * a * 2.0, -a * 40.0), (a * a * 10.0, 0.0)):
        capa = CapaF(f, zeta=1.0, mu=0.0, dLam=dl, d2Lam=d2)
        for w in (-50.0, -400.0, -3000.0, 200.0 + 50j):
            _, r_tt = vector_garganta(f, capa, w)
            res.append(float(r_tt))
    out["bianchi_res_tt_max"] = max(res)
    return out


def signo_crecimiento(f, capa, w, separacion=40.0):
    """Signo del coeficiente de la solucion creciente para |w| muy grande, donde D(w) se anula por el factor e^{-k l1}.
    Se integra el dato de menores (renormalizando solo con factores positivos) hasta l1 con int k dl ~ separacion, donde
    la solucion creciente domina por e^{2 separacion}, y se toma el signo de q' + (k + r'/r) q con k local."""
    from scipy.integrate import solve_ivp
    from estabilidad_F import derivadas, vector_menores
    k0 = math.sqrt(-w)
    l1 = min(1.0, separacion / k0)
    y = np.real(vector_menores(f, capa, w)).astype(float)
    y = y / np.max(np.abs(y))
    l = 0.0

    def grande(x, yy):
        return np.max(np.abs(yy)) - 1e80
    grande.terminal = True
    while l < l1 * (1 - 1e-14):
        s = solve_ivp(lambda x, yy: derivadas(f, x, yy, w), (l, l1), y, method="DOP853", rtol=1e-10,
                      atol=1e-14 * np.max(np.abs(y)), events=[grande])
        y, l = s.y[:, -1], s.t[-1]
        y = y / np.max(np.abs(y))
    r, rp, Phi = f.sol.sol(l1)[:3]
    k = k0 * math.exp(-Phi)
    return float(np.sign(y[4] + (k + rp / r) * y[3]))


def raices_grandes(f, capa, a, Cs):
    """Raices reales con C = -w/a^2 en la malla Cs (creciente), por cambio de signo de signo_crecimiento y biseccion."""
    sg = [signo_crecimiento(f, capa, -a * a * C) for C in Cs]
    out = []
    for i in range(len(Cs) - 1):
        if sg[i] != sg[i + 1]:
            lo, hi, slo = Cs[i], Cs[i + 1], sg[i]
            for _ in range(40):
                mid = math.sqrt(lo * hi)
                sm = signo_crecimiento(f, capa, -a * a * mid)
                lo, hi, slo = (mid, hi, sm) if sm == slo else (lo, mid, slo)
            out.append(math.sqrt(lo * hi))
    return out


def grande(args):
    nom, z, t, s2, tipo, al, p0, p1, Cs = args
    f = FondoF(tipo, al, p0, p1).resolver()
    a = f.y0()[3]
    capa = CapaF(f, zeta=z, mu=0.0, dLam=a * a * t / p1, d2Lam=a * s2)
    return dict(nombre=nom, zeta=z, t=t, s2=s2, a=a, C=raices_grandes(f, capa, a, list(Cs)))


def caso_bajo(args):
    """Raices reales con C = -w/a^2 entre 1e-4 y 10 (metodo de menores y D renormalizado)."""
    nom, z, t, s2, tipo, al, p0, p1 = args
    f = FondoF(tipo, al, p0, p1).resolver()
    a = f.y0()[3]
    capa = CapaF(f, zeta=z, mu=0.0, dLam=a * a * t / p1, d2Lam=a * s2)
    ws = -a * a * np.logspace(-4, 1, 81)
    return [-x / a**2 for x in raices_menores(f, capa, ws) if np.isfinite(x)]


def imprimir(o):
    print(f"\n {o['nombre']}, zeta = {o['zeta']:g}, t = {o['t']:.4g}, s2 = {o['s2']:g} (a r0 = {o['a']:.3g}):")
    print("   predicho:  " + ", ".join(f"{w[0]:.5g}" + (f"{w[1]:+.4g}i" if abs(w[1]) > 1e-9 * abs(w[0]) else "")
                                     for w in o["W_pred"]))
    print("   completo reales: " + (", ".join(f"{x['C']:.5g} (|Dc| {x['Dc']:.1e} / {x['Dc_vecino']:.1e})"
                                             for x in o["reales"]) or "ninguna"))
    for x in o["complejas"]:
        print(f"   completo complejo: W = {x['W'][0]:+.5g}{x['W'][1]:+.4g}i desde {x['W_pred'][0]:+.5g}{x['W_pred'][1]:+.4g}i"
              f"  (|Dc| {x['Dc']:.1e} / {x['Dc_vecino']:.1e})")


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    parametros()
    if "barrido_t" in sys.argv[1:]:
        # zeta = 0, s2 = 0 y pendientes t > 2, donde el desplazamiento rigido ya no crece: todas las raices reales con
        # C de 1e-4 a 10 (raices_menores) y de 10 a 1e9 (signo_crecimiento)
        lista = [("KKK referencia", t) for t in (3.0, 4.0, 10.0, 20.0, 50.0)] + \
                [("cuadratico 0.3465", t) for t in (3.0, 4.0, 10.0, 20.0, 50.0)]
        args = [(nom, 0.0, t, 0.0, *FONDOS[nom], np.logspace(1, 9, 161)) for nom, t in lista]
        with ProcessPoolExecutor(max_workers=8) as ex:
            res = list(ex.map(grande, args))
            bajos = list(ex.map(caso_bajo, [(nom, 0.0, t, 0.0, *FONDOS[nom]) for nom, t in lista]))
        for o, b in zip(res, bajos):
            o["C_bajos"] = b
            print(f" {o['nombre']}, zeta = 0, t = {o['t']:g}, s2 = 0: C en [1e-4, 10] = {[round(c, 5) for c in b]};"
                  f"  C en [10, 1e9] = {[round(c, 2) for c in o['C']]}")
        json.dump(res, open("resultados/tension_completa_barrido_t.json", "w", encoding="utf-8"), indent=1)
        raise SystemExit
    if "grande" in sys.argv[1:]:
        # raices muy rapidas (C hasta 1e8), fuera del barrido general: control con KKK (t = 6, raiz conocida en
        # C = 422.5) y los dos casos cuadraticos con zeta = 0 cuya raiz predicha mas rapida queda fuera de |w| <= 1e6
        Cs = np.logspace(1, 8, 141)
        lista = [("KKK referencia", 0.0, 6.0, 0.0), ("KKK referencia", 0.0, 2.0, 0.0),
                 ("cuadratico 0.3465", 0.0, 2.0, 0.0), ("cuadratico 0.3465", 0.0, 6.0, 0.0)]
        cf = polinomio()
        with ProcessPoolExecutor(max_workers=4) as ex:
            res = list(ex.map(grande, [(nom, z, t, s2, *FONDOS[nom], Cs) for nom, z, t, s2 in lista]))
        for o in res:
            tipo, al, p0, p1 = FONDOS[o["nombre"]]
            F0, _, F2, _ = FondoF(tipo, al, p0, p1).Fd(p0)
            Wp = prediccion(cf, p1, F0, F2, o["zeta"], o["t"], o["s2"])
            o["W_pred"] = [[w.real, w.imag] for w in Wp]
            print(f" {o['nombre']}, zeta = {o['zeta']:g}, t = {o['t']:g}, s2 = {o['s2']:g}: C predichas = "
                  f"{[round(w.real, 4) for w in Wp if abs(w.imag) < 1e-9 * abs(w)]};  C completas (10 a 1e8) = "
                  f"{[round(c, 4) for c in o['C']]}")
        json.dump(res, open("resultados/tension_completa_grande.json", "w", encoding="utf-8"), indent=1)
        raise SystemExit
    if "sobre_zc" in sys.argv[1:] or "extremos" in sys.argv[1:]:
        lista, salida = ((CASOS_SOBRE_ZC, "sobre_zc") if "sobre_zc" in sys.argv[1:] else (CASOS_EXTREMOS, "extremos"))
        args = [(nom, z, t, s2, *FONDOS[nom]) for nom, z, t, s2 in lista]
        with ProcessPoolExecutor(max_workers=4) as ex:
            res = list(ex.map(caso, args))
        for o in res:
            imprimir(o)
        json.dump(res, open(f"resultados/tension_completa_{salida}.json", "w", encoding="utf-8"), indent=1)
        raise SystemExit
    val = validacion()
    print("Validacion:")
    print(f"   tension constante, KKK zeta = 1: Dc(w*(1 -+ 1e-6)) = {val['compat_Dc'][0]:+.3e}, {val['compat_Dc'][1]:+.3e}"
          "  (cambio de signo: la raiz se conserva)")
    print(f"   Lambda_s' y Lambda_s'' != 0: residuo relativo del empalme tt <= {val['bianchi_res_tt_max']:.2e} (Bianchi)")
    args = [(nom, z, t, s2, *FONDOS[nom]) for nom, z, t, s2 in CASOS]
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(caso, args))
    print("\nCasos: C = -w/a^2 predicho (zona cercana, raices con Re s > 0) y del calculo completo")
    for o in res:
        imprimir(o)
    json.dump(dict(validacion=val, casos=res), open("resultados/tension_completa.json", "w", encoding="utf-8"), indent=1)
    print("\nGuardado: resultados/tension_completa.json")
