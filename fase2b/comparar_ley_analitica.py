# -*- coding: utf-8 -*-
"""
Contraste de la relacion de dispersion de zona cercana (ley_analitica.py, resultados/ley_analitica.json) con todas las
raices numericas de la Parte II y de la Parte I:
  - modos reales con zeta = 0 y zeta = 1 de los 73 fondos de ley_aceleracion.json (KKK, cuadraticos, dobles, gaussianos);
  - trayectorias en zeta de trayectorias_robustas.json;
  - modos complejos de complejos_F.json y modos_complejos_F.json;
  - prediccion y verificacion: modos complejos que la relacion predice (KKK referencia con zeta = 1.8, cuadratico
    0.3087 con zeta = 2, y los cinco fondos que con zeta = 1 estan en la fase compleja: KKK f0 = 1.06, 1.22, 1.55 y
    cuadraticos 0.4921 y 0.4779), por secante sobre Dc desde w = -a^2 W.
La rama fisica es Re s > 0 (s^2 = W = -w/a^2).
Salida: resultados/comparar_ley_analitica.json.   Uso: python comparar_ley_analitica.py
"""
import sys
import json
import math
import numpy as np
import sympy as sp
from fondo_F import FondoF
from estabilidad_F import CapaF, secante, Dc

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ley = json.load(open("resultados/ley_analitica.json", encoding="utf-8"))
u_, F0_, F2_, z_ = sp.symbols("u F0 F2 zeta", real=True)
s_ = sp.symbols("s", positive=True)
coefs = sp.lambdify((u_, F0_, F2_, z_), sp.Poly(sp.sympify(ley["polinomio"], locals=dict(s=s_, u=u_, F0=F0_, F2=F2_,
                                                                                        zeta=z_)), s_).all_coeffs())


def prediccion(u, F0, F2, z):
    """W = s^2 de las raices con Re s > 0 (s = 0 excluida)."""
    rts = np.roots(np.array(coefs(u, F0, F2, z), dtype=float))
    out = [complex(x * x) for x in rts if x.real > 1e-9 and abs(x) > 1e-9]
    reales = sorted(w.real for w in out if abs(w.imag) <= 1e-9 * max(1.0, abs(w)))
    complejas = sorted((w for w in out if abs(w.imag) > 1e-9 * max(1.0, abs(w))), key=lambda w: (w.real, w.imag))
    return reales, complejas


def parametros():
    """nombre -> (tipo, alpha, phi0, phi1, kw)."""
    par = {}
    for o in json.load(open("resultados/estabilidad_robusta.json", encoding="utf-8")):
        par[o["nombre"]] = (o["tipo"], o["alpha"], o["phi0"], o["phi1"], {})
    for o in json.load(open("resultados/estabilidad_gauss.json", encoding="utf-8")):
        par[f"gauss {o['alpha_eta']:.4f}"] = ("gauss", o["alpha"], o["phi0"], o["phi1"], {"gamma": 1.0})
    for o in json.load(open("resultados/trayectorias_robustas.json", encoding="utf-8")):
        par.setdefault(o["nombre"], (o["tipo"], o["alpha"], o["phi0"], o["phi1"], {}))
    return par


def buscar(nom):
    if nom in par:
        return par[nom]
    claves = [k for k in par if k.startswith(nom + " ")]
    assert len(claves) == 1, (nom, claves)
    return par[claves[0]]


def fisicos(p):
    if isinstance(p, str):
        p = buscar(p)
    tipo, al, p0, p1, kw = p
    f = FondoF(tipo, al, p0, p1, kw=kw)
    F0, F1, F2, _ = f.Fd(p0)
    return f, p1, F0, F2, f.y0()[3]


par = parametros()
salida = dict(filas=[], trayectorias=[], complejos=[], prediccion=[])

# ---------------------------------------------------------------- modos reales, zeta = 0 y 1
print("Modos reales (C = -w/a^2): numerico frente a zona cercana")
filas = json.load(open("resultados/ley_aceleracion.json", encoding="utf-8"))["filas"]
err = {0.0: [], 1.0: []}
for fila in filas:
    f, u, F0, F2, a = fisicos(par[fila["nombre"]])
    reales, complejas = prediccion(u, F0, F2, fila["zeta"])
    pares = []
    for Cn in fila["C"]:
        if reales:
            Cp = min(reales, key=lambda c: abs(math.log(c / Cn)))
            pares.append((Cn, Cp, Cp / Cn - 1))
    salida["filas"].append(dict(nombre=fila["nombre"], familia=fila["familia"], zeta=fila["zeta"], a=a, u=u, F0=F0, F2=F2,
                                C_num=fila["C"], C_pred=reales, W_complejas=[[w.real, w.imag] for w in complejas],
                                pares=pares))
    for Cn, Cp, e in pares:
        err[fila["zeta"]].append((fila["familia"], a, Cn, Cp, e))
for z in (0.0, 1.0):
    E = err[z]
    rapidos = [x for x in E if x[2] > 0.3]
    lentos = [x for x in E if x[2] <= 0.3]
    print(f"\n zeta = {z:g}: {len(E)} modos reales numericos")
    for nombre, sub in (("rapidos (C > 0.3)", rapidos), ("lentos (C <= 0.3)", lentos)):
        if not sub:
            continue
        ae = np.array([abs(x[4]) for x in sub])
        print(f"   {nombre}: {len(sub)};  |C_pred/C_num - 1|: mediana {np.median(ae):.2e}, maximo {ae.max():.2e};"
              f"  <= 1 %: {int((ae <= 0.01).sum())},  <= 3 %: {int((ae <= 0.03).sum())},  <= 10 %: {int((ae <= 0.1).sum())}")
    for fam in ("KKK", "cuadratico", "doble", "gauss", "gauss doble"):
        sub = [x for x in rapidos if x[0] == fam]
        if sub:
            ae = np.array([abs(x[4]) for x in sub])
            peor = max(sub, key=lambda x: abs(x[4]))
            print(f"     {fam:12s}: {len(sub):3d} rapidos;  mediana {np.median(ae):.2e};  peor {peor[4]:+.3f}"
                  f" (C_num = {peor[2]:.4f}, C_pred = {peor[3]:.4f}, a r0 = {peor[1]:.3g})")
    # dependencia del error con a r0 (rapidos): ajuste |e| ~ k (a r0)^-p
    if len(rapidos) > 5:
        x = np.log([r[1] for r in rapidos])
        y = np.log([max(abs(r[4]), 1e-6) for r in rapidos])
        p, k = np.polyfit(x, y, 1)
        print(f"     ajuste |C_pred/C_num - 1| ~ {math.exp(k):.2f} (a r0)^({p:.2f})")

# rama de escalarizacion (F'' < 0, zeta = 0): C = 64 F''^2
print("\nRama de escalarizacion con zeta = 0 (F'' < 0): C = 64 F''(phi0)^2 r0^-4")
for fila in salida["filas"]:
    if fila["zeta"] == 0.0 and fila["F2"] < 0:
        Cesc = 64 * fila["F2"]**2
        cerca = min(fila["C_num"], key=lambda c: abs(c - Cesc)) if fila["C_num"] else None
        print(f"   {fila['nombre']:28s} F'' = {fila['F2']:+.4f}  64F''^2 = {Cesc:7.4f}   C numericos = "
              f"{[round(c, 4) for c in fila['C_num']]}   complejas predichas = "
              f"{[f'{w[0]:.4f}{w[1]:+.4f}i' for w in fila['W_complejas']]}")

# ---------------------------------------------------------------- trayectorias en zeta
print("\nTrayectorias en zeta (trayectorias_robustas.json): C numerico / C zona cercana")
for o in json.load(open("resultados/trayectorias_robustas.json", encoding="utf-8")):
    if o["nombre"].startswith("KKK 0.0029"):
        continue
    f, u, F0, F2, a = fisicos(par[o["nombre"]])
    reales, complejas = prediccion(u, F0, F2, o["zeta"])
    Cn = sorted(-w / a**2 for w in o["w"])
    salida["trayectorias"].append(dict(nombre=o["nombre"], zeta=o["zeta"], a=a, C_num=Cn, C_pred=reales,
                                       W_complejas=[[w.real, w.imag] for w in complejas]))
    print(f"   {o['nombre']:20s} zeta={o['zeta']:+g}:  num {[round(c, 4) for c in Cn]}   pred {[round(c, 4) for c in reales]}"
          + (f"   complejas {[f'{w.real:.3f}{w.imag:+.3f}i' for w in complejas]}" if complejas else ""))
zc = {}
for nom in ("KKK referencia", "cuadratico 0.3087", "cuadratico 0.3465"):
    f, u, F0, F2, a = fisicos(nom)
    zc[nom] = 2 * (1 - F0) / (1 + u * u / 4)
    print(f"   zeta_c = 2(1 - F0)/(1 + u^2/4) para {nom}: {zc[nom]:.4f}")
salida["zeta_c"] = zc

# ---------------------------------------------------------------- modos complejos conocidos
print("\nModos complejos conocidos: W = -w/a^2 numerico frente a zona cercana")
comp = []
for o in json.load(open("resultados/complejos_F.json", encoding="utf-8")):
    nom, z = o["caso"].rsplit(" zeta=", 1)
    for r_ in o["raices"]:
        comp.append((nom, float(z), complex(-r_["w_a2"][0], -r_["w_a2"][1])))
for o in json.load(open("resultados/modos_complejos_F.json", encoding="utf-8")):
    if o["convergido"]:
        comp.append((o["nombre"], o["zeta"], complex(-o["w_a2"][0], -o["w_a2"][1])))
vistos = set()
for nom, z, Wn in comp:
    clave = (nom, z, round(Wn.real, 6), round(abs(Wn.imag), 6))
    if clave in vistos:
        continue
    vistos.add(clave)
    f, u, F0, F2, a = fisicos(nom)
    reales, complejas = prediccion(u, F0, F2, z)
    cands = complejas + [complex(c) for c in reales]
    Wp = min(cands, key=lambda w: abs(w - Wn) if Wn.imag * w.imag >= 0 or w.imag == 0 else abs(w.conjugate() - Wn))
    Wp = Wp if Wn.imag * Wp.imag >= 0 else Wp.conjugate()
    salida["complejos"].append(dict(nombre=nom, zeta=z, a=a, W_num=[Wn.real, Wn.imag], W_pred=[Wp.real, Wp.imag]))
    print(f"   {nom:22s} zeta={z:g} (a r0 = {a:.3g}):  num {Wn.real:+.4f}{Wn.imag:+.4f}i   pred {Wp.real:+.4f}{Wp.imag:+.4f}i"
          f"   |dW|/|W| = {abs(Wp - Wn) / abs(Wn):.3f}")

# ---------------------------------------------------------------- prediccion y verificacion
print("\nPrediccion y verificacion (secante sobre Dc desde la prediccion):")
for nom, z in (("KKK referencia", 1.8), ("cuadratico 0.3087", 2.0), ("KKK 0.0035 1.0571", 1.0),
               ("KKK 0.0141 1.2186", 1.0), ("KKK 0.0333 1.5501", 1.0), ("cuadratico 0.4921", 1.0),
               ("cuadratico 0.4779", 1.0)):
    f, u, F0, F2, a = fisicos(nom)
    f.resolver()
    reales, complejas = prediccion(u, F0, F2, z)
    capa = CapaF(f, zeta=z)
    for Wp in complejas:
        if Wp.imag > 0:
            continue
        w0 = -a * a * Wp
        res = secante(f, capa, w0)
        w = complex(res[0]) if isinstance(res, tuple) else complex(res)
        Dw = abs(Dc(f, capa, w))
        Dv = abs(Dc(f, capa, w * (1 + 1e-3)))
        salida["prediccion"].append(dict(nombre=nom, zeta=z, a=a, W_pred=[Wp.real, Wp.imag],
                                         W_num=[(-w / a**2).real, (-w / a**2).imag], Dc=Dw, Dc_vecino=Dv))
        print(f"   {nom} zeta={z}: predicho W = {Wp.real:+.4f}{Wp.imag:+.4f}i  ->  secante W = {(-w / a**2).real:+.4f}"
              f"{(-w / a**2).imag:+.4f}i   |Dc| = {Dw:.2e} (vecino {Dv:.2e});"
              f"  omega_inf r0 = {complex(np.sqrt(w + 0j) * math.exp(-f.Phi_inf)):.4f}")

json.dump(salida, open("resultados/comparar_ley_analitica.json", "w", encoding="utf-8"), indent=1)
print("\nGuardado: resultados/comparar_ley_analitica.json")
