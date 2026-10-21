# -*- coding: utf-8 -*-
"""
Continuacion de la resonancia inestable compleja desde el fondo de referencia (alpha = 0.03, phi'(0+) = -0.7, zeta = 2,
w = 1149 + 171 i) hasta un fondo en el que la capa de KKK (zeta = 1) no tiene raiz real: alpha = 0.003, phi'(0+) = -3
(alpha/r0^2 = 0.0035, f0 = 1.057).
Camino s en [0, 2]: s <= 1 mueve el fondo con zeta = 2 (alpha = 0.03 * 0.1^s, phi' = -0.7 - 2.3 s); s > 1 baja zeta de 2 a 1.
Paso adaptativo (se divide por 2 si la secante no converge) con prediccion lineal desde las dos ultimas raices.
Salida: resultados/continuacion_compleja.json.   Uso: python continuacion_compleja.py
"""
import sys
import json
import math
import cmath
from kkk import FondoKKK
from estabilidad_kkk import CapaKKK
from modos_complejos import secante

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def punto(s):
    if s <= 1.0:
        return 0.03 * 0.1**s, -0.7 - 2.3 * s, 2.0
    return 0.003, -3.0, 2.0 - (s - 1.0)


def resolver(s, w0):
    al, u0, z = punto(s)
    f = FondoKKK(alpha=al, phi1=u0).resolver()
    c = CapaKKK(f, zeta=z)
    w, paso, dD = secante(f, c, w0)
    ok = paso < 1e-8 * abs(w) and w.imag > 0 and abs(w - w0) < 0.5 * abs(w0)
    return ok, w, f, (al, u0, z)


s, ds = 0.0, 0.125
hist = [(0.0, complex(1149.06, 170.71))]
res = []
while s < 2.0 - 1e-12 and ds > 1e-3:
    s_new = min(2.0, s + ds)
    if len(hist) >= 2:
        (s1, w1), (s2, w2) = hist[-2], hist[-1]
        pred = w2 + (w2 - w1) * (s_new - s2) / (s2 - s1)
    else:
        pred = hist[-1][1]
    ok, w, f, (al, u0, z) = resolver(s_new, pred)
    if not ok:
        ds /= 2
        print(f"  s={s_new:.4f}: sin convergencia, paso -> {ds:.4f}", flush=True)
        continue
    s = s_new
    hist.append((s, w))
    lam = math.exp(f.Phi_inf)
    om = cmath.sqrt(w)
    om = om if om.imag < 0 else -om
    o = dict(s=s, alpha=al, phi1=u0, zeta=z, alpha_kkk=f.alpha_kkk, f0=f.f0, w=[w.real, w.imag],
             omega_inf=[om.real / lam, om.imag / lam])
    res.append(o)
    print(f"s={s:.4f} alpha_kkk={f.alpha_kkk:.4f} f0={f.f0:.4f} zeta={z:.3f}: w = {w.real:.6g}{w.imag:+.6g}j; "
          f"omega_inf r0 = {om.real / lam:+.4f}{om.imag / lam:+.4f}j", flush=True)
    ds = min(0.25, ds * 1.5)
json.dump(res, open("resultados/continuacion_compleja.json", "w", encoding="utf-8"), indent=1)
print("final s =", s)
