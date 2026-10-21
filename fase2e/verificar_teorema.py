# -*- coding: utf-8 -*-
"""
Comprobacion numerica del teorema del coeficiente cinetico con lapso variable (teorema_chi_lapso.py):
  1. Fondos con potencial (kink, F = alpha phi^2, doble pozo): K y G de la forma cerrada frente a los coeficientes del
     polinomio de dispersion de las ecuaciones generadas (potencial.dispersion); en cada cero de
     chi = r r' + 2 F_l (1 - 3 r'^2): num ~ 0 (cruce regular), K = -16 F'^2 (1 - r'^2)^2 chi'/r'^2, y lambda < 0 en los
     ceros ascendentes.
  2. El agujero exacto de Canate et al. (fase2d/canate.py, Maxwell, Phi' = 0): la misma identidad en su cero de chi
     (es el caso sigma_m = 0 de la extension con materia).
Salida: resultados/verificar_teorema.json.   Uso: python verificar_teorema.py
"""
import os
import sys
import json
import math
import numpy as np
from scipy.optimize import brentq

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
sys.path.append(os.path.join(AQUI, "..", "fase2d"))
import potencial as P          # noqa: E402
import mapa_potencial as MP    # noqa: E402
import canate as C             # noqa: E402

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
out = {"kinks": [], "canate": []}

# ------------------------------------------------------------------ 1. fondos con potencial
for alpha, v, u in ((-0.5, 1.0, 1.0), (-0.3, 1.0, 1.5), (-1.0, 0.5, 2.0), (-2.0, 2.0, 3.0)):
    Fd = P.acoplamiento("cuadratico", alpha)
    amp = (2 + u * u / 2) / v**4
    L = P.LadoV(Fd, P.potencial("doble_pozo", amp, v=v), 0.0, u, 0.0, phi_v=v, l_max=300.0, rtol=1e-12).resolver()
    f = P.DobleV(L, L)
    err = 0.0
    for x in np.linspace(1e-3, 0.999 * L.l_end, 60):
        m = MP.magnitudes(L, x)
        c0, c1, c2 = P.dispersion(f, x)
        e2 = m["N"] ** -2
        Kd, Gd = -c1 * m["r"]**5 / (2 * e2), c0 * m["r"]**5 / 2
        err = max(err, abs(m["K"] - Kd) / (1 + abs(Kd)), abs(m["G"] - Gd) / (1 + abs(Gd)), abs(c2) / (1 + abs(c1)))
    zs = MP.ceros_chi(L, n=3000)
    ls = np.linspace(1e-3, 0.999 * L.l_end, 3000)
    lam = np.array([MP.magnitudes(L, x)["G"] * MP.magnitudes(L, x)["N"] ** 2 / MP.magnitudes(L, x)["K"] for x in ls])
    neg = ls[lam < 0]
    out["kinks"].append(dict(alpha=alpha, v=v, u=u, fin=L.fin, l_end=L.l_end, err_KG=err, ceros=zs,
                             lam_neg=[float(neg.min()), float(neg.max())] if len(neg) else None))
    print(f"1. kink alpha={alpha} v={v} u={u}: final '{L.fin}' en l={L.l_end:.3f}; K, G cerrados frente a la "
          f"dispersion: {err:.1e}")
    for z in zs:
        print(f"     cero de chi en l={z['l']:.6f} ({'sube' if z['sube'] else 'baja'}): r'={z['rp']:+.4f} "
              f"num={z['num']:+.1e} K={z['K']:+.6e} K_teo={z['K_teo']:+.6e} G={z['G']:+.4e} lambda={z['lam']:+.3e}")
    if len(neg):
        print(f"     lambda < 0 en l en [{neg.min():.4f}, {neg.max():.4f}]")

# ------------------------------------------------------------------ 2. Canate et al.


def magn_c(f, l):
    y = f.sol.sol(l)
    rpp, Ppp, ppp = f.exactas_segundas(l)
    r, rp, Phi, Pp, p, pp = y[:6]
    F0, F1, F2, F3 = f.Fd(p)
    Fl, Fll = F1 * pp, F2 * pp**2 + F1 * ppp
    A, B = r - 4 * Fl * rp, 1 - rp * rp
    chi = r * rp + 2 * Fl * (1 - 3 * rp * rp)
    dchi = rp * rp + r * rpp + 2 * Fll * (1 - 3 * rp * rp) - 12 * Fl * rp * rpp
    K = A**2 * r**2 - 16 * F1**2 * B * (4 * A * rpp + B * (1 - 4 * Fll))
    return chi, dchi, K, F1, B, rp


for Qs in (0.5, 1.0, 2.0, 5.0):
    f = C.FondoCanate(Qs)
    ls = -np.geomspace(1e-3, 50.0, 4000)                 # lado l < 0, donde F_l r' < 0 (c = 0)
    ch = [magn_c(f, l)[0] for l in ls]
    for i in range(len(ls) - 1):
        if ch[i] * ch[i + 1] < 0:
            z = brentq(lambda l: magn_c(f, l)[0], ls[i], ls[i + 1], xtol=1e-14)
            chi, dchi, K, F1, B, rp = magn_c(f, z)
            Kt = -16 * F1**2 * B**2 * dchi / rp**2
            out["canate"].append(dict(Qs=Qs, l=z, K=K, K_teo=Kt))
            print(f"2. Canate Qs/Q={Qs}: cero de chi en l={z:.6f}: K = {K:+.10e}, -16F'^2B^2chi'/r'^2 = {Kt:+.10e}")
json.dump(out, open(os.path.join(AQUI, "resultados", "verificar_teorema.json"), "w", encoding="utf-8"), indent=1)
