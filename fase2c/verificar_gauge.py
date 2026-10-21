# -*- coding: utf-8 -*-
"""
Comprobacion del modo gauge residual del gauge de longitud propia, xi^l = f(t), xi^t = f'(t) int e^{-2 Phi}:
  g(l) = (Phi'' - w e^{-2 Phi}, r', r'', phi', phi'')
debe ser solucion exacta del sistema lineal (estabilidad_F.derivadas) y cumplir las dos ligaduras para todo w.
Se comprueba sobre una mitad regular (lado atractivo de un fondo cuadratico de Antoniou), integrando g desde l1 y
comparando con la expresion analitica en l2, para w real y complejo.   Uso: python verificar_gauge.py
"""
import os
import sys
import json
import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fondo_doble import Lado                                  # noqa: E402
from estabilidad_lisa import derivadas, ligaduras, fondo_en  # noqa: E402


class Mitad:
    """Envoltorio de una mitad para estabilidad_F (solo l >= 0)."""
    def __init__(self, lado):
        self.Fd, self.sol = lado.Fd, lado.sol


def g(f, l, w):
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = fondo_en(f, l)
    return np.array([P2 - w * e2, r1, r2, p1, p2], dtype=complex if isinstance(w, complex) else float)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    lado = Lado("cuadratico", 0.1, 0.1, 0.32989, l_max=300, rtol=1e-12).resolver()
    f = Mitad(lado)
    res = []
    for w in (-3.7, 12.5, -40.0 + 7.0j, 0.0):
        for l1, l2 in ((1e-3, 0.5), (0.3, 5.0), (2.0, 40.0)):
            y1 = g(f, l1, w)
            s = solve_ivp(lambda x, y: derivadas(f, x, y, w), (l1, l2), y1, method="DOP853", rtol=1e-12, atol=1e-14)
            y2 = s.y[:, -1]
            ex = g(f, l2, w)
            err = np.max(np.abs(y2 - ex)) / np.max(np.abs(ex))
            lig = max(np.max(np.abs(ligaduras(f, l, g(f, l, w), w))) for l in np.linspace(l1, l2, 9))
            res.append(dict(w=str(w), l1=l1, l2=l2, err_rel=float(err), ligaduras=float(lig)))
            print(f"w = {w}: l {l1} -> {l2}: |y - g|/|g| = {err:.2e}, max |ligaduras(g)| = {lig:.2e}")
    os.makedirs("resultados", exist_ok=True)
    json.dump(res, open("resultados/verificar_gauge.json", "w", encoding="utf-8"), indent=1)
