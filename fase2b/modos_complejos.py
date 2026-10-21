# -*- coding: utf-8 -*-
"""
Modos inestables complejos (sobreestabilidad) de los agujeros de KKK con capa dinamica.
Con kappa = sqrt(-w e^{-2 Phi_inf}) en la rama principal (Re kappa > 0, corte en w > 0), un cero de D(w) fuera del
corte es un modo que crece en el tiempo (omega = -i kappa e^{Phi_inf}, Im omega < 0), real o complejo.  Los modos
estables (cuasinormales) estan en la otra hoja.  D(w̄) = conj D(w): basta el semiplano superior.
Dato inicial analitico y sin polos: el vector nulo de las 4 condiciones en la garganta dado por sus menores 4x4 con
signo (producto vectorial generalizado), entero en w.
Recuento: principio del argumento en rectangulos con muestreo uniforme y biseccion (umbral 0.4 rad).
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from estabilidad_kkk import matriz_garganta
from perturb_kkk import derivadas


def vector_menores(f, capa, w):
    K = matriz_garganta_c(f, capa, w)[[0, 1, 3, 4]]
    return np.array([(-1) ** k * np.linalg.det(np.delete(K, k, axis=1)) for k in range(5)])


def matriz_garganta_c(f, capa, w):
    """matriz_garganta admite w complejo si se construye con columnas complejas."""
    return matriz_garganta(f, capa, w)


def Dc(f, capa, w, l_extra=30.0, rtol=1e-10):
    w = complex(w)
    u = vector_menores(f, capa, w).astype(complex)
    kap = np.sqrt(-w * math.exp(-2 * f.Phi_inf))
    if kap.real < 0:
        kap = -kap
    l1 = min(10.0 * f.r0 + l_extra / max(abs(kap), 1e-3), 0.9 * f.l_end)
    sol = solve_ivp(lambda l, y: derivadas(f, l, y, w), (0.0, l1), u, method="DOP853", rtol=rtol, atol=1e-14)
    a1, d0, d1, q0, q1 = sol.y[:, -1]
    r, s = f.sol.sol(l1)[:2]
    return (q1 + (kap + s / r) * q0) * r * np.exp(-kap * l1)


def contar(f, capa, a, b, c0, c1, n=40, salto=0.4, rondas=20, pool=None):
    """Numero de ceros de Dc en [a, b] x [c0, c1] (principio del argumento)."""
    def borde(t):
        k, u = int(t), t - int(t)
        if k == 0:
            return complex(a + (b - a) * u, c0)
        if k == 1:
            return complex(b, c0 + (c1 - c0) * u)
        if k == 2:
            return complex(b + (a - b) * u, c1)
        return complex(a, c1 + (c0 - c1) * u)
    ev = (lambda ws: list(pool.map(_eval, [(f.alpha, f.phi1, capa.zeta, capa.mu, w) for w in ws]))) if pool else \
        (lambda ws: [Dc(f, capa, w) for w in ws])
    D, pend = {}, [float(t) for t in np.linspace(0, 4, 4 * n, endpoint=False)]
    for _ in range(rondas):
        if not pend:
            break
        for t, v in zip(pend, ev([borde(t) for t in pend])):
            D[t] = complex(v)
        ts = sorted(D)
        pend = [0.5 * (ta + tb) % 4.0 for ta, tb in zip(ts, ts[1:] + [ts[0] + 4.0])
                if abs(np.angle(D[tb % 4.0] / D[ta])) > salto]
    ts = sorted(D)
    dphi = [np.angle(D[tb % 4.0] / D[ta]) for ta, tb in zip(ts, ts[1:] + [ts[0] + 4.0])]
    return sum(dphi) / (2 * math.pi), len(ts), not pend


_cache = {}


def _eval(args):
    from kkk import FondoKKK
    from estabilidad_kkk import CapaKKK
    al, u0, zeta, mu, w = args
    key = (al, u0, zeta, mu)
    if key not in _cache:
        f = FondoKKK(alpha=al, phi1=u0).resolver()
        _cache[key] = (f, CapaKKK(f, zeta=zeta, mu=mu))
    f, c = _cache[key]
    return Dc(f, c, w)


def secante(f, capa, w0, it=40):
    """Raiz compleja de Dc por secante desde w0 (paso limitado al 30 % de |w|)."""
    w0 = complex(w0)
    w1 = w0 * (1 + 1e-3) + 1e-2j
    D0, D1 = Dc(f, capa, w0), Dc(f, capa, w1)
    paso = 0.0
    for _ in range(it):
        if D1 == D0:
            break
        paso = -D1 * (w1 - w0) / (D1 - D0)
        if abs(paso) > 0.3 * abs(w1):
            paso *= 0.3 * abs(w1) / abs(paso)
        w0, D0 = w1, D1
        w1 = w1 + paso
        D1 = Dc(f, capa, w1)
        if abs(paso) < 1e-10 * abs(w1):
            break
    return w1, abs(paso), abs(D1)


if __name__ == "__main__":
    # Uso: python modos_complejos.py  -> resultados/modos_complejos.json
    import sys
    import json
    import cmath
    from kkk import FondoKKK
    from estabilidad_kkk import CapaKKK
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    f = FondoKKK(alpha=0.03, phi1=-0.7).resolver()
    res = []
    for zeta, inicios in ((2.0, (1200 + 20j, 1200 + 200j)), (1.8, (3000 + 50j, 6000 + 500j))):
        c = CapaKKK(f, zeta=zeta)
        for w0 in inicios:
            w, p, d = secante(f, c, w0)
            ok = p < 1e-8 * abs(w)
            om = -cmath.sqrt(w) if (-cmath.sqrt(w)).imag < 0 else cmath.sqrt(w)     # rama que crece: Im omega < 0
            lam = math.exp(f.Phi_inf)
            o = dict(zeta=zeta, inicio=[w0.real, w0.imag], w=[w.real, w.imag], convergido=bool(ok), D=d,
                     omega_garganta=[om.real, om.imag], omega_inf=[om.real / lam, om.imag / lam])
            res.append(o)
            print(f"zeta={zeta}: w = {w.real:.6g}{w.imag:+.6g}j (convergido {ok}); omega_inf r0 = "
                  f"{om.real / lam:+.4f}{om.imag / lam:+.4f}j  (tasa de crecimiento {-om.imag / lam:.4f} c/r0)", flush=True)
    json.dump(res, open("resultados/modos_complejos.json", "w", encoding="utf-8"), indent=1)
