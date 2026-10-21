# -*- coding: utf-8 -*-
"""
Principio del argumento para el Modelo L (busqueda de autovalores complejos: sobreestabilidad).
D(w) es entera en w: los coeficientes del sistema son lineales en w y el intervalo es finito.  El
reescalado del metodo compuesto multiplica D por un factor real positivo, que no altera arg D.  Por tanto
    N = (1/2 pi) * variacion de arg D a lo largo del borde de un rectangulo
es el numero de autovalores encerrados (con multiplicidad).  Si N coincide con el numero de raices reales
de la misma region (cambios de signo de D sobre el eje real), no hay autovalores complejos en ella.

Region (w interno, reloj de la garganta):  Re w en [-1.5 |w_g|, B],  |Im w| < C = min(|w_g|, C_MAX),
con w_g el modo de la garganta (resultados/modo_rapido.json) y B entre el modo fundamental del fluido y el
primer sobretono.  Se divide en tres rectangulos:
    garganta  [-A, -d] x [-C, C],   lenta  [-d, B] x [-d, d],   franja  [-d, B] x [d, C],
con d entre la escala lenta y |w_g|.  La franja conjugada tiene el mismo numero de ceros, porque
D(conj w) = conj D(w); el total es N_garganta + N_lenta + 2 N_franja.
Por que no un unico rectangulo: el lado derecho "ve" el cumulo de sobretonos del fluido en el eje real
positivo y la fase gira muy deprisa cerca de las esquinas.  Con un muestreo agrupado hacia el eje, un salto
de casi 2 pi entre muestras pasaba por pequeno y el recuento salia +2 (por parejas, por la simetria).

Muestreo inicial hibrido.  arg D sigue, con un error < 1 %, a la fase de
    f_m(w) = sinh(T s)/(T s),   s = sqrt(-w),   T = int_0^ls e^-Phi dl / c_s
(T: tiempo de cruce sonico en el reloj de la garganta; los ceros de f_m, (n pi/T)^2, son el cumulo de
sobretonos).  Con xi0 = 0.1, 1e14 (T = 178) el lado derecho del rectangulo de la garganta gira -335 rad y
cada lado horizontal +143 rad, tanto mas deprisa cuanto mas cerca del origen: con 48 puntos uniformes por
lado quedaban saltos de casi 2 pi entre muestras que pasaban por pequenos (se perdian 7 vueltas en cada
lado horizontal y se ganaban 8 en el derecho: N_garganta = -5).  Por eso cada lado empieza con la union de
N_LADO puntos uniformes y de los puntos que avanzan DFASE rad en arg f_m.
Muestreo adaptativo: se biseca cada tramo cuyo salto de arg D supera SALTO rad.
Robustez: cada caso se repite con el doble de puntos iniciales (2 N_LADO, DFASE/2; los valores ya
calculados se reutilizan); el recuento de cada rectangulo debe ser el mismo entero no negativo.
Salida: resultados/contorno.json.
Uso: python contorno.py [XI0 RHO] [--cache ARCHIVO]   (tras test_h1.py y modo_rapido.py).  Con XI0 RHO se
recalcula solo ese caso y se sustituye en el json; con --cache los valores de D se guardan en ARCHIVO
(pickle) y un calculo interrumpido se reanuda.
"""
import os
import sys
import json
import math
import time
import pickle
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.integrate import quad
from hibrido import Politropo, Fondo, RHO_GEOM

eos = Politropo(k=0.1, gamma=2.0)
SEL = [(0.6, 1e15), (0.6, 1e14), (0.1, 1e15), (0.1, 1e14)]
N_LADO, DFASE, SALTO, RONDAS, C_MAX = 48, 0.3, 0.4, 24, 2.0
_fondos = {}


def fondo(xi0, rho):
    if (xi0, rho) not in _fondos:
        rb = rho * RHO_GEOM
        _fondos[(xi0, rho)] = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    return _fondos[(xi0, rho)]


def evaluar(args):
    xi0, rho, w = args
    return fondo(xi0, rho).det("ligado", w)


def borde(t, a, b, c0, c1):
    """t en [0, 4): rectangulo [a, b] x [c0, c1] en sentido antihorario; t = k + u recorre el lado k."""
    k, u = int(t), t - int(t)
    if k == 0:
        return complex(a + (b - a) * u, c0)
    if k == 1:
        return complex(b, c0 + (c1 - c0) * u)
    if k == 2:
        return complex(b + (a - b) * u, c1)
    return complex(a, c1 + (c0 - c1) * u)


def region(xi0, rho, h1, rap):
    k = math.exp(2.0 * fondo(xi0, rho).Phi_inf)
    wg = k * min(next(o for o in rap if o["modelo"] == "ligado" and o["xi0"] == xi0 and o["rho"] == rho)["raices"])
    lig = [k * w for w in next(o for o in h1 if o.get("xi0") == xi0 and o.get("rho") == rho)["ligado"]]
    pos = sorted(w for w in lig if w > 0)
    B = math.sqrt(pos[0] * pos[1]) if len(pos) >= 2 else (2.0 * pos[0] if pos else 1e-3)
    lentos = [abs(w) for w in lig if abs(w) < 0.05 * abs(wg)]
    d = math.sqrt(abs(wg) * max(lentos + [B]))
    return dict(xi0=xi0, rho=rho, A=1.5 * abs(wg), B=B, C=min(abs(wg), C_MAX), d=d, k=k, w_g=wg)


def tiempo_sonico(f):
    """T = int_0^ls e^-Phi dl / c_s, con c_s^2 = Gamma_1 P / (eps + P) (reloj de la garganta).
    c_s ~ (ls - l)^1/2 en la superficie: se integra en y, con l = ls - y^2."""
    def g(l):
        Phi = f.y_bg(l)[2]
        eps, P, _, _ = eos.termo(max(f.theta(Phi), 1e-300), f.rho_bc, f.sig)
        return math.exp(-Phi) * math.sqrt((eps + P) / (eos.Gam * P))
    return quad(lambda y: 2.0 * y * g(f.ls - y * y), 0.0, math.sqrt(f.ls), limit=400)[0]


def fase_modelo(T, a, b, c0, c1, m=20000):
    """Para cada lado, u en [0, 1] y d(arg f_m)/du; d log f_m/dw = -(T coth(T s) - 1/s)/(2 s)."""
    esq = [complex(a, c0), complex(b, c0), complex(b, c1), complex(a, c1), complex(a, c0)]
    u = np.linspace(0.0, 1.0, m + 1)
    lados = []
    for k in range(4):
        dw = esq[k + 1] - esq[k]
        s = np.sqrt(-(esq[k] + dw * u))                  # rama principal: Re s >= 0, |e| <= 1
        e = np.exp(-2.0 * T * s)
        lados.append((u, (-(T * (1.0 + e) / (1.0 - e) - 1.0 / s) / (2.0 * s) * dw).imag))
    return lados


class Cache(dict):
    """D(w) ya calculados, por (rectangulo, t).  Con ruta se guardan en disco y se recuperan al reanudar."""

    def __init__(self, ruta=None):
        super().__init__()
        self.ruta = ruta
        if ruta and os.path.exists(ruta):
            self.update(pickle.load(open(ruta, "rb")))

    def guardar(self):
        if self.ruta:
            pickle.dump(dict(self), open(self.ruta + ".tmp", "wb"))
            os.replace(self.ruta + ".tmp", self.ruta)


def contar(ex, xi0, rho, rect, T, n_lado, dfase, cache, bloque=400):
    lados = fase_modelo(T, *rect)
    ini = set(np.linspace(0.0, 4.0, 4 * n_lado, endpoint=False).tolist())   # incluye las esquinas
    for k, (u, v) in enumerate(lados):
        V = np.concatenate([[0.0], np.cumsum(0.5 * (np.abs(v[1:]) + np.abs(v[:-1])) * np.diff(u))])
        ini |= set((k + np.interp(dfase * np.arange(1, math.ceil(V[-1] / dfase)), V, u)).tolist())
    D, pend = {}, sorted(ini)
    for _ in range(RONDAS):
        if not pend:
            break
        nuevos = [t for t in pend if (rect, t) not in cache]
        for i in range(0, len(nuevos), bloque):
            for t, val in zip(nuevos[i:i + bloque], ex.map(evaluar, [(xi0, rho, borde(t, *rect))
                                                                     for t in nuevos[i:i + bloque]])):
                cache[(rect, t)] = complex(val)
            cache.guardar()
        D.update((t, cache[(rect, t)]) for t in pend)
        ts = sorted(D)
        pend = [0.5 * (ta + tb) % 4.0 for ta, tb in zip(ts, ts[1:] + [ts[0] + 4.0])
                if abs(np.angle(D[tb % 4.0] / D[ta])) > SALTO]
    ts = sorted(D)
    dphi = [np.angle(D[tb % 4.0] / D[ta]) for ta, tb in zip(ts, ts[1:] + [ts[0] + 4.0])]
    fase = [[sum(x for ta, x in zip(ts, dphi) if int(ta) == k), float(np.sum(0.5 * (v[1:] + v[:-1]) * np.diff(u)))]
            for k, (u, v) in enumerate(lados)]                 # variacion de arg D y de arg f_m por lado
    return dict(N=sum(dphi) / (2.0 * math.pi), puntos=len(ts), iniciales=len(ini),
                salto_max=max(abs(x) for x in dphi), convergido=not pend, fase_lados=fase)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    t0 = time.time()
    h1 = json.load(open("resultados/h1.json", encoding="utf-8"))
    rap = json.load(open("resultados/modo_rapido.json", encoding="utf-8"))
    args, ruta = sys.argv[1:], None
    if "--cache" in args:
        i = args.index("--cache")
        ruta, args[i:i + 2] = args[i + 1], []
    casos, res, cache = SEL, [], Cache(ruta)
    if len(args) == 2:
        casos = [(float(args[0]), float(args[1]))]
        try:
            res = json.load(open("resultados/contorno.json", encoding="utf-8"))
        except FileNotFoundError:
            pass
    with ProcessPoolExecutor(max_workers=8) as ex:
        for x, r in casos:
            c = region(x, r, h1, rap)
            A, B, C, d = c["A"], c["B"], c["C"], c["d"]
            g = np.concatenate([-np.logspace(math.log10(A), -6.0, 120), np.linspace(B / 60, B, 60)])
            v = np.sign(list(ex.map(evaluar, [(x, r, float(w)) for w in g])))
            c["reales_en"] = [float(g[j]) for j in np.nonzero(v[:-1] * v[1:] < 0)[0]]
            c["reales"] = len(c["reales_en"])
            c["T"] = T = tiempo_sonico(fondo(x, r))
            c["muestreo"] = dict(n_lado=N_LADO, dfase=DFASE, salto=SALTO)
            rects = dict(garganta=(-A, -d, -C, C), lenta=(-d, B, -d, d), franja=(-d, B, d, C))
            doble = {}
            for o, n, df in ((c, N_LADO, DFASE), (doble, 2 * N_LADO, 0.5 * DFASE)):
                for s, R in rects.items():
                    o[s] = contar(ex, x, r, R, T, n, df, cache)
                    print(f"   {'doble' if o is doble else 'base '} {s}: N = {o[s]['N']:.3f}, puntos = {o[s]['puntos']}"
                          f" (iniciales {o[s]['iniciales']}), salto_max = {o[s]['salto_max']:.3f}  "
                          f"({time.time() - t0:.0f} s)", flush=True)
                o["N"] = o["garganta"]["N"] + o["lenta"]["N"] + 2.0 * o["franja"]["N"]
                o["puntos"] = sum(o[s]["puntos"] for s in rects)
                o["convergido"] = all(o[s]["convergido"] for s in rects)
            c["doble"] = doble
            c["robusto"] = all(round(c[s]["N"]) == round(doble[s]["N"]) >= 0 and
                               max(abs(o[s]["N"] - round(o[s]["N"])) for o in (c, doble)) < 0.01 for s in rects)
            i = next((j for j, o in enumerate(res) if (o["xi0"], o["rho"]) == (x, r)), len(res))
            res[i:i + 1] = [c]
            print(f"xi0={x} rho={r:.0e}: Re w en [{-A:.4g}, {B:.3e}], |Im w| < {C:.4g}, d = {d:.3g}, T = {T:.1f} "
                  f"(w interno) -> N = {c['garganta']['N']:.3f} + {c['lenta']['N']:.3f} + 2 x {c['franja']['N']:.3f}"
                  f" = {c['N']:.3f}; doble: {doble['garganta']['N']:.3f} + {doble['lenta']['N']:.3f} + 2 x "
                  f"{doble['franja']['N']:.3f} = {doble['N']:.3f}; raices reales = {c['reales']} en "
                  f"{[f'{y:.3e}' for y in c['reales_en']]}; puntos = {c['puntos']} / {doble['puntos']}; "
                  f"convergido = {c['convergido'] and doble['convergido']}; robusto = {c['robusto']}  "
                  f"({time.time() - t0:.0f} s)", flush=True)
            for s in rects:
                print(f"   {s}: arg D / arg f_m por lado (rad) = "
                      f"{', '.join(f'{a:.2f} / {m:.2f}' for a, m in c[s]['fase_lados'])}", flush=True)
            json.dump(res, open("resultados/contorno.json", "w", encoding="utf-8"), indent=1)
    print(f"tiempo: {time.time() - t0:.0f} s")
