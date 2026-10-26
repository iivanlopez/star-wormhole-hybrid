# -*- coding: utf-8 -*-
"""
Fase 2h: validacion del espectro de elementos finitos (espectro2.py).
  1. Garganta de referencia de la Fase 2f (sin fluido, su propia reconstruccion): formulacion de vacio de los elementos
     finitos en todo el dominio frente a los autovalores pares de Sturm-Liouville de la Fase 2f (-0.11750 y -0.00485).
  2. Hibrido con fluido: convergencia del espectro completo con la malla (par e impar).
Uso: python validar_espectro.py
"""
import os
import sys
import json
import time
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2f"), os.path.join(AQUI, "..", "fase2g"), AQUI):
    sys.path.insert(0, p_)
import espectro as ES                 # noqa: E402  (Fase 2f, Sturm-Liouville)
import fondo as FH                    # noqa: E402
import teoria as TE                   # noqa: E402
import espectro2 as E2                # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


def KPU(T, ls, trozo=400):
    """K, P, U del monopolo sin fluido (T3) sobre los coeficientes de la teoria reconstruida."""
    K, P, U = [], [], []
    for i in range(0, ls.size, trozo):
        c, V = T.coeficientes(ls[i:i + trozo])
        N, r, r1, Gc = V["N"], V["r0"], V["r1"], V["Gc"]
        import simbolico_num as SN
        import completo as CO
        gam, al1 = CO.COEF["gam"](V), CO.COEF["al1"](V)
        kc, e0, e1 = CO.COEF["k_c"](V), CO.COEF["e0"](V), CO.COEF["e1"](V)
        A1 = al1 / gam
        A0 = -V["Phi1"] * A1
        u_d = -(2 / r) * (N * r1 * Gc).d()
        p_d = -2 * N * Gc
        Pj = p_d + A1 * e1 - (kc / 2) * A1 * A1
        Uj = u_d - (kc / 2) * A0 * A0 + A0 * e0 + 0.5 * (kc * A1 * A0 - A1 * e0 - A0 * e1).d()
        K.append(V["K"].c[0])
        P.append(Pj.c[0])
        U.append(Uj.c[0])
    return np.concatenate(K), np.concatenate(P), np.concatenate(U)


class AdaptadorF2f:
    """Coeficientes de la formulacion de elementos finitos sobre una garganta reconstruida de la Fase 2f (sin fluido)."""
    def __init__(self, g):
        self.g = g

    def coeficientes(self, l):
        import completo as CO
        import jets as J
        g = self.g
        V = g.V(l)
        cero = 0 * V["N"]
        V["Pf0"], V["Pf1"] = cero, cero
        m = g._mono(V, que=("K",))
        N, r, r1, Gc = V["N"], V["r0"], V["r1"], V["Gc"]
        gam, al1 = CO.COEF["gam"](V), CO.COEF["al1"](V)
        a2 = 4 * N * N * r * Gc
        out = dict(N=N, Phi1=V["Phi1"], r=r, r1=r1, Gc=Gc, G4=V["G4"], gam=gam, al1=al1, kc=CO.COEF["k_c"](V),
                   e0=CO.COEF["e0"](V), e1=CO.COEF["e1"](V), ett=CO.COEF["e_tt"](V), a2=a2, a2p=a2.d(),
                   Ng=N * N * gam, Ngp=(N * N * gam).d(), rp=cero, K=m["K"], kct=CO.COEF["k_c"](V))
        out = {k: (v.c[0] if isinstance(v, J.Jet) else v) for k, v in out.items()}
        out["cs2"] = 0 * out["N"]
        out["sig"] = out["r1"] / out["r"] - out["Phi1"]
        out["Np"] = out["N"] * out["Phi1"]
        return out, V


if __name__ == "__main__":
    x = json.load(open(os.path.join(AQUI, "..", "fase2g", "resultados", "criba_optimiza_2.json")))["x"]
    # ------------------------------------------------------------- 1. vacio: referencia de la Fase 2f
    import busqueda as BQ
    BASE = dict(n=2, b=4.0, frac=0.9, wdb=0.5, M=2.0, a=2.0, wK=10.0, wPK=1.5, cm2=0.8)
    ad = AdaptadorF2f(BQ.construir(BASE))
    log("referencia de la Fase 2f: autovalores pares de Sturm-Liouville -0.11750, -0.00485")
    for n_out, L in ((150, 400.0), (300, 400.0), (300, 1000.0)):
        o = E2.espectro(ad, "par", 0.0, n_in=1, L=L, n_out=n_out, k=3)
        log(f"   elementos finitos (vacio, {n_out} elementos, L = {L}): Tmin = {o['Tmin']:+.3e}, w = {o.get('w')}")
    # ------------------------------------------------------------- 2. hibrido con fluido
    F = FH.Fondo(**x, kappa=1.0, phi1=1.0, Lmax=4000.0)
    T = TE.Teoria(F)
    lc = E2.lc_fluido(T)
    log(f"hibrido: lc = {lc:.4f} (ls = {F.ls:.4f}), B = {T.B}")
    for par in ("par", "impar"):
        for n_in, n_out in ((60, 150), (120, 250), (240, 400)):
            o = E2.espectro(T, par, lc, n_in=n_in, L=1500.0, n_out=n_out, k=4)
            log(f"   {par} ({n_in}+{n_out}): Tmin = {o['Tmin']:+.3e}, w = {o.get('w')}, fantasmas {o.get('n_T_neg')}")
