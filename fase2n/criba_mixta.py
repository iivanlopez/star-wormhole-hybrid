# -*- coding: utf-8 -*-
"""
Fase 2n: criba de la clase que deja el lema de la cola (fase2n_resultados.md, §3): nucleo de Ellis (densidad
efectiva negativa en la garganta) y cola de energia positiva con densidad decreciente.

Masa de Misner-Sharp del corte, m(r) = 1/(2r) + M (1 - r^{-n})^2 (r0 = 1):  m(1) = 1/2, m_r(1) = -1/2 (Ellis en la
garganta, r''(0) = 1), y lejos m_r ~ 2 M n r^{-n-1} > 0 con densidad e = m_r/(4 pi r^2) decreciente si n < 1.
Flujo, A4, lapso y criterios: los de fase2k/criba_masa.py (se reutiliza su evalua cambiando la funcion de masa).
Uso: python criba_mixta.py [n_fondos]
"""
import os
import sys
import json
import time
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "..", "fase2k"))
import criba_masa as cm                                  # noqa: E402


def masa_mixta(M, n):
    def m(r):
        return 0.5 / r + M * (1 - r ** (-n)) ** 2

    def mr(r):
        return -0.5 / r ** 2 + 2 * M * n * (1 - r ** (-n)) * r ** (-n - 1)

    def mrr(r):
        return 1 / r ** 3 + 2 * M * n * (n * r ** (-2 * n - 2) - (n + 1) * (1 - r ** (-n)) * r ** (-n - 2))
    return m, mr, mrr


def geometria_ok(M, n):
    m = masa_mixta(M, n)[0]
    rr = np.geomspace(1.0001, 1e4, 4000)
    return bool(np.all(1 - 2 * m(rr) / rr > 0))


cm.masa = masa_mixta

if __name__ == "__main__":
    nf = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    rng = np.random.default_rng(31)
    RES = []
    t0 = time.time()
    while len(RES) < nf:
        M = float(np.exp(rng.uniform(np.log(0.05), np.log(2.0))))
        nn = float(rng.uniform(0.2, 0.95))
        if not geometria_ok(M, nn):
            continue
        prm = dict(M=M, n=nn, MN=float(max(M, 0.05) * rng.uniform(1.0, 3.0)),
                   w=float(np.exp(rng.uniform(np.log(0.8), np.log(5)))),
                   a=float(np.exp(rng.uniform(np.log(0.3), np.log(300)))),
                   wq=float(np.exp(rng.uniform(np.log(0.1), np.log(3)))),
                   wn=0.0,
                   w2=float(np.exp(rng.uniform(np.log(0.08), np.log(0.6)))), lam=float(rng.uniform(0, 0.3)))
        prm["wn"] = float(prm["MN"] / rng.uniform(0.1, 0.9))    # Nc(0) en (0.1, 0.9)
        try:
            c = cm.evalua(prm)
        except Exception as ex_:                                # noqa: BLE001
            c = dict(prm, B=False, error=str(ex_))
        RES.append(c)
        if c.get("B"):
            print(f"[{time.time() - t0:6.0f} s] " + str({k_: (float(f"{v_:.4g}") if isinstance(v_, float) else v_)
                                                       for k_, v_ in c.items()}), flush=True)
    nB = sum(1 for c in RES if c.get("B"))
    print(f"{nf} fondos: beta > 0 en {nB}; hiperbolicos {sum(1 for c in RES if c.get('H'))}; cinetica positiva en "
          f"todo rho {sum(1 for c in RES if c.get('cin_pos_todo'))}; cinetica positiva lejos "
          f"{sum(1 for c in RES if c.get('lejos_cin'))}; (Q r')' > 0 lejos {sum(1 for c in RES if c.get('QRp_lejos'))}; "
          f"todo a la vez {sum(1 for c in RES if c.get('H') and c.get('cin_pos_todo') and c.get('lapso_ok'))}",
          flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", "criba_mixta.json"), "w"), indent=1)
