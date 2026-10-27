# -*- coding: utf-8 -*-
"""
Fase 2k: la prueba d = N f de las Fases 2g-2j en las variables ADM del corte phi = cte (adm_t.py).
En el fondo estatico (qt = 0) la prueba es c = a = beta = dphi = 0, d = N f, que deja X fijo.  Con phi = qt t + psi
la misma configuracion es n = beta~ = h = 0, d = N f (n = 0 <=> dX = 0).  Aqui se evaluan, a frecuencia cero y en
chorros exactos de la capa de masas, las ecuaciones de los auxiliares E_n y E_beta~ (y E_h) sobre la prueba:
si no se anulan, la prueba no cumple las ligaduras y hay que mover el lapso y el desplazamiento.
Uso: python prueba_dN_t.py [G5]
"""
import os
import sys
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_t import Phi, rr, PP, Ee, Xs, qt, D, chorro_en_capa   # noqa: E402
import adm_t as A                                                      # noqa: E402

f = sp.symbols("f0:9", real=True)
N = sp.sqrt(Ee)
dsub = {}
e = N * f[0]
for k in range(A.NM + 3):
    dsub[A.d_[k]] = e
    e = D(e, (f,))
cero = {s: 0 for fam in (A.n_, A.B_, A.h_) for s in fam}
RES = {}
for sem in range(3):
    ch = A.chorro_en_capa(A.fe, A.J, 950 + sem, orden=5, cero=() if A.CON_G5 else ("G5",))
    En = A.ecuaciones_en(ch)
    for k in range(9):
        ch[f[k]] = 0
    ch[f[0]] = 1
    val = {k: sp.expand(v.subs(A.iw, 0).xreplace(cero).xreplace({s: ch.ev(x) for s, x in dsub.items()}))
           for k, v in En.items()}
    sig = ch.ev(rr[1] / rr[0] - Phi[1])
    RES[sem] = {k: str(v) for k, v in val.items()}
    print(f"chorro {sem}: d = N (f = 1), n = B = h = 0, w = 0:  " +
          ", ".join(f"E_{k} = {float(v):.4g}" for k, v in val.items()) + f";  sigma = {float(sig):.4g}", flush=True)
json.dump(RES, open(os.path.join(AQUI, "resultados", f"prueba_dN_{A.NOMBRE}.json"), "w"), indent=1)
