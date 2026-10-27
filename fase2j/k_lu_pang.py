# -*- coding: utf-8 -*-
"""
Fase 2j: el coeficiente cinetico K del monopolo (reduccion_q.py) en chorros genericos de la capa de masas de la teoria
quintica (para descartar que se anule por construccion) y sobre el agujero negro de Lu-Pang con 40 cifras (mpmath sobre
las expresiones exactas de bck_semilla.py), junto con P.
Uso: python k_lu_pang.py
"""
import os
import sys
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_q import Phi, rr, PP, Ee, Xs, D, chorro_en_capa      # noqa: E402

R = pickle.load(open(os.path.join(AQUI, "resultados", "reduccion_q.pkl"), "rb"))
q = pickle.load(open(os.path.join(AQUI, "resultados", "q_simbolico.pkl"), "rb"))
horn = {s: 0 for e in q["fondo_eq"].values() for s in e.free_symbols if s.name.startswith(("F4_", "F5_"))}
fe = {k: sp.expand(v.xreplace(horn)) for k, v in q["fondo_eq"].items()}
g_, a1_, a0_ = R["extra"]
RES = {}

# 1. chorros genericos
gen = []
for sem in range(4):
    v = chorro_en_capa(fe, 500 + sem, orden=5)
    sub = {}
    for simb, expr in ((g_, R["gam"]), (a1_, R["al1"]), (a0_, R["al0"])):
        e = expr
        for k in range(4):
            sub[simb[k]] = v.ev(e)
            e = D(e)
    gen.append(float(v.ev(R["K"].xreplace(sub))))
RES["K_generico"] = gen
print("K en chorros genericos de la capa de masas:", gen, flush=True)

# 2. Lu-Pang con 40 cifras (mpmath sobre la expresion de bck_semilla.py: K del problema reducido en el fondo exacto)
import mpmath as mp            # noqa: E402
mp.mp.dps = 40
src = open(os.path.join(AQUI, "bck_semilla.py"), encoding="utf-8").read().split("RES = []")[0]
ns = {"__file__": os.path.join(AQUI, "bck_semilla.py"), "__name__": "bck_aux"}
exec(compile(src, "bck_semilla_setup", "exec"), ns)
lp = []
for a_, m_, r_ in (("0.1", "1", "2.7"), ("0.5", "1", "4.1"), ("0.9", "1", "7.7"), ("0.5", "2", "9")):
    Kv = ns["f_K"](mp.mpf(r_), mp.mpf(a_), mp.mpf(m_))
    Pv = ns["f_P"](mp.mpf(r_), mp.mpf(a_), mp.mpf(m_))
    lp.append(dict(alpha=float(a_), M=float(m_), r=float(r_), K=float(Kv), P=float(Pv)))
    print(f"Lu-Pang al {a_} M {m_} r {r_} (40 cifras): K = {mp.nstr(Kv, 5)}, P = {mp.nstr(Pv, 12)}", flush=True)
RES["K_lu_pang_40"] = lp
json.dump(RES, open(os.path.join(AQUI, "resultados", "k_lu_pang.json"), "w"), indent=1)
