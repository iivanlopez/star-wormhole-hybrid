# -*- coding: utf-8 -*-
"""
Fase 2m: comprobaciones exactas de la maquinaria de Proca (resultados/proca_simbolico.pkl).
  1. Limite Einstein-Maxwell (G4 = 1, G2 = F): Reissner-Nordstrom en longitud propia cumple las ecuaciones del fondo
     (simbolico en r, M, Q), con A0 = 2Q/r (la normalizacion -F_mn^2/4 frente a -F_mn^2 de las unidades geometricas).
  2. Correspondencia con Horndeski (A0 = 0): E_A, E_B, E_C y -D_l E_w coinciden con E_A, E_B, E_C, E_p de la Fase 2f
     (fase2f/resultados/bH_simbolico.pkl) con phi' = A1, G2^H(X) = G2(X, 0), G3^H = -G3, G4^H = G4; y la accion
     cuadratica con delta A = 0 (al = w = 0) coincide con la de la Fase 2f en el gauge unitario (q = 0).
  3. Paridad temporal (A0 = 0): E_u = 0 identicamente, y en L2 ningun termino acopla {al, beta} con {a, c, d, w} con un
     numero par de derivadas temporales.  Igual para A1 = 0, G3 = 0 con la paridad TC: {w, beta} impares.
  4. Modos gauge: L_xi del fondo, xi = (T, Z) e^{i w t}, resuelve las seis ecuaciones lineales en chorros exactos de la
     capa de masas (casos general, R y E).
Uso: python verificar_proca.py
"""
import os
import sys
import json
import time
import pickle
import random
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import chorros_proca as CP                              # noqa: E402
from chorros_proca import NN, rr, UU, VV, gsym          # noqa: E402

T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


d = pickle.load(open(os.path.join(AQUI, "resultados", "proca_simbolico.pkl"), "rb"))
fe, lin, L2, iw, amp, pj = (d[k] for k in ("fondo_eq", "lineales", "L2", "iw", "amp", "pj"))
CHK = {}
ORD = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)]

# ------------------------------------------------------------------ 1. Reissner-Nordstrom
R, M, Q = sp.symbols("R M Q", positive=True)
f = 1 - 2 * M / R + Q ** 2 / R ** 2


def Dl(e):
    return sp.sqrt(f) * sp.diff(e, R)


rn = {}
eN, er, eU = sp.sqrt(f), R, 2 * Q / R
for k in range(6):
    rn[NN[k]], rn[rr[k]], rn[UU[k]], rn[VV[k]] = eN, er, eU, 0
    eN, er, eU = Dl(eN), Dl(er), Dl(eU)
Xrn = sp.simplify((UU[0] ** 2 / (2 * NN[0] ** 2)).xreplace(rn))
Frn = sp.simplify((UU[1] ** 2 / (2 * NN[0] ** 2)).xreplace(rn))
teo = {}
for s in set().union(*[e.free_symbols for e in fe.values()]):
    if CP.es_G(s):
        teo[s] = 0
teo[gsym("G4", 0, 0)] = 1
teo[gsym("G2", 0, 0)] = Frn
teo[gsym("G2", 0, 1)] = 1
res = {}
for n, e in fe.items():
    res[n] = sp.simplify(e.xreplace(teo).xreplace(rn))
CHK["RN"] = {n: str(v) for n, v in res.items()}
assert all(v == 0 for v in res.values()), res
log("Reissner-Nordstrom (G4 = 1, G2 = F, A0 = 2Q/r): las seis ecuaciones del fondo se anulan (simbolico)")

# ------------------------------------------------------------------ 2. correspondencia con Horndeski (fase 2f)
F2F = os.path.join(AQUI, "..", "fase2f")
h = pickle.load(open(os.path.join(F2F, "resultados", "bH_simbolico.pkl"), "rb"))
feH, L2H, PhiH, rH, PH, EeH = h["fondo_eq"], h["L2"], h["Phi"], h["rr"], h["PP"], h["Ee"]
XsH = h["Xs"]
rnd = random.Random(7)


def rq():
    return sp.Rational(rnd.randint(-900, 900), 997)


val = {NN[0]: sp.Rational(rnd.randint(500, 1500), 1000), rr[0]: sp.Rational(rnd.randint(600, 2000), 1000)}
for k in range(1, CP.NMAXB):
    val[NN[k]] = rq()
    val[rr[k]] = rq()
for k in range(CP.NMAXB):
    val[VV[k]] = rq()
    val[UU[k]] = 0
# funciones: G2(X, F) -> G2_i0 libres, G2_ij (j > 0) libres (no deben aparecer con U = 0 salvo en E_u, L2 con al)
Gv = {}
for (nom, i, j), s in d["GS"].items():
    Gv[s] = rq()
for i in range(8):
    for j in range(4):
        for nom in ("G2", "G3", "G4"):
            Gv.setdefault(gsym(nom, i, j), rq())
val.update(Gv)
# mapa a la Fase 2f: Phi = ln N, P_k = V_(k-1), G_0j^H
lnN = sp.log(NN[0])
valH = {EeH: val[NN[0]] ** 2, PH[0]: sp.Rational(3, 7)}
e = lnN
for k in range(6):
    if k >= 1:
        valH[PhiH[k]] = sp.nsimplify(sp.expand(e).xreplace(val))
    e = CP._D0(e)
for k in range(6):
    valH[rH[k]] = val[rr[k]]
for k in range(1, 7):
    valH[PH[k]] = val[VV[k - 1]]
valH[XsH] = -val[VV[0]] ** 2 / 2
for s in set().union(*[x.free_symbols for x in feH.values()]) | L2H.free_symbols:
    m = __import__("re").match(r"^(G2|G3|G4|F4)_(\d)(\d)$", s.name)
    if m:
        nom, i, j = m.group(1), int(m.group(2)), int(m.group(3))
        if i > 0 or nom == "F4":
            valH[s] = 0
        else:
            sg = -1 if nom == "G3" else 1
            valH[s] = sg * val[gsym(nom, j, 0)]
comp = {}
for n in ("A", "B", "C"):
    a_, b_ = sp.expand(fe[n].xreplace(val)), sp.expand(feH[n].xreplace(valH))
    comp[n] = float(abs(a_ - b_))
    assert a_ == b_, (n, a_, b_)
ep = sp.expand(-CP.D(fe["w"]).xreplace(val))
comp["p"] = float(abs(ep - sp.expand(feH["p"].xreplace(valH))))
assert ep == sp.expand(feH["p"].xreplace(valH))
comp["E_u(U=0)"] = str(sp.expand(fe["u"].xreplace({UU[k]: 0 for k in range(CP.NMAXB)})))
assert comp["E_u(U=0)"] == "0"
CHK["horndeski_fondo"] = comp
log(f"fondo con A0 = 0 frente a la Fase 2f (E_A, E_B, E_C, E_p = -D E_w): identicos; E_u = 0 identicamente")
# accion cuadratica con delta A = 0
cero_v = {pj[(v, o)]: 0 for v in ("al", "w") for o in ORD}
pjH = {s.name: s for s in L2H.free_symbols if "_" in s.name and s.name.split("_")[0] in ("a", "be", "c", "d", "q")}
L2P = sp.expand(L2.xreplace(cero_v).xreplace(val))
L2Hq = sp.expand(L2H.xreplace({s: 0 for n, s in pjH.items() if n.startswith("q_")}).xreplace(valH))
L2Hq = L2Hq.xreplace({s: pj[(n.split("_")[0], (int(n.split("_")[1][0]), int(n.split("_")[1][1])))]
                      for n, s in pjH.items() if not n.startswith("q_")})
dif = sp.expand(L2P - L2Hq)
CHK["horndeski_L2_dif"] = str(dif)
assert dif == 0, dif
log("accion cuadratica con delta A = 0 (A0 = 0) = accion de la Fase 2f en el gauge unitario: identicas")

# ------------------------------------------------------------------ 3. paridad temporal
IMPAR_R, IMPAR_E = ("al", "be"), ("w", "be")


def paridad(L2s, impares):
    malos = 0
    for tm in sp.Add.make_args(L2s):
        par = 0
        for s, p in tm.as_powers_dict().items():
            if isinstance(s, sp.Symbol) and "_" in s.name and s.name.split("_")[0] in ("a", "be", "c", "d", "al", "w"):
                v, o = s.name.split("_")
                par += p * ((v in impares) + int(o[0]))
        malos += par % 2
    return malos


L2R = sp.expand(L2.xreplace({UU[k]: 0 for k in range(CP.NMAXB)}))
L2E = sp.expand(L2.xreplace({VV[k]: 0 for k in range(CP.NMAXB)}).xreplace(
    {s: 0 for s in L2.free_symbols if s.name.startswith("G3_")}))
CHK["paridad_R_terminos_malos"] = paridad(L2R, IMPAR_R)
CHK["paridad_E_terminos_malos"] = paridad(L2E, IMPAR_E)
CHK["paridad_general_terminos_que_la_rompen"] = paridad(L2, IMPAR_R)
log(f"paridad T (A0 = 0): {CHK['paridad_R_terminos_malos']} terminos la violan; paridad TC (A1 = 0, G3 = 0): "
    f"{CHK['paridad_E_terminos_malos']}; fondo general (A0, A1 != 0): {CHK['paridad_general_terminos_que_la_rompen']}")
assert CHK["paridad_R_terminos_malos"] == 0 and CHK["paridad_E_terminos_malos"] == 0

# ------------------------------------------------------------------ 4. modos gauge en la capa de masas
Ta = sp.symbols("T0:9", real=True)
Za = sp.symbols("Z0:9", real=True)
gauge_ = {"a": Za[0] * NN[1] / NN[0] + iw * Ta[0], "be": -NN[0] ** 2 * Ta[1] + iw * Za[0], "c": Za[1],
          "d": rr[1] * Za[0], "al": Za[0] * UU[1] + UU[0] * iw * Ta[0] + VV[0] * iw * Za[0],
          "w": Za[0] * VV[1] + UU[0] * Ta[1] + VV[0] * Za[1]}
gsub = {}
for v, e in gauge_.items():
    for k in range(6):
        gsub[amp[(v, k)]] = e
        e = CP.D(e, [Ta, Za])
res_g = {}
for caso in ("general", "R", "E"):
    for sem in range(2):
        v = CP.chorro_en_capa(fe, 100 * sem + {"general": 1, "R": 2, "E": 3}[caso], orden=5, caso=caso)
        if caso != "E":
            for k in range(5, CP.NMAXB):
                v.setdefault(VV[k], v.r())          # datos libres de A1 (los G2_k0 que fijan no aparecen)
        aparecen = set().union(*[e.free_symbols for e in lin.values()])
        assert not any(s.name == "G2_50" for s in aparecen)
        for k in range(9):
            v[Ta[k]] = v.r()
            v[Za[k]] = v.r()
        w2 = sp.Rational(3, 11)
        mx = 0
        for nom, e in lin.items():
            x = v.ev(sp.expand(e.xreplace(gsub)))
            x = sp.expand(x.subs(iw ** 2, -w2))
            assert not x.free_symbols - {iw}, x.free_symbols
            mx = max(mx, abs(float(x.coeff(iw, 0))) + abs(float(x.coeff(iw, 1))))
        res_g[f"{caso}_{sem}"] = float(mx)
        log(f"modo gauge, chorro {caso} {sem}: max |E_v| = {mx:.1e}")
        assert mx == 0
CHK["modos_gauge"] = res_g
json.dump(CHK, open(os.path.join(AQUI, "resultados", "verificar_proca.json"), "w"), indent=1, default=str)
log("todas las comprobaciones superadas")
