# -*- coding: utf-8 -*-
"""
Fase 2k: comprobacion con la maquinaria covariante (generar_t.py) de la estructura ADM del caso phi = qt t (psi' = 0,
v = 0) en la GLPV cuartica sin termino cubico (G3 = 0; ver el informe, §6).

Prediccion ADM (corte phi = cte = corte estatico, rho = l, Nc = N, n = N^2 a = N nu, w = beta, gauge h = c = 0):
    L2 = L_est,2(nu, d) - r (Q r')' w^2 + 2 r (Q d_t)' w + Q d_t^2,      Q = 2 A4/N,  A4 = -G4 + 2X G4_X + 4X^2 F4,
    dE_Nc(nu = 0, d) = -4 beta_N D[d],  beta_N = G4 - 2X G4_X,  D[d] = (r d')' + (1 - r'^2 - r r'') d/r.
Con a = n/N^2 (E_n = E_a/N^2), E_B = E_beta, se comprueba en chorros exactos de la capa de masas con psi' = 0:
  (1) E_A = E_B = E_C = E_p = E_b = J = 0;
  (2) E_n no tiene terminos con iw; su parte en d es -4 beta_N D[d]/N;
  (3) E_B = -2 r (Q r')' B + 2 iw r (Q d' + Q' d), sin n;
  (4) en E_d: coeficiente de iw^2 d = -2Q; de iw B' = 2 Q r; de iw B = 2 Q r';
  (5) la ecuacion de beta = N G4 (nogo_beta.py): r^2 beta'' - r r' beta' + (1 + r r'' - r'^2) beta = 0.
Las ecuaciones del fondo se restringen a psi^(k) = 0, G3 = G5 = 0 ANTES de derivar (valido: cada termino de D_l
proporcional a psi^(k+1) o a G3_0(j+1) se anula).
Uso: python caso_qt.py
"""
import os
import sys
import json
import pickle
import sympy as sp

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
from chorros_t import Phi, rr, PP, Ee, Xs, qt, D, Valores, gsym, GPAT   # noqa: E402

dat = pickle.load(open(os.path.join(AQUI, "resultados", "t_simbolico.pkl"), "rb"))
lin, iw, amp = dat["lineales"], dat["iw"], dat["amp"]
cero = {s: 0 for s in PP[1:]}


def restringe(e):
    e = e.xreplace(cero)
    e = e.xreplace({s: 0 for s in e.free_symbols if GPAT.match(s.name) and s.name[:2] in ("G3", "G5")})
    return sp.expand(e)


fe0 = {k: restringe(v) for k, v in dat["fondo_eq"].items()}
_CACHE = []


def chorro_qt(semilla, orden=5):
    """psi' = 0: nivel 2, (r_2, G2_00) de (E_A, E_B); nivel k >= 3, (r_k, Phi_{k-1}) de (D^{k-2} E_A, D^{k-2} E_B)
    (con psi' = 0, E_B es la ligadura de primer orden H = 0 y sus derivadas fijan el lapso); resto de G2_0j al azar."""
    v = Valores(semilla, ("G3", "G5"))
    v[Ee] = sp.Rational(121, 100)
    v[qt] = sp.Rational(v.rnd.randint(200, 900), 1000)
    v[rr[0]] = sp.Rational(v.rnd.randint(600, 2000), 1000)
    v[rr[1]] = v.r()
    for k in range(1, 10):
        v[PP[k]] = 0
    v[Phi[1]] = v.r()
    v[Xs] = v[qt] ** 2 / (2 * v[Ee])
    if not _CACHE:
        _CACHE.append([fe0["A"], fe0["B"]])
    for k in range(2, orden + 2):
        while len(_CACHE) <= k - 2:
            _CACHE.append([D(e) for e in _CACHE[-1]])
        eqs = _CACHE[k - 2]
        g = gsym("G2", 0) if k == 2 else Phi[k - 1]
        inc = [rr[k], g]
        ev = [v.ev(e, salvo=(g,)) for e in eqs]
        Mk = sp.Matrix([[x.coeff(u) for u in inc] for x in ev])
        bk = sp.Matrix([sp.expand(x - sum(x.coeff(u) * u for u in inc)) for x in ev])
        assert all(not x.free_symbols for x in list(bk) + list(Mk)), (k, [x.free_symbols for x in bk])
        assert Mk.det() != 0, (k, Mk)
        for u, s_ in zip(inc, Mk.LUsolve(-bk)):
            v[u] = sp.Rational(s_)
    return v


NM = 5
n_ = sp.symbols(f"n0:{NM + 3}", real=True)
B_ = sp.symbols(f"B0:{NM + 3}", real=True)
d_ = sp.symbols(f"dd0:{NM + 3}", real=True)
AMPS = (n_, B_, d_)
base = {"a": n_[0] / Ee, "be": B_[0], "c": sp.Integer(0), "d": d_[0], "q": sp.Integer(0)}
sub = {}
for vn, e in base.items():
    for k in range(NM + 1):
        sub[amp[(vn, k)]] = e
        e = D(e, AMPS)
N = sp.sqrt(Ee)
A4 = -gsym("G4", 0) + 2 * Xs * gsym("G4", 1) + 4 * Xs ** 2 * gsym("F4", 0)
Qe = 2 * A4 / N
betaN = gsym("G4", 0) - 2 * Xs * gsym("G4", 1)
DDd = sp.expand(D(rr[0] * D(d_[0], (d_,)), (d_,)) + (1 - rr[1] ** 2 - rr[0] * rr[2]) * d_[0] / rr[0])
beta = N * gsym("G4", 0)
b1 = D(beta)
ecb = rr[0] ** 2 * D(b1) - rr[0] * rr[1] * b1 + (1 + rr[0] * rr[2] - rr[1] ** 2) * beta

if __name__ == "__main__":
    RES = {}
    for sem in range(3):
        ch = chorro_qt(1200 + sem, orden=5)
        res = {n2: ch.ev(e) for n2, e in dat["fondo_eq"].items()}
        res["J"] = ch.ev(dat["J"])
        assert all(x == 0 for x in res.values()), res
        assert ch.ev(D(fe0["A"])) == 0 and ch.ev(D(fe0["B"])) == 0
        subn = {k: ch.ev(e) for k, e in sub.items()}
        Eu = {vn: sp.expand(ch.ev(lin[vn]).xreplace(subn)) for vn in ("a", "be", "d")}
        eN, eB, eD = sp.expand(Eu["a"] / ch[Ee]), Eu["be"], Eu["d"]
        ok2a = all(eN.coeff(iw, p) == 0 for p in (1, 2, 3))
        parte_d = sp.expand(eN.subs(iw, 0).xreplace({s: 0 for s in n_ + B_}))
        ok2 = ok2a and sp.expand(parte_d - (-4 * ch.ev(betaN) / ch.ev(N) * ch.ev(DDd))) == 0
        Qv, Q1 = ch.ev(Qe), ch.ev(D(Qe))
        k_w = ch.ev(-2 * rr[0] * D(Qe * rr[1]))
        ok3 = sp.expand(eB - (k_w * B_[0] + 2 * iw * ch[rr[0]] * (Qv * d_[1] + Q1 * d_[0]))) == 0
        ok4 = (sp.expand(eD.coeff(iw, 2).coeff(d_[0])) == -2 * Qv and
               sp.expand(eD.coeff(iw, 1).coeff(B_[1])) == 2 * Qv * ch[rr[0]] and
               sp.expand(eD.coeff(iw, 1).coeff(B_[0])) == 2 * Qv * ch[rr[1]])
        ok5 = ch.ev(ecb) == 0
        RES[sem] = dict(fondo=True, E_n=ok2, E_B=ok3, E_d_cinetico=ok4, ecuacion_beta=ok5)
        print(f"chorro {sem}: (1) fondo, J, E_b = 0; (2) E_n sin iw y parte en d = -4 beta_N D[d]/N: {ok2}; "
              f"(3) E_B: {ok3}; (4) E_d cinetico: {ok4}; (5) ecuacion de beta: {ok5}", flush=True)
        if not (ok2 and ok3 and ok4 and ok5):
            print("   E_n parte d:", parte_d, " esperado:", sp.expand(-4 * ch.ev(betaN) / ch.ev(N) * ch.ev(DDd)))
            print("   E_B:", eB, " esperado k_w =", k_w, " Q =", Qv, " Q' =", Q1)
            print("   E_d iw^2 d0:", eD.coeff(iw, 2).coeff(d_[0]), " iw B1:", eD.coeff(iw, 1).coeff(B_[1]),
                  " iw B0:", eD.coeff(iw, 1).coeff(B_[0]))
    json.dump(RES, open(os.path.join(AQUI, "resultados", "caso_qt.json"), "w"), indent=1)
