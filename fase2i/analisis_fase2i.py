# -*- coding: utf-8 -*-
"""
Fase 2i: datos del informe (resultados/analisis_fase2i.json y resultados/perfiles_fase2i.npz).
  A. Validacion: criba c = 0 con la forma lagrangiana frente a la euleriana de la Fase 2h (candidato rejilla_pozo).
  B. El candidato de la Fase 2h como solido: margenes frente a c_T^2 (c_L^2 = 1) y modo par (c_T^2 = 0 y 0.2).
  C. La traslacion casi gauge en la estrella con nucleo hueco (la = 3, c_T^2 = 0.3): desviacion de la relacion de RG
     nu = r w/(4 Gcal) y densidad de energia (vacio y material) de la prueba h = r'/N.
  D. Cascaras de membrana (p~_r = 0): margenes frente a r_b/M, y el cociente analitico Q(x) de la cascara delgada.
  E. f_m (anclaje del material) en la estrella con nucleo hueco.
Uso: python analisis_fase2i.py
"""
import os
import sys
import json
import time
import warnings
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2g"), os.path.join(AQUI, "..", "fase2h"), AQUI):
    sys.path.insert(0, p_)
import criba as C                    # noqa: E402
import fondo_nu as FN                # noqa: E402
import criba_solido as CS            # noqa: E402
import fondo_cascara as FC           # noqa: E402
import fondo_solido as FS            # noqa: E402

warnings.filterwarnings("ignore")
T0 = time.time()
RES = {}
PERF = {}


def log(m):
    print(f"[{time.time() - T0:7.1f} s] {m}", flush=True)


x = json.load(open(os.path.join(AQUI, "..", "fase2h", "resultados", "rejilla_pozo.json")))["x"]
F = FN.FondoNu(**x)
BASE = {k: x[k] for k in ("wt", "nu0", "phi1", "b", "n", "d0", "wd", "lam", "lT")}


class PerfilSolido(CS.PerfilFluido):
    """La estrella fluida de la Fase 2h con las mismas tensiones (isotropas) y modulos de solido (c_L, c_T)."""

    def __init__(self, F, cL2, cT2):
        super().__init__(F)
        self.cL2, self.cT2 = cL2, cT2

    def perfil(self, l):
        p = super().perfil(l)
        rho = -p["G0"]
        pr = p["w"] - rho
        p.update(CS.modulos(rho, pr, pr, self.cL2, self.cT2))
        return p


def resumen(r):
    return {k: float(v) for k, v in r.items() if not k.startswith("modo")}


# ------------------------------------------------------------------ A. validacion
RES["A"] = dict(lagrangiana=resumen(CS.criba(CS.PerfilFluido(F), n=160)),
                euleriana=resumen(C.criba_normalizada(F, n=160, libre=True)), E1_euleriana=C.prueba_sin_compresion(F)["E1"])
log(f"A: Epar lagrangiana {RES['A']['lagrangiana']['Epar']:+.5f}, euleriana {RES['A']['euleriana']['Epar']:+.5f}")

# ------------------------------------------------------------------ B. cizalla en la estrella de la Fase 2h
RES["B"] = []
for cT2 in (0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.37):
    r = CS.criba(PerfilSolido(F, 1.0, cT2), n=160, devolver=cT2 in (0.0, 0.2))
    RES["B"].append(dict(cT2=cT2, **resumen(r)))
    if "modo_par" in r:
        nod, h, hp = r["modo_par"]
        k = int(np.argmax(np.abs(h)))
        PERF[f"B_l_{int(100 * cT2)}"], PERF[f"B_h_{int(100 * cT2)}"] = nod, h / h[k]
    log(f"B: cT2 = {cT2}: E1n {r['E1n']:+.4f}, Epar {r['Epar']:+.4f}, Eimpar {r['Eimpar']:+.4f}")

# ------------------------------------------------------------------ C. traslacion casi gauge
G = FS.FondoSolido(la=3.0, lchi=2.0, A=x["A"], Ns=x["Ns"], cL2=1.0, cT2=0.3, **BASE)
l = np.linspace(G.la, G.lb * (1 - 1e-9), 8001)
p = G.perfil(l)
N, r, r1, Gc, w, sig, G4 = (p[k] for k in ("N", "r", "r1", "Gc", "w", "sig", "G4"))
ph = G.phi1 * G.b ** 2 / (l ** 2 + G.b ** 2)
nu = -G.nu0 * np.exp(-(l / G.wt) ** 2) + G.f_fl * r * w / (4 * G4) + 0.5 * G.kappa * r * ph ** 2
h = r1 / N
hp = np.gradient(h, l)
hpp = np.gradient(hp, l)
q = 4 * N * Gc / (r * w)
qp = q * (p["Phi1"] + p["Gcp"] / Gc - r1 / r - p["wp"] / w)
xi = q * hp
V = np.stack([xi, qp * hp + q * hpp, N * h, N * p["Phi1"] * h + N * hp])
em = -(r ** 2) * N * np.einsum("in,ijn,jn->n", V, CS.forma_lag(p), V)
ev = -2 * N ** 3 * Gc * (hp - sig * h) ** 2
PERF.update(C_l=l, C_desv=nu / (r * w / (4 * Gc)) - 1, C_xi=xi, C_ev=ev, C_em=em)
RES["C"] = dict(la=G.la, lb=G.lb, M=G.M_adm, E_out_h2=float(G.E_out() * h[-1] ** 2), regiones=[])
for a, b in ((3, 10), (10, 20), (20, 40), (40, 60), (60, G.lb)):
    m = (l >= a) & (l <= b)
    RES["C"]["regiones"].append(dict(a=a, b=float(b), xi=float(xi[m].mean()), vacio=float(np.trapezoid(ev[m], l[m])),
                                     materia=float(np.trapezoid(em[m], l[m]))))
log("C: traslacion: " + "; ".join(f"[{e['a']},{e['b']:.0f}] {e['vacio'] + e['materia']:+.2e}" for e in RES["C"]["regiones"]))

# ------------------------------------------------------------------ D. cascaras de membrana
RES["D"] = []
base_c = dict(wt=2.0, nu0=2.5, phi1=0.5, b=3.0, n=2, d0=0.5, wd=1.0, lam=1.0, lT=30.0)
for la, lb, W0 in ((60, 80, 1.0e-3), (60, 80, 1.5e-3), (60, 80, 2e-3), (60, 70, 3.5e-3), (100, 110, 1.5e-3),
                   (100, 110, 2e-3), (30, 40, 5e-3), (30, 40, 7e-3)):
    H = FC.FondoCascara(la=la, lb=lb, W0=W0, s_w=0.0, pi_p=0.0, cL2=1.0, cT2=0.1, **base_c)
    if not H.ok:
        continue
    r = CS.criba(H, n=120)
    pa = H.perfil(np.array([lb]))
    fis = H.fisica()
    RES["D"].append(dict(la=la, lb=lb, W0=W0, M=H.M_adm, rb_M=float(pa["r"][0] / H.M_adm), masa=fis["masa"],
                         pt_max=fis["pt_max"], **resumen(r)))
    log(f"D: [{la},{lb}] W0 {W0}: r_b/M {pa['r'][0] / H.M_adm:.2f}, E1n {r['E1n']:+.3f}, Epar {r['Epar']:+.3f}")
xs = np.linspace(0.05, 0.45, 81)
Nx = np.sqrt(1 - 2 * xs)
PERF.update(D_x=xs, D_Q=4 * xs * (0.9 * Nx ** 3 + Nx * xs / 2) / (1 - 3 * xs + 3 * xs ** 2))

# ------------------------------------------------------------------ E. f_m
RES["E"] = []
for fm in (0.4, 0.7, 1.0, 1.3, 1.6, 2.0):
    G = FS.FondoSolido(la=3.0, lchi=2.0, A=x["A"], Ns=x["Ns"], cL2=1.0, cT2=0.2, f_m=fm, **BASE)
    if not G.ok:
        continue
    r = CS.criba(G, n=160)
    RES["E"].append(dict(f_m=fm, M=G.M_adm, ls=G.ls, **resumen(r)))
    log(f"E: f_m {fm}: E1n {r['E1n']:+.4f}, Epar {r['Epar']:+.4f}, Eimpar {r['Eimpar']:+.4f}")

json.dump(RES, open(os.path.join(AQUI, "resultados", "analisis_fase2i.json"), "w"), indent=1, default=float)
np.savez(os.path.join(AQUI, "resultados", "perfiles_fase2i.npz"), **PERF)
log("escrito resultados/analisis_fase2i.json y perfiles_fase2i.npz")
