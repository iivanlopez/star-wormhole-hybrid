# -*- coding: utf-8 -*-
"""
Modo impar inestable del Modelo L (flujo de fluido a traves de la garganta):
 (a) escala con el tamano de la garganta (xi0 = 1e-3 ... 0.1, rho_bc = 1e14);
 (b) tasa local frente a 4 pi (eps0 + 3 P0) a lo largo de las secuencias (de resultados/h1_impar.json);
 (c) autofuncion: perfil xi/xi(0) frente a (r0/r)^2, deslizamiento C del andamio y campo Psi(0) en la garganta.
Salida: resultados/impar_garganta.json.   Uso: python impar_garganta.py   (tras test_h1_impar.py)
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq
from hibrido import Politropo, Fondo, RHO_GEOM, PI4
from autofuncion import autofuncion_gs

eos = Politropo(k=0.1, gamma=2.0)


def fondo(xi0, rho):
    rb = rho * RHO_GEOM
    f = Fondo(r0=xi0 * eos.inv_A(rb), eos=eos, rho_bc=rb).resolver()
    f.paridad = "impar"
    return f


def escala(xi0, rho=1e14):
    f = fondo(xi0, rho)
    k = math.exp(2 * f.Phi_inf)
    ws = np.sort(np.concatenate([-np.logspace(-1.5, -6.5, 36) * k, np.linspace(2e-6, 6e-4, 25) * k]))
    D = [f.det("ligado", w) for w in ws]
    r = [f.w_fisico(brentq(lambda w: f.det("ligado", w), ws[i], ws[i + 1], rtol=1e-10))
         for i in range(len(ws) - 1) if D[i] * D[i + 1] < 0]
    return dict(xi0=xi0, rho=rho, r0=f.r0, raices=r)


def autof(xi0, rho=1e14):
    f = fondo(xi0, rho)
    k = math.exp(2 * f.Phi_inf)
    w = brentq(lambda x: f.det("ligado", x), -1.2e-3 * k, -0.6e-3 * k, xtol=1e-16, rtol=1e-12)
    l, X, sv = autofuncion_gs(f, "ligado", w)
    xi, dr, Ps = X[0] / X[0][0], X[2] / X[0][0], X[3] / X[0][0]
    r = np.array([f.y_bg(x)[0] for x in l])
    eps0, P0 = f.bg(f.y0(), True)[12:14]
    perfil = []
    for x in (f.r0, 3 * f.r0, 10 * f.r0):
        i = int(np.argmin(abs(l - x)))
        perfil.append(dict(l=float(l[i]), r=float(r[i]), xi=float(xi[i]), r0_r2=float((f.r0 / r[i]) ** 2)))
    return dict(xi0=xi0, rho=rho, r0=f.r0, w_fis=f.w_fisico(w), w_int=w,
                C=float(dr[-1] * math.exp(-f.y_bg(l[-1])[2])), Psi0=float(Ps[0]), Phi2_0=PI4 * (eps0 + 3 * P0),
                perfil=perfil, sv=[float(x) for x in sv])


def local(o):
    f = fondo(o["xi0"], o["rho"])
    eps0, P0 = f.bg(f.y0(), True)[12:14]
    w1 = o["impar"][0] * math.exp(2 * f.Phi_inf)
    return dict(o, eps0=eps0, P0=P0, Phi_inf=f.Phi_inf, w1_int=w1, cociente=abs(w1) / (PI4 * (eps0 + 3 * P0)))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    sec = json.load(open("resultados/h1_impar.json", encoding="utf-8"))
    with ProcessPoolExecutor(max_workers=3) as ex:
        fe = [ex.submit(escala, x) for x in (0.001, 0.003, 0.01, 0.03, 0.1)]
        fa = [ex.submit(autof, x) for x in (0.1, 0.01)]
        fl = [ex.submit(local, o) for o in sec]
        out = dict(escala=[x.result() for x in fe], autofuncion=[x.result() for x in fa],
                   secuencia=[x.result() for x in fl])
    for o in out["escala"]:
        print(f"escala: xi0={o['xi0']:<6} r0={o['r0']:.4f} km  omega^2 impar = {[f'{x:.5e}' for x in o['raices']]}")
    for o in out["autofuncion"]:
        print(f"autofuncion xi0={o['xi0']}: omega^2={o['w_fis']:.5e}  C/xi(0)={o['C']:.3e}  Psi(0)/xi(0)={o['Psi0']:.4e}"
              f"  Phi''(0)={o['Phi2_0']:.4e}  perfil xi vs (r0/r)^2: "
              + ", ".join(f"{p['xi']:.4f}/{p['r0_r2']:.4f}" for p in o["perfil"]))
    for o in out["secuencia"]:
        print(f"secuencia xi0={o['xi0']:<4} rho={o['rho']:.0e}: w1 (garganta) = {o['w1_int']:.4e}  "
              f"|w1| / 4pi(eps0+3P0) = {o['cociente']:.3f}")
    json.dump(out, open("resultados/impar_garganta.json", "w", encoding="utf-8"), indent=1)
