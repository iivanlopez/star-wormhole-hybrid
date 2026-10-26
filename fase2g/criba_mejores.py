# -*- coding: utf-8 -*-
"""
Fase 2g, parte III: las mejores configuraciones de la optimizacion de la criba (criba_optimiza.py, semillas 1 y 2):
convergencia de los margenes en la malla, valores sin normalizar, presupuesto de fantasma, sector impar l = 2, 3
(sin funciones libres, impar_fluido.py), cantidades fisicas (densidad, presion, compacidad) y perfiles para las figuras.
Escribe resultados/criba_mejores.json y resultados/perfiles_criba.npz.
Uso: python criba_mejores.py
"""
import os
import sys
import json
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import criba as C          # noqa: E402


def fisica(F):
    l = np.linspace(0, F.ls, 4001)
    f = F.campos(l)
    rho_t = f["rp"] - f["p"]                   # rho~ = (rho~ + p~) - p~
    # masa de la materia (en r) y fraccion de la masa ADM; rho~ = 16 pi rho
    masa_mat = np.trapezoid(4 * np.pi * f["r"] ** 2 * rho_t / (16 * np.pi) * f["r1"], l)
    return dict(rho0=float(rho_t[0] / (16 * np.pi)), p0=float(f["p"][0] / (16 * np.pi)), cs2_0=float(f["cs2"][0]),
                r_s=float(f["r"][-1]), masa_materia=float(masa_mat), M=F.M, N0=F.N0, ls=F.ls,
                G4_0=float(f["G4"][0]), Gc_0=float(f["Gc"][0]))


if __name__ == "__main__":
    out = {}
    perf = {}
    for sem in (1, 2):
        d = json.load(open(os.path.join(AQUI, "resultados", f"criba_optimiza_{sem}.json")))
        x = d["x"]
        F = C.FondoFluido(**x)
        conv = {n: C.criba_normalizada(F, n=n) for n in (120, 240, 480)}
        crudo = C.criba(F, n=240)
        pf = C.presupuesto_fantasma(F)
        imp = {LL: C.impar(F, LL) for LL in (2, 3)}
        out[f"semilla_{sem}"] = dict(x=x, convergencia={str(k): v for k, v in conv.items()}, crudo=crudo, fantasma=pf,
                                     impar={str(k): v for k, v in imp.items()}, fisica=fisica(F))
        print(f"semilla {sem}: x = { {k: round(v, 4) if isinstance(v, float) else v for k, v in x.items()} }")
        for n, v in conv.items():
            print(f"   n = {n}: " + ", ".join(f"{k} = {val:+.4f}" for k, val in v.items()))
        print("   sin normalizar: " + ", ".join(f"{k} = {float(v):+.4f}" for k, v in crudo.items()))
        print(f"   fantasma: gasto = {pf['gasto']:.4f} frente a 3M = {pf['presupuesto']:.4f}")
        print(f"   impar: {imp}")
        print(f"   fisica: {out[f'semilla_{sem}']['fisica']}")
        l = np.concatenate([np.linspace(1e-4, F.ls * 1.5, 1500), np.geomspace(F.ls * 1.5, 400, 500)[1:]])
        f = F.campos(l)
        for k in ("N", "r", "r1", "Gc", "G4", "rp", "cs2", "p", "sig"):
            perf[f"s{sem}_{k}"] = f[k]
        perf[f"s{sem}_l"] = l
        perf[f"s{sem}_ls"] = np.array([F.ls])
    json.dump(out, open(os.path.join(AQUI, "resultados", "criba_mejores.json"), "w"), indent=1, default=float)
    np.savez(os.path.join(AQUI, "resultados", "perfiles_criba.npz"), **perf)
