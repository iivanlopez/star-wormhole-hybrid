# -*- coding: utf-8 -*-
"""
Barrido robusto de modos inestables reales (estabilidad_F.raices_menores: dato inicial por menores 4x4, analitico, sin
seguimiento de signo; integracion renormalizada).  Sustituye a los barridos con el vector nulo unitario, cuyo seguimiento
de signo por continuidad puede perder raices cuando el vector gira deprisa entre puntos de la malla.
Fondos: los 19 de KKK (estabilidad_dominio), los hiperbolicos de la familia cuadratica (familia_cuadratica +
hiperbolicidad) y las dos dobles gargantas exponenciales hiperbolicas; zeta = 1 y 0.
Salida: resultados/estabilidad_robusta.json.   Uso: python estabilidad_robusta.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from fondo_F import FondoF
from estabilidad_F import CapaF, raices_menores

WS = -np.logspace(6, -4, 121)


def tarea(args):
    nombre, tipo, al, p0, u, z = args
    f = FondoF(tipo, al, p0, u).resolver()
    rs = raices_menores(f, CapaF(f, zeta=z), WS)
    return dict(nombre=nombre, tipo=tipo, alpha=al, phi0=p0, phi1=u, zeta=z, ecuador=bool(f.ecuador), f0=f.f0,
                alpha_eta=f.alpha_eta, alpha_kkk=f.alpha * math.exp(-f.phi_inf) if tipo == "exp" else None,
                e2Phi=math.exp(2 * f.Phi_inf), w=rs, winf=[x * math.exp(-2 * f.Phi_inf) for x in rs])


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    casos = []
    for o in json.load(open("resultados/estabilidad_dominio.json", encoding="utf-8")):
        if "w_zeta1" in o:
            casos += [(f"KKK {o['alpha_kkk']:.4f} {o['f0']:.4f}", "exp", o["alpha"], 0.0, o["phi1"], z) for z in (1.0, 0.0)]
    hip = {o["nombre"]: o["hiperbolico"] for o in json.load(open("resultados/hiperbolicidad.json", encoding="utf-8"))}
    for o in json.load(open("resultados/familia_cuadratica.json", encoding="utf-8")):
        nom = f"cuadratico {o['alpha_eta']:.4f}"
        if hip[nom]:
            casos += [(nom, "cuadratico", o["alpha"], o["phi0"], o["phi1"], z) for z in (1.0, 0.0)]
    for a, u in ((0.02, -12.340), (0.04, -4.534)):
        casos += [(f"doble exp a={a} u={u}", "exp", a, 0.0, u, z) for z in (1.0, 0.0)]
    print(f"{len(casos)} barridos", flush=True)
    res = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        for o in ex.map(tarea, casos):
            res.append(o)
            print(f"{o['nombre']:28s} zeta={o['zeta']:g}: omega_inf^2 r0^2 = {[f'{x:.4g}' for x in o['winf']]}", flush=True)
            json.dump(res, open("resultados/estabilidad_robusta.json", "w", encoding="utf-8"), indent=1)
