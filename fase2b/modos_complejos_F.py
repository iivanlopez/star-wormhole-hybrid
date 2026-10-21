# -*- coding: utf-8 -*-
"""
Resonancias inestables complejas por encima de zeta_c para un fondo cuadratico de Antoniou et al.
(alpha/eta0^2 = 0.3087) y, como control, para el fondo de referencia de KKK (zeta = 2: w = 1149.06 + 170.71 i).
Secante sobre la funcion entera Dc desde varios arranques escalados con a^2 = Phi'(0+)^2 (en KKK, w/a^2 = 6.43 + 0.96 i).
Un cero con Im w > 0 en la hoja principal es un modo que crece (oscilacion creciente).
Salida: resultados/modos_complejos_F.json.   Uso: python modos_complejos_F.py
"""
import sys
import json
import math
from concurrent.futures import ProcessPoolExecutor
from fondo_F import FondoF
from estabilidad_F import CapaF, Dc, secante

ARRANQUES = [6.43 + 0.96j, 5.0 + 1.5j, 8.0 + 0.5j, 4.0 + 0.3j]


def tarea(args):
    nombre, par, z, x0 = args
    f = FondoF(*par).resolver()
    c = CapaF(f, zeta=z)
    a2 = f.y0()[3] ** 2
    w, paso, dD = secante(f, c, x0 * a2)
    # control: |Dc| en la raiz frente a puntos cercanos
    vecino = max(abs(Dc(f, c, w * (1 + 1e-3 * e))) for e in (1, -1, 1j, -1j))
    om = -1j * (-(w * math.exp(-2 * f.Phi_inf))) ** 0.5     # omega_inf r0 con Im < 0 (rama principal de kappa)
    return dict(nombre=nombre, zeta=z, arranque=[x0.real, x0.imag], w=[w.real, w.imag], w_a2=[(w / a2).real, (w / a2).imag],
                paso=paso, Dc=dD, Dc_vecino=vecino, convergido=bool(paso < 1e-8 * abs(w) and dD < 1e-6 * vecino),
                omega_inf=[om.real, om.imag], Phi_inf=f.Phi_inf)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    fam = {f"{o['alpha_eta']:.4f}": ("cuadratico", o["alpha"], o["phi0"], o["phi1"])
           for o in json.load(open("resultados/familia_cuadratica.json", encoding="utf-8"))}
    zetas = [float(x) for x in sys.argv[1:]] or [2.0]
    casos = [("KKK referencia", ("exp", 0.03, 0.0, -0.7), 2.0, ARRANQUES[0])]
    casos += [("cuadratico 0.3087", fam["0.3087"], z, x0) for z in zetas for x0 in ARRANQUES]
    with ProcessPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(tarea, casos))
    for o in res:
        print(f"{o['nombre']:18s} zeta={o['zeta']:g} arranque w/a^2={complex(*o['arranque']):.2f}: w = {o['w'][0]:.6g} "
              f"{o['w'][1]:+.6g} i  (w/a^2 = {o['w_a2'][0]:.4f} {o['w_a2'][1]:+.4f} i)  |Dc| = {o['Dc']:.1e} frente a "
              f"{o['Dc_vecino']:.1e}  convergido={o['convergido']}  omega_inf r0 = {o['omega_inf'][0]:+.4f} "
              f"{o['omega_inf'][1]:+.4f} i")
    json.dump(res, open("resultados/modos_complejos_F.json", "w", encoding="utf-8"), indent=1)
