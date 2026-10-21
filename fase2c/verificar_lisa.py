# -*- coding: utf-8 -*-
"""
Pruebas unitarias de la maquinaria sin empalmes (estabilidad_lisa.py) sobre una mitad regular AF (lado atractivo de una
garganta lisa de Antoniou, F = alpha phi^2).  No es un calculo de modos: sin fondos lisos AF no hay problema fisico.
  1. el dato asintotico decreciente cumple las ligaduras en l = L, y la solucion las conserva hasta l = 0 (Bianchi);
  2. invariancia gauge: anadir c g(L) al dato asintotico no cambia (delta r_t, delta phi_t) en l = 0;
  3. la integracion renormalizada conserva la fase y el cociente de componentes (frente a la integracion sin
     renormalizar), con w real y complejo;
  4. prueba de humo del calculo completo de Delta(w) sobre un fondo FICTICIO (la misma mitad AF a los dos lados, sin
     significado fisico): Delta renormalizado y sin renormalizar con la misma fase, busqueda de raices e
     hiperbolicidad ejecutadas de punta a punta.
Salida: resultados/verificar_lisa.json.   Uso: python verificar_lisa.py
"""
import os
import sys
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fondo_doble import Lado                                            # noqa: E402
import estabilidad_lisa as E                                            # noqa: E402


class Mitad:
    def __init__(self, lado):
        self.Fd, self.sol = lado.Fd, lado.sol
        self.rpp0, self.phi1 = lado.rpp0, lado.phi1
        y0 = lado.y0()
        self.phi1 = y0[5]


def g_en(f, l, w):
    (F1, F2, F3), (r, r1, r2, P1, P2, p1, p2, e2) = E.fondo_en(f, l)
    return np.array([P2 - w * e2, r1, r2, p1, p2], dtype=complex if isinstance(w, complex) else float)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    lado = Lado("cuadratico", 0.1, 0.1, 0.32989, l_max=400, rtol=1e-12).resolver()
    f = Mitad(lado)
    out = []
    for w in (-5.0, -400.0, complex(-5.0, 2.0), complex(-400.0, 150.0)):
        kap = abs(np.sqrt(-complex(w))) * np.exp(-lado.Phi_inf)
        L = min(10.0 + 30.0 / kap, 0.9 * lado.l_end)
        y0 = E.dato_asintotico(f, L, w)
        lig_L = np.max(np.abs(E.ligaduras(f, L, y0, w))) / np.max(np.abs(y0))
        yr, logs = E.integrar_hacia_garganta(f, L, w, rtol=1e-11, renorm=True)
        yn, _ = E.integrar_hacia_garganta(f, L, w, rtol=1e-11, renorm=False)
        lig_0 = np.max(np.abs(E.ligaduras(f, 0.0, yr, w))) / np.max(np.abs(yr))
        # 3. fase y cocientes: yn = yr * e^{logs} (factor positivo)
        coc = np.max(np.abs(yn / np.max(np.abs(yn)) - yr / np.max(np.abs(yr))))
        # 2. invariancia gauge
        c = 0.37 * np.max(np.abs(y0))
        from scipy.integrate import solve_ivp
        tipo = complex if isinstance(w, complex) else float
        s = solve_ivp(lambda x, yy: E.derivadas(f, x, yy, w), (L, 0.0), (y0 + c * g_en(f, L, w)).astype(tipo),
                      method="DOP853", rtol=1e-11, atol=1e-14 * np.max(np.abs(y0)))
        yg = s.y[:, -1]
        I0 = np.array(E.invariantes_garganta(f, yn))
        Ig = np.array(E.invariantes_garganta(f, yg))
        dinv = np.max(np.abs(Ig - I0)) / np.max(np.abs(I0))
        dvec = np.max(np.abs(yg - yn)) / np.max(np.abs(yn))
        o = dict(w=str(w), L=L, lig_L=float(lig_L), lig_0=float(lig_0), cociente_renorm=float(coc),
                 cambio_vector=float(dvec), cambio_invariantes=float(dinv), log_escala=float(logs))
        out.append(o)
        print(f"w = {w}: L = {L:.3g}; ligaduras en L {lig_L:.1e}, en 0 {lig_0:.1e}; |renorm - sin renorm| = {coc:.1e}; "
              f"+c g(L): cambio de y(0) {dvec:.2e}, de (delta r_t, delta phi_t) {dinv:.1e}")
    # 4. prueba de humo sobre un fondo ficticio (sin significado fisico)
    from fondo_doble import Doble, _Densa
    fic = Doble("cuadratico", 0.1, 0.1, 0.32989, l_max=400, rtol=1e-11)
    fic.izq = Lado("cuadratico", 0.1, 0.1, 0.32989, l_max=400, rtol=1e-11).resolver()
    fic.sol = _Densa(fic.der, fic.izq)
    humo = {}
    for w in (-50.0, complex(-50.0, 10.0)):
        D = E.Delta(fic, w)
        Dn, logD = E.Delta_sin_renormalizar(fic, w)
        dfase = float(abs(np.angle(D) - np.angle(Dn)))
        humo[str(w)] = dfase
        print(f"humo, w = {w}: diferencia de fase entre Delta renormalizado y sin renormalizar = {dfase:.1e} rad")
    raices, _ = E.raices_reales(fic, -np.logspace(0, 3, 12))
    lam = E.hiperbolicidad(fic, [-0.5, 0.0, 0.5])
    print(f"humo: raices_reales y hiperbolicidad ejecutadas ({len(raices)} cambios de signo; lambda_min = {lam})")
    os.makedirs("resultados", exist_ok=True)
    json.dump(dict(pruebas=out, humo_fase=humo), open("resultados/verificar_lisa.json", "w", encoding="utf-8"), indent=1)
