# -*- coding: utf-8 -*-
"""
Naturaleza del final de las mitades con a = Phi'(0) < 0 (lapso decreciente hacia fuera).
Para una muestra de mitades se integra hasta |Phi'| = 1e7 y, cerca del final l_s:
  - posicion l_s y exponente p de N = e^Phi ~ (l_s - l)^p: 1/Phi' es lineal en l con pendiente 1/p (ajuste en la cola);
  - r(l_s), r'(l_s) (r' != 0 => las esferas no son extremas: no es un horizonte regular);
  - Kretschmann K = 4 (N''/N)^2 + 8 (r' N'/(r N))^2 + 8 (r''/r)^2 + 4 ((1 - r'^2)/r^2)^2 y su divergencia K (l_s - l)^2;
  - comparacion de l_s con la distancia de Rindler 1/|a|;
  - relacion de horizonte: si Phi' -> -inf con r, r', phi, phi' finitos, la ligadura ll (lineal en Phi') exige que se anule
    el coeficiente de Phi':  r_s r'_s = 2 F'(phi_s) phi'_s (3 r'_s^2 - 1)  (con F' phi' < 0: 0 < r'_s < 1/sqrt 3).
Ademas, para los finales 'phi' y 'r0' (patada escalar): masa de Misner-Sharp minima, |r'| maximo y N minimo.
Salida: resultados/singularidad.json.   Uso: python singularidad.py
"""
import os
import sys
import json
import math
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fondo_doble import Lado, segundas, H     # noqa: E402


def ligadura_rel(L, l):
    """residuo relativo de la ligadura ll (suma de |terminos| como escala)."""
    r, rp, Phi, P, p, u = L.sol.sol(l)[:6]
    F1 = L.Fd(p)[1]
    t = [(rp * rp - 1) / r**2, 2 * P * rp / r, -u * u / 4, 4 * F1 * P * u * (1 - 3 * rp * rp) / r**2]
    return abs(H.est_ligadura(F1, r, rp, P, u)) / sum(abs(x) for x in t)


def kretschmann(L, l):
    y = L.sol.sol(l)
    x, _ = segundas(y, L.Fd)
    r, rp, Phi, P, p, u = y[:6]
    rpp, Ppp = x[0], x[1]
    NppN = Ppp + P * P
    return 4 * NppN**2 + 8 * (rp * P / r) ** 2 + 8 * (rpp / r) ** 2 + 4 * ((1 - rp * rp) / r**2) ** 2


def analiza(tipo, alpha, phi0, phi1):
    L = Lado(tipo, alpha, phi0, phi1, l_max=400, rtol=1e-12).resolver()
    out = dict(tipo=tipo, alpha=alpha, phi0=phi0, phi1=phi1, a=L.a, fin=L.fin, l_end=L.l_end)
    if L.fin != "horizonte":
        return out
    le = L.l_end
    # cola: puntos donde |Phi'| > 30 |a| (ya en el regimen singular)
    ls = np.linspace(0.0, le, 4001)[:-1]
    ys = np.array([L.sol.sol(x) for x in ls])
    m = np.abs(ys[:, 3]) > 30 * abs(L.a)
    if m.sum() < 8:
        ls = le - np.logspace(-9, math.log10(le), 3000)[::-1]
        ys = np.array([L.sol.sol(x) for x in ls])
        m = np.abs(ys[:, 3]) > 30 * abs(L.a)
    lt, invP = ls[m], 1.0 / ys[m, 3]
    A = np.vstack([lt, np.ones_like(lt)]).T
    (pend, ord_), *_ = np.linalg.lstsq(A, invP, rcond=None)
    p_exp = 1.0 / pend
    l_s = -ord_ / pend
    yend = L.sol.y[:, -1]
    Ks = [kretschmann(L, l) * (l_s - l) ** 2 for l in lt[-5:]]
    r, rp, u = yend[0], yend[1], yend[5]
    A = 4 * L.Fd(yend[4])[1] * u
    raices = [x.real for x in np.roots([3 * A, -2 * r, -A]) if abs(x.imag) < 1e-12 and x.real > 0] if A else [0.0]
    out.update(l_s=float(l_s), p=float(p_exp), r_s=float(r), rp_s=float(rp), phi_s=float(yend[4]),
               phip_end=float(u), Phi_end=float(yend[2]), K_x2=[float(k) for k in Ks],
               ls_a=float(l_s * abs(L.a)), Fp_phip=float(A / 4), rp_relacion=float(raices[0]) if raices else None,
               ligadura_max=float(max(ligadura_rel(L, l) for l in np.linspace(0, le, 400)[1:])))
    return out


def patada(tipo, alpha, phi0, phi1):
    """Finales 'phi' y 'r0' de la mitad repulsiva: masa de Misner-Sharp minima, |r'| maximo, N minimo, ecuador."""
    L = Lado(tipo, alpha, phi0, phi1, l_max=400, rtol=1e-12).resolver()
    ls = np.linspace(0.0, L.l_end, 20001)
    ys = np.array([L.sol.sol(x) for x in ls])
    r, rp, N = ys[:, 0], ys[:, 1], np.exp(ys[:, 2])
    m = r * (1 - rp**2) / 2
    k = int(np.argmin(N))
    return dict(tipo=tipo, alpha=alpha, phi0=phi0, phi1=phi1, a=float(L.a), fin=L.fin, l_end=float(L.l_end),
                m_min=float(m.min()), rp_abs_max=float(np.abs(rp).max()), N_min=float(N.min()), l_Nmin=float(ls[k]),
                ecuador=[float(x) for x in L.cruces], r_fin=float(r[-1]), phi_fin=float(ys[-1, 4]))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    pat = []
    for c in [("cuadratico", 0.02, 0.01, -20.0), ("cuadratico", 0.02, 0.1, -5.0), ("cuadratico", 0.5, 0.1, -2.0),
              ("cuadratico", 0.5, 0.03, -5.0), ("cuadratico", 2.0, 0.02, -1.0), ("cuadratico", 2.0, 0.1, -3.0),
              ("lineal", 10.0, 0.0, -5.0), ("lineal", 1.0, 0.0, -20.0)]:
        o = patada(*c)
        pat.append(o)
        print(f"patada {c}: a={o['a']:+.4g} fin={o['fin']} l_fin={o['l_end']:.4g}  m_min={o['m_min']:+.4g}  "
              f"|r'|max={o['rp_abs_max']:.3f}  N_min={o['N_min']:.3e} en l={o['l_Nmin']:.3g}  ecuador en {o['ecuador'][:2]}  "
              f"r_fin={o['r_fin']:.3g} phi_fin={o['phi_fin']:.3g}")
    casos = [("cuadratico", 0.1, 0.1, -0.3299), ("cuadratico", 0.1, 0.01, -0.3551), ("cuadratico", 0.09, 0.3, -0.07226),
             ("cuadratico", 0.5, 1.0, -1.0), ("cuadratico", -0.5, 1.0, 1.0), ("cuadratico", -2.0, 0.3, 3.0),
             ("exp", 0.5, 0.0, 1.0), ("exp", 2.0, 0.0, 2.0), ("exp", -1.0, 0.0, -1.0), ("exp", 0.03, 0.0, 0.7),
             ("exp", -30.0, 0.0, -20.0), ("gauss", 0.2, 0.1, 0.5), ("gauss", 1.0, 0.5, 2.0), ("lineal", 0.3, 0.0, -1.0),
             ("lineal", 0.01, 0.0, -0.5)]
    res = []
    for c in casos:
        o = analiza(*c)
        res.append(o)
        if o["fin"] == "horizonte":
            print(f"{c}: a={o['a']:+.4g} l_s={o['l_s']:.5g} (l_s|a| = {o['ls_a']:.3f})  N ~ (l_s-l)^{o['p']:.3f}  "
                  f"r_s={o['r_s']:.4g} r'_s={o['rp_s']:+.6f} (relacion: {o['rp_relacion']:.6f})  F'phi'={o['Fp_phip']:+.3g}  "
                  f"phi'_fin={o['phip_end']:+.3g}  K (l_s-l)^2 -> {', '.join(f'{k:.3g}' for k in o['K_x2'][-3:])}  "
                  f"ligadura ll <= {o['ligadura_max']:.1e}")
        else:
            print(f"{c}: a={o['a']:+.4g} fin={o['fin']} l_end={o['l_end']:.4g}")
    os.makedirs("resultados", exist_ok=True)
    json.dump(dict(horizonte=res, patada=pat), open("resultados/singularidad.json", "w", encoding="utf-8"), indent=1)
