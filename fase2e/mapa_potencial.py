# -*- coding: utf-8 -*-
"""
Mapa de existencia de agujeros lisos SIMETRICOS de EsGB con potencial (potencial.LadoV), con a = Phi'(0) = 0:
  clase A (escalar par):    phi'(0) = 0, garganta en phi0 = x, V(phi0) = 2/r0^2
  clase B (kink, impar):    phi0 = 0, phi'(0) = u = x, F y V pares, V(0) = 2/r0^2 + u^2/2
La amplitud del potencial la fija la ligadura (potencial.amp_garganta).  Por la simetria, el agujero existe si y solo
si la mitad l > 0 es AF: el escalar tiene que llegar al vacio phi_v (masivo: un disparo, la rama creciente e^{+m l}
debe anularse).

Para cada familia y cada acoplamiento se barre x y se clasifica el primer suceso del escalar (LadoV.tiro):
  +1 cruza el vacio (se pasa), -1 se vuelve antes (no llega), 0 final singular antes de decidir.
Cada cambio +1/-1 entre vecinos se biseca hasta |dx| < 3e-15 x: es la separatriz, la unica candidata a mitad AF.
Se registra el estado con que llega al vacio (r', N, chi, det M, masa de Misner-Sharp) y cuanto aguanta.
Ademas, en cada mitad integrada se localizan los ceros regulares de chi = r r' + 2 F_l (1 - 3 r'^2) y se comprueba
K = -16 F'^2 (1 - r'^2)^2 chi'/r'^2 (teorema_chi_lapso.py) y el signo de lambda = N^2 G/K.
Familias (r0 = 1):
  B-cuad   V = lam (phi^2 - v^2)^2, F = alpha phi^2       v en {0.3, 0.5, 1, 2}
  B-gauss  V = lam (phi^2 - v^2)^2, F = alpha e^{-phi^2}  v en {0.3, 0.5, 1, 2}
  A-lin    V = lam (phi^2 - v^2)^2, F = alpha phi         v en {0.3, 1}, phi0 en (0, v) y (v, 3v)
  A-exp    V = lam (phi^2 - v^2)^2, F = alpha e^{-phi}    v en {0.3, 1}
  A-cuad   V = lam (phi^2 - v^2)^2, F = alpha phi^2       v en {0.3, 1}
  A-masivo V = m^2 phi^2/2, F en {alpha phi^2, alpha e^{-phi}, alpha phi}
Salida: resultados/mapa_potencial.json.   Uso: python mapa_potencial.py   (4 procesos, ~1 h)
"""
import os
import sys
import json
import math
import time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.optimize import brentq

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import potencial as P      # noqa: E402

L_MAX = 300.0
NX = 40


def construir(fam, tipoF, alpha, tipoV, v, x, rtol=1e-11):
    """LadoV de la familia con el dato de garganta x; None si no hay potencial positivo o vacio en esa direccion."""
    Fd = P.acoplamiento(tipoF, alpha)
    Vkw = {"v": v} if tipoV == "doble_pozo" else {}
    Ud = P.forma(tipoV, **Vkw)
    phi0, phi1 = (0.0, x) if fam == "B" else (x, 0.0)
    if Ud(phi0)[0] <= 0:
        return None
    amp = P.amp_garganta(Ud, Fd, phi0, phi1, 0.0)
    if not amp > 0:
        return None
    Vd = P.potencial(tipoV, amp, **Vkw)
    x0, _ = P.segundas_V(np.array([1.0, 0.0, 0.0, 0.0, phi0, phi1]), Fd, Vd)
    direc = math.copysign(1.0, phi1 if phi1 != 0 else x0[2])
    vacios = (-v, v) if tipoV == "doble_pozo" else (0.0,)
    cands = [c for c in vacios if (c - phi0) * direc > 0]
    if not cands:
        return None
    phi_v = min(cands, key=lambda c: abs(c - phi0))
    L = P.LadoV(Fd, Vd, phi0, phi1, 0.0, phi_v=phi_v, l_max=L_MAX, rtol=rtol).resolver()
    L.amp = amp
    return L


def lambda_garganta(L):
    """velocidad radial^2 en la garganta (G/K con N = 1)."""
    f = P.DobleV(L, L)
    c0, c1, _ = P.dispersion(f, 0.0)
    return -c0 / c1


def magnitudes(L, l):
    y = L.sol.sol(l)
    x, det = P.segundas_V(y, L.Fd, L.Vd)
    r, rp, Phi, Pp, p, pp = y[:6]
    rpp, Ppp, ppp = x
    F0, F1, F2, F3 = L.Fd(p)
    Fl, Fll = F1 * pp, F2 * pp**2 + F1 * ppp
    A, B = r - 4 * Fl * rp, 1 - rp * rp
    chi = r * rp + 2 * Fl * (1 - 3 * rp * rp)
    dchi = rp * rp + r * rpp + 2 * Fll * (1 - 3 * rp * rp) - 12 * Fl * rp * rpp
    K = A**2 * r**2 - 16 * F1**2 * B * (4 * A * rpp + B * (1 - 4 * Fll))
    G = A**2 * r**2 - 16 * F1**2 * B * (4 * A * Pp * rp + B * (1 - 4 * Fl * Pp))
    num = B / r**2 + pp**2 / 4 - L.Vd(p)[0] / 2
    return dict(chi=chi, dchi=dchi, K=K, G=G, det=det, rp=rp, F1=F1, B=B, num=num, N=math.exp(Phi), r=r)


def ceros_chi(L, n=1500):
    """Ceros de chi a lo largo de la mitad (excluida la garganta) con K, -a chi' y G en cada uno."""
    if L.l_end < 1e-3:
        return []
    ls = np.linspace(L.l_end * 1e-3, L.l_end * (1 - 1e-4), n)
    ch = [magnitudes(L, x)["chi"] for x in ls]
    out = []
    for i in range(n - 1):
        if ch[i] * ch[i + 1] < 0:
            z = brentq(lambda x: magnitudes(L, x)["chi"], ls[i], ls[i + 1], xtol=1e-13)
            m = magnitudes(L, z)
            if abs(m["rp"]) < 1e-8:
                continue
            ka = -16 * m["F1"]**2 * m["B"]**2 * m["dchi"] / m["rp"]**2
            out.append(dict(l=z, sube=bool(m["dchi"] > 0), rp=m["rp"], num=m["num"], K=m["K"], K_teo=ka, G=m["G"],
                            lam=m["G"] * m["N"]**2 / m["K"] if m["K"] != 0 else float("nan")))
    return out


def llegada(L):
    """Estado en el punto de maxima aproximacion al vacio (phi cerca de phi_v y phi' pequeno)."""
    ls = np.linspace(1e-4, L.l_end, 3000)
    Y = L.sol.sol(ls)
    escala = abs(L.phi_v - L.phi0) + 1e-12
    i = int(np.argmin(np.abs(Y[4] - L.phi_v) / escala + np.abs(Y[5])))
    r, rp, Phi, Pp, p, pp, I = Y[:, i]
    m = magnitudes(L, ls[i])
    return dict(l=float(ls[i]), dphi=float(abs(p - L.phi_v) / escala), r=float(r), rp=float(rp),
                N=float(math.exp(Phi)), chi=float(m["chi"]), det=float(m["det"]), m=float(r * (1 - rp * rp) / 2),
                Pp=float(Pp))


def tarea(args):
    fam, tipoF, alpha, tipoV, v, xs = args
    puntos, tiros, Ls = [], [], []
    for x in xs:
        L = construir(fam, tipoF, alpha, tipoV, v, x)
        Ls.append(L)
        if L is None:
            tiros.append(9)
            # sin potencial positivo en la garganta o el escalar arranca alejandose de todos los vacios
            puntos.append(dict(x=x, valido=False))
            continue
        t = L.tiro()
        tiros.append(t)
        zs = ceros_chi(L)
        puntos.append(dict(x=x, valido=True, fin=L.fin, tiro=t, l_end=L.l_end, rpp0=L.rpp0, phi_v=L.phi_v,
                           lam0=lambda_garganta(L), amp=L.amp, ceros_chi=zs))
    seps = []
    for i in range(len(xs) - 1):
        if {tiros[i], tiros[i + 1]} != {-1, 1}:
            continue
        lo, hi, slo = xs[i], xs[i + 1], tiros[i]
        ok = True
        for _ in range(60):
            m = 0.5 * (lo + hi)
            Lm = construir(fam, tipoF, alpha, tipoV, v, m)
            s = Lm.tiro() if Lm is not None else 9
            if s == slo:
                lo = m
            elif s == -slo:
                hi = m
            else:
                ok = False
                break
            if abs(hi - lo) < 3e-15 * abs(m):
                break
        Llo = construir(fam, tipoF, alpha, tipoV, v, lo, rtol=1e-12)
        Lhi = construir(fam, tipoF, alpha, tipoV, v, hi, rtol=1e-12)
        ll, lh = llegada(Llo), llegada(Lhi)
        sana = (ll["rp"] > 0 and lh["rp"] > 0 and ll["chi"] > 0 and ll["det"] > 0 and ll["N"] > 1e-3
                and ll["dphi"] < 1e-3)
        seps.append(dict(x=lo, dx=hi - lo, biseccion_completa=ok, fin_lo=Llo.fin, fin_hi=Lhi.fin, l_lo=Llo.l_end,
                         l_hi=Lhi.l_end, llegada_lo=ll, llegada_hi=lh, sana=bool(sana), AF=Llo.fin == "AF" or Lhi.fin == "AF",
                         ceros_chi=ceros_chi(Llo)))
    return dict(fam=fam, F=tipoF, alpha=alpha, V=tipoV, v=v, puntos=puntos, separatrices=seps)


def tareas():
    out = []
    xB = [float(x) for x in np.geomspace(0.02, 20.0, NX)]
    for v in (0.3, 0.5, 1.0, 2.0):
        for al in (-0.07, -0.1, -0.2, -0.5, -1.0, -2.0, -5.0, 0.05, 0.2, 1.0):
            out.append(("B", "cuadratico", al, "doble_pozo", v, xB))
        for al in (0.07, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, -0.1, -0.5):
            out.append(("B", "gauss", al, "doble_pozo", v, xB))
    for v in (0.3, 1.0):
        xA = [float(x) for x in np.r_[np.linspace(0.005, 0.995, 30), np.linspace(1.01, 3.0, 20)] * v]
        for al in (-0.13, -0.1, -0.05, -0.02, 0.02, 0.05, 0.1, 0.13, 0.3, 1.0):
            out.append(("A", "lineal", al, "doble_pozo", v, xA))
        for al in (-1.0, -0.3, -0.1, -0.03, 0.03, 0.1, 0.3, 1.0):
            out.append(("A", "exp", al, "doble_pozo", v, xA))
            out.append(("A", "cuadratico", al, "doble_pozo", v, xA))
    xM = [float(x) for x in np.geomspace(0.05, 20.0, NX)]
    for tipoF, als in (("cuadratico", (-0.1, -0.03, 0.01, 0.02, 0.04, 0.06, 0.1, 0.3)),
                       ("exp", (-1.0, -0.3, -0.1, 0.03, 0.1, 0.3, 1.0)),
                       ("lineal", (-0.13, -0.05, 0.05, 0.1, 0.13, 0.3, 1.0))):
        for al in als:
            out.append(("A", tipoF, al, "masivo", 0.0, xM))
    return out


def resumen(res):
    print(f"{len(res)} barridos; {sum(len(o['puntos']) for o in res)} mitades")
    por = {}
    for o in res:
        k = f"{o['fam']}-{o['F']}-{o['V']}"
        d = por.setdefault(k, dict(mitades=0, fines={}, seps=0, sanas=0, AF=0, rp_llegada=[]))
        for p in o["puntos"]:
            if p["valido"]:
                d["mitades"] += 1
                d["fines"][p["fin"]] = d["fines"].get(p["fin"], 0) + 1
        for s in o["separatrices"]:
            d["seps"] += 1
            d["sanas"] += int(s["sana"])
            d["AF"] += int(s["AF"])
            if s["llegada_lo"]["dphi"] < 1e-3:
                d["rp_llegada"].append(s["llegada_lo"]["rp"])
    for o in res:
        k = f"{o['fam']}-{o['F']}-{o['V']}"
        por[k]["sin_vacio"] = por[k].get("sin_vacio", 0) + sum(not p["valido"] for p in o["puntos"])
    for k, d in por.items():
        rp = d.pop("rp_llegada")
        extra = f"; r' al llegar al vacio en [{min(rp):+.3f}, {max(rp):+.3f}] ({len(rp)} llegan)" if rp else ""
        print(f"  {k}: {d['mitades']} mitades (+{d['sin_vacio']} sin vacio hacia el que avanzar), finales "
              f"{dict(sorted(d['fines'].items()))}; separatrices {d['seps']}, sanas {d['sanas']}, AF {d['AF']}{extra}")
    zs = [z for o in res for p in o["puntos"] if p.get("valido") for z in p["ceros_chi"]]
    zs += [z for o in res for s in o["separatrices"] for z in s["ceros_chi"]]
    sub = [z for z in zs if z["sube"]]
    err = max((abs(z["K"] - z["K_teo"]) / (abs(z["K"]) + abs(z["K_teo"]) + 1e-300) for z in zs), default=0.0)
    print(f"  ceros regulares de chi: {len(zs)} ({len(sub)} ascendentes); max |K - K_teo|/(|K|+|K_teo|) = {err:.1e}; "
          f"max |num| = {max((abs(z['num']) for z in zs), default=0):.1e}")
    print(f"  en los ascendentes: K < 0 en {sum(z['K'] < 0 for z in sub)}, G > 0 en {sum(z['G'] > 0 for z in sub)}, "
          f"lambda < 0 en {sum(z['lam'] < 0 for z in sub)}")


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    t0 = time.time()
    ts = tareas()
    with ProcessPoolExecutor(max_workers=4) as ex:
        res = list(ex.map(tarea, ts, chunksize=1))
    json.dump(res, open(os.path.join(AQUI, "resultados", "mapa_potencial.json"), "w", encoding="utf-8"), indent=0)
    resumen(res)
    print(f"tiempo: {time.time() - t0:.0f} s")
