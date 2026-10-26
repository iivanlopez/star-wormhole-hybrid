# -*- coding: utf-8 -*-
"""
Fase 2h: espectro completo del monopolo del hibrido (fluido + escalar/metrica) con elementos finitos de Hermite cubicos.

Variables: h (d = delta r = N h) en todo el dominio y c (= delta g_ll/2) en la region con fluido [0, lc]; fuera la
ligadura de momento sin fluido da W = 0, c = chi1 h' con chi1 = 4 r Gcal/gam.  Formas (Lagrangiano sin integrar por
partes; reduccion_fluido.py y coeficientes_fluido.py de la Fase 2g):
    dentro:  E = Q^2/(4 kap_a) - 2 N^3 Gcal (h' - sig h)^2 - k_c c^2/2 - c (e0 N h + e1 (N h)')
             T = (al1 N h' + gam c)^2/(4 kap_b) - 2 G4 N h^2 - e_tt N c h
             Q = -(4 N^2 r Gcal h')' - 2 N^2 r (rho~+p~) h + (N^2 gam c)' - N (rho~+p~) r^2 c
             1/(4 kap_a) = c_s^2/(2 N r^2 (rho~+p~)),   1/(4 kap_b) = N/(2 r^2 (rho~+p~))
    fuera:   E = -2 N^3 Gcal (h' - sig h)^2 - k_c (chi1 h')^2/2 - chi1 h' (e0 N h + e1 (N h)')
             T = -2 G4 N h^2 - e_tt N chi1 h' h
Condiciones: paridad en la garganta (par: h' = c' = 0; impar: h = c = 0), W(lc) = 0 (c(lc) = chi1 h'(lc); la masa por
encima de lc es despreciable) y h(L) = 0 (caja grande).  Problema generalizado E v = w T v: T definida positiva (sin
fantasmas) y los autovalores negativos son los taquiones.
Modos de comprobacion:
    formulacion "vacio" en todo el dominio (sin fluido), para comparar con Sturm-Liouville (K, P, U de T3);
    restriccion c = 0 dentro (la familia de la Fase 2g).
"""
import os
import sys
import math
import numpy as np
from scipy.linalg import eigh

AQUI = os.path.dirname(os.path.abspath(__file__))
for p_ in (os.path.join(AQUI, "..", "fase2f"), os.path.join(AQUI, "..", "fase2g"), AQUI):
    sys.path.insert(0, p_)
import criba as CR                    # noqa: E402

GX, GW = np.polynomial.legendre.leggauss(8)


def lc_fluido(T, umbral=1e-6):
    """Ultimo punto donde (rho~+p~) > umbral (rho~+p~)(0)."""
    F = T.F
    lt = np.linspace(0, F.ls, 40001)
    rp = F.fluido(F.lapso(lt)[0])[0]
    dentro = np.nonzero(rp > umbral * rp[0])[0]
    return float(lt[dentro[-1]])


def mallas(lc, n_in, L, n_out, h0=None):
    a = np.linspace(0.0, lc, n_in + 1) if lc > 0 else np.array([0.0])
    h0 = h0 if h0 is not None else (lc / n_in if lc > 0 else 0.02)
    # geometrica desde lc con primer paso h0 hasta L
    q = 1.0
    if n_out > 0:
        from scipy.optimize import brentq
        f = lambda q: h0 * (q ** n_out - 1) / (q - 1) - (L - lc) if q != 1 else h0 * n_out - (L - lc)  # noqa: E731
        q = brentq(f, 1.0 + 1e-12, 2.0)
    pasos = h0 * q ** np.arange(n_out)
    b = lc + np.concatenate([[0.0], np.cumsum(pasos)])
    b[-1] = L
    return a, b


class Coefs:
    """Coeficientes en los puntos de Gauss de todos los elementos (evaluacion por trozos)."""
    def __init__(self, T, nodos, trozo=400):
        t = 0.5 * (GX + 1)
        H = np.diff(nodos)
        self.lg = (nodos[:-1, None] + t[None, :] * H[:, None]).ravel()
        acc = {}
        for i in range(0, self.lg.size, trozo):
            c, _ = T.coeficientes(self.lg[i:i + trozo])
            for k, v in c.items():
                acc.setdefault(k, []).append(np.asarray(v))
        self.c = {k: np.concatenate(v) for k, v in acc.items()}
        self.ng = len(GX)

    def elem(self, e):
        s = slice(e * self.ng, (e + 1) * self.ng)
        return {k: v[s] for k, v in self.c.items()}


def ensamblar(T, paridad, lc, n_in=120, L=1500.0, n_out=200, vacio=False, c_cero=False, coefs=None):
    a, b = mallas(lc, n_in, L, n_out)
    nodos = np.concatenate([a, b[1:]]) if lc > 0 else b
    Kn = len(nodos)
    Ki = len(a) if (lc > 0 and not vacio) else 1
    nh = 2 * Kn
    nd = nh + 2 * Ki
    E = np.zeros((nd, nd))
    Tm = np.zeros((nd, nd))
    t = 0.5 * (GX + 1)
    cf = coefs if coefs is not None else Coefs(T, nodos)
    for e in range(Kn - 1):
        x0, x1 = nodos[e], nodos[e + 1]
        H = x1 - x0
        w = 0.5 * GW * H
        f = cf.elem(e)
        v, d1, d2 = CR.hermite(t, H)
        ih = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
        N, r, Gc, sig = f["N"], f["r"], f["Gc"], f["sig"]
        Dm = d1 - sig * v
        Nh_p = f["Np"] * v + N * d1
        if e < Ki - 1:
            ic = [nh + 2 * e, nh + 2 * e + 1, nh + 2 * e + 2, nh + 2 * e + 3]
            idx = ih + ic
            rp = np.maximum(f["rp"], 1e-300)
            AQ = f["cs2"] / (2 * N * r * r * rp)
            AB = N / (2 * r * r * rp)
            Qh = -f["a2"] * d2 - f["a2p"] * d1 - 2 * N * N * r * f["rp"] * v
            Qc = f["Ng"] * d1 + f["Ngp"] * v - N * f["rp"] * r * r * v
            Qv = np.vstack([Qh, Qc])
            Dv = np.vstack([Dm, 0 * Dm])
            cv = np.vstack([0 * v, v])
            Eloc = (Qv * AQ * w) @ Qv.T - (Dv * 2 * N ** 3 * Gc * w) @ Dv.T - (cv * f["kc"] / 2 * w) @ cv.T
            acop = np.vstack([f["e0"] * N * v + f["e1"] * Nh_p, 0 * v])
            Eloc -= 0.5 * ((cv * w) @ acop.T + (acop * w) @ cv.T)
            Wv = np.vstack([f["al1"] * N * d1, f["gam"] * v])
            hv = np.vstack([v, 0 * v])
            Tloc = (Wv * AB * w) @ Wv.T - (hv * 2 * f["G4"] * N * w) @ hv.T
            cr = (cv * (-f["ett"] * N / 2) * w) @ hv.T
            Tloc += cr + cr.T
            E[np.ix_(idx, idx)] += Eloc
            Tm[np.ix_(idx, idx)] += Tloc
        else:
            chi1 = 4 * r * Gc / f["gam"]
            cvec = chi1 * d1
            Eloc = -(Dm * 2 * N ** 3 * Gc * w) @ Dm.T - (cvec * f["kc"] / 2 * w) @ cvec.T
            acop = f["e0"] * N * v + f["e1"] * Nh_p
            Eloc -= 0.5 * ((cvec * w) @ acop.T + (acop * w) @ cvec.T)
            Tloc = -(v * 2 * f["G4"] * N * w) @ v.T
            cr = (cvec * (-f["ett"] * N / 2) * w) @ v.T
            Tloc += cr + cr.T
            E[np.ix_(ih, ih)] += Eloc
            Tm[np.ix_(ih, ih)] += Tloc
    # transformacion de grados de libertad: c(lc) = chi1(lc) h'(lc)
    quitar = set()
    Tr = np.eye(nd)
    if Ki > 1:
        fl = T.coeficientes(np.array([lc]))[0]
        chi_c = float(4 * fl["r"][0] * fl["Gc"][0] / fl["gam"][0])
        jc = nh + 2 * (Ki - 1)
        jh1 = 2 * (Ki - 1) + 1
        Tr[jc, jc] = 0.0
        Tr[jc, jh1] = chi_c
        quitar.add(jc)
        if c_cero:
            for j in range(nh, nd):
                quitar.add(j)
    else:
        quitar.update(range(nh, nd))
    if paridad == "par":
        quitar.update({1, nh + 1})
    else:
        quitar.update({0, nh})
    quitar.add(2 * (Kn - 1))
    keep = [j for j in range(nd) if j not in quitar]
    Tr = Tr[:, keep]
    E2 = Tr.T @ E @ Tr
    T2 = Tr.T @ Tm @ Tr
    return 0.5 * (E2 + E2.T), 0.5 * (T2 + T2.T), nodos, cf


def espectro(T, paridad, lc, n_in=120, L=1500.0, n_out=200, k=6, **kw):
    E, Tm, nodos, cf = ensamblar(T, paridad, lc, n_in, L, n_out, **kw)
    tv = np.linalg.eigvalsh(Tm)
    out = dict(Tmin=float(tv[0]), Tmax=float(tv[-1]), ndof=E.shape[0])
    if tv[0] > 0:
        w, vec = eigh(E, Tm, subset_by_index=[0, k - 1])
        out["w"] = [float(x) for x in w]
    else:
        # fantasmas: autovalores negativos de T (en la norma de la masa)
        out["n_T_neg"] = int(np.sum(tv < 0))
    return out
