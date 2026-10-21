# -*- coding: utf-8 -*-
"""
Autofuncion de un modo con crecimiento evanescente fuerte (p. ej. el modo rapido del Modelo L):
ortonormalizacion de Godunov por tramos + sustitucion hacia atras c_k = R_{k+1}^{-1} c_{k+1}.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp


def autofuncion_gs(f, modelo, w, n_por_tramo=60, rtol=1e-10):
    a, b = (np.array(v, dtype=float) for v in f.bases(modelo, w))
    Q = np.stack([a, b], axis=1)
    nrm = np.linalg.norm(Q, axis=0)
    Q = Q / nrm
    c0_escala = np.diag(1.0 / nrm)                     # X = [a b] diag(1/nrm) -> coeficientes en la base original

    def rhs(l, y, fluido):
        A = f.matriz_A(f.bg(f.y_bg(l), fluido), w, modelo, fluido)
        return (A @ y.reshape(6, 2)).ravel()

    def grande(l, y, fluido):
        return np.max(np.abs(y)) - 1e6
    grande.terminal = True

    tramos, Rs = [], []

    def avanzar(l, l1, Q, fluido):
        while l < l1 * (1 - 1e-14):
            sol = solve_ivp(rhs, (l, l1), Q.ravel(), args=(fluido,), events=[grande], method="DOP853",
                            rtol=rtol, atol=1e-12, dense_output=True)
            M = sol.y[:, -1].reshape(6, 2)
            q, r = np.linalg.qr(M)
            s = np.sign(np.diag(r))
            q, r = q * s, (r.T * s).T
            tramos.append((sol.t[0], sol.t[-1], sol.sol))
            Rs.append(r)
            Q, l = q, sol.t[-1]
        return Q

    Q = avanzar(0.0, f.ls, Q, True)
    l_end = f.l_max
    Q = avanzar(f.ls, l_end, Q, False)
    r, s, Phi, q, ff, psi = f.y_bg(l_end)[:6]
    F1 = Q[1]
    F2 = math.exp(Phi) * r * r * Q[4] - w * psi * psi * r**3 * math.exp(-Phi) * Q[2]
    B = np.vstack([F1, F2])
    c = np.linalg.svd(B)[2][-1]                        # vector nulo (aprox.) de las condiciones de contorno
    coefs = [None] * len(tramos)
    ck = c
    for k in range(len(tramos) - 1, -1, -1):           # c_k = R_{k+1}^{-1} c_{k+1}
        ck = np.linalg.solve(Rs[k], ck)
        coefs[k] = ck
    ls, X = [], []
    for (l0, l1, sol), ck in zip(tramos, coefs):
        x = np.linspace(l0, l1, n_por_tramo)
        Y = sol(x).reshape(6, 2, -1)
        X.append(np.einsum("ijk,j->ik", Y, ck))
        ls.append(x)
    return np.concatenate(ls), np.concatenate(X, axis=1), np.linalg.svd(B)[1]
