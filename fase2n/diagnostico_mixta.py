# -*- coding: utf-8 -*-
"""
Fase 2n: donde falla el monopolo en un fondo de la clase mixta.  Reutiliza fase2k/criba_masa.evalua (con la masa de
criba_mixta.py), devolviendo ademas los perfiles s2, s0 y el discriminante en la malla.
Uso: python diagnostico_mixta.py  (lee los fondos con beta > 0 de resultados/criba_mixta.json o del log)
"""
import os
import sys
import ast
import json
import inspect
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import criba_mixta as cx                                 # noqa: E402
cm = cx.cm
src = inspect.getsource(cm.evalua).replace("def evalua(", "def evalua_perfiles(")
cab, ult = src.rsplit("\n    return out\n", 1)
src = cab + "\n    out.update(xs=xs, s2=s2, s0=s0, disc=disc, b=b, v=V0 * v, rK=rK, QR=QR)\n    return out\n" + ult
exec(src, cm.__dict__)


def intervalos(xs, mask):
    """intervalos [a, b] de xs donde mask es cierto."""
    out, ini = [], None
    for x_, m_ in zip(xs, mask):
        if m_ and ini is None:
            ini = x_
        if not m_ and ini is not None:
            out.append((round(float(ini), 3), round(float(x_), 3)))
            ini = None
    if ini is not None:
        out.append((round(float(ini), 3), round(float(xs[-1]), 3)))
    return out


def diagnostica(prm):
    c = cm.evalua_perfiles(prm)
    if not c.get("B"):
        return dict(B=False)
    xs = c["xs"]
    return dict(beta0=c["beta0"], vmax=c["vmax"], s2_neg=intervalos(xs, c["s2"] <= 0),
                s0_pos=intervalos(xs, c["s0"] >= 0), no_hip=intervalos(xs, c["disc"] < -1e-8),
                rho_vmax=float(xs[np.argmax(np.abs(c["v"]))]))


if __name__ == "__main__":
    lineas = [l_ for l_ in open(os.path.join(AQUI, "criba_mixta_log.txt"), encoding="utf-8") if l_.startswith("[")]
    claves = ("M", "n", "MN", "w", "a", "wq", "wn", "w2", "lam")
    RES = []
    for l_ in lineas:
        d = ast.literal_eval(l_.split("] ", 1)[1])
        prm = {k: float(d[k]) for k in claves}
        r = diagnostica(prm)
        RES.append(dict(prm=prm, **r))
        print(prm, "\n   ", r, flush=True)
    json.dump(RES, open(os.path.join(AQUI, "resultados", "diagnostico_mixta.json"), "w"), indent=1, default=float)
