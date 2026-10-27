# -*- coding: utf-8 -*-
"""Comprobacion estructural del LaTeX sin compilador: llaves, entornos, cite/bibitem, ref/label, $ pares, marcadores."""
import re
import sys

p = sys.argv[1] if len(sys.argv) > 1 else "gargantas_monopolo.tex"
s = open(p, encoding="utf-8").read()
t = re.sub(r"(?<!\\)%.*", "", s)
err = []
nivel = 0
for i, ch in enumerate(re.sub(r"\\[{}]", "", t)):
    nivel += (ch == "{") - (ch == "}")
    if nivel < 0:
        err.append(f"llave de cierre sin abrir cerca del caracter {i}")
        nivel = 0
if nivel:
    err.append(f"{nivel} llaves sin cerrar")
pila = []
for m in re.finditer(r"\\(begin|end)\{([^}]*)\}", t):
    if m.group(1) == "begin":
        pila.append(m.group(2))
    elif not pila or pila.pop() != m.group(2):
        err.append(f"entorno mal cerrado: {m.group(2)}")
if pila:
    err.append(f"entornos sin cerrar: {pila}")
citas = {c.strip() for g in re.findall(r"\\cite\{([^}]*)\}", t) for c in g.split(",")}
bib = set(re.findall(r"\\bibitem\{([^}]*)\}", t))
if citas - bib:
    err.append(f"citas sin bibitem: {citas - bib}")
if bib - citas:
    print("aviso: bibitem sin citar:", bib - citas)
refs = set(re.findall(r"\\(?:eq)?ref\{([^}]*)\}", t))
labs = set(re.findall(r"\\label\{([^}]*)\}", t))
if refs - labs:
    err.append(f"ref sin label: {refs - labs}")
if re.sub(r"\\\$", "", t).count("$") % 2:
    err.append("numero impar de $")
for pal in ("NMIX", "RESUMENMIX", "TABLA_", "RESUMEN_"):
    if pal in t:
        err.append(f"marcador sin rellenar: {pal}")
print("\n".join(err) if err else f"{p}: estructura correcta ({len(bib)} referencias, {len(labs)} etiquetas)")
sys.exit(1 if err else 0)
