# -*- coding: utf-8 -*-
"""
Verificacion bibliografica con la API de arXiv (export.arxiv.org/api/query): metadatos de los trabajos citados y
busquedas de agujeros de gusano EsGB asimetricos / lisos / sin cascara.
Salida: resultados/arxiv_busqueda.json (y listado por pantalla).   Uso: python arxiv_busqueda.py
"""
import sys
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

NS = {"a": "http://www.w3.org/2005/Atom", "x": "http://arxiv.org/schemas/atom"}
API = "https://export.arxiv.org/api/query?"

CITADOS = ["1904.13091", "1108.3003", "1111.4049", "1804.11170", "1610.06070", "1503.01409", "2006.13008", "1708.02324",
           "1907.09463", "2204.08232", "2609.18475", "2412.20296", "1711.04152", "1509.08808", "1812.07022"]

BUSQUEDAS = [
    'abs:wormhole AND abs:"Gauss-Bonnet" AND abs:asymmetric',
    'abs:wormholes AND abs:"Gauss-Bonnet" AND abs:asymmetric',
    'abs:wormhole AND abs:"Gauss-Bonnet" AND abs:"thin shell"',
    'abs:wormhole AND abs:"Gauss-Bonnet" AND abs:smooth',
    'abs:wormhole AND abs:"Gauss-Bonnet" AND abs:"without matter"',
    'abs:wormhole AND abs:"scalar-Gauss-Bonnet"',
    'abs:wormholes AND abs:"scalar-Gauss-Bonnet"',
    'abs:wormhole AND abs:"dilatonic" AND abs:"Gauss-Bonnet"',
    'au:Kleihaus AND abs:wormhole AND abs:"Gauss-Bonnet"',
    'au:Kanti AND abs:wormhole',
    'au:Antoniou AND abs:wormhole',
    'au:Bakopoulos AND abs:wormhole',
    'abs:wormhole AND abs:"Gauss-Bonnet" AND abs:stability AND abs:radial',
    'abs:wormholes AND abs:Horndeski AND abs:"no-go"',
    'au:Evseev AND abs:wormholes',
    'au:Bronnikov AND au:Baleevskikh',
]


def consulta(params):
    url = API + urllib.parse.urlencode(params)
    for intento in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return ET.fromstring(r.read())
        except Exception as e:                       # la API limita la frecuencia: reintentar
            time.sleep(5 * (intento + 1))
            err = e
    raise err


def entradas(raiz):
    out = []
    for e in raiz.findall("a:entry", NS):
        j = e.find("x:journal_ref", NS)
        out.append(dict(id=e.find("a:id", NS).text.strip().split("/abs/")[-1],
                        titulo=" ".join(e.find("a:title", NS).text.split()),
                        autores=[a.find("a:name", NS).text for a in e.findall("a:author", NS)],
                        revista=j.text if j is not None else None,
                        fecha=e.find("a:published", NS).text[:10],
                        resumen=" ".join(e.find("a:summary", NS).text.split())))
    return out


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass
    res = {"fecha_consulta": time.strftime("%Y-%m-%d"), "citados": [], "busquedas": {}}
    res["citados"] = entradas(consulta({"id_list": ",".join(CITADOS), "max_results": 50}))
    for e in res["citados"]:
        print(f"[{e['id']}] {e['titulo']} — {', '.join(e['autores'])} — {e['revista']} ({e['fecha']})")
    for q in BUSQUEDAS:
        time.sleep(3.5)
        lst = entradas(consulta({"search_query": q, "max_results": 100, "sortBy": "submittedDate",
                                 "sortOrder": "descending"}))
        res["busquedas"][q] = lst
        print(f"\n### {q}: {len(lst)} resultados")
        for e in lst:
            print(f"  [{e['id']}] {e['fecha']} {e['titulo']} — {', '.join(e['autores'][:4])}")
    json.dump(res, open("resultados/arxiv_busqueda.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
