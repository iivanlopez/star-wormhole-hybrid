# -*- coding: utf-8 -*-
"""Verifica con la API de arXiv las referencias del articulo (titulo, autores, revista). Uso: python arxiv_bib.py"""
import re
import json
import time
import urllib.request

IDS = ["0806.0608", "1108.3003", "1202.4893", "1312.3204", "1402.6740", "1402.7026", "1404.6495", "1509.08808",
       "1510.07400", "1603.05806", "1705.09662", "1711.04152", "1712.04398", "1801.01787", "1803.11444", "1811.05481",
       "1811.05832", "1812.07022", "1901.07183", "1907.00699", "1907.09463", "2003.11552", "2104.08049", "2111.09857",
       "2212.05969", "2404.06297", "2409.15606", "2412.19408", "gr-qc/0211026"]
url = "https://export.arxiv.org/api/query?max_results=100&id_list=" + ",".join(IDS)
t = urllib.request.urlopen(url, timeout=90).read().decode()
out = {}
for e in t.split("<entry>")[1:]:
    i = re.search(r"<id>http://arxiv.org/abs/(\S+?)</id>", e).group(1)
    ti = " ".join(re.search(r"<title>(.*?)</title>", e, re.S).group(1).split())
    au = re.findall(r"<name>(.*?)</name>", e)
    jr = re.search(r"<arxiv:journal_ref[^>]*>(.*?)</arxiv:journal_ref>", e, re.S)
    yr = re.search(r"<published>(\d{4})", e).group(1)
    out[i] = dict(titulo=ti, autores=au, revista=" ".join(jr.group(1).split()) if jr else "", anio=yr)
    print(i, "|", yr, "|", ti, "|", ", ".join(au[:6]), "|", out[i]["revista"])
json.dump(out, open("resultados/arxiv_bib.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
