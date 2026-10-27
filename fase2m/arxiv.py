# -*- coding: utf-8 -*-
"""Consulta la API de arXiv (verificacion de citas).  Uso: python arxiv.py "consulta" [n]   o   python arxiv.py id:XXXX.XXXXX"""
import re
import sys
import urllib.parse
import urllib.request

q = sys.argv[1]
n = int(sys.argv[2]) if len(sys.argv) > 2 else 20
if q.startswith("id:"):
    url = "https://export.arxiv.org/api/query?id_list=" + q[3:]
else:
    url = "https://export.arxiv.org/api/query?search_query=" + urllib.parse.quote(q.replace(" ", "+"), safe=":+()") + f"&max_results={n}"
t = urllib.request.urlopen(url, timeout=60).read().decode()
for e in t.split("<entry>")[1:]:
    i = re.search(r"<id>http://arxiv.org/abs/(\S+?)</id>", e).group(1)
    ti = " ".join(re.search(r"<title>(.*?)</title>", e, re.S).group(1).split())
    au = re.findall(r"<name>(.*?)</name>", e)
    jr = re.search(r"<arxiv:journal_ref[^>]*>(.*?)</arxiv:journal_ref>", e, re.S)
    print(i, "|", ti, "|", ", ".join(au[:5]), "|", jr.group(1).strip() if jr else "")
    if q.startswith("id:"):
        print("   ", " ".join(re.search(r"<summary>(.*?)</summary>", e, re.S).group(1).split())[:1500])
