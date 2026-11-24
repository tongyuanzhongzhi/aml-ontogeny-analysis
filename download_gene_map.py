#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""分页查 genes，构建 entrezGeneId -> hugoGeneSymbol 映射"""
import urllib.request, json, time

OUT = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/原始数据"

id2symbol = {}
page = 0
while True:
    url = f"https://www.cbioportal.org/api/genes?pageSize=2000&pageNumber={page}&projection=SUMMARY"
    genes = None
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                chunks = []
                while True:
                    c = r.read(65536)
                    if not c:
                        break
                    chunks.append(c)
                genes = json.loads(b"".join(chunks))
            break
        except Exception as e:
            if attempt == 2:
                raise
            time.sleep(3)
    if not genes:
        break
    for g in genes:
        id2symbol[g["entrezGeneId"]] = g.get("hugoGeneSymbol", "")
    print(f"page {page}: 累计 {len(id2symbol)}", flush=True)
    if len(genes) < 2000:
        break
    page += 1

with open(f"{OUT}/gene_symbol_map.json", "w") as f:
    json.dump(id2symbol, f)
print(f"映射完成: {len(id2symbol)}", flush=True)
print("DONE", flush=True)
