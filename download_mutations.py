#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""稳健下载 Beat AML 全量突变数据（分页 genes + 分块读取 + 重试）"""
import urllib.request, json, time

OUT = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/原始数据"

def fetch(url, data=None, headers=None, timeout=120, retries=3):
    h = {"User-Agent": "Mozilla/5.0"}
    if headers:
        h.update(headers)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers=h)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                chunks = []
                while True:
                    c = r.read(65536)
                    if not c:
                        break
                    chunks.append(c)
                return b"".join(chunks)
        except Exception as e:
            print(f"  尝试 {attempt+1}/{retries} 失败: {e}", flush=True)
            if attempt < retries - 1:
                time.sleep(4)
            else:
                raise

# 1. 分页查全部基因 entrezGeneId
all_genes = []
page = 0
while True:
    url = f"https://www.cbioportal.org/api/genes?pageSize=2000&pageNumber={page}&projection=ID"
    data = fetch(url, timeout=60)
    genes = json.loads(data)
    if not genes:
        break
    all_genes.extend(genes)
    print(f"genes 第 {page} 页: {len(genes)} 个, 累计 {len(all_genes)}", flush=True)
    if len(genes) < 2000:
        break
    page += 1

entrez = [g["entrezGeneId"] for g in all_genes]
print(f"总基因数: {len(entrez)}", flush=True)
with open(f"{OUT}/entrez_gene_ids.json", "w") as f:
    json.dump(entrez, f)

# 2. 查 sequenced 样本 ID
data = fetch("https://www.cbioportal.org/api/sample-lists/aml_ohsu_2022_sequenced/sample-ids", timeout=30)
sample_ids = json.loads(data)
print(f"样本数: {len(sample_ids)}", flush=True)

# 3. fetch 全量突变
url = "https://www.cbioportal.org/api/molecular-profiles/aml_ohsu_2022_mutations/mutations/fetch"
body = json.dumps({"sampleIds": sample_ids, "entrezGeneIds": entrez}).encode()
data = fetch(url, data=body, headers={"Content-Type": "application/json"}, timeout=900, retries=3)
print(f"mutations: {len(data)} bytes", flush=True)
with open(f"{OUT}/mutations.json", "wb") as f:
    f.write(data)
print(f"已保存: {OUT}/mutations.json", flush=True)
print("DONE", flush=True)
