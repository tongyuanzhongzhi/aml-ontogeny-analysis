#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""上传代码包到 Zenodo，生成永久 DOI"""
import os, sys, json, urllib.request, urllib.parse

TOKEN = os.environ.get("ZENODO_TOKEN")
if not TOKEN:
    print("错误：缺少 ZENODO_TOKEN 环境变量")
    sys.exit(1)

ZIP = "/Users/daizhensheng/Desktop/学术文章发表/已投文章/继发性AML三分类预后_AMJ/aml_ontogeny_analysis_code_20261124.zip"
BASE = "https://zenodo.org/api"

def api(method, path, data=None, headers_extra=None):
    url = f"{BASE}{path}?access_token={TOKEN}"
    headers = {"Content-Type": "application/json"}
    if headers_extra:
        headers.update(headers_extra)
    body = json.dumps(data).encode() if data is not None else None
    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        print(f"HTTP {e.code}: {e.read().decode()[:500]}")
        raise

# 1. 创建 deposition
dep = api("POST", "/deposit/depositions", data={})
dep_id = dep["id"]
bucket_url = dep["links"]["bucket"]
print(f"Deposition 创建: id={dep_id}")

# 2. 填元数据
metadata = {
    "title": "Analysis code: Clinical ontogeny and MR gene mutations in acute myeloid leukemia",
    "upload_type": "software",
    "publication_type": "other",
    "description": (
        "Reproducible analysis scripts for the manuscript 'Clinical ontogeny and MR gene "
        "mutations are independently associated with ex vivo drug resistance and poor "
        "prognosis in acute myeloid leukemia'. Includes sequentially adjusted and stratified "
        "Cox models, four-factor type II ANOVA of ex vivo drug sensitivity, initial-diagnosis "
        "sensitivity analyses, ELN 2022 re-derivation, mutation-spectrum analysis, and figure "
        "generation. Data sources: Beat AML (aml_ohsu_2022, cBioPortal) and TCGA LAML "
        "(laml_tcga_pub, cBioPortal)."
    ),
    "creators": [
        {"name": "Dai, Zhensheng", "affiliation": "Shanghai Pudong Hospital, Fudan University Pudong Medical Center", "orcid": "0009-0009-7476-8163"},
        {"name": "Jin, Yueling", "affiliation": "Shanghai Science and Technology Museum"},
    ],
    "keywords": ["acute myeloid leukemia", "secondary AML", "clinical ontogeny",
                 "myelodysplasia-related mutations", "ex vivo drug resistance", "venetoclax"],
    "license": "MIT",
    "access_right": "open",
}
dep = api("PUT", f"/deposit/depositions/{dep_id}", data={"metadata": metadata})
print("元数据已填写")

# 3. 上传 zip 文件
filename = os.path.basename(ZIP)
with open(ZIP, "rb") as f:
    zip_bytes = f.read()
upload_url = f"{bucket_url}/{urllib.parse.quote(filename)}?access_token={TOKEN}"
req = urllib.request.Request(upload_url, data=zip_bytes, method="PUT")
req.add_header("Content-Type", "application/octet-stream")
with urllib.request.urlopen(req) as r:
    upload_resp = json.load(r)
print(f"文件已上传: {filename} ({len(zip_bytes)} bytes)")

# 4. 发布
dep = api("POST", f"/deposit/depositions/{dep_id}/actions/publish")
doi = dep.get("doi")
print(f"\n=== 发布成功 ===")
print(f"DOI: {doi}")
print(f"URL: https://doi.org/{doi}")
print(f"Concept DOI: {dep.get('conceptdoi')}")
