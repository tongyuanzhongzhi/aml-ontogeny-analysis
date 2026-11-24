#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 Table 3（突变谱，补 SRSF2 等 MR 基因）+ Table 4（药敏）"""
import json
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact

DATA = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/其他文件/原始数据"
OUT = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/投稿包"

with open(f"{DATA}/gene_symbol_map.json") as f:
    id2symbol = {int(k): v for k, v in json.load(f).items()}
with open(f"{DATA}/mutations.json") as f:
    mutations = json.load(f)
sample_genes = defaultdict(set)
for m in mutations:
    sym = id2symbol.get(m["entrezGeneId"], "")
    if sym:
        sample_genes[m["sampleId"]].add(sym)

with open(f"{DATA}/clinical_sample.json") as f:
    raw_sample = json.load(f)
sample_info = defaultdict(dict)
for r in raw_sample:
    sid = r.get("sampleId"); a = r.get("clinicalAttributeId"); v = r.get("value")
    if sid and a and v not in (None, "", "NA", "N/A"):
        sample_info[sid][a] = v

with open(f"{DATA}/clinical_patient.json") as f:
    rawp = json.load(f)
patient = defaultdict(dict)
for r in rawp:
    pid = r.get("patientId"); a = r.get("clinicalAttributeId"); v = r.get("value")
    if pid and a and (a not in patient[pid] or patient[pid][a] in (None, "", "NA", "N/A")):
        patient[pid][a] = v

def yes(d, k): return str(d.get(k, "")).strip().lower() == "yes"
def classify(d):
    if yes(d, "PRIOR_MALIGNANCY_NON_MYELOID"): return "therapy-related"
    if yes(d, "PRIOR_MDS") or yes(d, "PRIOR_MDS_MPN") or yes(d, "PRIOR_MPN"): return "secondary"
    return "de novo"

pid2samples = defaultdict(list)
for r in raw_sample:
    sid = r.get("sampleId"); pid = r.get("patientId")
    if sid and pid:
        pid2samples[pid].append(sid)

def baseline_samples(pid):
    samples = pid2samples.get(pid, [])
    initial = [s for s in samples
               if sample_info[s].get("DISEASE_STAGE_AT_SPECIMEN_COLLECTION", "") == "Initial Diagnosis"]
    return initial if initial else samples

def patient_mut(pid, gene):
    return any(gene in sample_genes.get(s, set()) for s in baseline_samples(pid))

# 患者三分类
groups = defaultdict(list)
for pid, d in patient.items():
    groups[classify(d)].append(pid)

# Table 3：MR 基因 + 其他
genes = ["SRSF2", "SF3B1", "U2AF1", "ZRSR2", "ASXL1", "EZH2", "STAG2", "RUNX1", "BCOR",
         "TP53", "NPM1", "FLT3"]
rows = []
for g in genes:
    freqs = {}
    for grp in ["de novo", "therapy-related", "secondary"]:
        pids = groups[grp]
        n_mut = sum(1 for pid in pids if patient_mut(pid, g))
        freqs[grp] = (n_mut, len(pids))
    # chi2 de novo vs secondary
    dn_mut, dn_n = freqs["de novo"]
    sec_mut, sec_n = freqs["secondary"]
    table = [[sec_mut, sec_n - sec_mut], [dn_mut, dn_n - dn_mut]]
    try:
        odds, p = fisher_exact(table)
    except:
        p = np.nan
    rows.append({
        "Gene": g,
        "de novo": f"{dn_mut}/{dn_n} ({dn_mut/dn_n*100:.1f}%)",
        "therapy-related": f"{freqs['therapy-related'][0]}/{freqs['therapy-related'][1]} ({freqs['therapy-related'][0]/freqs['therapy-related'][1]*100:.1f}%)",
        "secondary": f"{sec_mut}/{sec_n} ({sec_mut/sec_n*100:.1f}%)",
        "P (sec vs de novo)": f"{p:.2e}" if not np.isnan(p) else "—",
    })

table3 = pd.DataFrame(rows)
with pd.ExcelWriter(f"{OUT}/Table3.xlsx", engine="openpyxl") as w:
    table3.to_excel(w, index=False, sheet_name="Table 3")
print("Table3.xlsx 完成")
print(table3.to_string(index=False))

# Table 4：药敏 FDR 显著药
drug = pd.read_excel(f"{OUT}/TableS2_药敏.xlsx")
print("\n=== Table 4 药敏 ===")
print(drug.head(15).to_string(index=False))
# 重命名 Table4
drug.to_excel(f"{OUT}/Table4.xlsx", index=False, sheet_name="Table 4")
print("\nTable4.xlsx 完成")
