#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
三分类突变谱差异分析（Beat AML）
比较 de novo / therapy-related / secondary 三组的核心基因突变率
基因来源：clinical_sample.json 的临床检测字段（NPM1/FLT3-ITD/TP53/RUNX1/ASXL1/CEBPA）
"""
import json
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

DATA = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/原始数据"
OUT = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/结果"

# ---------- 读 patient，三分类 ----------
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

# ---------- 读 sample，分子字段映射到 patient（任意样本有突变=有突变）----------
with open(f"{DATA}/clinical_sample.json") as f:
    raws = json.load(f)

# 突变判定规则
def is_mut(v):
    s = str(v).strip()
    return s not in ("", "NA", "N/A", "Negative", "Wildtype", "None", "WT", "wt")

gene_rules = {
    "NPM1": lambda v: str(v).strip().lower() == "positive",
    "FLT3_ITD": lambda v: str(v).strip().lower() == "positive",
    "TP53": is_mut,
    "RUNX1": is_mut,
    "ASXL1": is_mut,
    "CEBPA_BIALLELIC": lambda v: str(v).strip().lower() == "bi",
}

patient_mut = defaultdict(set)  # pid -> 突变基因集合
for r in raws:
    pid = r.get("patientId")
    attr = r.get("clinicalAttributeId")
    val = r.get("value")
    if not pid or attr not in gene_rules:
        continue
    if gene_rules[attr](val):
        patient_mut[pid].add(attr)

# ---------- 汇总三组突变率 ----------
order = ["de novo", "therapy-related", "secondary"]
groups = defaultdict(list)
for pid, d in patient.items():
    groups[classify(d)].append(pid)

print("=== 三分类核心基因突变率 ===")
results = []
for gene in ["NPM1", "FLT3_ITD", "TP53", "RUNX1", "ASXL1", "CEBPA_BIALLELIC"]:
    row = {"gene": gene}
    counts = []
    for g in order:
        pids = groups[g]
        n_mut = sum(1 for p in pids if gene in patient_mut[p])
        n = len(pids)
        counts.append(n_mut)
        row[g] = f"{n_mut}/{n} ({n_mut/n*100:.1f}%)"
    # 卡方检验（三组）
    table = [[counts[i], len(groups[order[i]]) - counts[i]] for i in range(3)]
    try:
        chi2, p, dof, _ = chi2_contingency(table)
        row["chi2_p"] = f"{p:.4g}"
    except ValueError:
        from scipy.stats import fisher_exact
        ps = []
        for gi in [1, 2]:
            tbl = [[counts[0], len(groups[order[0]]) - counts[0]],
                   [counts[gi], len(groups[order[gi]]) - counts[gi]]]
            _, p2 = fisher_exact(tbl)
            ps.append(p2)
        row["chi2_p"] = f"Fisher dn-vs-th {ps[0]:.3g} / dn-vs-sec {ps[1]:.3g}"
    results.append(row)
    print(f"  {gene:10s} | de novo {row['de novo']:18s} | therapy {row['therapy-related']:18s} | secondary {row['secondary']:18s} | {row['chi2_p']}")

# ---------- 存结果 ----------
res_df = pd.DataFrame(results)
res_df.to_csv(f"{OUT}/突变谱差异.csv", index=False, encoding="utf-8-sig")
print(f"\n结果已存: {OUT}/突变谱差异.csv")
print("完成。")
