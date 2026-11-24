#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""原创点3：红系偏倚 → BCL2/BCL-XL 转换 → Venetoclax 耐药 机制闭环验证"""
import json
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

DATA = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/原始数据"
OUT = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/结果"

# 三分类
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
subj2group = {}
for pid, d in patient.items():
    subj = pid.split("_")[-1]
    if subj.isdigit():
        subj2group[int(subj)] = classify(d)

mapping = pd.read_excel(f"{DATA}/BeatAML2_harmonized/beataml_waves1to4_sample_mapping.xlsx")
mapping = mapping[~mapping["rna_control"].astype(str).str.contains("Healthy", na=False)]
rna2subject = {}
for _, row in mapping.iterrows():
    rna = row["dbgap_rnaseq_sample"]; subj = row["dbgap_subject_id"]
    if pd.notna(rna) and pd.notna(subj):
        rna2subject[rna] = int(subj)

print("读表达矩阵（281MB）...")
exp = pd.read_csv(f"{DATA}/BeatAML2_harmonized/beataml_waves1to4_norm_exp_dbgap.txt", sep="\t", index_col=0, low_memory=False)
meta_cols = ["display_label", "description", "biotype"]
gene_labels = exp["display_label"].copy()
sample_cols = [c for c in exp.columns if c not in meta_cols]

# 样本 -> group
sample_group = {}
for c in sample_cols:
    subj = rna2subject.get(c)
    if subj is not None:
        sample_group[c] = subj2group.get(subj, None)
mapped = {c: g for c, g in sample_group.items() if g is not None}
print(f"映射样本: {len(mapped)}, 分布: {dict(Counter(mapped.values()))}")

# 目标基因：BCL2 家族 + 红系/巨核系标志
targets = {
    "BCL2": "anti-apoptotic (venetoclax target)",
    "BCL2L1": "BCL-XL anti-apoptotic",
    "MCL1": "MCL1 anti-apoptotic",
    "BCL2L2": "BCL-W anti-apoptotic",
    "GATA1": "erythroid TF",
    "TAL1": "erythroid TF",
    "EPOR": "erythropoietin receptor",
    "KLF1": "erythroid TF",
    "GYPA": "glycophorin A (erythroid marker)",
    "ITGA2B": "CD41 (megakaryocytic marker)",
    "GP9": "megakaryocytic marker",
    "HBG1": "fetal hemoglobin",
    "HBA1": "alpha globin",
    "BCL2A1": "BCL2A1 anti-apoptotic",
}
sym2idx = {}
for i, s in enumerate(gene_labels.values):
    if str(s) in targets and str(s) not in sym2idx:
        sym2idx[str(s)] = i

exp_mat = exp[sample_cols].apply(pd.to_numeric, errors="coerce")
groups = ["de novo", "therapy-related", "secondary"]

print("\n=== BCL2 家族 + 红系标志表达（三分类中位数）===")
print(f"{'基因':10s} {'de novo':>9s} {'therapy':>9s} {'secondary':>9s} {'sec-vs-dn P':>12s}")
for sym in targets:
    if sym not in sym2idx:
        print(f"{sym:10s} [未找到]")
        continue
    i = sym2idx[sym]
    meds = {}
    for g in groups:
        cols = [c for c, gg in mapped.items() if gg == g]
        meds[g] = np.nanmedian(exp_mat.iloc[i][cols].values)
    dn = exp_mat.iloc[i][[c for c, gg in mapped.items() if gg == "de novo"]].values
    sec = exp_mat.iloc[i][[c for c, gg in mapped.items() if gg == "secondary"]].values
    dn = dn[~np.isnan(dn)]; sec = sec[~np.isnan(sec)]
    _, p = mannwhitneyu(sec, dn, alternative="two-sided")
    print(f"{sym:10s} {meds['de novo']:9.1f} {meds['therapy-related']:9.1f} {meds['secondary']:9.1f} {p:12.3g}")

# BCL2L1/BCL2 和 MCL1/BCL2 比值（per sample）
print("\n=== BCL2L1/BCL2 和 MCL1/BCL2 比值（三分类）===")
for num, den, label in [("BCL2L1", "BCL2", "BCL2L1/BCL2"), ("MCL1", "BCL2", "MCL1/BCL2")]:
    if num not in sym2idx or den not in sym2idx:
        continue
    i_num, i_den = sym2idx[num], sym2idx[den]
    print(f"\n{label}:")
    for g in groups:
        cols = [c for c, gg in mapped.items() if gg == g]
        ratios = []
        for c in cols:
            n = exp_mat.iloc[i_num][c]; d = exp_mat.iloc[i_den][c]
            if not np.isnan(n) and not np.isnan(d) and d > 0:
                ratios.append(n / d)
        print(f"  {g:15s} 中位比值 {np.median(ratios):.3f} (n={len(ratios)})")
    dn_r = [exp_mat.iloc[i_num][c]/exp_mat.iloc[i_den][c] for c, gg in mapped.items()
            if gg == "de novo" and exp_mat.iloc[i_den][c] > 0]
    sec_r = [exp_mat.iloc[i_num][c]/exp_mat.iloc[i_den][c] for c, gg in mapped.items()
             if gg == "secondary" and exp_mat.iloc[i_den][c] > 0]
    _, p = mannwhitneyu(sec_r, dn_r, alternative="two-sided")
    print(f"  secondary vs de novo: P={p:.3g}")

# ABT-737 药敏结果
print("\n=== BCL-XL 抑制剂 ABT-737 药敏 ===")
drug = pd.read_csv(f"{OUT}/药物敏感性差异_FDR校正.csv")
for name in ["ABT-737", "Venetoclax", "Navitoclax"]:
    r = drug[drug["inhibitor"].str.contains(name, case=False, na=False)]
    if len(r):
        for _, row in r.iterrows():
            print(f"  {row['inhibitor']}: de_novo {row['de_novo_med']:.1f} vs secondary {row['secondary_med']:.1f}, "
                  f"P={row['p']:.3g}, FDR q={row['q_fdr']:.3g}")
    else:
        print(f"  {name}: 未找到")
print("完成。")
