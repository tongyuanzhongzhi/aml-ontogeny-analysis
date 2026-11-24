#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重新生成 Fig 4：红系基因柱状图 + 显著性星号（Mann-Whitney U 实测 P 值）"""
import json
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

DATA = "/Users/daizhensheng/Desktop/学术文章发表/已投文章/继发性AML三分类预后_AMJ/其他文件/原始数据"
OUT = "/Users/daizhensheng/Desktop/学术文章发表/已投文章/继发性AML三分类预后_AMJ/投稿包_CMJ"
BLUE, ORANGE = "#0072B2", "#D55E00"

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

exp = pd.read_csv(f"{DATA}/BeatAML2_harmonized/beataml_waves1to4_norm_exp_dbgap.txt", sep="\t", index_col=0, low_memory=False)
meta_cols = ["display_label", "description", "biotype"]
gene_labels = exp["display_label"].copy()
sample_cols = [c for c in exp.columns if c not in meta_cols]
sample_group = {}
for c in sample_cols:
    subj = rna2subject.get(c)
    if subj is not None:
        sample_group[c] = subj2group.get(subj, None)
mapped = {c: g for c, g in sample_group.items() if g is not None}
sym2idx = {}
for i, s in enumerate(gene_labels.values):
    if str(s) not in sym2idx:
        sym2idx[str(s)] = i
exp_mat = exp[sample_cols].apply(pd.to_numeric, errors="coerce")

groups = ["de novo", "therapy-related", "secondary"]
ery_genes = ["TAL1", "GATA1", "GP9", "ITGA2B", "EPOR"]
ery_dn, ery_sec, pvals = [], [], []
for g in ery_genes:
    if g in sym2idx:
        i = sym2idx[g]
        dn = exp_mat.iloc[i][[c for c, gg in mapped.items() if gg == "de novo"]].values
        sec = exp_mat.iloc[i][[c for c, gg in mapped.items() if gg == "secondary"]].values
        dn = dn[~np.isnan(dn)]; sec = sec[~np.isnan(sec)]
        ery_dn.append(np.median(dn)); ery_sec.append(np.median(sec))
        _, p = mannwhitneyu(dn, sec, alternative="two-sided")
        pvals.append(p)

def stars(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "ns"

# BCL2L1/BCL2, MCL1/BCL2 比值
ratio_data = {}
for label, num, den in [("BCL2L1/BCL2", "BCL2L1", "BCL2"), ("MCL1/BCL2", "MCL1", "BCL2")]:
    if num in sym2idx and den in sym2idx:
        i_num, i_den = sym2idx[num], sym2idx[den]
        ratio_data[label] = {}
        for g in groups:
            cols = [c for c, gg in mapped.items() if gg == g]
            vals = [exp_mat.iloc[i_num][c]/exp_mat.iloc[i_den][c] for c in cols
                    if exp_mat.iloc[i_den][c] > 0 and not np.isnan(exp_mat.iloc[i_num][c])]
            ratio_data[label][g] = np.median(vals)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
x = np.arange(len(ery_genes)); w = 0.35
ax1.bar(x - w/2, ery_dn, w, label="de novo", color=BLUE, alpha=0.85)
ax1.bar(x + w/2, ery_sec, w, label="secondary", color=ORANGE, alpha=0.85)
# 星号
for xi, (dn, sec, p) in enumerate(zip(ery_dn, ery_sec, pvals)):
    top = max(dn, sec)
    ax1.text(xi, top + 0.03 * top, stars(p), ha="center", va="bottom", fontsize=10, fontweight="bold")
ax1.set_xticks(x); ax1.set_xticklabels(ery_genes)
ax1.set_ylabel("Median expression")
ax1.set_title("Erythroid/megakaryocytic skewing", fontsize=10)
ax1.legend(); ax1.grid(axis="y", alpha=0.3)
ax1.set_ylim(0, max(max(ery_dn), max(ery_sec)) * 1.15)

rlabels = list(ratio_data.keys())
dn_r = [ratio_data[l]["de novo"] for l in rlabels]
sec_r = [ratio_data[l]["secondary"] for l in rlabels]
x2 = np.arange(len(rlabels))
ax2.bar(x2 - w/2, dn_r, w, label="de novo", color=BLUE, alpha=0.85)
ax2.bar(x2 + w/2, sec_r, w, label="secondary", color=ORANGE, alpha=0.85)
ax2.set_xticks(x2); ax2.set_xticklabels(rlabels)
ax2.set_ylabel("Ratio (median)")
ax2.set_title("Anti-apoptotic switch (BCL2 → BCL-XL/MCL1)", fontsize=10)
ax2.legend(); ax2.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig4_BCL2_switch.jpg", dpi=300, bbox_inches="tight")
plt.close(fig)

print("Fig4 改进版完成（星号）")
for g, p in zip(ery_genes, pvals):
    print(f"  {g}: P = {p:.2e} → {stars(p)}")
