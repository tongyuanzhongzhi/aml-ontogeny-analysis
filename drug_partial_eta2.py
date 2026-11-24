#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P0-1 修正：药敏拆解加入 TP53 + 复杂核型，四因素方差分解比较 partial eta²
检验 ontogeny 是否仍主导（而非显著个数 58 vs 1）
"""
import json, re
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.stats import rankdata
import statsmodels.api as sm
from statsmodels.formula.api import ols

DATA = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/其他文件/原始数据"
OUT = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/其他文件/结果"

with open(f"{DATA}/gene_symbol_map.json") as f:
    id2symbol = {int(k): v for k, v in json.load(f).items()}
with open(f"{DATA}/mutations.json") as f:
    mutations = json.load(f)
sample_genes = defaultdict(set)
for m in mutations:
    sym = id2symbol.get(m["entrezGeneId"], "")
    if sym:
        sample_genes[m["sampleId"]].add(sym)

def tp53_vaf_positive(sid):
    for m in mutations:
        if m["sampleId"] == sid and id2symbol.get(m["entrezGeneId"]) == "TP53":
            alt = m.get("tumorAltCount") or 0; ref = m.get("tumorRefCount") or 0
            if alt + ref > 0 and alt / (alt + ref) >= 0.10:
                return True
    return False

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

ELN2022_SEC = ["ASXL1", "BCOR", "EZH2", "RUNX1", "SF3B1", "SRSF2", "STAG2", "U2AF1", "ZRSR2"]

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

def patient_tp53(pid):
    return any(tp53_vaf_positive(s) for s in baseline_samples(pid))

def patient_complex(pid):
    """复杂核型：≥3 染色体异常"""
    for s in baseline_samples(pid):
        k = str(sample_info[s].get("KARYOTYPE", "")).upper()
        marks = (re.findall(r'DEL\(', k) + re.findall(r'ADD\(', k) + re.findall(r'INV\(', k)
                 + re.findall(r'T\(', k) + re.findall(r'DER\(', k) + re.findall(r'\bMAR\b', k)
                 + re.findall(r'[+-]\d+', k))
        if len(marks) >= 3:
            return True
    return False

subj2info = {}
for pid, d in patient.items():
    g = classify(d)
    if g == "therapy-related":
        continue
    subj = pid.split("_")[-1]
    if subj.isdigit():
        subj2info[int(subj)] = {
            "ontogeny": 1 if g == "secondary" else 0,
            "mr_mut": int(any(patient_mut(pid, gene) for gene in ELN2022_SEC)),
            "tp53": int(patient_tp53(pid)),
            "complex": int(patient_complex(pid)),
        }

# 药敏 AUC
drug = pd.read_csv(f"{DATA}/BeatAML2_harmonized/beataml_probit_curve_fits_v4_dbgap.txt", sep="\t")
drug["auc"] = pd.to_numeric(drug["auc"], errors="coerce")
drug = drug.dropna(subset=["auc"])
for k in ["ontogeny", "mr_mut", "tp53", "complex"]:
    drug[k] = drug["dbgap_subject_id"].map(lambda s: subj2info.get(int(s), {}).get(k))
d = drug[drug["ontogeny"].notna() & drug["mr_mut"].notna()].copy()

# 对每个药物，四因素 ANOVA（秩转换），计算 partial eta²
rows = []
for inhibitor, sub in d.groupby("inhibitor"):
    cnt = sub.groupby(["ontogeny", "mr_mut"]).size()
    if cnt.min() < 3 or len(cnt) < 4:
        continue
    sub = sub.copy()
    sub["auc_rank"] = rankdata(sub["auc"])
    try:
        m = ols("auc_rank ~ C(ontogeny) + C(mr_mut) + C(tp53) + C(complex)", data=sub).fit()
        anova = sm.stats.anova_lm(m, typ=2)
        ss_total = anova["sum_sq"].sum()
        ss_error = anova.loc["Residual", "sum_sq"]
        row = {"inhibitor": inhibitor}
        for factor in ["C(ontogeny)", "C(mr_mut)", "C(tp53)", "C(complex)"]:
            if factor in anova.index:
                ss = anova.loc[factor, "sum_sq"]
                p = anova.loc[factor, "PR(>F)"]
                eta2 = ss / (ss + ss_error)  # partial eta²
                row[factor + "_eta2"] = eta2
                row[factor + "_p"] = p
        rows.append(row)
    except Exception:
        continue

res = pd.DataFrame(rows)
print(f"成功分析药物数: {len(res)}")

# 比较 partial eta² 均值
print("\n=== partial eta² 均值（秩转换，四因素）===")
factors = {"ontogeny": "C(ontogeny)", "MR基因": "C(mr_mut)", "TP53": "C(tp53)", "复杂核型": "C(complex)"}
for name, col in factors.items():
    eta_col = col + "_eta2"
    if eta_col in res.columns:
        print(f"  {name}: 中位 partial eta² = {res[eta_col].median():.4f}, 均值 = {res[eta_col].mean():.4f}")

# ontogeny 主导的药物数（partial eta² 最大）
print("\n=== 各因素 partial eta² 最大（主导）的药物数 ===")
dom_cols = [c + "_eta2" for c in factors.values() if c + "_eta2" in res.columns]
if dom_cols:
    dom = res[dom_cols].idxmax(axis=1)
    from collections import Counter
    print("  " + str(Counter(dom)))

# ontogeny 仍是显著（p<0.05）且 eta² 最大的药物
onto_eta = "C(ontogeny)_eta2"; onto_p = "C(ontogeny)_p"
if onto_eta in res.columns:
    n_onto_sig = (res[onto_p] < 0.05).sum()
    n_onto_dom = (res[dom_cols].idxmax(axis=1) == onto_eta).sum()
    print(f"\nontogeny 显著(p<0.05)的药物: {n_onto_sig}")
    print(f"ontogeny partial eta² 最大的药物: {n_onto_dom}")

# 关键药
print("\n=== 关键药 partial eta² ===")
for name in ["Venetoclax", "Cytarabine", "Sorafenib", "Midostaurin"]:
    r = res[res["inhibitor"].str.contains(name, case=False, na=False)]
    if len(r):
        row = r.iloc[0]
        print(f"  {name}: ontogeny eta²={row.get('C(ontogeny)_eta2', 0):.3f}, "
              f"MR eta²={row.get('C(mr_mut)_eta2', 0):.3f}, "
              f"TP53 eta²={row.get('C(tp53)_eta2', 0):.3f}, "
              f"复杂核型 eta²={row.get('C(complex)_eta2', 0):.3f}")

res.to_csv(f"{OUT}/药敏_四因素partial_eta2.csv", index=False, encoding="utf-8-sig")
print("\n完成。")
