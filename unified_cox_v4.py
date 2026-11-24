#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P0-5 修正：TP53 按 VAF≥10% 判不良（ELN 2022 口径），而非任意突变
重跑完整 Cox，检验 HR 1.40 是否上移
"""
import json, re
from collections import defaultdict
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter

DATA = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/其他文件/原始数据"
OUT = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/其他文件/结果"

with open(f"{DATA}/gene_symbol_map.json") as f:
    id2symbol = {int(k): v for k, v in json.load(f).items()}
with open(f"{DATA}/mutations.json") as f:
    mutations = json.load(f)

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

# sample -> 突变基因（全 MAF）
sample_genes = defaultdict(set)
for m in mutations:
    sym = id2symbol.get(m["entrezGeneId"], "")
    if sym:
        sample_genes[m["sampleId"]].add(sym)

# TP53 突变 VAF 判据：VAF >= 10% 才判为 TP53 阳性（ELN 2022 口径）
def tp53_vaf_positive(sid):
    for m in mutations:
        if m["sampleId"] == sid and id2symbol.get(m["entrezGeneId"]) == "TP53":
            alt = m.get("tumorAltCount") or 0
            ref = m.get("tumorRefCount") or 0
            if alt + ref > 0 and alt / (alt + ref) >= 0.10:
                return True
    return False

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

def patient_age(pid):
    ages = [float(sample_info[s]["AGE_AT_DIAGNOSIS"]) for s in baseline_samples(pid)
            if "AGE_AT_DIAGNOSIS" in sample_info[s]]
    return np.median(ages) if ages else np.nan

# CEBPA bZIP
def cebpa_bzip(sid):
    for m in mutations:
        if m["sampleId"] == sid and id2symbol.get(m["entrezGeneId"]) == "CEBPA":
            pos = m.get("proteinPosStart"); change = str(m.get("proteinChange", ""))
            if pos is not None and pos >= 278 and "fs" not in change and "*" not in change:
                return True
    return False
def patient_cebpa_bzip(pid):
    return any(cebpa_bzip(s) for s in baseline_samples(pid))

# 细胞遗传学
def cyto_risk(karyo, fish):
    k = str(karyo).upper(); f = str(fish).upper()
    if "T(15;17)" in k: return "apl"
    if any(x in k for x in ["INV(16)", "T(16;16)", "T(8;21)"]): return "favorable"
    marks = (re.findall(r'DEL\(', k) + re.findall(r'ADD\(', k) + re.findall(r'INV\(', k)
             + re.findall(r'T\(', k) + re.findall(r'DER\(', k) + re.findall(r'\bMAR\b', k)
             + re.findall(r'[+-]\d+', k))
    is_complex = len(marks) >= 3
    is_adverse_specific = any(x in k for x in [
        "INV(3", "T(3;3", "3Q26", "ADD(3", "T(6;9", "T(9;22", "11Q23",
        "-5", "DEL(5", "-7", "DEL(7", "-17", "17P", "I(17", "12P"])
    fish_abnormal = (any(x in f for x in ["DELETION", "TRISOMY", "REARRANGEMENT",
                                          "MONOSOMY", "LOSS OF", "GAIN OF"])
                     and "NOT DETECTED" not in f and "NORMAL" not in f)
    if is_complex or is_adverse_specific or fish_abnormal: return "adverse"
    if k.startswith(("46,XX", "46,XY")): return "normal"
    return "other"

def patient_cyto(pid):
    risks = [cyto_risk(sample_info[s].get("KARYOTYPE", ""), sample_info[s].get("OTHER_CYTOGENETICS", ""))
             for s in baseline_samples(pid)]
    for prio in ["apl", "adverse", "favorable", "normal", "other"]:
        if prio in risks: return prio
    return "other"

def eln2022(pid):
    cyto = patient_cyto(pid)
    has_sec = any(patient_mut(pid, g) for g in ELN2022_SEC)
    has_tp53 = patient_tp53(pid)  # VAF >= 10%
    has_npm1 = patient_mut(pid, "NPM1")
    has_flt3 = patient_mut(pid, "FLT3")
    cebpa_bzip_mut = patient_cebpa_bzip(pid)
    if cyto == "apl": return "APL"
    if has_tp53 or cyto == "adverse" or has_sec: return "Adverse"
    if cebpa_bzip_mut or (has_npm1 and not has_flt3) or cyto == "favorable": return "Favorable"
    return "Intermediate"

# 构建数据
rows = []
for pid, d in patient.items():
    try: osm = float(d.get("OS_MONTHS"))
    except: osm = np.nan
    rows.append({
        "pid": pid, "group": classify(d), "os_months": osm,
        "event": 1 if str(d.get("OS_STATUS", "")).startswith("1:") else 0,
        "age": patient_age(pid),
        "sex_male": 1 if str(d.get("SEX", "")).strip() == "Male" else 0,
        "chemo": 1 if yes(d, "CUMULATIVE_CHEMO") else 0,
        "cyto": patient_cyto(pid), "eln2022": eln2022(pid),
        "eln_adv": 1 if eln2022(pid) == "Adverse" else 0,
        "eln_fav": 1 if eln2022(pid) == "Favorable" else 0,
        "has_sec_mut": int(any(patient_mut(pid, g) for g in ELN2022_SEC)),
        "has_tp53": int(patient_tp53(pid)),
    })

df = pd.DataFrame(rows).dropna(subset=["os_months"])
df = df[df["os_months"] >= 0]
df["is_therapy"] = (df["group"] == "therapy-related").astype(int)
df["is_secondary"] = (df["group"] == "secondary").astype(int)

# TP53 判据对比
print("=== TP53 判据变化（VAF≥10% vs 任意突变）===")
tp53_any = sum(1 for pid in patient if any("TP53" in sample_genes.get(s, set()) for s in baseline_samples(pid)))
tp53_vaf = df["has_tp53"].sum()
print(f"任意 TP53 突变患者: {tp53_any}")
print(f"VAF≥10% TP53 患者: {tp53_vaf}")
print(f"被降级的亚克隆 TP53 患者: {tp53_any - tp53_vaf}")

df_main = df[df["cyto"] != "apl"].copy()

def run_cox(d, cols, label):
    m = d[["os_months", "event"] + cols].dropna()
    cph = CoxPHFitter()
    cph.fit(m, duration_col="os_months", event_col="event")
    print(f"\n{label}（n={len(m)}）:")
    for c in ["is_secondary", "is_therapy", "chemo", "eln_adv", "eln_fav", "has_sec_mut", "has_tp53", "age", "sex_male"]:
        if c in cols:
            hr = np.exp(cph.params_[c])
            lo = np.exp(cph.confidence_intervals_.loc[c, "95% lower-bound"])
            hi = np.exp(cph.confidence_intervals_.loc[c, "95% upper-bound"])
            p = cph.summary.loc[c, "p"]
            print(f"  {c:12s} HR={hr:.3f} (95% CI {lo:.2f}-{hi:.2f})  p={p:.4g}")

print("\n========== 统一 Cox（TP53 VAF≥10% 判据）==========")
run_cox(df_main, ["is_therapy", "is_secondary"], "模型1：三分类")
run_cox(df_main, ["is_therapy", "is_secondary", "age", "sex_male"], "模型2：+年龄性别")
run_cox(df_main, ["is_therapy", "is_secondary", "age", "sex_male", "chemo"], "模型3：+化疗")
run_cox(df_main, ["is_therapy", "is_secondary", "age", "sex_male", "chemo", "eln_adv", "eln_fav"], "模型4：+完整ELN2022(VAF版)")
run_cox(df_main, ["is_therapy", "is_secondary", "age", "sex_male", "chemo", "has_sec_mut", "has_tp53"], "模型5：+继发突变+TP53")

# 保存
with open(f"{OUT}/Cox_TP53_VAF判据.txt", "w") as f:
    f.write("=== Cox（TP53 按 VAF≥10% 判不良，ELN 2022 口径）===\n")
    f.write(f"任意 TP53: {tp53_any}, VAF≥10%: {tp53_vaf}, 降级 {tp53_any-tp53_vaf}\n\n")
print("\n完成。")
