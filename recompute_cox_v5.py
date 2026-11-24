#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""审稿意见响应：重跑 Cox + stratified Cox + 仅初诊敏感性 + 共线性诊断 + 森林图数据"""
import json, re
from collections import defaultdict
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

DATA = "/Users/daizhensheng/Desktop/学术文章发表/已投文章/继发性AML三分类预后_AMJ/其他文件/原始数据"
OUT = "/Users/daizhensheng/Desktop/学术文章发表/已投文章/继发性AML三分类预后_AMJ/其他文件/结果"

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

def has_initial_specimen(pid):
    return any(sample_info[s].get("DISEASE_STAGE_AT_SPECIMEN_COLLECTION", "") == "Initial Diagnosis"
               for s in pid2samples.get(pid, []))

def patient_mut(pid, gene):
    return any(gene in sample_genes.get(s, set()) for s in baseline_samples(pid))
def patient_tp53(pid):
    return any(tp53_vaf_positive(s) for s in baseline_samples(pid))
def patient_age(pid):
    ages = [float(sample_info[s]["AGE_AT_DIAGNOSIS"]) for s in baseline_samples(pid)
            if "AGE_AT_DIAGNOSIS" in sample_info[s]]
    return np.median(ages) if ages else np.nan
def cebpa_bzip(sid):
    for m in mutations:
        if m["sampleId"] == sid and id2symbol.get(m["entrezGeneId"]) == "CEBPA":
            pos = m.get("proteinPosStart"); change = str(m.get("proteinChange", ""))
            if pos is not None and pos >= 278 and "fs" not in change and "*" not in change:
                return True
    return False
def patient_cebpa_bzip(pid):
    return any(cebpa_bzip(s) for s in baseline_samples(pid))

def cyto_risk(karyo, fish):
    k = str(karyo).upper(); f = str(fish).upper()
    if "T(15;17)" in k: return "apl"
    if any(x in k for x in ["INV(16)", "T(16;16)", "T(8;21)"]): return "favorable"
    marks = (re.findall(r'DEL\(', k) + re.findall(r'ADD\(', k) + re.findall(r'INV\(', k)
             + re.findall(r'T\(', k) + re.findall(r'DER\(', k) + re.findall(r'\bMAR\b', k)
             + re.findall(r'[+-]\d+', k))
    is_complex = len(marks) >= 3
    is_adverse_specific = any(x in k for x in ["INV(3", "T(3;3", "3Q26", "ADD(3", "T(6;9", "T(9;22", "11Q23",
                                               "-5", "DEL(5", "-7", "DEL(7", "-17", "17P", "I(17", "12P"])
    fish_abnormal = (any(x in f for x in ["DELETION", "TRISOMY", "REARRANGEMENT", "MONOSOMY", "LOSS OF", "GAIN OF"])
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
    has_tp53 = patient_tp53(pid)
    has_npm1 = patient_mut(pid, "NPM1")
    has_flt3 = patient_mut(pid, "FLT3")
    if cyto == "apl": return "APL"
    if has_tp53 or cyto == "adverse" or has_sec: return "Adverse"
    if patient_cebpa_bzip(pid) or (has_npm1 and not has_flt3) or cyto == "favorable": return "Favorable"
    return "Intermediate"

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
        "eln2022": eln2022(pid),
        "eln_adv": 1 if eln2022(pid) == "Adverse" else 0,
        "eln_fav": 1 if eln2022(pid) == "Favorable" else 0,
        "has_sec_mut": int(any(patient_mut(pid, g) for g in ELN2022_SEC)),
        "has_tp53": int(patient_tp53(pid)),
        "has_initial": has_initial_specimen(pid),
    })

df = pd.DataFrame(rows).dropna(subset=["os_months"])
df = df[df["os_months"] >= 0]
df["is_therapy"] = (df["group"] == "therapy-related").astype(int)
df["is_secondary"] = (df["group"] == "secondary").astype(int)
df["age_group"] = pd.cut(df["age"], bins=[0, 60, 70, 200], labels=["<60", "60-69", "≥70"])

df_main = df[df["eln2022"] != "APL"].copy()
print(f"=== 样本量 ===\n总(有OS): {len(df)}, 排除APL: {len(df_main)}, 仅初诊标本: {df_main['has_initial'].sum()}")

def run_cox(d, cols, label, strata=None):
    m = d[["os_months", "event"] + cols + (strata or [])].dropna()
    cph = CoxPHFitter()
    if strata:
        cph.fit(m, duration_col="os_months", event_col="event", strata=strata)
    else:
        cph.fit(m, duration_col="os_months", event_col="event")
    print(f"\n[{label}] n={len(m)}")
    for c in cols:
        if c in cph.params_:
            hr = np.exp(cph.params_[c])
            lo = np.exp(cph.confidence_intervals_.loc[c, "95% lower-bound"])
            hi = np.exp(cph.confidence_intervals_.loc[c, "95% upper-bound"])
            p = cph.summary.loc[c, "p"]
            print(f"  {c:12s} HR={hr:.3f} (95% CI {lo:.2f}-{hi:.2f}) p={p:.4g}")
    return cph

print("\n========== 主 Cox（Model 4 完整 ELN）==========")
run_cox(df_main, ["is_secondary", "is_therapy", "age", "sex_male", "chemo", "eln_adv", "eln_fav"], "Model 4")
print("\n========== 主 Cox（Model 5 TP53+MR）==========")
run_cox(df_main, ["is_secondary", "is_therapy", "age", "sex_male", "chemo", "has_sec_mut", "has_tp53"], "Model 5")

print("\n========== Stratified Cox（对违背PH的变量分层）==========")
run_cox(df_main, ["is_secondary", "is_therapy", "age", "sex_male"], "Stratified by ELN+chemo+age_group",
        strata=["eln2022", "chemo", "age_group"])

print("\n========== 仅初诊标本敏感性（n=580）==========")
df_init = df_main[df_main["has_initial"]].copy()
run_cox(df_init, ["is_secondary", "is_therapy", "age", "sex_male", "chemo", "eln_adv", "eln_fav"], "Model 4 (仅初诊)")
run_cox(df_init, ["is_secondary", "is_therapy", "age", "sex_male", "chemo", "has_sec_mut", "has_tp53"], "Model 5 (仅初诊)")

print("\n========== 共线性诊断：ontogeny × ELN 列联表 ==========")
ct = pd.crosstab(df_main["group"], df_main["eln2022"], normalize="index") * 100
print(ct.round(1))
print("\nELN adverse 在三组中的占比:")
print((df_main.groupby("group")["eln_adv"].mean() * 100).round(1))

# 森林图数据输出（n=784 最新 Model 5）
print("\n========== 森林图数据（Model 5, n=784 最新）==========")
model5 = run_cox(df_main, ["is_secondary", "is_therapy", "chemo", "has_sec_mut", "has_tp53", "age", "sex_male"], "Model 5 森林图")
