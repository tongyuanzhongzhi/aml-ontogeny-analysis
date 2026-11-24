#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重算 Table 1（基线）和 Table 4（突变谱），确保与正文定稿值一致"""
import json, re
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact

DATA = "/Users/daizhensheng/Desktop/学术文章发表/已投文章/继发性AML三分类预后_AMJ/其他文件/原始数据"

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

# sample -> genes
sample_genes = defaultdict(set)
for m in mutations:
    sym = id2symbol.get(m["entrezGeneId"], "")
    if sym:
        sample_genes[m["sampleId"]].add(sym)

# TP53 VAF >= 10%
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

groups = defaultdict(list)
for pid, d in patient.items():
    groups[classify(d)].append(pid)

for g, pids in groups.items():
    print(f"{g}: n={len(pids)}")

# Table 4 突变谱（患者级，TP53 VAF≥10%）
genes = ["SRSF2", "SF3B1", "U2AF1", "ZRSR2", "ASXL1", "EZH2", "STAG2", "RUNX1", "BCOR",
         "TP53", "NPM1", "FLT3"]
print("\n=== Table 4 突变谱（患者级；TP53 用 VAF≥10%）===")
print("Gene\tde novo\ttherapy\tsecondary\tP(sec vs de novo)")
for g in genes:
    freqs = {}
    for grp in ["de novo", "therapy-related", "secondary"]:
        pids = groups[grp]
        n_mut = sum(1 for pid in pids if (patient_tp53(pid) if g == "TP53" else patient_mut(pid, g)))
        freqs[grp] = (n_mut, len(pids))
    dn_mut, dn_n = freqs["de novo"]
    sec_mut, sec_n = freqs["secondary"]
    # 注意：正文 TP53 是 tAML vs de novo，其他基因是 sec vs de novo
    tbl = [[sec_mut, sec_n - sec_mut], [dn_mut, dn_n - dn_mut]]
    odds, p = fisher_exact(tbl)
    print(f"{g}\t{freqs['de novo'][0]}/{freqs['de novo'][1]} ({freqs['de novo'][0]/freqs['de novo'][1]*100:.1f}%)\t"
          f"{freqs['therapy-related'][0]}/{freqs['therapy-related'][1]} ({freqs['therapy-related'][0]/freqs['therapy-related'][1]*100:.1f}%)\t"
          f"{sec_mut}/{sec_n} ({sec_mut/sec_n*100:.1f}%)\t{p:.2e}")

# TP53 tAML vs de novo 的 Fisher P（正文口径）
t_mut, t_n = freqs['therapy-related']
dn_mut2, dn_n2 = freqs['de novo']
tbl2 = [[t_mut, t_n - t_mut], [dn_mut2, dn_n2 - dn_mut2]]
odds2, p2 = fisher_exact(tbl2)
print(f"\nTP53 tAML vs de novo: {t_mut}/{t_n} ({t_mut/t_n*100:.1f}%) vs {dn_mut2}/{dn_n2} ({dn_mut2/dn_n2*100:.1f}%), Fisher P={p2:.2e}")
