#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重新生成 Supplementary.pdf：Table S1（补 n=759 + 澄清 ANOVA）+ Fig S1（流向图改 validation/ANOVA/n）+ Fig S2（森林图）+ Fig S3（原 Fig5 签名）"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.image as mpimg
import textwrap

OUT = "/Users/daizhensheng/Desktop/学术文章发表/已投文章/继发性AML三分类预后_AMJ/投稿包_CMJ"
RES = "/Users/daizhensheng/Desktop/学术文章发表/已投文章/继发性AML三分类预后_AMJ/其他文件/结果"

# ---------- 1. 生成修正版患者流向图（临时保存） ----------
fig, ax = plt.subplots(figsize=(8, 9))
ax.set_xlim(0, 10); ax.set_ylim(0, 14); ax.axis("off")
def box(x, y, w, h, text, color="#EAF2F8", fontsize=10, bold=False):
    b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02",
                       linewidth=1.2, edgecolor="#2C3E50", facecolor=color)
    ax.add_patch(b)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fontsize, fontweight="bold" if bold else "normal")
def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=14, color="#2C3E50", linewidth=1.2))

box(2.5, 12.5, 5, 1.2, "Beat AML cohort (aml_ohsu_2022)\n805 patients, 942 specimens", "#D5E8D4", bold=True)
box(2.5, 10.6, 5, 1.2, "804 patients with OS data", "#EAF2F8")
box(2.5, 8.8, 5, 1.2, "784 patients (exclude 20 APL)\n→ Cox (complete-case n=759)", "#EAF2F8")
arrow(5, 12.5, 5, 11.8); arrow(5, 10.6, 5, 10.0)
box(0.3, 6.8, 4.4, 1.4, "Ontogeny classification\nby clinical history\nde novo 539 / tAML 143 / sAML 123", "#FDEBD0")
box(5.3, 6.8, 4.4, 1.4, "Molecular annotation\nbaseline specimens\n(580 Initial Diagnosis)", "#EAF2F8")
arrow(5, 8.8, 2.5, 8.2); arrow(5, 8.8, 7.5, 8.2)
box(0.3, 4.6, 4.4, 1.4, "2×2 analysis (n=661)\nde novo + secondary\n380/158/39/84", "#FDEBD0")
box(5.3, 4.6, 4.4, 1.4, "Drug sensitivity\n144 inhibitors (four-factor\ntype II ANOVA)", "#EAF2F8")
arrow(2.5, 6.8, 2.5, 6.0); arrow(7.5, 6.8, 7.5, 6.0)
box(2.5, 2.4, 5, 1.4, "Transcriptomic analysis\n613 patients (RNA-seq mapped)\n10,476 DEGs + signature", "#E8DAEF")
arrow(7.5, 4.6, 5, 3.8)
box(2.5, 0.4, 5, 1.2, "TCGA LAML (n=200)\nexternal reference comparison", "#F5B7B1")
arrow(5, 2.4, 5, 1.6)
plt.tight_layout()
flowchart_path = f"{RES}/FigS1_流向图_v2.png"
fig.savefig(flowchart_path, dpi=300, bbox_inches="tight")
plt.close(fig)
print("患者流向图 v2 完成")

# ---------- 2. 合并 Supplementary.pdf（4 页） ----------
with PdfPages(f"{OUT}/Supplementary.pdf") as pdf:
    # 第 1 页：Supplementary Table 1
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.5, 0.95, "Supplementary Table 1. ELN 2022 re-derivation rules",
             ha="center", fontsize=13, fontweight="bold")
    text = (
        "Cytogenetic risk (from karyotype and FISH):\n"
        "  APL: t(15;17) on karyotype\n"
        "  Favorable: inv(16)/t(16;16); t(8;21)\n"
        "  Adverse: complex karyotype (≥3 abnormalities); −5/del(5q); −7; −17/17p; "
        "inv(3)/t(3;3)/3q26; t(6;9); t(9;22); 11q23/KMT2A; 12p\n"
        "  Normal: 46,XX / 46,XY without other abnormalities\n"
        "  Other: all remaining karyotypes\n\n"
        "FISH was considered abnormal only when reporting DELETION/TRISOMY/REARRANGEMENT/"
        "MONOSOMY/LOSS OF/GAIN OF and not 'not detected'/'normal' (to avoid counting routine "
        "PML/MLL probes as true abnormalities).\n\n"
        "Mutation-based risk:\n"
        "  Adverse: TP53 mutation (variant allele fraction ≥10%); any of the nine myelodysplasia-related (MR) genes "
        "ASXL1, BCOR, EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1, ZRSR2\n"
        "  Favorable: in-frame mutation in the CEBPA bZIP domain (protein position ≥278, "
        "excluding frameshift/nonsense); or NPM1 mutation without FLT3-ITD\n\n"
        "Priority order (highest wins): APL > Adverse > Favorable > Intermediate\n\n"
        "Patient flow:\n"
        "  805 patients (942 specimens) in Beat AML (aml_ohsu_2022)\n"
        "  → 804 with overall survival data\n"
        "  → 784 after excluding 20 APL cases → stepwise Cox models (complete-case n = 759)\n"
        "  → molecular/cytogenetic annotation: baseline (Initial Diagnosis) specimens prioritized "
        "(580 initial-diagnosis specimens; remainder annotated from earliest specimen)\n"
        "  → drug sensitivity: 144 inhibitors → four-factor type II ANOVA (ontogeny, MR mutations, "
        "TP53, complex karyotype); complementary two-way ANOVA (ontogeny × MR) in de novo and "
        "secondary cases with ≥3 samples per cell\n"
        "  → 2×2 analysis: 661 patients (de novo + secondary)\n"
        "  → transcriptomic signature: 613 patients with RNA-seq mapped to ontogeny\n"
    )
    wrapped_lines = []
    for line in text.split("\n"):
        wrapped_lines.extend(textwrap.wrap(line, width=88) if len(line) > 88 else [line])
    fig.text(0.08, 0.88, "\n".join(wrapped_lines), ha="left", va="top", fontsize=9.5)
    pdf.savefig(fig); plt.close(fig)

    # 第 2 页：Fig S1 患者流向图（修正版）
    fig = plt.figure(figsize=(8.5, 11))
    img = mpimg.imread(flowchart_path)
    plt.imshow(img); plt.axis("off")
    plt.title("Supplementary Fig S1. Patient flow diagram", fontsize=12, fontweight="bold")
    pdf.savefig(fig); plt.close(fig)

    # 第 3 页：Fig S2 Cox 森林图（n=759 正确值）
    fig = plt.figure(figsize=(8.5, 11))
    img = mpimg.imread(f"{RES}/FigS1_Cox森林图.png")
    plt.imshow(img); plt.axis("off")
    plt.title("Supplementary Fig S2. Cox forest plot (n=759)", fontsize=12, fontweight="bold")
    pdf.savefig(fig); plt.close(fig)

    # 第 4 页：Fig S3 转录组签名（原正文 Fig 5，年龄混杂阴性结果）
    fig = plt.figure(figsize=(8.5, 11))
    img = mpimg.imread(f"{OUT}/Supplementary_FigS3_transcriptomic.jpg")
    plt.imshow(img); plt.axis("off")
    plt.title("Supplementary Fig S3. Transcriptomic secondary-like signature "
              "(age-confounded negative result)", fontsize=12, fontweight="bold")
    pdf.savefig(fig); plt.close(fig)

print("Supplementary.pdf 完成（4 页：Table S1 + Fig S1/S2/S3）")
