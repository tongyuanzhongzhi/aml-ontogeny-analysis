#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""重绘 Supplementary Fig S3（转录组签名），内部标题改为中性表述，去掉旧编号与阳性口吻"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

BLUE, ORANGE, GREEN = "#0072B2", "#D55E00", "#009E73"
RES = "/Users/daizhensheng/Desktop/学术文章发表/已投文章/继发性AML三分类预后_AMJ/其他文件/结果"

sc = pd.read_csv(f"{RES}/ontogeny_score_v2.csv")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
box_data = [sc[sc["is_secondary"] == 1]["score"],
            sc[sc["is_therapy"] == 1]["score"],
            sc[(sc["is_secondary"] == 0) & (sc["is_therapy"] == 0)]["score"]]
bp = ax1.boxplot(box_data, tick_labels=["secondary", "therapy-\nrelated", "de novo"], patch_artist=True)
for patch, color in zip(bp["boxes"], [ORANGE, GREEN, BLUE]):
    patch.set_facecolor(color); patch.set_alpha(0.7)
ax1.set_ylabel("Secondary-like signature score")
ax1.set_title("Signature score by ontogeny group", fontsize=10)
ax1.grid(axis="y", alpha=0.3)

dn = sc[sc["is_secondary"] == 0].copy()
dn["hi"] = (dn["score"] >= dn["score"].median()).astype(int)
kmf = KaplanMeierFitter()
for label, sub in [("de novo, high signature", dn[dn["hi"] == 1]),
                   ("de novo, low signature", dn[dn["hi"] == 0])]:
    kmf.fit(sub["os_months"], sub["event"], label=label)
    color = ORANGE if "high" in label else BLUE
    kmf.plot_survival_function(ax=ax2, color=color, ci_show=False)
ax2.set_xlabel("Months"); ax2.set_ylabel("Overall survival probability")
ax2.set_title("Signature score within de novo AML", fontsize=10)
ax2.set_ylim(0, 1); ax2.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(f"{RES}/Fig5_转录组签名.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Fig S3 重绘完成（中性标题，无旧编号）")
