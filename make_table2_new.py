#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成新的 Table2.xlsx（Cox 表，含 TP53 VAF 判据 + 补充分析）"""
import pandas as pd

OUT = "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/投稿包"

# 主 Cox 表（sequentially adjusted models）
table = pd.DataFrame({
    "Covariate": ["Secondary AML", "Therapy-related AML", "Chemotherapy", "ELN adverse", "ELN favorable"],
    "Model 1 (ontogeny) HR (95% CI)": ["2.24 (1.76–2.86)", "1.73 (1.38–2.18)", "—", "—", "—"],
    "Model 2 (+age/sex) HR (95% CI)": ["1.59 (1.24–2.04)", "1.28 (1.02–1.62)", "—", "—", "—"],
    "Model 3 (+chemo) HR (95% CI)": ["1.55 (1.20–1.99)", "1.28 (1.02–1.62)", "0.48 (0.34–0.66)", "—", "—"],
    "Model 4 (+ELN2022) HR (95% CI)": ["1.40 (1.09–1.81)", "1.26 (0.99–1.59)", "0.47 (0.34–0.65)", "1.03 (0.83–1.27)", "0.45 (0.31–0.65)"],
    "Model 4 P": ["0.0089", "0.057", "4.5e-6", "0.82", "1.4e-5"],
})
with pd.ExcelWriter(f"{OUT}/Table2.xlsx", engine="openpyxl") as w:
    table.to_excel(w, index=False, sheet_name="Table 2")
    # 补充分析 sheet
    supp = pd.DataFrame({
        "Supplementary analyses": [
            "TP53 re-judged by VAF≥10% (vs any mutation): secondary HR unchanged (1.40); 1 of 75 TP53-mutant patients reclassified",
            "Transplant sensitivity (excluding 242 transplanted): secondary HR 1.36 (1.03–1.80), P=0.033",
            "E-value for HR 1.40 = 2.15 (unmeasured confounder HR required to nullify)",
            "RMST at 24 months: secondary 10.3 vs de novo 16.6 months (loss 6.3 months)",
            "Schoenfeld: secondary-AML coefficient P=0.33 (PH satisfied); age/chemo/ELN-adverse violated PH",
        ],
    })
    supp.to_excel(w, index=False, sheet_name="Supplementary analyses")
print("Table2.xlsx 完成")
