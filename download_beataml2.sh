#!/bin/bash
# 下载 BeatAML2 官方 harmonized 数据（走代理）
set -e
cd "/Users/daizhensheng/Desktop/继发性AML三分类预后_BeatAML/原始数据"
mkdir -p BeatAML2_harmonized
cd BeatAML2_harmonized
PROXY="http://127.0.0.1:7897"
BASE="https://github.com/biodev/beataml2.0_data/raw/main"

echo "=== 下载 BeatAML2 官方 harmonized 数据 ==="
curl -sS -L -x "$PROXY" "$BASE/beataml_wv1to4_clinical.xlsx" -o beataml_wv1to4_clinical.xlsx -w "clinical: HTTP %{http_code}, %{size_download} bytes, %{time_total}s\n"
curl -sS -L -x "$PROXY" "$BASE/beataml_wes_wv1to4_mutations_dbgap.txt" -o beataml_wes_wv1to4_mutations_dbgap.txt -w "mutations: HTTP %{http_code}, %{size_download} bytes, %{time_total}s\n"
curl -sS -L -x "$PROXY" "$BASE/beataml_probit_curve_fits_v4_dbgap.txt" -o beataml_probit_curve_fits_v4_dbgap.txt -w "drug_AUC: HTTP %{http_code}, %{size_download} bytes, %{time_total}s\n"
curl -sS -L -x "$PROXY" "$BASE/beataml_drug_families.xlsx" -o beataml_drug_families.xlsx -w "drug_families: HTTP %{http_code}, %{size_download} bytes, %{time_total}s\n"
curl -sS -L -x "$PROXY" "$BASE/beataml_waves1to4_sample_mapping.xlsx" -o beataml_waves1to4_sample_mapping.xlsx -w "sample_mapping: HTTP %{http_code}, %{size_download} bytes, %{time_total}s\n"
curl -sS -L -x "$PROXY" "$BASE/beataml_waves1to4_norm_exp_dbgap.txt" -o beataml_waves1to4_norm_exp_dbgap.txt -w "expression: HTTP %{http_code}, %{size_download} bytes, %{time_total}s\n"

echo ""
echo "=== 下载完成，文件清单 ==="
ls -la
