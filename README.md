# Analysis code: Clinical ontogeny and MR gene mutations in AML

Reproducible analysis scripts for the manuscript:

**"Clinical ontogeny and MR gene mutations are independently associated with ex vivo drug resistance and poor prognosis in acute myeloid leukemia"**

## Data sources

All data are public and re-downloadable; no new patient-level data were collected.

1. **Beat AML cohort** (`aml_ohsu_2022`, Bottomly et al., *Cancer Cell* 2022; 805 patients / 942 specimens)
   - cBioPortal REST API (https://www.cbioportal.org/api/)
   - Files: `clinical_patient.json`, `clinical_sample.json`, `clinical_attributes.json`, `mutations.json`, `entrez_gene_ids.json`, `gene_symbol_map.json`
   - Example: `/api/studies/aml_ohsu_2022/clinical-data?clinicalDataType=PATIENT`

2. **BeatAML2 harmonized data** (drug AUC, RNA-seq TPM, sample mapping)
   - GitHub: https://github.com/biodev/beataml2.0_data (raw)
   - Files: `beataml_wv1to4_clinical.xlsx`, `beataml_wes_wv1to4_mutations_dbgap.txt`, `beataml_probit_curve_fits_v4_dbgap.txt`, `beataml_drug_families.xlsx`, `beataml_waves1to4_sample_mapping.xlsx`, `beataml_waves1to4_norm_exp_dbgap.txt`
   - `download_beataml2.sh` automates these downloads (uses a local proxy; remove `-x "$PROXY"` if not needed)

3. **TCGA LAML cohort** (`laml_tcga_pub`) — used only for de novo mutation-spectrum reference comparison
   - cBioPortal REST API

## Environment

- Python 3.11+
- Packages: `lifelines`, `scipy`, `numpy`, `pandas`, `statsmodels`, `matplotlib`, `openpyxl`, `requests`

```
pip install lifelines scipy numpy pandas statsmodels matplotlib openpyxl
```

## Script inventory

| Script | Purpose |
|---|---|
| `download_beataml2.sh` | Download BeatAML2 harmonized data from GitHub |
| `download_mutations.py` / `download_gene_map.py` | Download cBioPortal mutation + gene-map JSON |
| `unified_cox_v4.py` | Sequentially adjusted Cox models (Models 1–5), TP53 VAF≥10% |
| `recompute_cox_v5.py` | Reproduces all Cox models, stratified Cox (PH), initial-diagnosis sensitivity, collinearity diagnostic |
| `make_table2_new.py` | Table 2 (Cox models) + supplementary analyses |
| `make_table3_4.py` | Table 3 (mutation spectrum) + Table 4 (drug resistance) |
| `recompute_tables_cmj.py` | Re-derive Table 1/4 from raw data (patient-level, TP53 VAF≥10%) |
| `mutation_spectrum.py` | Mutation spectrum by ontogeny |
| `drug_partial_eta2.py` | Four-factor type II ANOVA (partial η²) across 144 inhibitors |
| `bcl2_bclxl.py` | BCL2L1/BCL2 and MCL1/BCL2 ratio analysis |
| `make_fig1_improved.py` | Fig 1 (Kaplan–Meier + number-at-risk + P + median OS + HR) |
| `make_fig4_improved.py` | Fig 4 (erythroid markers + significance asterisks) |
| `make_supplementary_v2.py` | Supplementary PDF (Table S1 + Fig S1–S3) |
| `remake_figs3.py` | Re-draw Supplementary Fig S3 (transcriptomic signature, neutral titles) |

## Key analytical decisions

- **TP53** judged adverse at VAF ≥10% (ELN 2022 convention), not any mutation.
- **Complete-case Cox**: n = 759 after excluding APL (n = 784) and dropping 25 cases with missing covariates.
- **Sensitivity analyses**: stratified Cox (strata = ELN + chemo + age group, for PH violations); initial-diagnosis-specimen subset (n = 543); transplant-excluded.
- **Drug-resistance partitioning**: four-factor type II ANOVA (ontogeny, MR mutations, TP53, complex karyotype) on rank-transformed AUC; complementary two-way ANOVA (ontogeny × MR) in the 2×2 analysis.

## Citation

If you reuse this code, please cite the manuscript and the primary data source (Bottomly D, et al. *Cancer Cell*. 2022;40(8):850-864). A permanent archived version of this code is available on Zenodo (https://doi.org/10.5281/zenodo.24514059).

## License

MIT
