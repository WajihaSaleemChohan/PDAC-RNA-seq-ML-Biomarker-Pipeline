<div align="center">

# 🧬 TCGA-PAAD RNA-seq Biomarker Pipeline

### *Pancreatic Cancer · Differential Expression · Machine Learning · Survival Analysis*

<br>

![R](https://img.shields.io/badge/R-4.5.0-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.10-3776AB?style=for-the-badge&logo=python&logoColor=white)
![License](https://img.shields.io/badge/License-MIT-yellow?style=for-the-badge)
![Status](https://img.shields.io/badge/Status-Complete-success?style=for-the-badge)

<br>

![Samples](https://img.shields.io/badge/Samples-177%20TCGA--PAAD-blueviolet?style=flat-square)
![AUC](https://img.shields.io/badge/SVM%20AUC-0.72-orange?style=flat-square)
![pval](https://img.shields.io/badge/Log--rank%20p-1.03e--7-red?style=flat-square)
![Gene](https://img.shields.io/badge/Top%20Gene-SYT3-brightgreen?style=flat-square)
![DESeq2](https://img.shields.io/badge/Tool-DESeq2-blue?style=flat-square)

<br>

> **End-to-end ML/AI-enabled bioinformatics pipeline for PDAC biomarker discovery**  
> *Transforming raw RNA-seq counts into clinically actionable survival insights*

</div>

---

## 🔬 Overview

This project builds a **production-grade, reproducible bioinformatics pipeline** analyzing RNA-seq data from **177 TCGA-PAAD** (Pancreatic Ductal Adenocarcinoma) samples. The pipeline covers the full workflow from raw counts to clinical biomarker discovery.

> 💡 Pancreatic cancer has a 5-year survival rate of under 12%. Identifying transcriptomic biomarkers for patient stratification is critical for improving clinical outcomes.

---

## 📊 Key Results

| Metric | Value |
|--------|-------|
| 🧫 Dataset | TCGA-PAAD — 177 matched samples |
| 🔬 Genes Analyzed | 20,531 |
| 🤖 Best ML Model | SVM (Radial Kernel) |
| 📈 SVM Test AUC | **0.72** |
| 📉 Log-rank p-value | **1.03e-7** |
| ⏱️ Median OS — Low Risk | **71.7 months** |
| ⏱️ Median OS — High Risk | **15.5 months** |
| 🏆 Top Biomarker Gene | **SYT3** |

---

## 📈 Output Figures

> After pushing to GitHub, go to each figure → click **Raw** → copy URL → paste below as `![name](raw_url)`

**Volcano Plot — Differential Expression**
![Volcano](results/figures/volcano_DESeq2.png)

**ROC Curves — RF vs LASSO vs SVM**
![ROC](results/figures/roc_curves_v2.png)

**Feature Importance — Top 20 Biomarker Genes**
![Importance](results/figures/feature_importance_v2.png)

**Heatmap — Top 20 Genes Across 177 Samples**
![Heatmap](results/figures/heatmap_top20_v2.png)

**Kaplan-Meier — SVM Risk Stratification**
![KM SVM](results/figures/kaplan_meier_svm.png)

**Kaplan-Meier — SYT3 Expression**
![KM SYT3](results/figures/kaplan_meier_SYT3.png)

---

## 🔄 Pipeline

```
TCGA-PAAD Input (RNA-seq + Clinical)
         │
         ▼
  01_preprocess.R       →  Filter · Match · Label (177 samples)
         │
         ▼
  02_deseq2.R           →  DE Analysis · Volcano · Heatmap
         │
         ▼
  03_ml_classifier.R    →  RF · LASSO · SVM · ROC · Feature Importance
         │
         ▼
  04_survival.R         →  Kaplan-Meier · Log-rank · Risk Groups
         │
         ▼
  05_visualization.R →  6-panel publication figure (300 DPI)
```

---

## 📁 Project Structure

```
rnaseq-ml-biomarker-pdac/
│
├── data/
│   ├── raw/
│   │   ├── data_mrna_seq_v2_rsem.txt        # RNA-seq matrix (20,531 x 179)
│   │   └── data_clinical_patient.txt         # Clinical metadata (184 x 38)
│   └── processed/
│       ├── expr_matched.rds
│       └── clinical_matched.rds
│
├── scripts/
│   ├── 01_preprocess.R
│   ├── 02_deseq2.R
│   ├── 03_ml_classifier.R
│   ├── 04_survival.R
│   └── 05_visualization.R
│
├── results/
│   ├── deg/
│   │   ├── sig_degs.rds
│   │   └── survival_risk_scores.csv
│   ├── models/
│   │   ├── rf_model.rds
│   │   ├── lasso_model.rds
│   │   └── svm_model.rds
│   └── figures/
│       ├── volcano_DESeq2.png
│       ├── roc_curves_v2.png
│       ├── feature_importance_v2.png
│       ├── heatmap_top20_v2.png
│       ├── kaplan_meier_svm.png
│       ├── kaplan_meier_SYT3.png
│       └── visualization.png
│
├── README.md
└── .gitignore
```

---

## 🚀 Quick Start

### 1. Clone
```bash
git clone https://github.com/YOURUSERNAME/rnaseq-ml-biomarker-pdac.git
cd rnaseq-ml-biomarker-pdac
```

### 2. Download Data
From [cBioPortal TCGA-PAAD](https://www.cbioportal.org/study/summary?id=paad_tcga_pan_can_atlas_2018), place in `data/raw/`:
```
data_mrna_seq_v2_rsem.txt
data_clinical_patient.txt
```

### 3. Install R Packages
```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano"))
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap",
                   "caret", "randomForest", "glmnet", "kernlab",
                   "pROC", "survival", "survminer", "patchwork"))
```

### 4. Run
```r
setwd("path/to/rnaseq-ml-biomarker-pdac")
source("scripts/01_preprocess.R")
source("scripts/02_deseq2.R")
source("scripts/03_ml_classifier.R")
source("scripts/04_survival.R")
source("scripts/05_visualization.R")
```

---

## 🧪 Methods

**Preprocessing** — Matched 177 samples between RNA-seq matrix and clinical metadata. Removed duplicate/zero genes. Binary survival label: DECEASED=1, LIVING=0.

**DESeq2 DE Analysis** — Negative binomial model. Significance: padj < 0.05 and |log2FC| > 1. VST normalization for downstream ML features.

**Machine Learning** — Top 50 DEGs as features. 80/20 train-test split. 5-fold repeated cross-validation (10 repeats). Random Forest (mtry=7), LASSO (lambda=0.011), SVM-RBF. Best model: **SVM AUC = 0.72**.

**Survival Analysis** — Patients split into High/Low risk by SVM probability (median threshold). Kaplan-Meier + log-rank test. p = 1.03e-7. Median OS: Low Risk 71.7 months vs High Risk 15.5 months.

---

## 🛠️ Tech Stack

| Category | Tools |
|----------|-------|
| Language | R 4.5.0, Python 3.10 |
| DE Analysis | DESeq2, EnhancedVolcano |
| ML Framework | caret, randomForest, glmnet, kernlab |
| Survival | survival, survminer |
| Visualization | ggplot2, pheatmap, patchwork |
| Data | TCGA-PAAD via cBioPortal |

---

## 🔭 Future Work

- [ ] Multi-omics integration (mutations + methylation + expression)
- [ ] SHAP values for ML interpretability
- [ ] GSEA pathway enrichment analysis
- [ ] Docker containerization for full reproducibility
- [ ] scRNA-seq subtype deconvolution

---

## 📄 Data Citation

```
Cancer Genome Atlas Research Network. (2017).
Integrated Genomic Characterization of Pancreatic Ductal Adenocarcinoma.
Cancer Cell, 32(2), 185-203.
```

Data from [cBioPortal for Cancer Genomics](https://www.cbioportal.org).

---

## 👩‍💻 Author

**Wajiha Saleem**  
Bioinformatician | RNA-seq · Machine Learning · Clinical Genomics

[![GitHub](https://img.shields.io/badge/GitHub-Follow-black?style=for-the-badge&logo=github)](https://github.com/YOURUSERNAME)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Connect-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://linkedin.com/in/YOURPROFILE)

---

<div align="center">

⭐ **Star this repo if it helped you!** ⭐

*TCGA-PAAD · Bioinformatics · Clinical Genomics · 2026*

</div>
