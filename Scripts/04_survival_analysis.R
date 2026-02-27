install.packages("survminer")

# ============================================
# Script 04: Kaplan-Meier Survival Analysis
# ============================================
library(survival)
library(survminer)
library(dplyr)
library(tibble)
library(DESeq2)

setwd("F:/rnaseq-ml-biomarker-pdac")

# ---- Load Data ----
clinical_matched <- readRDS("data/processed/clinical_matched.rds")
sig_degs         <- readRDS("results/deg/sig_degs.rds")
dds              <- readRDS("results/deg/dds.rds")
svm_model        <- readRDS("results/models/svm_model.rds")
imp_df           <- readRDS("results/models/rf_importance.rds")

# ---- Rebuild full ML dataset (all 177 samples) ----
vst_mat <- assay(vst(dds, blind = FALSE))

top_genes <- sig_degs %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(gene)

feat_mat <- vst_mat[top_genes, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("PATIENT_ID")

ml_data_full <- feat_mat %>%
  left_join(
    clinical_matched %>% select(PATIENT_ID, survival_label),
    by = "PATIENT_ID"
  ) %>%
  mutate(survival_label = factor(survival_label,
                                 levels = c(0, 1),
                                 labels = c("LIVING", "DECEASED")))

# ---- Get SVM Predicted Probabilities for All Samples ----
probs_full <- predict(svm_model,
                      newdata = ml_data_full %>% select(-PATIENT_ID, -survival_label),
                      type = "prob")[, "DECEASED"]

# ---- Build Survival Dataset ----
surv_df <- clinical_matched %>%
  select(PATIENT_ID, OS_MONTHS, survival_label) %>%
  mutate(
    risk_score = probs_full,
    risk_group = ifelse(risk_score >= median(risk_score),
                        "High Risk", "Low Risk"),
    risk_group = factor(risk_group, levels = c("Low Risk", "High Risk")),
    OS_MONTHS  = as.numeric(OS_MONTHS),
    event      = as.integer(survival_label == 1)
  ) %>%
  filter(!is.na(OS_MONTHS))

cat("High Risk:", sum(surv_df$risk_group == "High Risk"), "\n")
cat("Low Risk: ", sum(surv_df$risk_group == "Low Risk"),  "\n")

# ---- Kaplan-Meier Fit ----
km_fit <- survfit(Surv(OS_MONTHS, event) ~ risk_group, data = surv_df)
print(km_fit)

# Log-rank test
lr_test <- survdiff(Surv(OS_MONTHS, event) ~ risk_group, data = surv_df)
pval    <- 1 - pchisq(lr_test$chisq, df = 1)
cat("\nLog-rank p-value:", signif(pval, 3), "\n")

# ---- KM Plot ----
km_plot <- ggsurvplot(
  km_fit,
  data          = surv_df,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.height = 0.25,
  palette       = c("#377EB8", "#E41A1C"),
  legend.labs   = c("Low Risk", "High Risk"),
  legend.title  = "SVM Risk Group",
  title         = "Overall Survival by SVM-Predicted Risk — TCGA-PAAD",
  xlab          = "Time (Months)",
  ylab          = "Survival Probability",
  ggtheme       = theme_bw(base_size = 13),
  font.main     = c(13, "bold"),
  surv.median.line = "hv"   # draw median survival lines
)

# Save
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
png("results/figures/kaplan_meier_svm.png",
    width = 9, height = 7, units = "in", res = 300)
print(km_plot)
dev.off()
cat("KM plot saved.\n")

# ---- Bonus: KM for Top Individual Gene (most important from RF) ----
top_gene <- imp_df$gene[1]
cat("\nTop biomarker gene:", top_gene, "\n")

gene_expr <- vst_mat[top_gene, ]

surv_df$gene_expr  <- gene_expr[surv_df$PATIENT_ID]
surv_df$gene_group <- ifelse(surv_df$gene_expr >= median(surv_df$gene_expr, na.rm = TRUE),
                             "High Expression", "Low Expression")
surv_df$gene_group <- factor(surv_df$gene_group,
                             levels = c("Low Expression", "High Expression"))

km_fit_gene <- survfit(Surv(OS_MONTHS, event) ~ gene_group, data = surv_df)

km_gene_plot <- ggsurvplot(
  km_fit_gene,
  data          = surv_df,
  pval          = TRUE,
  pval.method   = TRUE,
  conf.int      = TRUE,
  risk.table    = TRUE,
  risk.table.height = 0.25,
  palette       = c("#4DAF4A", "#984EA3"),
  legend.labs   = c("Low Expression", "High Expression"),
  legend.title  = top_gene,
  title         = paste0("Overall Survival by ", top_gene,
                         " Expression — TCGA-PAAD"),
  xlab          = "Time (Months)",
  ylab          = "Survival Probability",
  ggtheme       = theme_bw(base_size = 13),
  surv.median.line = "hv"
)

png(paste0("results/figures/kaplan_meier_", top_gene, ".png"),
    width = 9, height = 7, units = "in", res = 300)
print(km_gene_plot)
dev.off()
cat("Gene-level KM plot saved.\n")

# ---- Save survival dataframe ----
write.csv(surv_df, "results/deg/survival_risk_scores.csv", row.names = FALSE)

cat("\n✅ Script 04 complete.\n")

cat("Figures saved to results/figures/\n")