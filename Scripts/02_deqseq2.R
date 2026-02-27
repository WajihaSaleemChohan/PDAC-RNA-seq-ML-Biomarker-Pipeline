# ============================================
# Script 02: Differential Expression (DESeq2)
# ============================================
library(DESeq2)
library(dplyr)
library(ggplot2)

setwd("F:/rnaseq-ml-biomarker-pdac")

# ---- Load Processed Data ----
expr_matched     <- readRDS("data/processed/expr_matched.rds")
clinical_matched <- readRDS("data/processed/clinical_matched.rds")

# ---- Check Class Balance ----
cat("DECEASED:", sum(clinical_matched$survival_label == 1), "\n")
cat("LIVING:  ", sum(clinical_matched$survival_label == 0), "\n")

# ---- Build DESeq2 Object ----
# survival_label must be a factor
col_data <- clinical_matched %>%
  select(PATIENT_ID, survival_label) %>%
  mutate(survival_label = factor(survival_label,
                                 levels = c(0, 1),
                                 labels = c("LIVING", "DECEASED")))

rownames(col_data) <- col_data$PATIENT_ID

# Ensure column order matches
expr_matched <- expr_matched[ , rownames(col_data)]

dds <- DESeqDataSetFromMatrix(
  countData = expr_matched,
  colData   = col_data,
  design    = ~ survival_label
)

# ---- Pre-filter: Keep genes with >= 10 counts in at least 10 samples ----
keep <- rowSums(counts(dds) >= 10) >= 10
dds  <- dds[keep, ]
cat("Genes after pre-filtering:", nrow(dds), "\n")

# ---- Run DESeq2 ----
dds <- DESeq(dds)
res <- results(dds,
               contrast  = c("survival_label", "DECEASED", "LIVING"),
               alpha     = 0.05)

summary(res)

# ---- Extract Significant DEGs ----
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  filter(!is.na(padj))

# Significant: padj < 0.05 and |log2FC| > 1
sig_degs <- res_df %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
  arrange(padj)

cat("\nSignificant DEGs (padj<0.05, |LFC|>1):", nrow(sig_degs), "\n")
cat("Upregulated in DECEASED:  ",
    sum(sig_degs$log2FoldChange > 0), "\n")
cat("Downregulated in DECEASED:",
    sum(sig_degs$log2FoldChange < 0), "\n")

# ---- Volcano Plot ----
res_df <- res_df %>%
  mutate(significance = case_when(
    padj < 0.05 & log2FoldChange >  1 ~ "Up in DECEASED",
    padj < 0.05 & log2FoldChange < -1 ~ "Down in DECEASED",
    TRUE                               ~ "NS"
  ))

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj),
                              color = significance)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_color_manual(values = c(
    "Up in DECEASED"   = "#E41A1C",
    "Down in DECEASED" = "#377EB8",
    "NS"               = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "TCGA-PAAD: DECEASED vs LIVING",
       x = "log2 Fold Change", y = "-log10(adjusted p-value)",
       color = "") +
  theme_bw(base_size = 13)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("results/figures/volcano_DESeq2.png", volcano,
       width = 8, height = 6, dpi = 300)
cat("Volcano plot saved.\n")

# ---- Save DEG Results ----
dir.create("results/deg", showWarnings = FALSE, recursive = TRUE)
write.csv(res_df,    "results/deg/all_genes_deseq2.csv",  row.names = FALSE)
write.csv(sig_degs,  "results/deg/significant_degs.csv",  row.names = FALSE)
saveRDS(dds,         "results/deg/dds.rds")
saveRDS(sig_degs,    "results/deg/sig_degs.rds")

cat("\n✅ DESeq2 complete. Results saved to results/deg/\n")
