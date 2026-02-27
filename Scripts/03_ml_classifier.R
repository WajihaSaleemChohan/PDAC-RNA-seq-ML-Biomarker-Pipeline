# ============================================
# Script 03 (v2): Feature Selection + ML (Fixed)
# ============================================
library(dplyr)
library(tibble)
library(caret)
library(randomForest)
library(pROC)
library(ggplot2)
library(pheatmap)
library(DESeq2)

setwd("F:/rnaseq-ml-biomarker-pdac")
install.packages("kernlab")
# ---- Load Data ----
clinical_matched <- readRDS("data/processed/clinical_matched.rds")
sig_degs         <- readRDS("results/deg/sig_degs.rds")
dds              <- readRDS("results/deg/dds.rds")

# ---- VST Normalization ----
vst_mat <- assay(vst(dds, blind = FALSE))
cat("VST matrix:", nrow(vst_mat), "x", ncol(vst_mat), "\n")

# ---- Use Top 50 DEGs (more conservative) ----
top_genes <- sig_degs %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(gene)

# Build ML dataset
feat_mat <- vst_mat[top_genes, ] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("PATIENT_ID")

ml_data <- feat_mat %>%
  left_join(
    clinical_matched %>% select(PATIENT_ID, survival_label),
    by = "PATIENT_ID"
  ) %>%
  select(-PATIENT_ID) %>%
  mutate(survival_label = factor(survival_label,
                                 levels = c(0, 1),
                                 labels = c("LIVING", "DECEASED")))

cat("ML dataset:", nrow(ml_data), "samples x",
    ncol(ml_data) - 1, "features\n")

# ---- Train/Test Split (80/20 — more training data) ----
set.seed(42)
train_idx  <- createDataPartition(ml_data$survival_label,
                                  p = 0.8, list = FALSE)
train_data <- ml_data[ train_idx, ]
test_data  <- ml_data[-train_idx, ]
cat("Train:", nrow(train_data), "| Test:", nrow(test_data), "\n")

# ---- Cross-Validation: Repeated CV for stability ----
cv_ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 5,
  repeats         = 10,        # repeated CV much more stable with small N
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# ---- Model 1: Random Forest ----
cat("\nTraining Random Forest...\n")
set.seed(42)
rf_model <- train(
  survival_label ~ .,
  data     = train_data,
  method   = "rf",
  metric   = "ROC",
  trControl = cv_ctrl,
  tuneGrid  = expand.grid(mtry = c(5, 7, 10, 15)),
  ntree     = 1000
)
cat("Best mtry:", rf_model$bestTune$mtry, "\n")
cat("CV AUC:   ", round(max(rf_model$results$ROC), 4), "\n")

# ---- Model 2: LASSO ----
cat("\nTraining LASSO...\n")
set.seed(42)
lasso_model <- train(
  survival_label ~ .,
  data      = train_data,
  method    = "glmnet",
  metric    = "ROC",
  trControl = cv_ctrl,
  tuneGrid  = expand.grid(
    alpha  = 1,
    lambda = 10^seq(-4, 0, length = 50)
  )
)
cat("Best lambda:", round(lasso_model$bestTune$lambda, 6), "\n")
cat("CV AUC:     ", round(max(lasso_model$results$ROC), 4), "\n")
library(kernlab)

# ---- Model 3: SVM ----
cat("\nTraining SVM...\n")
set.seed(42)
svm_model <- train(
  survival_label ~ .,
  data       = train_data,
  method     = "svmRadial",
  metric     = "ROC",
  trControl  = cv_ctrl,
  preProcess = c("center", "scale"),
  tuneLength = 5
)
cat("CV AUC:", round(max(svm_model$results$ROC), 4), "\n")

# ---- Evaluate on Test Set ----
evaluate_model <- function(model, test_data, model_name) {
  preds   <- predict(model, test_data)
  probs   <- predict(model, test_data, type = "prob")[, "DECEASED"]
  roc_obj <- roc(test_data$survival_label, probs,
                 levels = c("LIVING", "DECEASED"), quiet = TRUE)
  auc_val <- auc(roc_obj)
  cm      <- confusionMatrix(preds, test_data$survival_label,
                             positive = "DECEASED")
  cat("\n---", model_name, "---\n")
  cat("AUC:        ", round(auc_val, 4), "\n")
  cat("Accuracy:   ", round(cm$overall["Accuracy"], 4), "\n")
  cat("Sensitivity:", round(cm$byClass["Sensitivity"], 4), "\n")
  cat("Specificity:", round(cm$byClass["Specificity"], 4), "\n")
  list(roc = roc_obj, auc = auc_val, cm = cm)
}

rf_eval    <- evaluate_model(rf_model,    test_data, "Random Forest")
lasso_eval <- evaluate_model(lasso_model, test_data, "LASSO")
svm_eval   <- evaluate_model(svm_model,   test_data, "SVM")

# ---- ROC Plot ----
roc_plot <- ggroc(list(
  "Random Forest" = rf_eval$roc,
  "LASSO"         = lasso_eval$roc,
  "SVM"           = svm_eval$roc
), linewidth = 1) +
  scale_color_manual(values = c(
    "Random Forest" = "#E41A1C",
    "LASSO"         = "#377EB8",
    "SVM"           = "#4DAF4A")) +
  geom_abline(slope = 1, intercept = 1,
              linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.3, y = 0.20,
           label = paste0("RF AUC = ", round(rf_eval$auc, 3)),
           color = "#E41A1C", size = 4) +
  annotate("text", x = 0.3, y = 0.13,
           label = paste0("LASSO AUC = ", round(lasso_eval$auc, 3)),
           color = "#377EB8", size = 4) +
  annotate("text", x = 0.3, y = 0.06,
           label = paste0("SVM AUC = ", round(svm_eval$auc, 3)),
           color = "#4DAF4A", size = 4) +
  labs(title = "TCGA-PAAD: ROC Curves on Test Set", color = "Model") +
  theme_bw(base_size = 13)

ggsave("results/figures/roc_curves_v2.png", roc_plot,
       width = 7, height = 6, dpi = 300)
cat("ROC plot saved.\n")

# ---- Feature Importance ----
imp_df <- varImp(rf_model)$importance %>%
  rownames_to_column("gene") %>%
  arrange(desc(Overall)) %>%
  slice_head(n = 20)

imp_plot <- ggplot(imp_df, aes(x = reorder(gene, Overall), y = Overall)) +
  geom_col(fill = "#E41A1C", alpha = 0.8) +
  coord_flip() +
  labs(title = "Top 20 Biomarker Genes — Random Forest",
       x = "", y = "Mean Decrease Gini") +
  theme_bw(base_size = 12)

ggsave("results/figures/feature_importance_v2.png", imp_plot,
       width = 7, height = 6, dpi = 300)
cat("Feature importance plot saved.\n")

# ---- Heatmap ----
top20_genes <- intersect(imp_df$gene, rownames(vst_mat))
cat("Genes in heatmap:", length(top20_genes), "\n")

heatmap_mat <- vst_mat[top20_genes, ] %>% t() %>% scale() %>% t()

annotation_col <- data.frame(
  Survival = factor(clinical_matched$survival_label,
                    levels = c(0,1), labels = c("LIVING","DECEASED"))
)
rownames(annotation_col) <- clinical_matched$PATIENT_ID

png("results/figures/heatmap_top20_v2.png",
    width = 10, height = 7, units = "in", res = 300)
pheatmap(heatmap_mat,
         annotation_col    = annotation_col,
         ann_colors        = list(Survival = c("LIVING"="#377EB8","DECEASED"="#E41A1C")),
         show_colnames     = FALSE,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
         main = "Top 20 Biomarker Genes — TCGA-PAAD")
dev.off()

# ---- Save ----
dir.create("results/models", showWarnings = FALSE, recursive = TRUE)
saveRDS(rf_model,    "results/models/rf_model.rds")
saveRDS(lasso_model, "results/models/lasso_model.rds")
saveRDS(svm_model,   "results/models/svm_model.rds")
saveRDS(imp_df,      "results/models/rf_importance.rds")

cat("\n✅ Script 03 v2 complete.\n")
