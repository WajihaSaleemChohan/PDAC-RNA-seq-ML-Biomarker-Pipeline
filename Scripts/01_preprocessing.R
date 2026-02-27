install.packages("BiocManager")
BiocManager::install(c("DESeq2"))
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"))
BiocManager::install("EnhancedVolcano")

# ============================================
# Script 01: Load and Preprocess TCGA-PAAD
# ============================================

library(dplyr)
library(tibble)

# ---- Set Working Directory ----
setwd("F:/")

# ---- Load Data ----
cat("Loading expression matrix...\n")
expr <- read.delim("data_mrna_seq_v2_rsem.txt", 
                   header=TRUE, stringsAsFactors=FALSE)

cat("Loading clinical data...\n")
clinical <- read.delim("data_clinical_patient.txt",
                       header=TRUE, stringsAsFactors=FALSE, 
                       comment.char="#")

# ---- Check what loaded ----
cat("Expression dimensions:", dim(expr), "\n")
cat("Clinical dimensions:", dim(clinical), "\n")
cat("Clinical column names:\n")
print(colnames(clinical))

# ============================================
# Script 01: Preprocess TCGA-PAAD Data
# ============================================

library(dplyr)
library(tibble)

setwd("F:/rnaseq-ml-biomarker-pdac")

# ---- Load Data ----
expr <- read.delim("data/raw/data_mrna_seq_v2_rsem.txt",
                   header=TRUE, stringsAsFactors=FALSE)

clinical <- read.delim("data/raw/data_clinical_patient.txt",
                       header=TRUE, stringsAsFactors=FALSE,
                       comment.char="#")

# ---- Clean Expression Matrix ----
# Remove duplicated genes
expr <- expr[!duplicated(expr$Hugo_Symbol), ]

# Set gene symbols as rownames
rownames(expr) <- expr$Hugo_Symbol
expr <- expr[ , !(names(expr) %in% c("Hugo_Symbol", "Entrez_Gene_Id"))]

# Round to integers for DESeq2
expr <- round(expr)

# Remove genes where all samples have zero counts
expr <- expr[rowSums(expr) > 0, ]

# Fix sample ID format: TCGA.XX.XXXX.01 -> TCGA-XX-XXXX
colnames(expr) <- gsub("\\.", "-", colnames(expr))
colnames(expr) <- substr(colnames(expr), 1, 12)

cat("Expression after cleaning:", nrow(expr), "genes x", ncol(expr), "samples\n")

# ---- Clean Clinical Data ----
clinical_clean <- clinical %>%
  select(PATIENT_ID, OS_STATUS, OS_MONTHS,
         AGE, SEX, PATH_T_STAGE) %>%
  filter(!is.na(OS_STATUS) & OS_STATUS != "")

# Create binary label: 1 = DECEASED, 0 = LIVING
clinical_clean$survival_label <- ifelse(
  grepl("DECEASED", clinical_clean$OS_STATUS), 1, 0)

cat("Clinical samples after cleaning:", nrow(clinical_clean), "\n")
cat("DECEASED:", sum(clinical_clean$survival_label == 1), "\n")
cat("LIVING:  ", sum(clinical_clean$survival_label == 0), "\n")

# ---- Match Samples Between Expression and Clinical ----
common_samples <- intersect(colnames(expr), clinical_clean$PATIENT_ID)
cat("Matched samples:", length(common_samples), "\n")

expr_matched <- expr[ , common_samples]
clinical_matched <- clinical_clean[
  match(common_samples, clinical_clean$PATIENT_ID), ]

# Final check
cat("\n✅ Final dimensions:\n")
cat("Expression matrix:", nrow(expr_matched), "genes x",
    ncol(expr_matched), "samples\n")
cat("Clinical table:   ", nrow(clinical_matched), "samples x",
    ncol(clinical_matched), "columns\n")

# ---- Save Processed Files ----
dir.create("data/processed", showWarnings=FALSE, recursive=TRUE)

saveRDS(expr_matched,     "data/processed/expr_matched.rds")
saveRDS(clinical_matched, "data/processed/clinical_matched.rds")
write.csv(clinical_matched, "data/processed/clinical_matched.csv",
          row.names=FALSE)

cat("\n✅ Saved to data/processed/\n")