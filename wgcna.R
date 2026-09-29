# ==========================================
# WGCNA + ML PIPELINE: SETUP
# ==========================================

library(WGCNA)
library(limma)
library(dynamicTreeCut)
library(clusterProfiler)
library(org.Hs.eg.db)
library(glmnet)
library(randomForest)
library(e1071)
library(pROC)
library(caret)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
set.seed(123)

# ---- Load processed data (from your existing dataset_list) ----
# GSE99039 = discovery (n=558), GSE6613 = validation (n=72)

disc <- dataset_list[["GSE99039"]]
val  <- dataset_list[["GSE6613"]]

# WGCNA requires samples in rows, genes in columns
# Transpose the expression matrices
datExpr_disc <- t(disc$expr)
datExpr_val  <- t(val$expr)

# Ensure common genes for validation
common_genes <- intersect(colnames(datExpr_disc), colnames(datExpr_val))
datExpr_disc <- datExpr_disc[, common_genes]
datExpr_val  <- datExpr_val[, common_genes]

cat(paste("Discovery:", nrow(datExpr_disc), "samples x", ncol(datExpr_disc), "genes\n"))
cat(paste("Validation:", nrow(datExpr_val), "samples x", ncol(datExpr_val), "genes\n"))

# ---- Trait data (PD vs Control) ----
traitData <- data.frame(
  PD = ifelse(disc$disease == "PD", 1, 0)
)
rownames(traitData) <- rownames(datExpr_disc)

cat(paste("\nDiscovery samples: PD =", sum(traitData$PD),
          "| Control =", sum(traitData$PD == 0), "\n"))
# ==========================================
# BLOCK 2: DIFFERENTIAL EXPRESSION ANALYSIS
# Corrected: genes × samples for limma
# ==========================================

disc <- dataset_list[["GSE99039"]]

# Verify alignment
stopifnot(ncol(disc$expr) == length(disc$disease))

# Design matrix
disease_factor <- factor(disc$disease, levels = c("Control", "PD"))
design <- model.matrix(~ 0 + disease_factor)
colnames(design) <- c("Control", "PD")

cat("Design matrix:", nrow(design), "x", ncol(design), "\n")

# ---- DE analysis: genes in rows, samples in columns ----
fit  <- lmFit(disc$expr, design)          # NOT t(disc$expr)
fit2 <- contrasts.fit(fit, makeContrasts(PD - Control, levels = design))
fit2 <- eBayes(fit2)
deg_results <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")

cat("\nTop 10 DEGs:\n")
print(head(deg_results, 10))

# ---- Relaxed significance for WGCNA intersection ----
# Strict criteria give only 10 genes — too few for intersection.
# Use relaxed thresholds for intersection with WGCNA modules.
deg_sig_strict <- deg_results[deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 0.3, ]
deg_sig_relaxed <- deg_results[deg_results$P.Value < 0.05 & abs(deg_results$logFC) > 0.1, ]

cat(paste("\nStrict DEGs (adj.P < 0.05, |logFC| > 0.3):", nrow(deg_sig_strict), "\n"))
cat(paste("Relaxed DEGs (P < 0.05, |logFC| > 0.1):", nrow(deg_sig_relaxed), "\n"))

# Save both
write.csv(deg_results, "PD_Biomarker_Project/tables/WGCNA_DEG_results.csv")
write.csv(deg_sig_relaxed, "PD_Biomarker_Project/tables/WGCNA_DEG_relaxed.csv")

cat("\n✅ BLOCK 2 COMPLETE\n")

# ==========================================
# BLOCK 3: WGCNA NETWORK CONSTRUCTION
# ==========================================

# ---- Quality control ----
gsg <- goodSamplesGenes(datExpr_disc, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", sum(!gsg$goodGenes)))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", sum(!gsg$goodSamples)))
  datExpr_disc <- datExpr_disc[gsg$goodSamples, gsg$goodGenes]
}

cat("After QC:", nrow(datExpr_disc), "samples x", ncol(datExpr_disc), "genes\n")

# ---- Sample clustering ----
sampleTree <- hclust(dist(datExpr_disc), method = "average")
pdf("PD_Biomarker_Project/figures/WGCNA_sample_clustering.pdf", width = 12, height = 9)
par(cex = 0.6, mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

# ---- Choose soft-thresholding power ----
powers <- c(1:10, seq(12, 20, by = 2))

sft <- pickSoftThreshold(datExpr_disc, powerVector = powers,
                         networkType = "signed", verbose = 5)

pdf("PD_Biomarker_Project/figures/WGCNA_soft_threshold.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.80, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5],
     labels = powers, cex = cex1, col = "red")
dev.off()

# Pick lowest power with R^2 > 0.80
chosen_power <- sft$fitIndices$Power[which(sft$fitIndices[,2] > 0.80)[1]]
if (is.na(chosen_power)) chosen_power <- 12
cat(paste("\nSelected soft-thresholding power:", chosen_power, "\n"))

# ---- Build network and detect modules ----
net <- blockwiseModules(
  datExpr_disc,
  power = chosen_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  deepSplit = 2,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
cat("\nModule sizes:\n")
print(table(moduleColors))

# ---- Dendrogram plot ----
pdf("PD_Biomarker_Project/figures/WGCNA_dendrogram.pdf", width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

cat("\n✅ BLOCK 3 COMPLETE\n")

# ==========================================
# BLOCK 4: MODULE-TRAIT CORRELATION
# ==========================================

MEs <- net$MEs

# Correlate module eigengenes with PD status
moduleTraitCor <- cor(MEs, traitData$PD, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr_disc))

moduleTraitResults <- data.frame(
  Module = rownames(moduleTraitCor),
  Correlation = round(moduleTraitCor, 3),
  Pvalue = signif(moduleTraitPvalue, 3)
)
moduleTraitResults <- moduleTraitResults[order(-abs(moduleTraitResults$Correlation)), ]
print(moduleTraitResults)

# Significant modules (|cor| > 0.2, p < 0.05)
sig_modules <- moduleTraitResults$Module[
  abs(moduleTraitResults$Correlation) > 0.2 &
    moduleTraitResults$Pvalue < 0.05
]
cat(paste("\nSignificant PD-associated modules:",
          paste(sig_modules, collapse = ", "), "\n"))

# ---- Module-trait heatmap ----
pdf("PD_Biomarker_Project/figures/WGCNA_module_trait_heatmap.pdf", width = 8, height = 6)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = t(moduleTraitCor),
               xLabels = rownames(moduleTraitCor),
               yLabels = "PD status",
               ySymbols = "PD",
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = round(moduleTraitCor, 2),
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

# ---- Extract genes from significant modules ----
module_genes <- list()
for (mod in sig_modules) {
  mod_color <- sub("ME", "", mod)
  module_genes[[mod_color]] <- colnames(datExpr_disc)[moduleColors == mod_color]
  cat(paste("Module", mod_color, ":", length(module_genes[[mod_color]]), "genes\n"))
}

all_module_genes <- unique(unlist(module_genes))
cat(paste("\nTotal genes in significant modules:", length(all_module_genes), "\n"))

cat("\n✅ BLOCK 4 COMPLETE\n")


# ==========================================
# BLOCK 5: INTERSECT DEGs WITH MODULE GENES
# ==========================================

# Use relaxed DEGs to get enough genes for intersection
overlapping_genes <- intersect(rownames(deg_sig_relaxed), all_module_genes)

cat(paste("Relaxed DEGs:", nrow(deg_sig_relaxed), "\n"))
cat(paste("Genes in significant modules:", length(all_module_genes), "\n"))
cat(paste("Overlapping genes:", length(overlapping_genes), "\n"))

# If too few, fall back to all module genes
if (length(overlapping_genes) < 20) {
  cat("\n⚠️ Too few overlapping genes. Using all significant module genes.\n")
  overlapping_genes <- all_module_genes
}

writeLines(overlapping_genes,
           "PD_Biomarker_Project/tables/WGCNA_overlapping_genes.txt")

cat(paste("\nFinal candidate gene set:", length(overlapping_genes), "genes\n"))

cat("\n✅ BLOCK 5 COMPLETE\n")

# ==========================================
# BLOCK 6: LASSO HUB GENE SELECTION
# ==========================================

# Prepare data
X <- t(disc$expr[overlapping_genes, ])
y <- as.numeric(disc$disease == "PD")

cat("Input to LASSO:", nrow(X), "samples x", ncol(X), "genes\n")

# LASSO with cross-validation
set.seed(123)
cv_lasso <- cv.glmnet(as.matrix(X), y,
                      family = "binomial",
                      alpha = 1,
                      nfolds = 10)

# Extract selected genes
lasso_coef <- coef(cv_lasso, s = "lambda.min")
selected_genes <- rownames(lasso_coef)[lasso_coef[,1] != 0][-1]

cat(paste("\nLASSO selected:", length(selected_genes), "hub genes\n"))
cat(paste("Genes:", paste(selected_genes, collapse = ", "), "\n"))

# LASSO CV plot
pdf("PD_Biomarker_Project/figures/WGCNA_LASSO_cv.pdf", width = 8, height = 6)
plot(cv_lasso)
dev.off()

writeLines(selected_genes,
           "PD_Biomarker_Project/tables/WGCNA_LASSO_hub_genes.txt")

cat("\n✅ BLOCK 6 COMPLETE\n")

# ==========================================
# BLOCK 7: MODEL CONSTRUCTION AND VALIDATION
# ==========================================

# ---- Logistic regression on discovery ----
X_selected <- as.data.frame(t(disc$expr[selected_genes, ]))
X_selected$PD <- ifelse(disc$disease == "PD", 1, 0)

model <- glm(PD ~ ., data = X_selected, family = "binomial")

# Predictions on discovery
pred_disc <- predict(model, type = "response")
roc_disc <- roc(disc$disease, pred_disc, quiet = TRUE)
cat(paste("Discovery (GSE99039) AUC:", round(auc(roc_disc), 3), "\n"))

# ---- Validate on GSE6613 ----
# Check which selected genes are available
available_val <- intersect(selected_genes, rownames(val$expr))
cat(paste("\nGenes available in GSE6613:", length(available_val), "/",
          length(selected_genes), "\n"))

if (length(available_val) >= 3) {
  X_val <- as.data.frame(t(val$expr[available_val, ]))
  # Fill any missing columns (should be none)
  pred_val <- predict(model, newdata = X_val, type = "response")
  roc_val <- roc(val$disease, pred_val, quiet = TRUE)
  
  cat(paste("Validation (GSE6613) AUC:", round(auc(roc_val), 3), "\n"))
  
  ci_val <- ci.auc(roc_val)
  cat(paste("Validation 95% CI:", round(ci_val[1], 3), "-",
            round(ci_val[3], 3), "\n"))
  
  threshold <- coords(roc_val, "best", method = "youden")$threshold
  pred_class <- ifelse(pred_val >= threshold, "PD", "Control")
  conf_mat <- table(Predicted = pred_class, Actual = val$disease)
  print(conf_mat)
  
  # ROC plot
  pdf("PD_Biomarker_Project/figures/WGCNA_validation_ROC.pdf", width = 8, height = 8)
  plot(roc_val, main = paste("WGCNA + LASSO Validation (GSE6613)\nAUC =",
                             round(auc(roc_val), 3)),
       col = "blue", lwd = 3)
  abline(0, 1, lty = 2, col = "gray")
  dev.off()
} else {
  cat("⚠️ Too few genes available in validation cohort.\n")
}

cat("\n✅ BLOCK 7 COMPLETE\n")

# ==========================================
# BLOCK 8: FUNCTIONAL ENRICHMENT
# ==========================================

library(clusterProfiler)
library(org.Hs.eg.db)

entrez <- bitr(selected_genes, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
cat("Mapped to Entrez:", nrow(entrez), "/", length(selected_genes), "\n")

# GO Biological Process
ego_bp <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, readable = TRUE)

# KEGG
ekegg <- enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Save and print
if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
  write.csv(as.data.frame(ego_bp),
            "PD_Biomarker_Project/tables/WGCNA_GO_BP.csv", row.names = FALSE)
  cat("\nTop GO-BP terms:\n")
  print(head(as.data.frame(ego_bp)[, c("Description", "p.adjust", "Count")], 10))
} else {
  cat("\nNo significant GO-BP terms\n")
}

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  write.csv(as.data.frame(ekegg),
            "PD_Biomarker_Project/tables/WGCNA_KEGG.csv", row.names = FALSE)
  cat("\nTop KEGG pathways:\n")
  print(head(as.data.frame(ekegg)[, c("Description", "p.adjust", "Count")], 10))
} else {
  cat("\nNo significant KEGG pathways\n")
}

cat("\n✅ BLOCK 8 COMPLETE\n")

# ==========================================
# BLOCK 9: FINAL SUMMARY
# ==========================================

cat("\n================================================================\n")
cat("WGCNA + ML FINAL SUMMARY\n")
cat("================================================================\n\n")

summary_table <- data.frame(
  Analysis = c("Discovery (GSE99039)", "Validation (GSE6613)"),
  AUC = c(round(auc(roc_disc), 3), round(auc(roc_val), 3)),
  CI_Lower = c(round(ci.auc(roc_disc)[1], 3), round(ci_val[1], 3)),
  CI_Upper = c(round(ci.auc(roc_disc)[3], 3), round(ci_val[3], 3))
)

print(summary_table)
write.csv(summary_table,
          "PD_Biomarker_Project/tables/WGCNA_final_summary.csv",
          row.names = FALSE)

cat("\nSelected hub genes:\n")
print(selected_genes)

# Save WGCNA results for comparison with ML approach
saveRDS(list(
  auc_disc  = as.numeric(auc(roc_disc)),
  auc_val   = as.numeric(auc(roc_val)),
  ci_val    = ci_val,
  hub_genes = selected_genes,
  conf_mat  = conf_mat
), "PD_Biomarker_Project/wgcna_results.rds")

cat("\n✅ WGCNA + ML ANALYSIS COMPLETE\n")
