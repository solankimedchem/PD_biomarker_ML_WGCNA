# ============================================================
# UNIVERSAL PD BIOMARKER PIPELINE
# All ML algorithms + cached downloads + auto version detection
# ============================================================

# ============================================================
# PART 0: ENVIRONMENT SETUP
# ============================================================

r_major_minor <- paste0(R.version$major, ".",
                        strsplit(R.version$minor, "\\.")[[1]][1])
user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R",
                      "win-library", r_major_minor)
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

cat("=========================================================\n")
cat("UNIVERSAL PD BIOMARKER PIPELINE\n")
cat("=========================================================\n")
cat("R version:", R.version.string, "\n")
cat("R major.minor:", r_major_minor, "\n")
cat("User library:", user_lib, "\n")
print(.libPaths())
cat("=========================================================\n\n")

# ---- Match Bioconductor version to R version ----
bioc_version_for_r <- function() {
  rv <- as.numeric(paste0(R.version$major, ".",
                          strsplit(R.version$minor, "\\.")[[1]][1]))
  if (rv >= 4.6) return("3.23")
  if (rv >= 4.5) return("3.22")
  if (rv >= 4.4) return("3.19")
  if (rv >= 4.3) return("3.18")
  return("3.17")
}
bioc_ver <- bioc_version_for_r()
cat("Selected Bioconductor version:", bioc_ver, "\n\n")

options(repos = c(
  CRAN     = "https://cloud.r-project.org/",
  BioCsoft = paste0("https://bioconductor.org/packages/", bioc_ver, "/bioc"),
  BioCann  = paste0("https://bioconductor.org/packages/", bioc_ver, "/data/annotation"),
  BioCexp  = paste0("https://bioconductor.org/packages/", bioc_ver, "/data/experiment")
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib, quiet = TRUE)
}

# ---- Install Matrix FIRST ----
cat("Ensuring Matrix is up to date...\n")
tryCatch({
  install.packages("Matrix", lib = user_lib, quiet = TRUE,
                   dependencies = TRUE, type = "binary")
  cat("✓ Matrix installed\n")
}, error = function(e) cat("⚠️ Matrix:", conditionMessage(e), "\n"))

# ---- Universal installer ----
install_pkg <- function(pkg, lib = user_lib) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "\n")
    tryCatch({
      BiocManager::install(pkg, lib = lib, update = FALSE,
                           ask = FALSE, quiet = TRUE)
    }, error = function(e) {
      tryCatch(
        install.packages(pkg, lib = lib, quiet = TRUE),
        error = function(e2) cat("  ✗ Failed:", pkg, "\n")
      )
    })
  }
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("  ✓", pkg, "\n")
    return(TRUE)
  } else {
    cat("  ✗", pkg, "(not installed)\n")
    return(FALSE)
  }
}

cran_pkgs <- c("glmnet", "randomForest", "e1071", "pROC", "caret",
               "ggplot2", "pheatmap", "reshape2", "RColorBrewer",
               "xgboost", "WGCNA", "dynamicTreeCut",
               "MASS", "nnet", "gbm", "kernlab")

bioc_pkgs <- c("GEOquery", "limma", "sva", "clusterProfiler",
               "org.Hs.eg.db", "ReactomePA", "DOSE", "GSVA",
               "AnnotationDbi")

ann_pkgs  <- c("hgu133a.db", "hgu133plus2.db", "illuminaHumanv4.db")

cat("\n📦 CRAN packages:\n")
for (p in cran_pkgs) install_pkg(p)

cat("\n📦 Bioconductor packages:\n")
for (p in bioc_pkgs) install_pkg(p)

cat("\n📦 Annotation packages:\n")
for (p in ann_pkgs) install_pkg(p)

options(stringsAsFactors = FALSE)
set.seed(123)

dir.create("PD_Biomarker_Project", showWarnings = FALSE)
dir.create("PD_Biomarker_Project/figures", showWarnings = FALSE)
dir.create("PD_Biomarker_Project/tables", showWarnings = FALSE)
dir.create("PD_Biomarker_Project/models", showWarnings = FALSE)
dir.create("GEO_cache", showWarnings = FALSE)  # Cache folder for downloads

cat("\n✅ PART 0 COMPLETE\n")

# ============================================================
# PART A: DATA LOADING (WITH CACHING) + MAPPING + PREPROCESSING
# ============================================================

cat("\n=========================================================\n")
cat("PART A: DATA LOADING (CACHED) + MAPPING + PREPROCESSING\n")
cat("=========================================================\n\n")

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(AnnotationDbi)
})

# ---- Cached GEO loader ----
# Uses getGEO() only if the file is not already in the cache folder.
# Skips re-downloading if you have already downloaded the dataset.
load_geo_cached <- function(gse_id) {
  cache_file <- file.path("GEO_cache", paste0(gse_id, "_eset.rds"))
  
  if (file.exists(cache_file)) {
    cat("  Loading from cache:", gse_id, "\n")
    eset <- readRDS(cache_file)
    return(eset)
  }
  
  cat("  Downloading:", gse_id, "(first time only)\n")
  gset <- getGEO(gse_id, GSEMatrix = TRUE)
  eset <- gset[[1]]
  saveRDS(eset, cache_file)
  cat("  ✓ Cached to:", cache_file, "\n")
  return(eset)
}

load_dataset_correct <- function(gse_id) {
  cat("Loading:", gse_id, "...\n")
  eset <- load_geo_cached(gse_id)
  expr <- exprs(eset)
  pheno <- pData(eset)
  ds <- rep(NA_character_, nrow(pheno))
  
  if (gse_id == "GSE6613") {
    for (i in seq_len(nrow(pheno))) {
      ch <- as.character(pheno[i, "characteristics_ch1"])
      if (grepl("Parkinson", ch, ignore.case = TRUE)) ds[i] <- "PD"
      else if (grepl("healthy", ch, ignore.case = TRUE)) ds[i] <- "Control"
    }
    k <- !is.na(ds); expr <- expr[, k]; ds <- ds[k]
  } else if (gse_id == "GSE54536") {
    for (i in seq_len(nrow(pheno))) {
      if ("diagnosis:ch1" %in% colnames(pheno)) {
        d <- as.character(pheno[i, "diagnosis:ch1"])
        if (grepl("Parkinson", d, ignore.case = TRUE)) ds[i] <- "PD"
        else if (grepl("healthy", d, ignore.case = TRUE)) ds[i] <- "Control"
      }
    }
    k <- !is.na(ds); expr <- expr[, k]; ds <- ds[k]
  } else if (gse_id == "GSE99039") {
    for (i in seq_len(nrow(pheno))) {
      if ("disease label:ch1" %in% colnames(pheno)) {
        l <- as.character(pheno[i, "disease label:ch1"])
        if (grepl("IPD", l, ignore.case = TRUE)) ds[i] <- "PD"
        else if (grepl("CONTROL", l, ignore.case = TRUE)) ds[i] <- "Control"
      }
    }
    k <- !is.na(ds); expr <- expr[, k]; ds <- ds[k]
  } else if (gse_id == "GSE7621") {
    for (i in seq_len(nrow(pheno))) {
      ch <- as.character(pheno[i, "characteristics_ch1"])
      if (grepl("Parkinson", ch, ignore.case = TRUE)) ds[i] <- "PD"
      else if (grepl("Old Control", ch, ignore.case = TRUE)) ds[i] <- "Control"
    }
    k <- !is.na(ds); expr <- expr[, k]; ds <- ds[k]
  }
  
  colnames(expr) <- paste(gse_id, seq_len(ncol(expr)), sep = "_")
  cat("  Samples:", ncol(expr), "| PD:", sum(ds == "PD"),
      "| Control:", sum(ds == "Control"), "\n")
  list(expr = expr, pheno = pheno, disease = ds, id = gse_id)
}

dataset_list <- list()
for (g in c("GSE6613", "GSE54536", "GSE99039", "GSE7621")) {
  dataset_list[[g]] <- load_dataset_correct(g)
}

# ---- Gene mapping ----
map_to_symbols_bioc <- function(data_list) {
  gse_id <- data_list$id
  expr <- data_list$expr
  cat("\nMapping", gse_id, "...\n")
  
  ann_pkg <- switch(gse_id,
                    "GSE6613"  = "hgu133a.db",
                    "GSE54536" = "illuminaHumanv4.db",
                    "GSE99039" = "hgu133plus2.db",
                    "GSE7621"  = "hgu133plus2.db")
  
  if (!requireNamespace(ann_pkg, quietly = TRUE)) {
    tryCatch(
      BiocManager::install(ann_pkg, lib = user_lib, update = FALSE,
                           ask = FALSE, quiet = TRUE),
      error = function(e) NULL
    )
  }
  
  if (!requireNamespace(ann_pkg, quietly = TRUE)) {
    cat("  ⚠️ Annotation package not available:", ann_pkg, "\n")
    data_list$genes <- rownames(expr)
    return(data_list)
  }
  
  suppressPackageStartupMessages(
    library(ann_pkg, character.only = TRUE)
  )
  ann_obj <- get(ann_pkg)
  
  probes <- rownames(expr)
  mapped <- tryCatch(
    suppressMessages(
      AnnotationDbi::select(ann_obj, keys = probes,
                            columns = "SYMBOL", keytype = "PROBEID")
    ),
    error = function(e) NULL
  )
  
  if (is.null(mapped) || nrow(mapped) == 0) {
    data_list$genes <- rownames(expr)
    return(data_list)
  }
  
  mapped <- mapped[!is.na(mapped$SYMBOL) & mapped$SYMBOL != "", ]
  mapped <- mapped[!duplicated(mapped$PROBEID), ]
  common <- intersect(rownames(expr), mapped$PROBEID)
  cat("  Probes mapped:", length(common), "/", length(probes), "\n")
  
  if (length(common) < 100) {
    data_list$genes <- rownames(expr)
    return(data_list)
  }
  
  expr_sub <- expr[common, , drop = FALSE]
  sym <- mapped$SYMBOL[match(common, mapped$PROBEID)]
  expr_df <- as.data.frame(expr_sub, stringsAsFactors = FALSE)
  expr_df$SYMBOL <- sym
  num_cols <- setdiff(colnames(expr_df), "SYMBOL")
  expr_agg <- aggregate(expr_df[, num_cols],
                        by = list(SYMBOL = expr_df$SYMBOL), FUN = mean)
  rownames(expr_agg) <- expr_agg$SYMBOL
  expr_agg <- expr_agg[, -1, drop = FALSE]
  
  data_list$expr <- as.matrix(expr_agg)
  data_list$genes <- rownames(expr_agg)
  cat("  ✓ Final:", nrow(expr_agg), "genes\n")
  data_list
}

for (g in names(dataset_list)) {
  dataset_list[[g]] <- map_to_symbols_bioc(dataset_list[[g]])
}

# ---- Preprocess ----
preprocess_data <- function(data_list) {
  expr <- data_list$expr
  if (max(expr, na.rm = TRUE) > 50) {
    expr[expr < 0] <- 0
    expr <- log2(expr + 1)
  }
  keep <- rowMeans(is.na(expr)) < 0.2
  expr <- expr[keep, , drop = FALSE]
  expr <- normalizeBetweenArrays(expr, method = "quantile")
  expr <- expr[rowMeans(expr > 0) > 0.5, , drop = FALSE]
  data_list$expr <- expr
  data_list
}

for (g in names(dataset_list)) {
  dataset_list[[g]] <- preprocess_data(dataset_list[[g]])
  cat(g, ":", nrow(dataset_list[[g]]$expr), "genes ×",
      ncol(dataset_list[[g]]$expr), "samples\n")
}

cat("\n✅ PART A COMPLETE\n")

# ============================================================
# PART B: ML PIPELINE FUNCTIONS — ALL ALGORITHMS
# ============================================================

cat("\n=========================================================\n")
cat("PART B: ML PIPELINE — ALL ALGORITHMS\n")
cat("=========================================================\n\n")

suppressPackageStartupMessages({
  library(randomForest); library(glmnet); library(e1071)
  library(xgboost); library(pROC); library(caret); library(sva); library(limma)
  library(MASS); library(nnet); library(gbm)
})

# ---- Safe xgboost ----
safe_xgb_train <- function(X, y, nrounds = 100, max_depth = 3, eta = 0.1) {
  y_int <- as.integer(y)
  tryCatch({
    xgboost(x = as.matrix(X), y = y_int,
            nrounds = nrounds, max_depth = max_depth,
            learning_rate = eta, objective = "binary:logistic",
            verbose = 0)
  }, error = function(e1) {
    tryCatch({
      xgboost(data = as.matrix(X), label = y_int,
              nrounds = nrounds, max_depth = max_depth,
              eta = eta, objective = "binary:logistic",
              verbose = 0)
    }, error = function(e2) NULL)
  })
}
safe_xgb_predict <- function(model, X) {
  if (is.null(model)) return(rep(0.5, nrow(X)))
  tryCatch(predict(model, as.matrix(X)),
           error = function(e) rep(0.5, nrow(X)))
}

# ---- Train all models on a training fold ----
train_all_models <- function(X_tr, y_tr) {
  y_num <- as.integer(y_tr) - 1L
  models <- list()
  
  # 1. Random Forest
  models$rf <- tryCatch(
    randomForest(x = X_tr, y = y_tr, ntree = 500),
    error = function(e) NULL)
  
  # 2. LASSO (alpha = 1)
  models$lasso <- tryCatch(
    cv.glmnet(as.matrix(X_tr), y_num, family = "binomial",
              alpha = 1, nfolds = min(5, max(2, sum(table(y_num)) - 1))),
    error = function(e) NULL)
  
  # 3. Elastic Net (alpha = 0.5)
  models$enet <- tryCatch(
    cv.glmnet(as.matrix(X_tr), y_num, family = "binomial",
              alpha = 0.5, nfolds = min(5, max(2, sum(table(y_num)) - 1))),
    error = function(e) NULL)
  
  # 4. SVM (RBF)
  models$svm <- tryCatch(
    svm(x = X_tr, y = y_tr, kernel = "radial", probability = TRUE),
    error = function(e) NULL)
  
  # 5. XGBoost
  models$xgb <- safe_xgb_train(X_tr, y_num)
  
  # 6. Naive Bayes
  models$nb <- tryCatch(
    naiveBayes(x = as.data.frame(X_tr), y = y_tr),
    error = function(e) NULL)
  
  # 7. kNN (k = 5)
  models$knn <- tryCatch(
    knn3(x = X_tr, y = y_tr, k = 5),
    error = function(e) NULL)
  
  # 8. LDA
  models$lda <- tryCatch(
    lda(x = X_tr, grouping = y_tr),
    error = function(e) NULL)
  
  # 9. Logistic Regression
  models$lr <- tryCatch(
    glm(y_tr ~ ., data = data.frame(y_tr, X_tr),
        family = "binomial"),
    error = function(e) NULL)
  
  # 10. Gradient Boosting (gbm)
  models$gbm <- tryCatch(
    gbm(y_tr ~ ., data = data.frame(y_tr, X_tr),
        distribution = "bernoulli", n.trees = 100,
        interaction.depth = 3, shrinkage = 0.1, verbose = FALSE),
    error = function(e) NULL)
  
  # 11. Neural Network (nnet)
  models$nnet <- tryCatch(
    nnet(x = X_tr, y = y_tr, size = 5, maxit = 200,
         decay = 0.1, trace = FALSE),
    error = function(e) NULL)
  
  models
}

# ---- Predict from each model ----
predict_all_models <- function(models, X_te) {
  n <- nrow(X_te)
  preds <- data.frame(
    RF      = rep(0.5, n),
    LASSO   = rep(0.5, n),
    ENET    = rep(0.5, n),
    SVM     = rep(0.5, n),
    XGB     = rep(0.5, n),
    NB      = rep(0.5, n),
    KNN     = rep(0.5, n),
    LDA     = rep(0.5, n),
    LR      = rep(0.5, n),
    GBM     = rep(0.5, n),
    NNET    = rep(0.5, n)
  )
  
  if (!is.null(models$rf)) {
    preds$RF <- tryCatch(
      predict(models$rf, X_te, type = "prob")[, "PD"],
      error = function(e) rep(0.5, n))
  }
  if (!is.null(models$lasso)) {
    preds$LASSO <- tryCatch(
      as.vector(predict(models$lasso, as.matrix(X_te), type = "response")),
      error = function(e) rep(0.5, n))
  }
  if (!is.null(models$enet)) {
    preds$ENET <- tryCatch(
      as.vector(predict(models$enet, as.matrix(X_te), type = "response")),
      error = function(e) rep(0.5, n))
  }
  if (!is.null(models$svm)) {
    preds$SVM <- tryCatch(
      attr(predict(models$svm, X_te, probability = TRUE),
           "probabilities")[, "PD"],
      error = function(e) rep(0.5, n))
  }
  preds$XGB <- safe_xgb_predict(models$xgb, X_te)
  
  if (!is.null(models$nb)) {
    preds$NB <- tryCatch(
      predict(models$nb, as.data.frame(X_te), type = "raw")[, "PD"],
      error = function(e) rep(0.5, n))
  }
  if (!is.null(models$knn)) {
    preds$KNN <- tryCatch(
      predict(models$knn, X_te, type = "prob")[, "PD"],
      error = function(e) rep(0.5, n))
  }
  if (!is.null(models$lda)) {
    preds$LDA <- tryCatch(
      predict(models$lda, X_te)$posterior[, "PD"],
      error = function(e) rep(0.5, n))
  }
  if (!is.null(models$lr)) {
    preds$LR <- tryCatch(
      predict(models$lr, newdata = data.frame(X_te), type = "response"),
      error = function(e) rep(0.5, n))
  }
  if (!is.null(models$gbm)) {
    preds$GBM <- tryCatch(
      predict(models$gbm, newdata = data.frame(X_te),
              n.trees = 100, type = "response"),
      error = function(e) rep(0.5, n))
  }
  if (!is.null(models$nnet)) {
    preds$NNET <- tryCatch(
      predict(models$nnet, X_te, type = "raw")[, 1],
      error = function(e) rep(0.5, n))
  }
  
  preds
}

# ---- Fit a locked model on a training set ----
fit_locked_model <- function(train_expr, train_labels, train_batch,
                             n_features = 100, seed = 123) {
  set.seed(seed)
  
  if (length(unique(train_batch)) > 1) {
    mod <- model.matrix(~ train_labels)
    train_expr <- ComBat(dat = train_expr, batch = as.factor(train_batch),
                         mod = mod, par.prior = TRUE)
  }
  
  train_labels <- factor(train_labels, levels = c("Control", "PD"))
  design <- model.matrix(~ 0 + train_labels)
  colnames(design) <- levels(train_labels)
  contrast <- makeContrasts(PD - Control, levels = design)
  fit <- lmFit(train_expr, design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  deg <- topTable(fit2, coef = 1, number = Inf, sort.by = "P")
  
  features <- rownames(deg)[1:min(n_features, nrow(deg))]
  X_tr <- t(train_expr[features, , drop = FALSE])
  train_means <- colMeans(X_tr)
  train_sds <- apply(X_tr, 2, sd)
  train_sds[train_sds == 0] <- 1
  X_tr_scaled <- scale(X_tr, center = train_means, scale = train_sds)
  X_tr_scaled[is.na(X_tr_scaled)] <- 0
  
  y_tr <- factor(train_labels, levels = c("Control", "PD"))
  
  # Train all models on the training fold
  models <- train_all_models(X_tr_scaled, y_tr)
  
  # ---- Ensemble weights via INNER split ----
  set.seed(seed + 1)
  inner_idx <- tryCatch(
    createDataPartition(y_tr, p = 0.7, list = FALSE),
    error = function(e) NULL)
  
  if (!is.null(inner_idx) && length(unique(y_tr[inner_idx])) == 2) {
    X_in <- X_tr_scaled[inner_idx, , drop = FALSE]
    y_in <- y_tr[inner_idx]
    X_ho <- X_tr_scaled[-inner_idx, , drop = FALSE]
    y_ho <- y_tr[-inner_idx]
    
    inner_models <- train_all_models(X_in, y_in)
    inner_preds <- predict_all_models(inner_models, X_ho)
    
    model_names <- colnames(inner_preds)
    aucs <- sapply(model_names, function(nm) {
      tryCatch(as.numeric(auc(roc(y_ho, inner_preds[[nm]], quiet = TRUE))),
               error = function(e) 0.5)
    })
    aucs[is.na(aucs) | aucs < 0.5] <- 0.5
    weights <- aucs / sum(aucs)
  } else {
    model_names <- c("RF", "LASSO", "ENET", "SVM", "XGB",
                     "NB", "KNN", "LDA", "LR", "GBM", "NNET")
    weights <- setNames(rep(1/length(model_names), length(model_names)),
                        model_names)
    aucs <- setNames(rep(0.5, length(model_names)), model_names)
  }
  
  list(features = features, train_means = train_means, train_sds = train_sds,
       models = models, weights = weights, train_aucs = aucs)
}

predict_locked <- function(locked, new_expr, new_labels = NULL) {
  feats <- locked$features
  missing_feats <- setdiff(feats, rownames(new_expr))
  if (length(missing_feats) > 0) {
    new_expr <- rbind(new_expr,
                      matrix(0, nrow = length(missing_feats), ncol = ncol(new_expr),
                             dimnames = list(missing_feats, colnames(new_expr))))
  }
  X_new <- t(new_expr[feats, , drop = FALSE])
  X_new <- scale(X_new, center = locked$train_means, scale = locked$train_sds)
  X_new[is.na(X_new)] <- 0
  
  preds <- predict_all_models(locked$models, X_new)
  
  # Weighted ensemble
  w <- locked$weights[colnames(preds)]
  w[is.na(w)] <- 0
  ens <- as.vector(as.matrix(preds) %*% w)
  
  out <- cbind(preds, Ensemble = ens)
  if (!is.null(new_labels)) out$TrueLabel <- new_labels
  out
}

compute_metrics <- function(y_true, y_prob, threshold = NULL) {
  y_true <- factor(y_true, levels = c("Control", "PD"))
  roc_obj <- roc(y_true, y_prob, quiet = TRUE)
  if (is.null(threshold))
    threshold <- as.numeric(coords(roc_obj, "best", method = "youden")$threshold)
  pred <- factor(ifelse(y_prob >= threshold, "PD", "Control"),
                 levels = c("Control", "PD"))
  cm <- table(Predicted = pred, Actual = y_true)
  tp <- if ("PD" %in% rownames(cm) && "PD" %in% colnames(cm)) cm["PD","PD"] else 0
  tn <- if ("Control" %in% rownames(cm) && "Control" %in% colnames(cm)) cm["Control","Control"] else 0
  fp <- if ("PD" %in% rownames(cm) && "Control" %in% colnames(cm)) cm["PD","Control"] else 0
  fn <- if ("Control" %in% rownames(cm) && "PD" %in% colnames(cm)) cm["Control","PD"] else 0
  list(AUC = as.numeric(auc(roc_obj)), Threshold = threshold,
       Accuracy = (tp + tn) / sum(cm),
       Sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA,
       Specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA,
       TP = tp, TN = tn, FP = fp, FN = fn, ROC = roc_obj)
}

cat("\n✅ PART B COMPLETE\n")

# ============================================================
# PART C: DISCOVERY ON GSE99039 + REPEATED 5-FOLD CV
# ============================================================

cat("\n=========================================================\n")
cat("PART C: DISCOVERY (GSE99039) + REPEATED 5-FOLD CV\n")
cat("=========================================================\n\n")

disc <- dataset_list[["GSE99039"]]
disc_expr <- disc$expr
disc_lab  <- disc$disease

cat("Discovery:", ncol(disc_expr), "samples |",
    sum(disc_lab == "PD"), "PD |", sum(disc_lab == "Control"), "Control\n\n")

# Track per-model AUCs across folds
model_names <- c("RF", "LASSO", "ENET", "SVM", "XGB",
                 "NB", "KNN", "LDA", "LR", "GBM", "NNET", "Ensemble")
per_model_aucs <- setNames(
  lapply(model_names, function(x) numeric(0)),
  model_names)
ens_aucs <- c(); ens_acc <- c(); ens_sens <- c(); ens_spec <- c()

for (rep in 1:10) {
  set.seed(200 + rep)
  folds <- createFolds(disc_lab, k = 5, list = TRUE, returnTrain = FALSE)
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(disc_lab), test_idx)
    tr_expr <- disc_expr[, train_idx, drop = FALSE]
    tr_lab  <- disc_lab[train_idx]
    te_expr <- disc_expr[, test_idx, drop = FALSE]
    te_lab  <- disc_lab[test_idx]
    
    locked <- fit_locked_model(tr_expr, tr_lab,
                               rep("GSE99039", length(train_idx)),
                               n_features = 100, seed = 5000 * rep + i)
    preds <- predict_locked(locked, te_expr, te_lab)
    
    for (nm in model_names) {
      if (nm %in% colnames(preds)) {
        a <- tryCatch(
          as.numeric(auc(roc(preds$TrueLabel, preds[[nm]], quiet = TRUE))),
          error = function(e) NA)
        if (!is.na(a)) per_model_aucs[[nm]] <- c(per_model_aucs[[nm]], a)
      }
    }
    
    m <- compute_metrics(preds$TrueLabel, preds$Ensemble)
    ens_aucs <- c(ens_aucs, m$AUC)
    ens_acc  <- c(ens_acc,  m$Accuracy)
    ens_sens <- c(ens_sens, m$Sensitivity)
    ens_spec <- c(ens_spec, m$Specificity)
  }
  cat("Repeat", rep, "done\n")
}

cat("\n--- REPEATED 5-FOLD CV (GSE99039) ---\n")
cat("Ensemble AUC:", round(mean(ens_aucs), 3),
    "±", round(sd(ens_aucs), 3), "\n")
cat("Ensemble Accuracy:", round(mean(ens_acc, na.rm = TRUE), 3),
    "±", round(sd(ens_acc, na.rm = TRUE), 3), "\n")
cat("Ensemble Sensitivity:", round(mean(ens_sens, na.rm = TRUE), 3),
    "±", round(sd(ens_sens, na.rm = TRUE), 3), "\n")
cat("Ensemble Specificity:", round(mean(ens_spec, na.rm = TRUE), 3),
    "±", round(sd(ens_spec, na.rm = TRUE), 3), "\n")

cat("\n--- Per-model AUC (mean ± SD) ---\n")
per_model_summary <- data.frame(
  Model = model_names,
  Mean_AUC = sapply(model_names, function(nm) round(mean(per_model_aucs[[nm]], na.rm = TRUE), 3)),
  SD_AUC = sapply(model_names, function(nm) round(sd(per_model_aucs[[nm]], na.rm = TRUE), 3)),
  N_folds = sapply(model_names, function(nm) length(per_model_aucs[[nm]]))
)
per_model_summary <- per_model_summary[order(-per_model_summary$Mean_AUC), ]
print(per_model_summary)

write.csv(per_model_summary,
          "PD_Biomarker_Project/tables/per_model_aucs.csv",
          row.names = FALSE)

saveRDS(list(
  ensemble = list(aucs = ens_aucs, acc = ens_acc,
                  sens = ens_sens, spec = ens_spec),
  per_model = per_model_aucs
), "PD_Biomarker_Project/repeated_cv_results.rds")

cat("\n✅ PART C COMPLETE\n")

# ============================================================
# PART D: FINAL MODEL + VALIDATION + CROSS-TISSUE
# ============================================================

cat("\n=========================================================\n")
cat("PART D: FINAL MODEL + VALIDATION\n")
cat("=========================================================\n\n")

final_locked <- fit_locked_model(disc_expr, disc_lab,
                                 rep("GSE99039", ncol(disc_expr)),
                                 n_features = 100, seed = 123)
saveRDS(final_locked, "PD_Biomarker_Project/final_locked_model.rds")
cat("Panel size:", length(final_locked$features), "genes\n")
cat("\nEnsemble weights (based on inner CV):\n")
print(round(final_locked$weights, 3))

# ---- Validation on GSE6613 ----
cat("\n--- Validation: GSE6613 (blood) ---\n")
val <- dataset_list[["GSE6613"]]
preds_val <- predict_locked(final_locked, val$expr, val$disease)
m_val <- compute_metrics(preds_val$TrueLabel, preds_val$Ensemble)
cat("Ensemble AUC:", round(m_val$AUC, 3),
    "| Acc:", round(m_val$Accuracy, 3),
    "| Sens:", round(m_val$Sensitivity, 3),
    "| Spec:", round(m_val$Specificity, 3), "\n")
print(table(Predicted = ifelse(preds_val$Ensemble >= m_val$Threshold, "PD", "Control"),
            Actual = preds_val$TrueLabel))
ci_val <- ci.auc(roc(preds_val$TrueLabel, preds_val$Ensemble, quiet = TRUE))
cat("95% CI AUC:", round(ci_val[1], 3), "-", round(ci_val[3], 3), "\n")

# Per-model AUCs on validation
cat("\nPer-model validation AUCs:\n")
val_aucs <- sapply(model_names[model_names != "Ensemble"], function(nm) {
  if (nm %in% colnames(preds_val))
    round(as.numeric(auc(roc(preds_val$TrueLabel, preds_val[[nm]], quiet = TRUE))), 3)
  else NA
})
print(sort(val_aucs, decreasing = TRUE, na.last = TRUE))

# ---- Cross-tissue: GSE7621 (brain) ----
cat("\n--- Cross-tissue: GSE7621 (brain, comparison only) ---\n")
brain <- dataset_list[["GSE7621"]]
preds_brain <- predict_locked(final_locked, brain$expr, brain$disease)
m_brain <- compute_metrics(preds_brain$TrueLabel, preds_brain$Ensemble)
cat("Ensemble AUC:", round(m_brain$AUC, 3),
    "| Acc:", round(m_brain$Accuracy, 3),
    "| Sens:", round(m_brain$Sensitivity, 3),
    "| Spec:", round(m_brain$Specificity, 3), "\n")
ci_brain <- ci.auc(roc(preds_brain$TrueLabel, preds_brain$Ensemble, quiet = TRUE))
cat("95% CI AUC:", round(ci_brain[1], 3), "-", round(ci_brain[3], 3), "\n")

# ---- Direction concordance ----
cat("\n--- Direction concordance (blood vs brain) ---\n")
design_b <- model.matrix(~ 0 + disc_lab)
colnames(design_b) <- c("Control", "PD")
fit_b <- eBayes(contrasts.fit(lmFit(disc_expr, design_b),
                              makeContrasts(PD - Control, levels = design_b)))
deg_blood <- topTable(fit_b, coef = 1, number = Inf)

brain_lab <- factor(brain$disease, levels = c("Control", "PD"))
design_br <- model.matrix(~ 0 + brain_lab)
colnames(design_br) <- c("Control", "PD")
fit_br <- eBayes(contrasts.fit(lmFit(brain$expr, design_br),
                               makeContrasts(PD - Control, levels = design_br)))
deg_brain <- topTable(fit_br, coef = 1, number = Inf)

panel <- final_locked$features
common_panel <- intersect(intersect(panel, rownames(deg_blood)),
                          rownames(deg_brain))
if (length(common_panel) > 0) {
  b  <- deg_blood[common_panel, "logFC"]
  br <- deg_brain[common_panel, "logFC"]
  conc <- sign(b) == sign(br)
  cat("Panel genes in both:", length(common_panel), "\n")
  cat("Concordance:", sum(conc), "/", length(conc), "=",
      round(mean(conc) * 100, 1), "%\n")
  write.csv(data.frame(Gene = common_panel, logFC_Blood = b,
                       logFC_Brain = br, Concordant = conc),
            "PD_Biomarker_Project/tables/cross_tissue_concordance.csv",
            row.names = FALSE)
}

# ---- GSE54536 small cohort ----
cat("\n--- GSE54536 (n=10, exploratory) ---\n")
gse54536 <- dataset_list[["GSE54536"]]
preds_54536 <- predict_locked(final_locked, gse54536$expr, gse54536$disease)
m_54536 <- compute_metrics(preds_54536$TrueLabel, preds_54536$Ensemble)
cat("Ensemble AUC:", round(m_54536$AUC, 3),
    "| Acc:", round(m_54536$Accuracy, 3), "\n")

cat("\n✅ PART D COMPLETE\n")

# ============================================================
# PART E: ENRICHMENT
# ============================================================

cat("\n=========================================================\n")
cat("PART E: FUNCTIONAL ENRICHMENT\n")
cat("=========================================================\n\n")

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(ReactomePA); library(DOSE)
})

panel <- final_locked$features
entrez <- tryCatch(bitr(panel, fromType = "SYMBOL", toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db),
                   error = function(e) NULL)

if (!is.null(entrez) && nrow(entrez) > 0) {
  cat("Mapped to Entrez:", nrow(entrez), "/", length(panel), "\n")
  
  for (ont in c("BP", "MF", "CC")) {
    res <- tryCatch(
      enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
               ont = ont, pAdjustMethod = "BH",
               pvalueCutoff = 0.05, readable = TRUE),
      error = function(e) NULL)
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
      write.csv(as.data.frame(res),
                paste0("PD_Biomarker_Project/tables/enrichment_GO_", ont, ".csv"),
                row.names = FALSE)
      cat("\nTop GO-", ont, " terms:\n", sep = "")
      print(head(as.data.frame(res)[, c("Description", "p.adjust", "Count")], 10))
    } else {
      cat("No significant GO-", ont, " terms\n", sep = "")
    }
  }
  
  ekegg <- tryCatch(enrichKEGG(gene = entrez$ENTREZID, organism = "hsa",
                               pAdjustMethod = "BH", pvalueCutoff = 0.05),
                    error = function(e) NULL)
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    write.csv(as.data.frame(ekegg),
              "PD_Biomarker_Project/tables/enrichment_KEGG.csv",
              row.names = FALSE)
    cat("\nTop KEGG pathways:\n")
    print(head(as.data.frame(ekegg)[, c("Description", "p.adjust", "Count")], 10))
  } else cat("No significant KEGG pathways\n")
}

cat("\n✅ PART E COMPLETE\n")

# ============================================================
# PART F: FINAL SUMMARY
# ============================================================

cat("\n=========================================================\n")
cat("FINAL SUMMARY\n")
cat("=========================================================\n\n")

rep_cv <- readRDS("PD_Biomarker_Project/repeated_cv_results.rds")

summary_table <- data.frame(
  Analysis = c(
    "Repeated 5-fold CV ensemble (GSE99039)",
    "Independent validation ensemble (GSE6613, blood)",
    "Cross-tissue ensemble (GSE7621, brain)",
    "Small cohort ensemble (GSE54536, n=10)"
  ),
  AUC = c(round(mean(rep_cv$ensemble$aucs), 3), round(m_val$AUC, 3),
          round(m_brain$AUC, 3), round(m_54536$AUC, 3)),
  AUC_SD = c(round(sd(rep_cv$ensemble$aucs), 3), NA, NA, NA),
  Accuracy = c(round(mean(rep_cv$ensemble$acc, na.rm = TRUE), 3),
               round(m_val$Accuracy, 3),
               round(m_brain$Accuracy, 3),
               round(m_54536$Accuracy, 3)),
  Sensitivity = c(round(mean(rep_cv$ensemble$sens, na.rm = TRUE), 3),
                  round(m_val$Sensitivity, 3),
                  round(m_brain$Sensitivity, 3),
                  round(m_54536$Sensitivity, 3)),
  Specificity = c(round(mean(rep_cv$ensemble$spec, na.rm = TRUE), 3),
                  round(m_val$Specificity, 3),
                  round(m_brain$Specificity, 3),
                  round(m_54536$Specificity, 3))
)

print(summary_table)
write.csv(summary_table,
          "PD_Biomarker_Project/tables/final_performance_summary.csv",
          row.names = FALSE)

cat("\nTop 20 panel genes:\n")
print(head(final_locked$features, 20))
writeLines(final_locked$features,
           "PD_Biomarker_Project/tables/final_panel_genes.txt")

cat("\n✅ ALL COMPLETE\n")
cat("=================================================================\n")