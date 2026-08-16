options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(limma)
  library(glmnet)
  library(e1071)
  library(randomForestSRC)
  library(gbm)
  library(MASS)
  library(mboost)
  library(xgboost)
  library(pROC)
})

set.seed(20260809)

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root_dir, "11.ML_reanalysis_20260809")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw_files <- list.files(
  root_dir,
  pattern = "^GSE(73661|75214|87466|107499|47908|13367)\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
raw_paths <- setNames(raw_files, tools::file_path_sans_ext(basename(raw_files)))
required_cohorts <- c("GSE73661", "GSE75214", "GSE87466", "GSE107499", "GSE47908", "GSE13367")
stopifnot(all(required_cohorts %in% names(raw_paths)))

candidate_file <- file.path(root_dir, "11.ML", "interGenes.txt")
candidate_genes <- unique(trimws(readLines(candidate_file, warn = FALSE)))
candidate_genes <- candidate_genes[nzchar(candidate_genes)]

normalize_header_file <- list.files(root_dir, pattern = "^normalize1\\.txt$", recursive = TRUE, full.names = TRUE)
stopifnot(length(normalize_header_file) == 1L)
header <- strsplit(readLines(normalize_header_file, n = 1L, warn = FALSE), "\\t", fixed = FALSE)[[1]][-1]
discovery_manifest <- data.frame(
  sample = sub("_(con|treat)$", "", header),
  label = ifelse(grepl("_treat$", header), 1L, 0L),
  stringsAsFactors = FALSE
)

valid_gene_names <- function(x) {
  !is.na(x) & nzchar(x) & x != "---" & grepl("^[A-Za-z0-9._-]+$", x)
}

gene_sets <- lapply(raw_paths, function(path) {
  first_col <- data.table::fread(path, select = 1L, data.table = FALSE, showProgress = FALSE)[[1]]
  unique(as.character(first_col[valid_gene_names(first_col)]))
})
common_genes <- Reduce(intersect, gene_sets)
candidate_genes <- intersect(candidate_genes, common_genes)
if (length(candidate_genes) < 5L) stop("Fewer than five candidate genes are shared by all six cohorts.")

sample_membership <- lapply(raw_paths[names(raw_paths) %in% required_cohorts[1:4]], function(path) {
  names(data.table::fread(path, nrows = 0L, data.table = FALSE, showProgress = FALSE))[-1]
})
discovery_manifest$cohort <- vapply(discovery_manifest$sample, function(s) {
  hit <- names(sample_membership)[vapply(sample_membership, function(ids) s %in% ids, logical(1))]
  if (length(hit) != 1L) stop("Discovery sample mapping failed for ", s)
  hit
}, character(1))
discovery_manifest$set <- "discovery"
discovery_manifest$subgroup <- ifelse(discovery_manifest$label == 1L, "UC", "Control")

gse47908_manifest <- data.frame(
  sample = c(sprintf("GSM%d", 1162227:1162241), sprintf("GSM%d", 1162248:1162286)),
  label = c(rep(0L, 15L), rep(1L, 39L)),
  cohort = "GSE47908",
  set = "external",
  subgroup = c(rep("Control", 15L), rep("UC_non_dysplasia", 39L)),
  stringsAsFactors = FALSE
)

gse13367_samples <- sprintf("GSM%d", 337490:337517)
gse13367_samples <- setdiff(gse13367_samples, "GSM337493")
gse13367_controls <- paste0("GSM", c(337492, 337498, 337501, 337502, 337504, 337505, 337509, 337511, 337512, 337516))
gse13367_inflamed <- paste0("GSM", c(337490, 337494, 337495, 337496, 337500, 337503, 337513, 337515))
gse13367_manifest <- data.frame(
  sample = gse13367_samples,
  label = as.integer(!gse13367_samples %in% gse13367_controls),
  cohort = "GSE13367",
  set = "external",
  subgroup = ifelse(gse13367_samples %in% gse13367_controls, "Control",
                    ifelse(gse13367_samples %in% gse13367_inflamed, "UC_inflamed", "UC_non_inflamed")),
  stringsAsFactors = FALSE
)

manifest <- rbind(discovery_manifest, gse47908_manifest, gse13367_manifest)
manifest <- manifest[, c("sample", "cohort", "set", "label", "subgroup")]
write.csv(manifest, file.path(out_dir, "sample_manifest.csv"), row.names = FALSE)

read_rank_candidate_matrix <- function(path, cohort, manifest, common_genes, candidate_genes) {
  tab <- data.table::fread(path, data.table = FALSE, check.names = FALSE, showProgress = FALSE)
  genes <- as.character(tab[[1]])
  keep <- valid_gene_names(genes)
  genes <- genes[keep]
  mat <- data.matrix(tab[keep, -1, drop = FALSE])
  rownames(mat) <- genes
  if (anyDuplicated(rownames(mat))) mat <- limma::avereps(mat)
  samples <- manifest$sample[manifest$cohort == cohort]
  missing_samples <- setdiff(samples, colnames(mat))
  if (length(missing_samples)) stop("Missing samples in ", cohort, ": ", paste(missing_samples, collapse = ","))
  mat <- mat[common_genes, samples, drop = FALSE]
  ranked <- apply(mat, 2L, function(x) {
    r <- rank(x, ties.method = "average", na.last = "keep")
    qnorm((r - 0.5) / sum(!is.na(r)))
  })
  rownames(ranked) <- common_genes
  ranked[candidate_genes, , drop = FALSE]
}

message("Reading six cohorts and applying per-sample rank-normal transformation...")
ranked_by_cohort <- lapply(required_cohorts, function(cohort) {
  read_rank_candidate_matrix(raw_paths[[cohort]], cohort, manifest, common_genes, candidate_genes)
})
names(ranked_by_cohort) <- required_cohorts
expression_matrix <- do.call(cbind, ranked_by_cohort)
expression_matrix <- expression_matrix[, manifest$sample, drop = FALSE]
X <- t(expression_matrix)
stopifnot(identical(rownames(X), manifest$sample), !anyNA(X))

write.csv(
  data.frame(sample = rownames(X), manifest[match(rownames(X), manifest$sample), -1], X, check.names = FALSE),
  file.path(out_dir, "rank_normalized_candidate_expression.csv"),
  row.names = FALSE
)

make_foldid <- function(y, k = 5L, seed = 1L) {
  set.seed(seed)
  foldid <- integer(length(y))
  for (lev in sort(unique(y))) {
    idx <- which(y == lev)
    foldid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  foldid
}

scale_fit <- function(x) {
  center <- colMeans(x)
  scale <- apply(x, 2L, sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(center = center, scale = scale)
}

scale_apply <- function(x, scaler) {
  z <- sweep(x, 2L, scaler$center, "-")
  sweep(z, 2L, scaler$scale, "/")
}

univariate_order <- function(x, y) {
  p <- apply(x, 2L, function(v) tryCatch(t.test(v[y == 1], v[y == 0])$p.value, error = function(e) 1))
  names(sort(p, decreasing = FALSE))
}

select_features <- function(selector, x, y, seed) {
  if (selector == "all") return(colnames(x))
  if (grepl("^univ_top", selector)) {
    n <- as.integer(sub("univ_top", "", selector))
    return(head(univariate_order(x, y), n))
  }
  if (grepl("^rf_top", selector)) {
    n <- as.integer(sub("rf_top", "", selector))
    set.seed(seed)
    d <- data.frame(Type = factor(y), x, check.names = FALSE)
    fit <- randomForestSRC::rfsrc(Type ~ ., data = d, ntree = 500, nodesize = 5, importance = "permute")
    imp <- fit$importance
    if (is.matrix(imp)) imp <- rowMeans(abs(imp), na.rm = TRUE)
    imp <- sort(setNames(as.numeric(imp), names(imp)), decreasing = TRUE)
    return(head(names(imp), n))
  }
  alpha <- switch(selector, lasso = 1, enet_025 = 0.25, enet_050 = 0.50, enet_075 = 0.75)
  foldid <- make_foldid(y, 5L, seed)
  cvfit <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = alpha, foldid = foldid,
                             type.measure = "auc", standardize = FALSE)
  get_nonzero <- function(s) {
    b <- as.matrix(coef(cvfit, s = s))
    setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
  }
  vars <- get_nonzero("lambda.1se")
  if (length(vars) < 2L) vars <- get_nonzero("lambda.min")
  if (length(vars) < 2L) vars <- unique(c(vars, head(univariate_order(x, y), 2L)))
  vars
}

fit_classifier <- function(classifier, x, y, seed) {
  set.seed(seed)
  if (grepl("^glmnet_", classifier)) {
    alpha <- as.numeric(sub("glmnet_", "", classifier))
    cvfit <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = alpha,
                               foldid = make_foldid(y, 5L, seed), type.measure = "auc", standardize = FALSE)
    fit <- glmnet::glmnet(x, y, family = "binomial", alpha = alpha,
                          lambda = cvfit$lambda.min, standardize = FALSE)
    b <- as.matrix(coef(fit))
    vars <- setdiff(rownames(b)[b[, 1] != 0], "(Intercept)")
    return(list(kind = "glmnet", fit = fit, vars = vars))
  }
  d <- data.frame(Type = factor(y), x, check.names = FALSE)
  if (classifier == "svm_linear") {
    fit <- e1071::svm(Type ~ ., data = d, kernel = "linear", probability = TRUE, scale = FALSE)
    return(list(kind = classifier, fit = fit, vars = colnames(x)))
  }
  if (classifier == "svm_radial") {
    fit <- e1071::svm(Type ~ ., data = d, kernel = "radial", probability = TRUE, scale = FALSE)
    return(list(kind = classifier, fit = fit, vars = colnames(x)))
  }
  if (classifier == "rf") {
    fit <- randomForestSRC::rfsrc(Type ~ ., data = d, ntree = 750, nodesize = 5, importance = "permute")
    return(list(kind = classifier, fit = fit, vars = colnames(x)))
  }
  if (classifier == "gbm") {
    gd <- data.frame(Type = y, x, check.names = FALSE)
    fit <- gbm::gbm(Type ~ ., data = gd, distribution = "bernoulli", n.trees = 1000,
                    interaction.depth = 2, n.minobsinnode = 10, shrinkage = 0.02,
                    bag.fraction = 0.8, cv.folds = 5, verbose = FALSE)
    best <- suppressWarnings(gbm::gbm.perf(fit, method = "cv", plot.it = FALSE))
    if (!is.finite(best) || best < 1) best <- 1000
    return(list(kind = classifier, fit = fit, vars = colnames(x), best = best))
  }
  if (classifier == "lda") {
    fit <- MASS::lda(Type ~ ., data = d)
    return(list(kind = classifier, fit = fit, vars = colnames(x)))
  }
  if (classifier == "naive_bayes") {
    fit <- e1071::naiveBayes(Type ~ ., data = d)
    return(list(kind = classifier, fit = fit, vars = colnames(x)))
  }
  if (classifier == "glmboost") {
    fit <- mboost::glmboost(Type ~ ., data = d,
                            family = mboost::Binomial(), control = mboost::boost_control(mstop = 200, trace = FALSE))
    cvfit <- mboost::cvrisk(fit, folds = mboost::cv(stats::model.weights(fit), type = "kfold", B = 5))
    mboost::mstop(fit) <- mboost::mstop(cvfit)
    b <- coef(fit)
    vars <- names(b)[abs(b) > 0]
    return(list(kind = classifier, fit = fit, vars = vars))
  }
  if (classifier == "xgboost") {
    dx <- xgboost::xgb.DMatrix(x, label = y)
    params <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 2,
                   eta = 0.05, subsample = 0.8, colsample_bytree = 0.8, nthread = 1)
    cvfit <- xgboost::xgb.cv(params = params, data = dx, nrounds = 250, nfold = 5,
                            stratified = TRUE, early_stopping_rounds = 20, verbose = 0)
    best <- cvfit$best_iteration
    if (is.null(best) || !is.finite(best)) best <- 100
    fit <- xgboost::xgb.train(params = params, data = dx, nrounds = best, verbose = 0)
    return(list(kind = classifier, fit = fit, vars = colnames(x)))
  }
  stop("Unknown classifier: ", classifier)
}

predict_classifier <- function(object, x) {
  kind <- object$kind
  if (kind == "glmnet") return(as.numeric(predict(object$fit, newx = x, type = "response")))
  nd <- data.frame(x, check.names = FALSE)
  if (kind %in% c("svm_linear", "svm_radial")) {
    pr <- predict(object$fit, nd, probability = TRUE)
    probs <- attr(pr, "probabilities")
    return(as.numeric(probs[, if ("1" %in% colnames(probs)) "1" else ncol(probs)]))
  }
  if (kind == "rf") {
    probs <- predict(object$fit, newdata = nd)$predicted
    return(as.numeric(probs[, if ("1" %in% colnames(probs)) "1" else ncol(probs)]))
  }
  if (kind == "gbm") return(as.numeric(predict(object$fit, nd, n.trees = object$best, type = "response")))
  if (kind == "lda") {
    probs <- predict(object$fit, nd)$posterior
    return(as.numeric(probs[, if ("1" %in% colnames(probs)) "1" else ncol(probs)]))
  }
  if (kind == "naive_bayes") {
    probs <- predict(object$fit, nd, type = "raw")
    return(as.numeric(probs[, if ("1" %in% colnames(probs)) "1" else ncol(probs)]))
  }
  if (kind == "glmboost") return(as.numeric(predict(object$fit, newdata = nd, type = "response")))
  if (kind == "xgboost") return(as.numeric(predict(object$fit, xgboost::xgb.DMatrix(x))))
  stop("Unknown fitted classifier: ", kind)
}

auc_fixed <- function(y, score) {
  if (length(unique(y)) < 2L || any(!is.finite(score))) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)))
}

selectors <- c("all", "lasso", "enet_025", "enet_050", "enet_075",
               "univ_top5", "univ_top10", "rf_top5", "rf_top10")
classifiers <- c("glmnet_0", "glmnet_0.25", "glmnet_0.5", "glmnet_0.75", "glmnet_1",
                 "svm_linear", "svm_radial", "rf", "gbm", "lda", "naive_bayes", "glmboost", "xgboost")
workflows <- expand.grid(selector = selectors, classifier = classifiers, stringsAsFactors = FALSE)
workflows$workflow <- paste(workflows$selector, workflows$classifier, sep = " + ")
stopifnot(nrow(workflows) == 117L)

discovery_idx <- which(manifest$set == "discovery")
discovery_cohorts <- required_cohorts[1:4]
auc_rows <- list()
feature_rows <- list()
selector_rows <- list()
row_i <- 0L

message("Running 117 workflows in four leave-one-cohort-out folds...")
for (fold_i in seq_along(discovery_cohorts)) {
  heldout <- discovery_cohorts[fold_i]
  train_idx <- discovery_idx[manifest$cohort[discovery_idx] != heldout]
  test_idx <- discovery_idx[manifest$cohort[discovery_idx] == heldout]
  scaler <- scale_fit(X[train_idx, , drop = FALSE])
  x_train <- scale_apply(X[train_idx, , drop = FALSE], scaler)
  x_test <- scale_apply(X[test_idx, , drop = FALSE], scaler)
  y_train <- manifest$label[train_idx]
  y_test <- manifest$label[test_idx]
  selected_by_method <- list()
  for (selector_i in seq_along(selectors)) {
    selector <- selectors[selector_i]
    seed <- 20260809L + fold_i * 1000L + selector_i * 10L
    selected_by_method[[selector]] <- tryCatch(
      select_features(selector, x_train, y_train, seed),
      error = function(e) {
        message("Selector failed [", heldout, ", ", selector, "]: ", conditionMessage(e))
        character()
      }
    )
    selector_rows[[length(selector_rows) + 1L]] <- data.frame(
      heldout_cohort = heldout,
      selector = selector,
      genes = paste(selected_by_method[[selector]], collapse = ";"),
      n_genes = length(selected_by_method[[selector]]),
      stringsAsFactors = FALSE
    )
  }
  for (workflow_i in seq_len(nrow(workflows))) {
    selector <- workflows$selector[workflow_i]
    classifier <- workflows$classifier[workflow_i]
    vars <- selected_by_method[[selector]]
    seed <- 20260809L + fold_i * 10000L + workflow_i
    result <- tryCatch({
      if (length(vars) < 1L) stop("No selected variables")
      fit <- fit_classifier(classifier, x_train[, vars, drop = FALSE], y_train, seed)
      pred <- predict_classifier(fit, x_test[, vars, drop = FALSE])
      list(auc = auc_fixed(y_test, pred), vars = fit$vars, error = "")
    }, error = function(e) list(auc = NA_real_, vars = character(), error = conditionMessage(e)))
    row_i <- row_i + 1L
    auc_rows[[row_i]] <- data.frame(
      workflow = workflows$workflow[workflow_i], selector = selector, classifier = classifier,
      heldout_cohort = heldout, auc = result$auc, n_selector_genes = length(vars),
      n_model_genes = length(result$vars), error = result$error, stringsAsFactors = FALSE
    )
    feature_rows[[row_i]] <- data.frame(
      workflow = workflows$workflow[workflow_i], heldout_cohort = heldout,
      selector_genes = paste(vars, collapse = ";"), model_genes = paste(result$vars, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  message("Completed held-out cohort: ", heldout)
}

auc_long <- do.call(rbind, auc_rows)
feature_long <- do.call(rbind, feature_rows)
selector_long <- do.call(rbind, selector_rows)
write.csv(auc_long, file.path(out_dir, "workflow_auc_by_outer_cohort.csv"), row.names = FALSE)
write.csv(feature_long, file.path(out_dir, "workflow_features_by_outer_cohort.csv"), row.names = FALSE)
write.csv(selector_long, file.path(out_dir, "selector_features_by_outer_cohort.csv"), row.names = FALSE)

summary_rows <- lapply(split(auc_long, auc_long$workflow), function(d) {
  ok <- is.finite(d$auc)
  data.frame(
    workflow = d$workflow[1], selector = d$selector[1], classifier = d$classifier[1],
    mean_auc = if (any(ok)) mean(d$auc[ok]) else NA_real_,
    median_auc = if (any(ok)) median(d$auc[ok]) else NA_real_,
    min_auc = if (any(ok)) min(d$auc[ok]) else NA_real_,
    sd_auc = if (sum(ok) > 1L) sd(d$auc[ok]) else NA_real_,
    successful_cohorts = sum(ok),
    stringsAsFactors = FALSE
  )
})
workflow_summary <- do.call(rbind, summary_rows)
workflow_summary <- workflow_summary[order(-workflow_summary$successful_cohorts, -workflow_summary$mean_auc,
                                           -workflow_summary$min_auc, workflow_summary$sd_auc), ]
rownames(workflow_summary) <- NULL
write.csv(workflow_summary, file.path(out_dir, "workflow_summary.csv"), row.names = FALSE)

sparse_selector_long <- selector_long[selector_long$selector != "all", ]
gene_frequency <- setNames(numeric(length(candidate_genes)), candidate_genes)
for (gene_string in sparse_selector_long$genes) {
  genes <- strsplit(gene_string, ";", fixed = TRUE)[[1]]
  genes <- genes[nzchar(genes)]
  gene_frequency[intersect(names(gene_frequency), genes)] <- gene_frequency[intersect(names(gene_frequency), genes)] + 1
}
gene_frequency <- gene_frequency / nrow(sparse_selector_long)

direction_table <- do.call(rbind, lapply(candidate_genes, function(gene) {
  diffs <- vapply(discovery_cohorts, function(cohort) {
    idx <- discovery_idx[manifest$cohort[discovery_idx] == cohort]
    mean(X[idx, gene][manifest$label[idx] == 1]) - mean(X[idx, gene][manifest$label[idx] == 0])
  }, numeric(1))
  data.frame(gene = gene, t(diffs), mean_abs_difference = mean(abs(diffs)),
             consistent_discovery_direction = length(unique(sign(diffs))) == 1L,
             discovery_direction = ifelse(mean(diffs) > 0, "UC_up", "UC_down"),
             stringsAsFactors = FALSE, check.names = FALSE)
}))
colnames(direction_table)[2:5] <- paste0("difference_", discovery_cohorts)

gene_stability <- merge(
  data.frame(gene = names(gene_frequency), selection_frequency = as.numeric(gene_frequency), stringsAsFactors = FALSE),
  direction_table,
  by = "gene", sort = FALSE
)
gene_stability <- gene_stability[order(-gene_stability$consistent_discovery_direction,
                                       -gene_stability$selection_frequency,
                                       -gene_stability$mean_abs_difference), ]
rownames(gene_stability) <- NULL

eligible <- gene_stability[gene_stability$consistent_discovery_direction, ]
stable_genes <- eligible$gene[eligible$selection_frequency >= 0.50]
stable_genes <- head(stable_genes, 10L)
if (length(stable_genes) < 5L) stable_genes <- head(eligible$gene, min(5L, nrow(eligible)))
gene_stability$stable_core <- gene_stability$gene %in% stable_genes
write.csv(gene_stability, file.path(out_dir, "gene_stability_summary.csv"), row.names = FALSE)
writeLines(stable_genes, file.path(out_dir, "stable_genes.txt"))

train_idx <- discovery_idx
external_idx <- which(manifest$set == "external")
scaler <- scale_fit(X[train_idx, stable_genes, drop = FALSE])
x_train <- scale_apply(X[train_idx, stable_genes, drop = FALSE], scaler)
x_external <- scale_apply(X[external_idx, stable_genes, drop = FALSE], scaler)
y_train <- manifest$label[train_idx]

final_cv <- glmnet::cv.glmnet(x_train, y_train, family = "binomial", alpha = 0,
                              foldid = make_foldid(y_train, 5L, 20260809L), type.measure = "auc",
                              standardize = FALSE)
final_fit <- glmnet::glmnet(x_train, y_train, family = "binomial", alpha = 0,
                            lambda = final_cv$lambda.min, standardize = FALSE)
external_score <- as.numeric(predict(final_fit, newx = x_external, type = "response"))
external_predictions <- data.frame(
  manifest[external_idx, ], score = external_score, stringsAsFactors = FALSE
)
write.csv(external_predictions, file.path(out_dir, "external_predictions.csv"), row.names = FALSE)

external_auc_rows <- lapply(unique(external_predictions$cohort), function(cohort) {
  d <- external_predictions[external_predictions$cohort == cohort, ]
  roc_obj <- pROC::roc(d$label, d$score, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  data.frame(cohort = cohort, n = nrow(d), controls = sum(d$label == 0), cases = sum(d$label == 1),
             auc = as.numeric(pROC::auc(roc_obj)), ci_low = ci[1], ci_high = ci[3], stringsAsFactors = FALSE)
})
external_auc <- do.call(rbind, external_auc_rows)
write.csv(external_auc, file.path(out_dir, "external_validation_auc.csv"), row.names = FALSE)

external_gene_auc <- do.call(rbind, lapply(stable_genes, function(gene) {
  direction <- gene_stability$discovery_direction[match(gene, gene_stability$gene)]
  do.call(rbind, lapply(unique(manifest$cohort[external_idx]), function(cohort) {
    idx <- external_idx[manifest$cohort[external_idx] == cohort]
    score <- X[idx, gene] * ifelse(direction == "UC_up", 1, -1)
    data.frame(gene = gene, cohort = cohort, auc_oriented = auc_fixed(manifest$label[idx], score),
               discovery_direction = direction, stringsAsFactors = FALSE)
  }))
}))
write.csv(external_gene_auc, file.path(out_dir, "external_single_gene_auc.csv"), row.names = FALSE)

saveRDS(
  list(model = final_fit, scaler = scaler, stable_genes = stable_genes,
       normalization = "within-sample rank inverse-normal over genes shared by all six cohorts",
       candidate_genes = candidate_genes, manifest = manifest),
  file.path(out_dir, "final_ridge_signature.rds")
)

source(file.path(out_dir, "plot_ml_figures.R"), local = TRUE)
write_ml_figures(out_dir)

pdf(file.path(out_dir, "external_roc.pdf"), width = 7, height = 6)
plot(0:1, 0:1, type = "n", xlab = "1 - Specificity", ylab = "Sensitivity", main = "Locked external validation")
abline(0, 1, lty = 2, col = "grey60")
cols <- c("#1B9E77", "#D95F02")
leg <- character()
for (i in seq_along(unique(external_predictions$cohort))) {
  cohort <- unique(external_predictions$cohort)[i]
  d <- external_predictions[external_predictions$cohort == cohort, ]
  ro <- pROC::roc(d$label, d$score, levels = c(0, 1), direction = "<", quiet = TRUE)
  lines(1 - ro$specificities, ro$sensitivities, col = cols[i], lwd = 2)
  leg[i] <- sprintf("%s AUC %.3f", cohort, as.numeric(pROC::auc(ro)))
}
legend("bottomright", legend = leg, col = cols, lwd = 2, bty = "n")
dev.off()

top10 <- head(workflow_summary[workflow_summary$successful_cohorts == 4L, ], 10L)
summary_lines <- c(
  "Leakage-controlled 117-workflow UC classifier reanalysis",
  paste0("Run date: ", Sys.Date()),
  paste0("Candidate panel (n=", length(candidate_genes), "): ", paste(candidate_genes, collapse = ", ")),
  paste0("Discovery samples: ", length(discovery_idx), " across ", length(discovery_cohorts), " cohorts"),
  paste0("External samples: ", length(external_idx), " across 2 locked cohorts"),
  paste0("Stable genes: ", paste(stable_genes, collapse = ", ")),
  "",
  "External validation:",
  apply(external_auc, 1, function(z) sprintf("%s: AUC %.3f (95%% CI %.3f-%.3f), n=%s", z[["cohort"]],
                                      as.numeric(z[["auc"]]), as.numeric(z[["ci_low"]]),
                                      as.numeric(z[["ci_high"]]), z[["n"]])),
  "",
  "Top 10 discovery LOCO workflows:",
  apply(top10, 1, function(z) sprintf("%s: mean AUC %.3f, minimum cohort AUC %.3f, SD %.3f",
                                     z[["workflow"]], as.numeric(z[["mean_auc"]]),
                                     as.numeric(z[["min_auc"]]), as.numeric(z[["sd_auc"]])))
)
writeLines(summary_lines, file.path(out_dir, "run_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
quit(save = "no", status = 0L)
