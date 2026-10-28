## UNIFIED DENGUE PATIENT CLASSIFIER - ELASTIC-NET LOGISTIC REGRESSION
## Tasks:  AP_RP  (acute vs recovery phase; all proteins; paired samples)
##         HP_SUB (hospitalized vs subclinical, acute phase)
##         DF_DHF (DF vs DHF/DSS, acute phase)

##CONFIG
CLASSIFICATION_TYPE <- "AP_RP"   # "AP_RP", "HP_SUB", "DF_DHF"
FEATURE_SET    <- "all"           # "panel": paper's selection rule, re-applied per fold
                                  #   (HP_SUB: BH p < 0.05, DF_DHF: p < 0.01, AP_RP: all proteins)
                                  # "all":   all proteins
USE_COVARIATES <- FALSE           # add age, sex, day of fever (not valid for AP_RP)

## command-line override (for SLURM):  Rscript script.R DF_DHF [all|panel] [covs]
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) CLASSIFICATION_TYPE <- .args[1]
if (length(.args) >= 2) FEATURE_SET         <- .args[2]
if (length(.args) >= 3) USE_COVARIATES      <- tolower(.args[3]) %in% c("covs", "true", "1")
FEATURE_SET <- tolower(trimws(FEATURE_SET))
if (!FEATURE_SET %in% c("all", "panel"))
  stop("FEATURE_SET must be 'all' or 'panel', got '", FEATURE_SET, "'")

RUN_HOLDOUT    <- TRUE            # also report a 75/25 patient-level hold-out split
HOLDOUT_FRAC   <- 0.25
ALPHA          <- 0.5             # used only when TUNE_ALPHA = FALSE
TUNE_ALPHA     <- TRUE            # choose alpha from ALPHA_GRID inside each fold (nested)
ALPHA_GRID     <- seq(0.1, 1, by = 0.1)
INNER_FOLDS    <- 5               # patient-grouped inner CV (fixed; not a model property)
CI_BOOT        <- 2000            # precision of the CI band only - does not change AUROC
N_CORES        <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",          # outer folds in parallel
                                        max(1, parallel::detectCores() - 1)))  # (always 1 on Windows)
SEED           <- 2026

LFQ_FILE      <- "20260120_LFQ_DENV_AP_RP.csv"
COHORT_FILE   <- "cohort_data_dengue_latest.csv"
OUTPUT_PREFIX <- paste0("LR_", CLASSIFICATION_TYPE, "_", FEATURE_SET, if (USE_COVARIATES) "_covs" else "")
OUT_DIR       <- "results_LR"
BASE_SIZE     <- 40               # 40 + 19 in matches the v1 figures; use 14 + 7 in for screen
FIG_IN        <- 19

##PACKAGES
pkgs <- c("dplyr", "tidyr", "ggplot2", "pROC", "glmnet", "gridExtra")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2)
                                 library(pROC); library(glmnet) })
set.seed(SEED)
dir.create(OUT_DIR, showWarnings = FALSE)

##Okabe-Ito palette (safe for protan/deutan/tritan vision) + line types
PAL_METHOD <- c("Train" = "#E69F00", "Test" = "#0072B2", "LOPO-CV" = "#000000")
LTY_METHOD <- c("Train" = "dashed", "Test" = "dotdash", "LOPO-CV" = "solid")
PAL_REF    <- "grey60"

##DATA
parse_day_of_fever <- function(x)                       # "d3" / "d10" -> 3 / 10
  suppressWarnings(as.numeric(sub("^d", "", trimws(tolower(as.character(x))))))

load_task <- function(type) {
  spec <- switch(type,
    AP_RP  = list(groups = c("AP", "RP"), col = "Group",             levels = c("AP", "RP")),
    HP_SUB = list(groups = "AP",          col = "Classification",    levels = c("hospitalized", "subclinical")),
    DF_DHF = list(groups = "AP",          col = "Classification_02", levels = c("DF", "DHF")),
    stop("CLASSIFICATION_TYPE must be AP_RP, HP_SUB or DF_DHF"))
  lfq <- read.csv(LFQ_FILE, stringsAsFactors = FALSE)
  lfq$Subject <- trimws(as.character(lfq$Subject))
  cohort <- read.csv(COHORT_FILE, sep = ";", colClasses = "character") %>%
    transmute(Subject = trimws(Subject), Age = as.numeric(Age),
              male = as.integer(trimws(Sex) == "M"), DoF = parse_day_of_fever(Day_of_fever))
  w <- lfq %>%
    filter(Group %in% spec$groups) %>%
    mutate(class = trimws(.data[[spec$col]])) %>%
    filter(class %in% spec$levels) %>%
    group_by(Subject, Group, class, Genes) %>%                  # average any technical replicates
    summarise(v = mean(log2_LFQ, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Genes, values_from = v) %>%
    left_join(cohort, by = "Subject")
  if (any(is.na(w$Age))) warning(sum(is.na(w$Age)), " samples have no cohort record")
  genes <- setdiff(names(w), c("Subject", "Group", "class", "Age", "male", "DoF"))
  n_dof_na <- sum(is.na(w$DoF))
  if (type != "AP_RP" && n_dof_na)               # used by covariate models and the baseline
    message("Day of fever missing for ", n_dof_na, " samples -> set to 0 ",
            "(asymptomatic cases: acute sample taken at inclusion). Confirm with clinical team.")
  list(X = as.matrix(w[, genes]), y = factor(w$class, levels = spec$levels),
       id = w$Subject, levels = spec$levels,
       C = as.matrix(data.frame(Age = w$Age, male = w$male,
                                DoF = ifelse(is.na(w$DoF), 0, w$DoF))))
}

##PREPROCESSING
impute_fill <- function(X)                               # per-protein fill value (log2 scale)
  apply(X, 2, function(v) if (all(is.na(v))) 0 else min(v, na.rm = TRUE) - 1)
apply_fill <- function(X, fill) { for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- fill[j]; X }

gene_p <- function(X, y) {                               # gene-wise linear model = equal-var t-test
  pos <- y == levels(y)[2]
  apply(X, 2, function(g) tryCatch(t.test(g[pos], g[!pos], var.equal = TRUE)$p.value,
                                   error = function(e) 1))
}

##Protein selection - called on TRAINING patients only
select_features <- function(X, y, method, type, min_genes = 5) {
  if (method == "all" || (method == "panel" && type == "AP_RP")) return(colnames(X))
  if (method == "panel") {
    p <- gene_p(X, y)
    s <- if (type == "HP_SUB") names(p)[p.adjust(p, "BH") < 0.05] else names(p)[p < 0.01]
    if (length(s) < min_genes) s <- names(sort(p))[seq_len(min_genes)]
    return(s)
  }
  stop("Unknown FEATURE_SET: ", method)
}

##INNER CV FOLDS (used by cv.glmnet)
##Patient-grouped inner folds, stratified by each patient's class pattern
inner_folds <- function(y, id, K) {
  pat <- tapply(as.character(y), id, function(v) paste(sort(unique(v)), collapse = "+"))
  f <- setNames(integer(length(pat)), names(pat))
  for (g in split(names(pat), pat)) { g <- g[sample.int(length(g))]; f[g] <- rep_len(seq_len(K), length(g)) }
  unname(f[id])
}
##Number of inner folds: at most INNER_FOLDS, and never more than the number of
##minority-class patients minus 1 (so every training part keeps >= 2 of each class)
safe_k <- function(y, id) {
  n_min <- min(tapply(id, y, function(v) length(unique(v))))
  max(2, min(INNER_FOLDS, n_min - 1, length(unique(id))))
}

##MODEL
##Elastic net. Features standardised on training data; covariates (if any)
##are appended last and left unpenalised so proteins must add to them.
fit_glmnet <- function(X, y, id, n_unpen = 0) {
  ctr <- colMeans(X); scl <- apply(X, 2, sd); scl[!is.finite(scl) | scl == 0] <- 1
  Xs <- scale(X, ctr, scl)
  wts <- as.numeric((max(table(y)) / table(y))[as.character(y)])      # balanced case weights
  fid <- inner_folds(y, id, safe_k(y, id))      # patient-grouped AND class-stratified
  pf <- c(rep(1, ncol(X) - n_unpen), rep(0, n_unpen))
  alphas <- if (TUNE_ALPHA) ALPHA_GRID else ALPHA
  fits <- lapply(alphas, function(a)
    withCallingHandlers(   # silence only glmnet's benign small-fold notice
      cv.glmnet(Xs, y, family = "binomial", alpha = a, weights = wts, foldid = fid,
                penalty.factor = pf, type.measure = "deviance"),
      warning = function(cond) if (grepl("dangerous ground", conditionMessage(cond)))
                                 invokeRestart("muffleWarning")))
  best <- which.min(sapply(fits, function(f) min(f$cvm)))            # same folds for every alpha
  list(fit = fits[[best]], alpha = alphas[best], center = ctr, scale = scl)
}

##Train on one set of patients; returns everything needed to predict new ones
train_model <- function(X, y, C = NULL, id = seq_len(nrow(X))) {
  fill  <- impute_fill(X); Xi <- apply_fill(X, fill)
  genes <- select_features(Xi, y, FEATURE_SET, CLASSIFICATION_TYPE)
  Xm <- cbind(Xi[, genes, drop = FALSE], C)
  list(fill = fill, genes = genes, use_C = !is.null(C),
       model = fit_glmnet(Xm, y, id, n_unpen = if (is.null(C)) 0 else ncol(C)))
}
predict_model <- function(obj, X, C = NULL) {
  Xm <- cbind(apply_fill(X, obj$fill)[, obj$genes, drop = FALSE], if (obj$use_C) C)
  Xs <- scale(Xm, obj$model$center, obj$model$scale)
  as.numeric(predict(obj$model$fit, Xs, s = "lambda.1se", type = "response"))   # P(second level)
}
n_nonzero <- function(obj) {                             # proteins kept by the elastic net
  b <- as.matrix(coef(obj$model$fit, s = "lambda.1se"))[-1, 1]
  sum(b[names(b) %in% obj$genes] != 0)
}
subC <- function(C, i) if (is.null(C)) NULL else C[i, , drop = FALSE]

##EVALUATION HELPERS
roc_fixed <- function(y, p, levels)
  roc(factor(y, levels = levels), p, levels = levels, direction = "<", quiet = TRUE)
youden_threshold <- function(y, p, levels)
  as.numeric(coords(roc_fixed(y, p, levels), "best", best.method = "youden",
                    ret = "threshold", transpose = FALSE)[1, 1])
summarise_perf <- function(y, p, levels, thr, label, id) {
  r <- roc_fixed(y, p, levels); ci <- as.numeric(ci.auc(r, method = "delong"))
  pos <- y == levels[2]; call <- p >= thr
  data.frame(set = label, n_patients = length(unique(id)), n_samples = length(y),
             n_positive = sum(pos), AUROC = ci[2], CI_low = ci[1], CI_high = ci[3],
             threshold = thr,
             sensitivity = mean(call[pos]), specificity = mean(!call[!pos]),
             PPV = if (any(call)) mean(pos[call]) else NA,
             NPV = if (any(!call)) mean(!pos[!call]) else NA)
}
empty_warning <- function(nz, label) {
  e <- sum(nz == 0)
  if (e > 0) warning(sprintf(paste0("%s: elastic net kept NO proteins in %d of %d folds (%.0f%%). ",
                     "Those patients received a constant prediction; the AUROC mostly reflects ",
                     "the absence of a model, not a weak model. Report as 'no stable model'."),
                     label, e, length(nz), 100 * e / length(nz)), call. = FALSE)
  e
}

##Covariates-only baseline: logistic regression on age, sex, day of fever (LOPO)
cov_only_lopo <- function(C, y, id) {
  d <- as.data.frame(C); pred <- rep(NA_real_, length(y))
  for (p in unique(id)) {
    te <- id == p; tr <- !te
    w <- as.numeric((max(table(y[tr])) / table(y[tr]))[as.character(y[tr])])
    f <- suppressWarnings(glm(yy ~ ., data = cbind(d[tr, , drop = FALSE], yy = y[tr]),
                              family = binomial, weights = w))
    pred[te] <- predict(f, d[te, , drop = FALSE], type = "response")
  }
  pred
}

##Nested leave-one-patient-out: everything (imputation, protein selection,
##alpha/lambda tuning) is redone in each fold on training patients only.
##Folds run in parallel; each fold has its own seed, so results do not depend
##on N_CORES.
nested_lopo <- function(X, y, id, C = NULL, verbose = TRUE) {
  ids <- unique(id)
  run_fold <- function(k) {
    set.seed(SEED + k)
    te <- id == ids[k]; tr <- !te
    obj <- train_model(X[tr, , drop = FALSE], droplevels(y[tr]), subC(C, tr), id[tr])
    list(te = which(te), p = predict_model(obj, X[te, , drop = FALSE], subC(C, te)),
         genes = obj$genes, nz = n_nonzero(obj), alpha = obj$model$alpha)
  }
  if (verbose) message("  ", length(ids), " outer folds on ", N_CORES, " core(s) ...")
  res <- if (N_CORES > 1 && .Platform$OS.type == "unix")
    parallel::mclapply(seq_along(ids), run_fold, mc.cores = N_CORES, mc.preschedule = FALSE)
  else lapply(seq_along(ids), run_fold)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("Outer fold(s) failed: ", paste(ids[bad], collapse = ", "), "\n", res[[which(bad)[1]]])
  pred <- rep(NA_real_, nrow(X))
  for (r in res) pred[r$te] <- r$p
  nz <- vapply(res, `[[`, numeric(1), "nz")
  list(pred = pred,
       chosen = setNames(lapply(res, `[[`, "genes"), ids),
       n_nonzero = setNames(nz, ids),
       choices = data.frame(patient = ids, alpha = vapply(res, `[[`, numeric(1), "alpha"),
                            n_proteins_selected = lengths(lapply(res, `[[`, "genes")),
                            n_proteins_kept = nz))
}

##MAIN
message("ELASTIC NET | ", CLASSIFICATION_TYPE, " | features: ", FEATURE_SET,
        " | covariates: ", USE_COVARIATES)
if (USE_COVARIATES && CLASSIFICATION_TYPE == "AP_RP")
  stop("Covariates are not meaningful for AP_RP: each patient contributes both phases, ",
       "and day of fever defines the phases.")
tk <- load_task(CLASSIFICATION_TYPE)
X <- tk$X; y <- tk$y; id <- tk$id; L <- tk$levels
C <- if (USE_COVARIATES) tk$C else NULL
message(length(unique(id)), " patients, ", nrow(X), " samples, ", ncol(X), " proteins; ",
        paste(names(table(y)), table(y), sep = " = ", collapse = ", "))

##Nested LOPO on all patients, pre-specified cutoff 0.5
message("Nested LOPO-CV ...")
lopo <- nested_lopo(X, y, id, C)
perf <- summarise_perf(y, lopo$pred, L, 0.5, "LOPO-CV", id)
perf$empty_model_folds <- empty_warning(lopo$n_nonzero, "LOPO-CV")
perf$median_proteins_kept <- median(lopo$n_nonzero)

##Hold-out split (patient level, stratified by class pattern)
hold <- NULL
if (RUN_HOLDOUT) {
  set.seed(SEED)                                   # same split as the RF script
  message("Hold-out split ...")
  pat <- tapply(as.character(y), id, function(v) paste(sort(unique(v)), collapse = "+"))
  test_ids <- unlist(lapply(split(names(pat), pat), function(s)
    sample(s, max(1, round(HOLDOUT_FRAC * length(s))))))
  te <- id %in% test_ids; tr <- !te
  Ctr <- subC(C, tr); Cte <- subC(C, te)
  inner <- nested_lopo(X[tr, ], droplevels(y[tr]), id[tr], Ctr, verbose = FALSE)   # training CV
  thr <- youden_threshold(y[tr], inner$pred, L)                                     # cutoff from training only
  set.seed(SEED)                                   # inner folds of the final fit: independent of N_CORES
  final_tr <- train_model(X[tr, ], droplevels(y[tr]), Ctr, id[tr])
  p_te <- predict_model(final_tr, X[te, ], Cte)
  nz_tr <- n_nonzero(final_tr)
  if (nz_tr == 0) warning("Hold-out: the model trained on the training split kept NO proteins - ",
                          "every test patient gets the same prediction (test AUROC = 0.5 by construction).",
                          call. = FALSE)
  pt <- summarise_perf(y[tr], inner$pred, L, thr, "Train", id[tr])
  pt$empty_model_folds <- empty_warning(inner$n_nonzero, "Train (inner CV)")
  pt$median_proteins_kept <- median(inner$n_nonzero)
  ps <- summarise_perf(y[te], p_te, L, thr, "Test", id[te])
  ps$empty_model_folds <- as.integer(nz_tr == 0); ps$median_proteins_kept <- nz_tr
  perf <- rbind(pt, ps, perf)
  hold <- list(tr = tr, te = te, p_tr = inner$pred, p_te = p_te)
}
num <- sapply(perf, is.numeric); perf[num] <- round(perf[num], 3)
print(perf, row.names = FALSE)
write.csv(perf, file.path(OUT_DIR, paste0(OUTPUT_PREFIX, "_performance.csv")), row.names = FALSE)
write.csv(data.frame(patient = id, sample_class = y, p_LOPO = lopo$pred),
          file.path(OUT_DIR, paste0(OUTPUT_PREFIX, "_lopo_predictions.csv")), row.names = FALSE)
if (CLASSIFICATION_TYPE != "AP_RP") {             # baseline for the supplementary overlay figure
  p_cov <- cov_only_lopo(tk$C, y, id)
  write.csv(data.frame(patient = id, sample_class = y, p_LOPO = p_cov),
            file.path(OUT_DIR, paste0(CLASSIFICATION_TYPE, "_covariates_only_lopo_predictions.csv")),
            row.names = FALSE)
}

##Panel stability across LOPO folds
stab <- sort(table(unlist(lopo$chosen)) / length(lopo$chosen), decreasing = TRUE)
stab <- data.frame(protein = names(stab), selected_fraction = round(as.numeric(stab), 3))
write.csv(stab, file.path(OUT_DIR, paste0(OUTPUT_PREFIX, "_selection_frequency.csv")), row.names = FALSE)
write.csv(lopo$choices, file.path(OUT_DIR, paste0(OUTPUT_PREFIX, "_model_choices.csv")), row.names = FALSE)
message("Feature set: ", FEATURE_SET,
        " | proteins kept per fold (median, range): ", median(lopo$n_nonzero),
        " (", min(lopo$n_nonzero), "-", max(lopo$n_nonzero), ")",
        " | alpha (median, range): ", median(lopo$choices$alpha),
        " (", min(lopo$choices$alpha), "-", max(lopo$choices$alpha), ")")

##Final model on ALL patients -> standardised coefficients (describes the
##model; performance comes from the nested LOPO above). Sign shows the
##direction: positive = higher in the second class.
set.seed(SEED)                                     # inner folds: independent of N_CORES
final <- train_model(X, y, C, id)
cf <- as.matrix(coef(final$model$fit, s = "lambda.1se"))[-1, 1]
imp <- data.frame(Feature = names(cf), Coefficient = cf) %>%
  filter(Coefficient != 0) %>% arrange(desc(abs(Coefficient)))
message("Final model: feature set = ", FEATURE_SET, ", alpha = ", final$model$alpha, ", ",
        nrow(imp), " non-zero coefficients")
write.csv(imp, file.path(OUT_DIR, paste0(OUTPUT_PREFIX, "_coefficients.csv")), row.names = FALSE)

##FIGURES
save_fig <- function(p, name, w = FIG_IN, h = FIG_IN) {
  base <- file.path(OUT_DIR, paste0(OUTPUT_PREFIX, "_", name))
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300, units = "in", bg = "white")
  ggsave(paste0(base, ".pdf"), p, width = w, height = h, units = "in", bg = "white")
  if (requireNamespace("svglite", quietly = TRUE))
    ggsave(paste0(base, ".svg"), p, width = w, height = h, units = "in", bg = "white")
}
roc_df <- function(y, p, method, band = FALSE) {
  r <- roc_fixed(y, p, L)
  out <- list(curve = data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, Method = method))
  if (band) {                                             # bootstrap CI band (LOPO-CV only)
    ci <- ci.se(r, specificities = seq(0, 1, 0.02), boot.n = CI_BOOT, progress = "none")
    out$band <- data.frame(FPR = 1 - as.numeric(rownames(ci)), lo = ci[, 1], hi = ci[, 3], Method = method)
  }
  out
}
sets <- list(`LOPO-CV` = list(y = y, p = lopo$pred))
if (!is.null(hold)) sets <- c(list(Train = list(y = y[hold$tr], p = hold$p_tr),
                                   Test  = list(y = y[hold$te], p = hold$p_te)), sets)
## shaded band = 95% CI of the nested LOPO-CV curve; train/test bands omitted (too few patients)
rd <- lapply(names(sets), function(m) roc_df(sets[[m]]$y, sets[[m]]$p, m, band = m == "LOPO-CV"))
curves <- bind_rows(lapply(rd, `[[`, "curve"))            # keep pROC point order
bands  <- bind_rows(lapply(rd, `[[`, "band"))
curves$Method <- factor(curves$Method, names(PAL_METHOD)); bands$Method <- factor(bands$Method, names(PAL_METHOD))

tab <- perf %>% transmute(Method = set, N = n_patients,
                          `AUROC (95% CI)` = sprintf("%.3f (%.2f-%.2f)", AUROC, CI_low, CI_high))
tg <- gridExtra::tableGrob(tab, rows = NULL, theme = gridExtra::ttheme_minimal(
  core = list(fg_params = list(fontsize = BASE_SIZE * 0.75,
                               col = matrix(c(PAL_METHOD[tab$Method], rep("black", 2 * nrow(tab))), nrow(tab)))),
  colhead = list(fg_params = list(fontsize = BASE_SIZE * 0.75, fontface = "bold"))))

p_roc <- ggplot(curves, aes(FPR, TPR, colour = Method, linetype = Method)) +
  geom_ribbon(data = bands, aes(x = FPR, ymin = lo, ymax = hi, fill = Method),
              inherit.aes = FALSE, alpha = 0.15) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = PAL_REF) +
  geom_path(linewidth = 1.6) +
  annotation_custom(tg, xmin = 0.40, xmax = 1, ymin = 0, ymax = 0.32) +
  scale_colour_manual(values = PAL_METHOD) + scale_fill_manual(values = PAL_METHOD, guide = "none") +
  scale_linetype_manual(values = LTY_METHOD) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)") +
  theme_minimal(base_size = BASE_SIZE) +
  theme(legend.position = "none", plot.margin = ggplot2::margin(10, BASE_SIZE, 10, 10))
save_fig(p_roc, "roc")

##Same style as the RF plot: importance = |standardised coefficient|, descending.
##(Gini importance exists only for tree models; for a penalised logistic model the
##magnitude of the standardised coefficient is the direct equivalent. Signs are
##kept in the *_coefficients.csv file.)
top <- imp %>% filter(!Feature %in% c("Age", "male", "DoF")) %>%       # proteins only
  mutate(Importance = abs(Coefficient)) %>% arrange(desc(Importance)) %>% head(20)
if (nrow(top) > 0) {
  top$Feature <- factor(top$Feature, rev(top$Feature))
  p_imp <- ggplot(top, aes(Feature, Importance, fill = Importance)) + geom_col() + coord_flip() +
    scale_fill_viridis_c(option = "plasma", begin = 0.05, end = 0.85, direction = -1,
                         name = "Importance\n(|Standardised\ncoefficient|)",
                         guide = guide_colourbar(barheight = grid::unit(0.45, "npc"),   # tall enough that
                                                 barwidth = grid::unit(BASE_SIZE * 0.6, "pt"))) +   # ticks never overlap
    labs(x = NULL, y = "Absolute standardised coefficient") + theme_minimal(base_size = BASE_SIZE)
  save_fig(p_imp, "feature_importance", w = FIG_IN, h = FIG_IN * 0.8)

  ##Lollipop version of the same ranking (reviewer request)
  p_lol <- ggplot(top, aes(x = Feature, y = Importance, colour = Importance)) +
    geom_segment(aes(xend = Feature, y = 0, yend = Importance), linewidth = BASE_SIZE / 16) +
    geom_point(size = BASE_SIZE / 4) + coord_flip() +
    scale_colour_viridis_c(option = "plasma", begin = 0.05, end = 0.85, direction = -1,
                           name = "Importance\n(|Standardised\ncoefficient|)",
                           guide = guide_colourbar(barheight = grid::unit(0.45, "npc"),
                                                   barwidth = grid::unit(BASE_SIZE * 0.6, "pt"))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Absolute standardised coefficient") + theme_minimal(base_size = BASE_SIZE) +
    theme(panel.grid.major.y = element_blank())
  save_fig(p_lol, "feature_importance_lollipop", w = FIG_IN, h = FIG_IN * 0.8)
} else message("No protein has a non-zero coefficient at lambda.1se - no importance plot.")

message("Done: outputs in ", OUT_DIR, "/")
