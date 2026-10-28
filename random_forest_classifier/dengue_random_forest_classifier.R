##UNIFIED DENGUE PATIENT CLASSIFIER:RANDOM FOREST

##Tasks: AP_RP  (acute vs recovery phase; all proteins; paired samples)
##       HP_SUB (hospitalized vs subclinical, acute phase)
##       DF_DHF (DF vs DHF/DSS, acute phase)

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
NTREE          <- 5000            # fixed: more trees never overfit, they only remove randomness
TUNE_RF        <- TRUE            # tune mtry/nodesize by OOB AUROC inside each fold (nested; slower)
CI_BOOT        <- 2000            # precision of the CI band only - does not change AUROC
N_CORES        <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",          # outer folds in parallel
                                        max(1, parallel::detectCores() - 1)))  # (always 1 on Windows)
SEED           <- 2026

LFQ_FILE      <- "20260120_LFQ_DENV_AP_RP.csv"
COHORT_FILE   <- "cohort_data_dengue_latest.csv"
OUTPUT_PREFIX <- paste0("RF_", CLASSIFICATION_TYPE, "_", FEATURE_SET, if (USE_COVARIATES) "_covs" else "")
OUT_DIR       <- "results_RF"
BASE_SIZE     <- 40               # 40 + 19 in matches the v1 figures; use 14 + 7 in for screen
FIG_IN        <- 19

##PACKAGES
pkgs <- c("dplyr", "tidyr", "ggplot2", "pROC", "randomForest", "gridExtra")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2)
                                 library(pROC); library(randomForest) })
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
    message("Day of fever missing for ", n_dof_na, " samples -> set to 0")
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

##MODEL
fit_rf <- function(X, y, ntree = NTREE, tune = TUNE_RF) {
  n_min <- min(table(y))
  grow <- function(mtry, nodesize) randomForest(x = X, y = y, ntree = ntree, mtry = mtry,
                                                nodesize = nodesize, strata = y,
                                                sampsize = rep(n_min, nlevels(y)))
  p <- ncol(X); mtry0 <- max(1, floor(sqrt(p)))
  if (!tune) return(grow(mtry0, 1))
  grid <- expand.grid(mtry = unique(pmax(1, c(floor(sqrt(p)), floor(log2(p)), floor(p / 3)))),
                      nodesize = c(1, 3, 5))
  oob <- apply(grid, 1, function(g) {                    # OOB AUROC: no extra data used
    f <- grow(g[["mtry"]], g[["nodesize"]])
    as.numeric(auc(roc(y, f$votes[, levels(y)[2]], levels = levels(y), direction = "<", quiet = TRUE)))
  })
  b <- grid[which.max(oob), ]
  grow(b$mtry, b$nodesize)
}

##Train on one set of patients; returns everything needed to predict new ones
train_model <- function(X, y, C = NULL) {
  fill  <- impute_fill(X); Xi <- apply_fill(X, fill)
  genes <- select_features(Xi, y, FEATURE_SET, CLASSIFICATION_TYPE)
  Xm <- cbind(Xi[, genes, drop = FALSE], C)
  list(fill = fill, genes = genes, model = fit_rf(Xm, y), use_C = !is.null(C))
}
predict_model <- function(obj, X, C = NULL) {
  Xm <- cbind(apply_fill(X, obj$fill)[, obj$genes, drop = FALSE], if (obj$use_C) C)
  as.numeric(predict(obj$model, Xm, type = "prob")[, 2])     # P(second level)
}

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

subC <- function(C, i) if (is.null(C)) NULL else C[i, , drop = FALSE]

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
##tuning) is redone in each fold on training patients only.
##Folds run in parallel; each fold has its own seed, so results do not depend
##on N_CORES.
nested_lopo <- function(X, y, id, C = NULL, verbose = TRUE) {
  ids <- unique(id)
  run_fold <- function(k) {
    set.seed(SEED + k)
    te <- id == ids[k]; tr <- !te
    obj <- train_model(X[tr, , drop = FALSE], droplevels(y[tr]), subC(C, tr))
    list(te = which(te), p = predict_model(obj, X[te, , drop = FALSE], subC(C, te)),
         genes = obj$genes)
  }
  if (verbose) message("  ", length(ids), " outer folds on ", N_CORES, " core(s) ...")
  res <- if (N_CORES > 1 && .Platform$OS.type == "unix")
    parallel::mclapply(seq_along(ids), run_fold, mc.cores = N_CORES, mc.preschedule = FALSE)
  else lapply(seq_along(ids), run_fold)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("Outer fold(s) failed: ", paste(ids[bad], collapse = ", "), "\n", res[[which(bad)[1]]])
  pred <- rep(NA_real_, nrow(X))
  for (r in res) pred[r$te] <- r$p
  chosen <- setNames(lapply(res, `[[`, "genes"), ids)
  list(pred = pred, chosen = chosen,
       n_proteins = data.frame(patient = ids, n_proteins = lengths(chosen)))
}

##MAIN
message("RANDOM FOREST | ", CLASSIFICATION_TYPE, " | features: ", FEATURE_SET,
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

##Hold-out split (patient level, stratified by class pattern)
hold <- NULL
if (RUN_HOLDOUT) {
  set.seed(SEED)
  message("Hold-out split ...")
  pat <- tapply(as.character(y), id, function(v) paste(sort(unique(v)), collapse = "+"))
  test_ids <- unlist(lapply(split(names(pat), pat), function(s)
    sample(s, max(1, round(HOLDOUT_FRAC * length(s))))))
  te <- id %in% test_ids; tr <- !te
  Ctr <- if (!is.null(C)) C[tr, , drop = FALSE]; Cte <- if (!is.null(C)) C[te, , drop = FALSE]
  inner <- nested_lopo(X[tr, ], droplevels(y[tr]), id[tr], Ctr, verbose = FALSE)   # training CV
  thr <- youden_threshold(y[tr], inner$pred, L)                                     # cutoff from training only
  set.seed(SEED)                                   # final RF fit: independent of N_CORES
  final_tr <- train_model(X[tr, ], droplevels(y[tr]), Ctr)
  p_te <- predict_model(final_tr, X[te, ], Cte)
  perf <- rbind(summarise_perf(y[tr], inner$pred, L, thr, "Train", id[tr]),
                summarise_perf(y[te], p_te, L, thr, "Test", id[te]), perf)
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
write.csv(lopo$n_proteins, file.path(OUT_DIR, paste0(OUTPUT_PREFIX, "_proteins_per_fold.csv")), row.names = FALSE)
n_p <- lopo$n_proteins$n_proteins
message("Feature set: ", FEATURE_SET, " | proteins per fold (median, range): ",
        median(n_p), " (", min(n_p), "-", max(n_p), ")")

##Final model on ALL patients -> feature importance (describes the model;
##performance comes from the nested LOPO above)
set.seed(SEED)                                     # independent of N_CORES
final <- train_model(X, y, C)
message("Final model: feature set = ", FEATURE_SET, ", ", length(final$genes), " proteins")
imp <- randomForest::importance(final$model)
imp <- data.frame(Feature = rownames(imp), MeanDecreaseGini = imp[, "MeanDecreaseGini"]) %>%
  arrange(desc(MeanDecreaseGini))
write.csv(imp, file.path(OUT_DIR, paste0(OUTPUT_PREFIX, "_feature_importance.csv")), row.names = FALSE)

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

top <- head(filter(imp, !Feature %in% c("Age", "male", "DoF")), 20)   # proteins only
top$Feature <- factor(top$Feature, rev(top$Feature))
p_imp <- ggplot(top, aes(Feature, MeanDecreaseGini, fill = MeanDecreaseGini)) + geom_col() + coord_flip() +
  scale_fill_viridis_c(option = "plasma", begin = 0.05, end = 0.85, direction = -1,
                       name = "Importance\n(Mean Decrease Gini)",
                       guide = guide_colourbar(barheight = grid::unit(0.45, "npc"),   # tall enough that
                                               barwidth = grid::unit(BASE_SIZE * 0.6, "pt"))) +   # ticks never overlap
  labs(x = NULL, y = "Mean decrease in Gini impurity") + theme_minimal(base_size = BASE_SIZE)
save_fig(p_imp, "feature_importance", w = FIG_IN, h = FIG_IN * 0.8)

##Lollipop version of the same ranking (reviewer request)
p_lol <- ggplot(top, aes(x = Feature, y = MeanDecreaseGini, colour = MeanDecreaseGini)) +
  geom_segment(aes(xend = Feature, y = 0, yend = MeanDecreaseGini), linewidth = BASE_SIZE / 16) +
  geom_point(size = BASE_SIZE / 4) + coord_flip() +
  scale_colour_viridis_c(option = "plasma", begin = 0.05, end = 0.85, direction = -1,
                         name = "Importance\n(Mean Decrease Gini)",
                         guide = guide_colourbar(barheight = grid::unit(0.45, "npc"),
                                                 barwidth = grid::unit(BASE_SIZE * 0.6, "pt"))) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Mean decrease in Gini impurity") + theme_minimal(base_size = BASE_SIZE) +
  theme(panel.grid.major.y = element_blank())
save_fig(p_lol, "feature_importance_lollipop", w = FIG_IN, h = FIG_IN * 0.8)

message("Analysis complete: Outputs in ", OUT_DIR, "/")
