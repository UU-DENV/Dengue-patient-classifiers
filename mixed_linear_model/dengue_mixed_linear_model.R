##Dengue pediatric proteomics: LFQ-per-protein vs clinical variables 

##Objectives:(1)Per-phase (AP, RP) association of each protein's log2 LFQ with 
##clinical variables -> strength, direction, FDR, driver proteins
##(2)Rank which clinical variable the proteome tracks most/least
##(3)Test whether an association HOLDS across phases (AP vs RP) via a clinical x phase interaction mixed model
##(4)Variance partitioning (variancePartition)

##CONFIG
LFQ_FILE       <- "20260120_LFQ_DENV_AP_RP.csv"
COHORT_FILE    <- "cohort_data_dengue_latest.csv"
OUT_DIR        <- "results_MLM"
CT_COLUMN      <- "Viral_load"    # column holding the Ct value ("Viral_load" = v1; or "CT")
ADJUST_AGE_SEX <- TRUE            # add Age + Sex to the association models
MIN_N          <- 20              # min patients with a value, per phase, to test a protein
N_BOOT_VP      <- 100             # bootstrap resamples for the variance-partition ranking
SEED           <- 1

##Packages

install.packages(c("tidyverse", "broom", "lme4", "lmerTest"))

BiocManager::install(c("variancePartition", "limma"))    

suppressPackageStartupMessages({ library(tidyverse); library(lme4); library(lmerTest) })

set.seed(SEED)

dir.create(OUT_DIR, showWarnings = FALSE)

##Okabe-Ito colour-blind-safe palette

okabe_ito <- c(black = "#000000", orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
               yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7")

grey_col <- "grey75"

sig_col  <- unname(okabe_ito["vermillion"])

pos_col  <- unname(okabe_ito["blue"])

neg_col  <- unname(okabe_ito["vermillion"])

save_all <- function(plot, base, width, height, dpi = 300) {
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, bg = "white")
  if (requireNamespace("svglite", quietly = TRUE))
    ggsave(paste0(base, ".svg"), plot, width = width, height = height, bg = "white")
}

##PART 1:CLEAN & RESHAPE

to_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", trimws(as.character(x)), fixed = TRUE)))
parse_day_of_fever <- function(x) suppressWarnings(as.numeric(sub("^d", "", trimws(tolower(as.character(x))))))

cohort <- read.csv(COHORT_FILE, sep = ";", colClasses = "character") %>%
  mutate(Subject = trimws(Subject))
lfq <- read.csv(LFQ_FILE, stringsAsFactors = FALSE) %>%
  mutate(Subject = trimws(as.character(Subject)))

##Ct consistency check (Viral_load vs CT columns)

ct_check <- cohort %>%
  transmute(Subject, Viral_load = to_num(Viral_load), CT = to_num(CT)) %>%
  filter(!is.na(Viral_load) | !is.na(CT), is.na(Viral_load) | is.na(CT) | abs(Viral_load - CT) > 0.01)

if (nrow(ct_check)) {
  message("Viral_load and CT disagree for ", nrow(ct_check), " subjects - verify against source records:")
  print(ct_check, row.names = FALSE)
  write_csv(ct_check, file.path(OUT_DIR, "Ct_column_mismatches.csv"))
}

df <- lfq %>%
  inner_join(cohort, by = "Subject") %>%
  mutate(
    log2_LFQ          = to_num(log2_LFQ),
    Hematocrit_D0     = to_num(Hematocrit_D0),
    Platelet_count_D0 = to_num(Platelet_count_D0),
    Day_of_fever      = parse_day_of_fever(Day_of_fever),
    Viral_burden      = -to_num(.data[[CT_COLUMN]]),     # -Ct: higher = more virus
    Age               = to_num(Age),
    Sex               = factor(trimws(Sex), levels = c("F", "M")),
    Genes             = trimws(as.character(Genes)),
    Group             = trimws(as.character(Group)),
    Vascular_leakage_binary = case_when(
      toupper(trimws(Vascular_leakage)) %in% c("Y", "YES", "1") ~ 1L,
      toupper(trimws(Vascular_leakage)) %in% c("N", "NO", "0")  ~ 0L,
      TRUE ~ NA_integer_),
    Classification    = trimws(Classification),
    Severity          = if_else(trimws(Classification_02) %in% c("DF", "DHF"),
                                trimws(Classification_02), NA_character_)) %>%  # DF vs DHF/DSS only
  filter(!is.na(Genes), Genes != "", !is.na(log2_LFQ))


clinical_vars <- c("Day_of_fever", "Platelet_count_D0", "Hematocrit_D0",
                   "Viral_burden", "Vascular_leakage_binary")
class_vars    <- c("Classification", "Severity")

adj_vars      <- if (ADJUST_AGE_SEX) c("Age", "Sex") else character(0)

pheno <- df %>% select(Patient_code, all_of(clinical_vars), all_of(class_vars), Age, Sex) %>% distinct()

stopifnot(nrow(pheno) == n_distinct(pheno$Patient_code))       # one clinical profile per patient

##Collapse technical replicates: one LFQ per Patient x Group x Protein

expr <- df %>%
  filter(Group %in% c("AP", "RP")) %>%
  group_by(Patient_code, Group, Genes) %>%
  summarise(log2_LFQ = mean(log2_LFQ, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group, Genes) %>% filter(n() >= MIN_N) %>% ungroup()    # filter WITHIN phase

adata <- expr %>% left_join(pheno, by = "Patient_code")

message(n_distinct(expr$Patient_code), " patients; ", n_distinct(expr$Genes), " proteins; ",
        "adjusted for age/sex: ", ADJUST_AGE_SEX, "; children only: ", PEDIATRIC_ONLY)

##PART 2: PER-PHASE, PER-PROTEIN ASSOCIATION
##Marginal: OLS slope + Spearman rho (as v1).
##Adjusted (if ADJUST_AGE_SEX): slope from log2_LFQ ~ var + Age + Sex, and
##partial Spearman rho (ranks residualised on Age + Sex, then correlated).

partial_spearman <- function(y, x, Z) {
  ry <- resid(lm(rank(y) ~ ., data = Z)); rx <- resid(lm(rank(x) ~ ., data = Z))
  r <- cor(ry, rx); dfree <- length(y) - 2 - ncol(Z)
  tt <- r * sqrt(dfree / (1 - r^2))
  c(rho = r, p = 2 * pt(-abs(tt), dfree))
}

fit_one <- function(d, v) {
  d <- d %>% filter(!is.na(log2_LFQ), !is.na(.data[[v]]), if_all(all_of(adj_vars), ~ !is.na(.x)))
  na <- tibble(clinical = v, n = nrow(d), slope = NA_real_, se = NA_real_, p_lm = NA_real_,
               rho = NA_real_, p_spear = NA_real_)
  if (nrow(d) < 10 || sd(d$log2_LFQ) == 0 || sd(d[[v]]) == 0) return(na)
  co <- summary(lm(reformulate(v, "log2_LFQ"), data = d))$coefficients
  sp <- suppressWarnings(cor.test(d$log2_LFQ, d[[v]], method = "spearman"))
  out <- tibble(clinical = v, n = nrow(d), slope = co[2, 1], se = co[2, 2], p_lm = co[2, 4],
                rho = unname(sp$estimate), p_spear = sp$p.value)
  if (ADJUST_AGE_SEX) {
    Z <- d %>% transmute(Age, Sex = as.integer(Sex == "M"))
    Z <- Z[, vapply(Z, function(z) sd(z) > 0, logical(1)), drop = FALSE]
    ca <- summary(lm(reformulate(c(v, names(Z)), "log2_LFQ"), data = bind_cols(d["log2_LFQ"], d[v], Z)))$coefficients
    ps <- partial_spearman(d$log2_LFQ, d[[v]], Z)
    out <- out %>% mutate(slope_adj = ca[v, 1], p_lm_adj = ca[v, 4],
                          rho_adj = ps[["rho"]], p_spear_adj = ps[["p"]])
  }
  out
}

assoc <- adata %>%
  group_by(Group, Genes) %>%
  group_modify(~ map_dfr(clinical_vars, fit_one, d = .x)) %>%
  ungroup() %>%
  group_by(Group, clinical) %>%                      # FDR within each phase x clinical variable
  mutate(fdr_lm = p.adjust(p_lm, "BH"), fdr_spear = p.adjust(p_spear, "BH"),
         across(any_of(c("p_lm_adj", "p_spear_adj")), ~ p.adjust(.x, "BH"),
                .names = "fdr_{sub('p_', '', .col)}")) %>%
  ungroup()
write_csv(assoc, file.path(OUT_DIR, "per_gene_association_AP_RP.csv"))

##Primary statistic for ranking/plots: adjusted if requested, else marginal

RHO <- if (ADJUST_AGE_SEX) "rho_adj" else "rho"
FDR <- if (ADJUST_AGE_SEX) "fdr_spear_adj" else "fdr_spear"
PV  <- if (ADJUST_AGE_SEX) "p_spear_adj" else "p_spear"

##PART 3:RANK VARIABLES & DRIVER PROTEINS, per phase

var_ranking <- assoc %>%
  group_by(Group, clinical) %>%
  summarise(n_sig_fdr05   = sum(.data[[FDR]] < 0.05, na.rm = TRUE),
            median_absrho = median(abs(.data[[RHO]]), na.rm = TRUE),
            max_absrho    = max(abs(.data[[RHO]]), na.rm = TRUE), .groups = "drop") %>%
  arrange(Group, desc(n_sig_fdr05), desc(median_absrho))

print(var_ranking)

write_csv(var_ranking, file.path(OUT_DIR, "variable_ranking.csv"))

drivers <- assoc %>%
  filter(.data[[FDR]] < 0.05) %>%
  group_by(Group, clinical) %>%
  arrange(desc(abs(.data[[RHO]])), .by_group = TRUE) %>%
  mutate(direction = if_else(.data[[RHO]] > 0, "up", "down")) %>%
  slice_head(n = 15) %>% ungroup()

write_csv(drivers, file.path(OUT_DIR, "driver_genes_top15.csv"))

##PART 4:CATEGORICAL OUTCOMES (two clean two-group contrasts per phase)

wilcox_one <- function(d, fct) {
  d <- d %>% filter(!is.na(.data[[fct]]), !is.na(log2_LFQ))
  if (n_distinct(d[[fct]]) < 2) return(tibble(factor = fct, n = nrow(d), p = NA_real_, delta = NA_real_))
  lv <- sort(unique(d[[fct]]))
  tibble(factor = fct, n = nrow(d),
         delta = median(d$log2_LFQ[d[[fct]] == lv[2]]) - median(d$log2_LFQ[d[[fct]] == lv[1]]),
         contrast = paste(lv[2], "-", lv[1]),
         p = suppressWarnings(wilcox.test(d$log2_LFQ ~ factor(d[[fct]]), exact = FALSE)$p.value))
}
class_res <- adata %>%
  group_by(Group, Genes) %>%
  group_modify(~ map_dfr(class_vars, wilcox_one, d = .x)) %>%
  ungroup() %>% group_by(Group, factor) %>% mutate(fdr = p.adjust(p, "BH")) %>% ungroup()

write_csv(class_res, file.path(OUT_DIR, "classification_diffexpr.csv"))

## AP vs RP, paired within patient

group_res <- expr %>%
  pivot_wider(names_from = Group, values_from = log2_LFQ) %>%
  filter(!is.na(AP), !is.na(RP)) %>%
  group_by(Genes) %>%
  summarise(n = n(), delta = mean(RP - AP),
            p = tryCatch(wilcox.test(AP, RP, paired = TRUE, exact = FALSE)$p.value,
                         error = function(e) NA_real_), .groups = "drop") %>%
  mutate(fdr = p.adjust(p, "BH")) %>% arrange(fdr)

write_csv(group_res, file.path(OUT_DIR, "AP_vs_RP_paired.csv"))

##PART 5:DOES THE ASSOCIATION HOLD ACROSS PHASE?
##log2_LFQ ~ var * Group (+ Age + Sex) + (1 | Patient_code), per protein.
##The var:GroupRP term tests whether the slope differs between AP and RP.

interaction_one <- function(d, v) {
  d <- d %>% filter(!is.na(log2_LFQ), !is.na(.data[[v]]), if_all(all_of(adj_vars), ~ !is.na(.x))) %>%
    mutate(Group = factor(Group, levels = c("AP", "RP")))
  na_row <- tibble(clinical = v, term = NA_character_, estimate = NA_real_,
                   std.error = NA_real_, statistic = NA_real_, p.value = NA_real_)
  if (n_distinct(d$Group) < 2 || sd(d[[v]]) == 0 || sum(duplicated(d$Patient_code)) < 3) return(na_row)
  f <- reformulate(c(paste0(v, "*Group"), adj_vars, "(1|Patient_code)"), "log2_LFQ")
  m <- tryCatch(suppressMessages(suppressWarnings(
         lmer(f, data = d, control = lmerControl(calc.derivs = FALSE)))), error = function(e) NULL)
  if (is.null(m)) return(na_row)
  co <- summary(m)$coefficients
  ix <- paste0(v, ":GroupRP")
  if (!ix %in% rownames(co)) return(na_row)
  tibble(clinical = v, term = ix, estimate = co[ix, "Estimate"], std.error = co[ix, "Std. Error"],
         statistic = co[ix, "t value"], p.value = co[ix, "Pr(>|t|)"])
}

interaction_res <- adata %>%
  group_by(Genes) %>%
  group_modify(~ map_dfr(clinical_vars, interaction_one, d = .x)) %>%
  ungroup() %>% group_by(clinical) %>% mutate(fdr = p.adjust(p.value, "BH")) %>% ungroup()
write_csv(interaction_res, file.path(OUT_DIR, "phase_interaction.csv"))
message("Interaction model: ", sum(interaction_res$fdr < 0.05, na.rm = TRUE),
        " protein x variable pairs with FDR < 0.05")

slope_wide <- assoc %>%
  transmute(Group, Genes, clinical, rho = .data[[RHO]], fdr = .data[[FDR]]) %>%
  pivot_wider(names_from = Group, values_from = c(rho, fdr)) %>%
  mutate(sign_concordant = sign(rho_AP) == sign(rho_RP),
         holds_both = fdr_AP < 0.05 & fdr_RP < 0.05 & sign_concordant,
         flips = (fdr_AP < 0.05 | fdr_RP < 0.05) & sign(rho_AP) != sign(rho_RP))
write_csv(slope_wide, file.path(OUT_DIR, "AP_vs_RP_slope_concordance.csv"))

##PART 6:PLOTS

rho_lab <- if (ADJUST_AGE_SEX) "Partial Spearman rho (adjusted for age and sex)" else "Spearman rho"

p_vol <- assoc %>% filter(!is.na(.data[[RHO]])) %>%
  ggplot(aes(.data[[RHO]], -log10(.data[[PV]]), colour = .data[[FDR]] < 0.05)) +
  geom_point(alpha = .5, size = .8) + facet_grid(Group ~ clinical) +
  scale_colour_manual(values = c(`FALSE` = grey_col, `TRUE` = sig_col)) +
  labs(x = rho_lab, y = expression(-log[10]~p), colour = "FDR < 0.05",
       caption = "Viral_burden = -Ct: positive rho = higher with more virus") + theme_bw()
save_all(p_vol, file.path(OUT_DIR, "volcano_by_group"), width = 12, height = 5)

p_conc <- slope_wide %>% filter(!is.na(rho_AP), !is.na(rho_RP)) %>%
  ggplot(aes(rho_AP, rho_RP)) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, linetype = 2, colour = "grey60") +
  geom_point(aes(colour = holds_both), alpha = .5, size = .8) + facet_wrap(~ clinical) +
  scale_colour_manual(values = c(`FALSE` = grey_col, `TRUE` = pos_col), na.value = "grey90") +
  labs(x = "rho in AP", y = "rho in RP", colour = "sig. & concordant\nin both phases") + theme_bw()
save_all(p_conc, file.path(OUT_DIR, "AP_vs_RP_concordance"), width = 9, height = 6)

sig_hits <- assoc %>% filter(.data[[FDR]] < 0.05) %>%
  transmute(Group, Genes, clinical, rho = .data[[RHO]], fdr = .data[[FDR]])
if (nrow(sig_hits) > 0) {
  sig_pts <- adata %>%
    filter(Genes %in% unique(sig_hits$Genes)) %>%
    pivot_longer(all_of(clinical_vars), names_to = "clinical", values_to = "xval") %>%
    transmute(Patient_code, Group, Genes, clinical, xval, yval = log2_LFQ) %>%
    inner_join(sig_hits, by = c("Group", "Genes", "clinical")) %>%
    filter(!is.na(xval), !is.na(yval)) %>%
    mutate(panel = paste0(Genes, "  (", Group, " \u00b7 ", clinical, ")"),
           dir = if_else(rho > 0, "positive", "negative"))
  write_csv(sig_pts, file.path(OUT_DIR, "FDR_significant_genes_points.csv"))
 
##  PART 7:VARIANCE PARTITION

if (requireNamespace("variancePartition", quietly = TRUE)) {
  suppressPackageStartupMessages(library(variancePartition))
  cont_vars <- c("Day_of_fever", "Platelet_count_D0", "Hematocrit_D0", "Viral_burden",
                 if (ADJUST_AGE_SEX) "Age")
  cat_vars  <- c("Vascular_leakage_binary", "Severity", if (ADJUST_AGE_SEX) "Sex")

  build_group <- function(grp) {
    M <- expr %>% filter(Group == grp) %>% select(Patient_code, Genes, log2_LFQ) %>%
      pivot_wider(names_from = Patient_code, values_from = log2_LFQ) %>%
      column_to_rownames("Genes") %>% as.matrix()
    M <- M[rowMeans(!is.na(M)) == 1, , drop = FALSE]
    meta <- pheno %>% filter(Patient_code %in% colnames(M)) %>% as.data.frame()
    rownames(meta) <- meta$Patient_code

    ## complete cases on every variable in the model -> same patients for every protein
    keep <- complete.cases(meta[, c(cont_vars, cat_vars)])
    message("  variancePartition ", grp, ": ", sum(keep), " of ", nrow(meta),
            " patients have complete clinical data")
    meta <- meta[keep, , drop = FALSE]; M <- M[, rownames(meta), drop = FALSE]
    for (v in cat_vars) meta[[v]] <- droplevels(as.factor(meta[[v]]))
    list(M = M, meta = meta)
  }
  
vp_formula <- function(meta) {
    cont <- cont_vars[vapply(cont_vars, function(v) sd(meta[[v]]) > 0, logical(1))]
    cat  <- cat_vars[vapply(cat_vars, function(v) nlevels(meta[[v]]) >= 2, logical(1))]
    reformulate(c(cont, sprintf("(1|%s)", cat)))
  }

  vp_by_group <- list()
  for (grp in c("AP", "RP")) {
    gg <- build_group(grp); form <- vp_formula(gg$meta)
    ## collinearity check: terms that are strongly related share variance arbitrarily
    cc <- canCorPairs(form, gg$meta)
    write.csv(round(cc, 2), file.path(OUT_DIR, paste0("variance_partition_collinearity_", grp, ".csv")))
    vp <- suppressMessages(fitExtractVarPartModel(gg$M, form, gg$meta))
    vp_df <- as.data.frame(vp) %>% rownames_to_column("Genes") %>% arrange(desc(Viral_burden))
    write_csv(vp_df, file.path(OUT_DIR, paste0("variance_partition_", grp, ".csv")))
    vp_by_group[[grp]] <- list(M = gg$M, meta = gg$meta, form = form)
    print(sort(colMeans(vp), decreasing = TRUE))
    save_all(plotVarPart(sortCols(vp)) + ggtitle(paste0("Variance partition - ", grp)),
             file.path(OUT_DIR, paste0("variance_partition_", grp)), width = 9, height = 5)
  }

  boot_varpart <- function(obj, var = "Viral_burden", topN = 20, B = N_BOOT_VP) {
    M <- obj$M; meta <- obj$meta; genes <- rownames(M); n <- ncol(M)
    top_count <- setNames(numeric(length(genes)), genes)
    ve <- matrix(NA_real_, length(genes), B, dimnames = list(genes, NULL))
    for (b in seq_len(B)) {
      idx <- sample(n, n, replace = TRUE)
      Mb <- M[, idx, drop = FALSE]; mb <- meta[idx, , drop = FALSE]
      colnames(Mb) <- rownames(mb) <- paste0("s", seq_len(n))
      vpb <- tryCatch(suppressMessages(suppressWarnings(fitExtractVarPartModel(Mb, obj$form, mb))),
                      error = function(e) NULL)
      if (is.null(vpb) || !var %in% colnames(vpb)) next
      v <- setNames(vpb[[var]], rownames(vpb)); ve[names(v), b] <- v
      tg <- names(sort(v, decreasing = TRUE))[seq_len(min(topN, length(v)))]
      top_count[tg] <- top_count[tg] + 1
    }
    tibble(Genes = genes, ve_mean = rowMeans(ve, na.rm = TRUE),
           ve_lo = apply(ve, 1, quantile, 0.025, na.rm = TRUE),
           ve_hi = apply(ve, 1, quantile, 0.975, na.rm = TRUE),
           topN_freq = top_count / B) %>% arrange(desc(topN_freq), desc(ve_mean))
  }
  write_csv(boot_varpart(vp_by_group[["AP"]]), file.path(OUT_DIR, "variance_partition_AP_bootstrap.csv"))
} else {
  message("variancePartition not installed")
}
message("Analysis complete: Outputs in '", OUT_DIR, "/'.")

