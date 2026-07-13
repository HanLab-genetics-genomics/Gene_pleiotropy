#1 load packages and set variables-------
setwd("/path/to/Gene_pleiotropy-main/Data")
outdir <- "/path/to/Gene_pleiotropy-main/Output/permutation"
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

sapply(c("data.table", "dplyr", "ggplot2", "corrplot", "ggpubr", "Gmisc", "openxlsx", "readxl", "paletteer", "MuMIn",
         "ggsci", "scales", "RColorBrewer", "gridExtra", "tidyr", "stringr", "colorspace", "cowplot",
         "ggbreak", "ggstatsplot", "viridis"), require, character.only = TRUE)
options(warn = -1)


# recheck the gene-level associations with GPS-----------------------
universe_pn <- fread("permutation/universe_pn.csv", na.strings = c("NA", ""))
universe_pm <- fread("permutation/universe_pm.csv", na.strings = c("NA", ""))
head(universe_pn)

# load matched matrices
matched_ctrl_mat_pn <- readRDS(file.path(outdir, "matched_ctrl_mat_pn.rds"))
matched_ctrl_mat_pm <- readRDS(file.path(outdir, "matched_ctrl_mat_pm.rds"))

# reconstruct high/control pools exactly as before
high_pn_genes <- universe_pn[high_pn == TRUE]
ctrl_pn_pool  <- universe_pn[high_pn == FALSE]
ctrl_pn_pool[, ctrl_id := .I]

high_pm_genes <- universe_pm[high_pm == TRUE]
ctrl_pm_pool  <- universe_pm[high_pm == FALSE]
ctrl_pm_pool[, ctrl_id := .I]

matched_feature_test <- function(
    var,
    high_genes,
    control_pool,
    matched_ctrl_mat,
    type = c("binary", "continuous"),
    positive_level = NULL,
    stat = c("mean", "median"),
    direction = c("greater", "less", "two.sided")
) {
  type <- match.arg(type)
  stat <- match.arg(stat)
  direction <- match.arg(direction)
  
  x_high <- high_genes[[var]]
  x_ctrl <- control_pool[[var]]
  
  if (type == "binary") {
    if (is.null(positive_level)) stop("For binary variables, provide positive_level.")
    
    obs_stat <- mean(x_high == positive_level, na.rm = TRUE)
    
    null_values <- matrix(
      x_ctrl[as.vector(matched_ctrl_mat)] == positive_level,
      nrow = nrow(matched_ctrl_mat),
      ncol = ncol(matched_ctrl_mat)
    )
    
    null_stats <- rowMeans(null_values, na.rm = TRUE)
    statistic_name <- "proportion"
  }
  
  if (type == "continuous") {
    if (stat == "median") {
      obs_stat <- median(x_high, na.rm = TRUE)
    } else {
      obs_stat <- mean(x_high, na.rm = TRUE)
    }
    
    null_values <- matrix(
      x_ctrl[as.vector(matched_ctrl_mat)],
      nrow = nrow(matched_ctrl_mat),
      ncol = ncol(matched_ctrl_mat)
    )
    
    if (stat == "median") {
      null_stats <- apply(null_values, 1, median, na.rm = TRUE)
    } else {
      null_stats <- rowMeans(null_values, na.rm = TRUE)
    }
    
    statistic_name <- stat
  }
  
  null_stats <- null_stats[!is.na(null_stats)]
  
  if (direction == "greater") {
    empirical_p <- (sum(null_stats >= obs_stat) + 1) / (length(null_stats) + 1)
  } else if (direction == "less") {
    empirical_p <- (sum(null_stats <= obs_stat) + 1) / (length(null_stats) + 1)
  } else {
    null_center <- mean(null_stats, na.rm = TRUE)
    empirical_p <- (
      sum(abs(null_stats - null_center) >= abs(obs_stat - null_center)) + 1
    ) / (length(null_stats) + 1)
  }
  
  data.table(
    feature = var,
    type = type,
    positive_level = ifelse(is.null(positive_level), NA, positive_level),
    statistic = statistic_name,
    direction = direction,
    n_high = length(x_high),
    n_high_nonmissing = sum(!is.na(x_high)),
    n_perm = length(null_stats),
    observed = obs_stat,
    null_mean = mean(null_stats, na.rm = TRUE),
    null_sd = sd(null_stats, na.rm = TRUE),
    observed_minus_null = obs_stat - mean(null_stats, na.rm = TRUE),
    empirical_p = empirical_p
  )
}

test_specs <- data.table(var = c("Euteleostomi", "Eutheria", "LOEUF_class", "ifhsd2", "ifspecific", "GES_class",
                                 "cis_egene", "cis_sgene", "iftf", "tf_sum", "ogee_connectivity"),
                         type = c("binary", "binary", "binary", "binary", "binary", "binary",
                                  "binary", "binary", "binary", "continuous", "continuous"),
                         positive_level = c("TRUE", "TRUE", "Intolerant", "Duplicates", "Specific", "Essentialty",
                                            "TRUE", "TRUE", "TRUE", NA, NA),
                         stat = c("mean", "mean", "mean", "mean", "mean", "mean",
                                  "mean", "mean", "mean", "median", "median"),
                         direction = c("greater", "less", "greater", "less", "less", "greater",
                                       "greater", "greater", "greater", "greater", "greater"))
test_specs

res_pn <- rbindlist(lapply(seq_len(nrow(test_specs)), function(i) {
  
  positive_level_i <- test_specs$positive_level[i]
  if (is.na(positive_level_i)) positive_level_i <- NULL
  
  matched_feature_test(
    var = test_specs$var[i],
    high_genes = high_pn_genes,
    control_pool = ctrl_pn_pool,
    matched_ctrl_mat = matched_ctrl_mat_pn,
    type = test_specs$type[i],
    positive_level = positive_level_i,
    stat = test_specs$stat[i],
    direction = test_specs$direction[i]
  )
}))

res_pn[, GPS_metric := "GPS-N"]
res_pn

res_pm <- rbindlist(lapply(seq_len(nrow(test_specs)), function(i) {
  
  positive_level_i <- test_specs$positive_level[i]
  if (is.na(positive_level_i)) positive_level_i <- NULL
  
  matched_feature_test(
    var = test_specs$var[i],
    high_genes = high_pm_genes,
    control_pool = ctrl_pm_pool,
    matched_ctrl_mat = matched_ctrl_mat_pm,
    type = test_specs$type[i],
    positive_level = positive_level_i,
    stat = test_specs$stat[i],
    direction = test_specs$direction[i]
  )
}))

res_pm[, GPS_metric := "GPS-M"]
res_pm


# output table------------------
res_all <- rbind(res_pn, res_pm)
head(res_all)
setDT(res_all)

fmt_num <- function(x) {sprintf("%.3f", x)}
fmt_p <- function(p) {ifelse(
  is.na(p), NA_character_,
  ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p)))}

res_all_fmt <- copy(res_all)[feature %in% c("Euteleostomi", "Eutheria", "LOEUF_class", "ifhsd2", "ifspecific", "GES_class",
                                            "cis_egene", "cis_sgene", "iftf", "tf_sum", "ogee_connectivity")]
head(res_all_fmt)

# format table
res_all_fmt <- res_all_fmt[, .(
  `GPS metric` = GPS_metric,
  Feature = rep(c("Euteleostomi", "Eutheria", "Intolerant", "Duplicates", "Tissue-specific", "Essentiality",
                  "cis_egene", "cis_sgene", "TF-encoding genes", "TF interactions", "PPIs"),2),
  Direction = direction,
  `Number of high-GPS genes` = n_high,
  `Number of high-GPS genes with non-missing feature` = n_high_nonmissing,
  `Observed statistic` = fmt_num(observed),
  `Matched-null mean` = fmt_num(null_mean),
  `Observed - null` = fmt_num(observed_minus_null),
  `Empirical P value` = fmt_p(empirical_p)
)]

res_all_fmt

fwrite(res_all_fmt, file = sprintf("%s/permutation_feature_table_for_supp.csv",outdir))
