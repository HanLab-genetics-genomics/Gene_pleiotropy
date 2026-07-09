# ============================================================
# Notes for running matched permutation analyses
# ============================================================

# This script was tested on a Linux server using 18 parallel workers.
# The permutation step can be computationally intensive, especially when
# using the full gene universe and 1,000 matched permutations.

# By default, the script is designed for Linux command-line execution. For
# Windows or RStudio environments, the parallel backend can be switched from
# future::multicore to future::multisession. However, the Windows/RStudio
# configuration has not been extensively tested in this study.

# We recommend running this script as a background job on a Linux server.
# For example:

# nohup Rscript /path/to/Gene_pleiotropy-main/Script/run_matched_permutation.R \
# > /path/to/Gene_pleiotropy-main/Output/permutation/run_matched_permutation.log 2>&1 &

# Please make sure that the output directory exists before running the
# command above.
# ============================================================


setwd("/path/to/Gene_pleiotropy-main/Data")
outdir <- "/path/to/Gene_pleiotropy-main/Data/permutation"
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

library(data.table)
library(future)
library(future.apply)


universe_pn <- fread("/permutation/universe_pn.csv", na.strings = c("NA", ""))
universe_pm <- fread("/permutation/universe_pm.csv", na.strings = c("NA", ""))


make_matched_matrix <- function(universe, high_col, outfile_prefix, n_perm = 1000, workers = 18, seed = 123,
                                parallel_strategy = c("multicore", "multisession", "auto")) {
  
  parallel_strategy <- match.arg(parallel_strategy)
  high_genes <- copy(universe[get(high_col) %in% TRUE])
  control_pool <- copy(universe[get(high_col) %in% FALSE])
  
  control_pool[, ctrl_id := .I]
  
  high_genes[, match_key_1 := paste(chr, length_bin, snp_bin, ld_bin, maf_bin, sep = "|")]
  control_pool[, match_key_1 := paste(chr, length_bin, snp_bin, ld_bin, maf_bin, sep = "|")]
  
  high_genes[, match_key_2 := paste(length_bin, snp_bin, ld_bin, maf_bin, sep = "|")]
  control_pool[, match_key_2 := paste(length_bin, snp_bin, ld_bin, maf_bin, sep = "|")]
  
  high_genes[, match_key_3 := paste(length_bin, snp_bin, ld_bin, sep = "|")]
  control_pool[, match_key_3 := paste(length_bin, snp_bin, ld_bin, sep = "|")]
  
  high_genes[, match_key_4 := paste(length_bin, snp_bin, sep = "|")]
  control_pool[, match_key_4 := paste(length_bin, snp_bin, sep = "|")]
  
  cand_1 <- split(control_pool$ctrl_id, control_pool$match_key_1)
  cand_2 <- split(control_pool$ctrl_id, control_pool$match_key_2)
  cand_3 <- split(control_pool$ctrl_id, control_pool$match_key_3)
  cand_4 <- split(control_pool$ctrl_id, control_pool$match_key_4)
  
  candidate_list <- lapply(seq_len(nrow(high_genes)), function(i) {
    x <- cand_1[[high_genes$match_key_1[i]]]
    level <- "strict_chr_length_snp_ld_maf"
    
    if (is.null(x) || length(x) == 0) {
      x <- cand_2[[high_genes$match_key_2[i]]]
      level <- "length_snp_ld_maf"
    }
    if (is.null(x) || length(x) == 0) {
      x <- cand_3[[high_genes$match_key_3[i]]]
      level <- "length_snp_ld"
    }
    if (is.null(x) || length(x) == 0) {
      x <- cand_4[[high_genes$match_key_4[i]]]
      level <- "length_snp"
    }
    if (is.null(x) || length(x) == 0) {
      warning("No matched controls for high gene: ", high_genes$gene[i],
              "; using full control pool as fallback.")
      x <- control_pool$ctrl_id
      level <- "full_control_pool"
    }
    
    list(ctrl_id = x, match_level = level)
  })
  
  candidate_ids <- lapply(candidate_list, `[[`, "ctrl_id")
  
  high_genes[, match_level := sapply(candidate_list, `[[`, "match_level")]
  high_genes[, n_candidates := sapply(candidate_ids, length)]
  
  fwrite(
    high_genes[, .(gene, match_level, n_candidates)],
    file.path(outdir, paste0(outfile_prefix, "_match_summary.csv"))
  )
  
  print(table(high_genes$match_level))
  print(summary(high_genes$n_candidates))
  
  #options(future.globals.maxSize = 8 * 1024^3)
  #plan(multicore, workers = workers)
  
  options(future.globals.maxSize = 8 * 1024^3)
  
  if (parallel_strategy == "auto") {
    if (.Platform$OS.type == "windows" || Sys.getenv("RSTUDIO") == "1") {
      parallel_strategy <- "multisession"
    } else {
      parallel_strategy <- "multicore"
    }
  }
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  if (parallel_strategy == "multicore") {
    future::plan(future::multicore, workers = workers)
  } else {
    future::plan(future::multisession, workers = workers)
  }
  
  set.seed(seed)
  
  matched_ctrl_mat <- future_sapply(
    candidate_ids,
    function(candidates) {
      sample(candidates, size = n_perm, replace = TRUE)
    },
    future.seed = TRUE
  )
  
  #plan(sequential)
  
  print(dim(matched_ctrl_mat))
  
  saveRDS(
    matched_ctrl_mat,
    file.path(outdir, paste0(outfile_prefix, ".rds"))
  )
  
  invisible(matched_ctrl_mat)
}

matched_ctrl_mat_pn <- make_matched_matrix(
  universe = universe_pn,
  high_col = "high_pn",
  outfile_prefix = "matched_ctrl_mat_pn",
  n_perm = 1000,
  workers = 18,
  seed = 123,
  parallel_strategy = "multicore"
)

matched_ctrl_mat_pm <- make_matched_matrix(
  universe = universe_pm,
  high_col = "high_pm",
  outfile_prefix = "matched_ctrl_mat_pm",
  n_perm = 1000,
  workers = 18,
  seed = 123,
  parallel_strategy = "multicore"
)
