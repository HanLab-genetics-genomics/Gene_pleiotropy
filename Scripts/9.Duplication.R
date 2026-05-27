#1 load packages and set variables-------
setwd("/path/to/Gene_pleiotropy-main/Data")
figure_file <- "/path/to/Gene_pleiotropy-main/Output/Duplication"
if (!dir.exists(figure_file)) {
  dir.create(figure_file, recursive = TRUE)
}

sapply(c("data.table", "dplyr", "ggplot2", "corrplot", "ggpubr", "Gmisc", "openxlsx", "readxl", "paletteer", "MuMIn",
         "ggsci", "scales", "RColorBrewer", "gridExtra", "tidyr", "stringr", "colorspace", "cowplot", 
         "ggbreak", "ggstatsplot", "viridis"), require, character.only = TRUE)
mycolor <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
pleio_color <- lighten(c("#80CDC1", "#85C6E0", "#223E84"), 0.45)
age_color <- lighten(c("#AB6191", "#DA7235", "#7FB65D", "#F2DD33"), 0.45)
age_color7 <- lighten(mycolor, 0.45)
age7_name <- c("Euteleostomi", "Tetrapoda", "Amniota", "Mammalia", "Theria", "Eutheria", "Primate")
x_text <-  c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria")
age_test <- c("Euteleostomi", "Eutheria")
options(warn = -1)

################Duplication (ancestral-higher vs ancestral-lower group)####################--------------------
#1 diff by pleio class------------------
load("Duplication/progenitor_data.RData")
progenitor_wide_pn <- merge(progenitor_data$pn_ld[progenitor == "Progenitor", .(use_score, group_id, gene, progenitor, pleio10_num = as.numeric(pleio10))],
                            progenitor_data$pn_ld[progenitor == "Copy", .(use_score, group_id, gene, progenitor, pleio10_num = as.numeric(pleio10))],
                            by = "group_id", all.x = T, all.y = T) %>%
  .[!is.na(use_score.x),] %>%
  .[, diff := use_score.x - use_score.y] %>%
  .[, diff_pleio10 := pleio10_num.x - pleio10_num.y] 


progenitor_wide_pm <- merge(progenitor_data$pm_ld[progenitor == "Progenitor", .(use_score, group_id, gene, progenitor, pleio10_num = as.numeric(pleio10))],
                            progenitor_data$pm_ld[progenitor == "Copy", .(use_score, group_id, gene, progenitor, pleio10_num = as.numeric(pleio10))],
                            by = "group_id", all.x = T, all.y = T) %>%
  .[!is.na(use_score.x),] %>%
  .[, diff := use_score.x - use_score.y] %>%
  .[, diff_pleio10 := pleio10_num.x - pleio10_num.y] 

head(progenitor_wide_pm)
nrow(progenitor_wide_pm) #1211 paralog pairs


#2 relative sizes of these two subsets------------------------
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggstatsplot)

plot_pairclass_bar <- function(dt,
                               metric_name = "GPS-N",
                               diff_col = "diff_pleio10",
                               use_percent = TRUE,
                               adjust_method = "BH") {
  
  dt2 <- copy(as.data.table(dt))
  
  ## define 3 groups
  dt2[, pair_class := fifelse(
    get(diff_col) > 0, "progenitor_higher",
    fifelse(get(diff_col) < 0, "progenitor_lower", "tie")
  )]
  
  class_levels <- c("progenitor_higher", "progenitor_lower", "tie")
  class_labels <- c("Ancestral\nhigher", "Ancestral\nlower", "Tie")
  
  dt2[, pair_class := factor(pair_class, levels = class_levels)]
  
  ## count + percentage
  plot_data <- dt2[, .N, by = pair_class][
    data.table(pair_class = class_levels), on = "pair_class"
  ]
  plot_data[is.na(N), N := 0]
  plot_data[, pct := 100 * N / sum(N)]
  plot_data[, label := sprintf("%.1f", pct)]
  
  ## y for plot
  if (use_percent) {
    plot_data[, y := pct]
    ylab_use <- "Percentage of paralog pairs (%)"
  } else {
    plot_data[, y := N]
    ylab_use <- "Count"
  }
  
  ## pairwise comparisons
  pair_list <- list(
    c("progenitor_higher", "progenitor_lower"),
    c("progenitor_higher", "tie"),
    c("progenitor_lower", "tie")
  )
  
  stat_test <- rbindlist(lapply(seq_along(pair_list), function(i) {
    g1 <- pair_list[[i]][1]
    g2 <- pair_list[[i]][2]
    
    n1 <- plot_data[pair_class == g1, N]
    n2 <- plot_data[pair_class == g2, N]
    n_pair <- n1 + n2
    
    p_raw <- if (n_pair > 0) {
      binom.test(n1, n_pair, p = 0.5)$p.value
    } else {
      NA_real_
    }
    
    data.table(
      group1 = g1,
      group2 = g2,
      xmin = g1,
      xmax = g2,
      p = p_raw
    )
  }))
  
  stat_test[, p_adj := p.adjust(p, method = adjust_method)]
  stat_test[, p_label := formatC(p_adj, format = "e", digits = 2)]
  
  ## set y positions
  ymax <- max(plot_data$y)
  step_increase <- ymax * 0.12
  if (step_increase == 0) step_increase <- 1
  
  stat_test[, y.position := ymax + c(1, 2, 3) * step_increase]
  
  ## optional overall test
  overall_p <- chisq.test(plot_data$N, p = rep(1/3, 3))$p.value
  
  ## plot
  p <- ggplot(plot_data, aes(x = pair_class, y = y, fill = pair_class)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = c("progenitor_higher" = "#4DBBD5FF", "progenitor_lower" = "#E64B35FF", "tie" =  "#8DA0CB")) +
    geom_text(aes(label = label), vjust = -0.3, size = 5, fontface = "bold") +
    stat_pvalue_manual(
      stat_test,
      label = "p_label",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      size = 4.5
    ) +
    scale_x_discrete(labels = class_labels) +
    labs(title = metric_name, x = "", y = ylab_use) +
    theme_bw() +
    theme(
      text = element_text(size = 13, face = "bold", color = "black"),
      axis.text = element_text(size = 12, face = "bold", color = "black"),
      axis.title = element_text(size = 15, face = "bold", color = "black"),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "grey40")
    )
  
  return(list(plot_data = plot_data,
              stat_test = stat_test,
              plot = p))
}


res_pn <- plot_pairclass_bar(dt = progenitor_wide_pn, metric_name = "GPS-N",
                             diff_col = "diff_pleio10", use_percent = TRUE)
res_pn$plot$layers <- res_pn$plot$layers[
  !sapply(res_pn$plot$layers, function(x) inherits(x$geom, "GeomBracket"))
]
res_pn$plot
res_pn$stat_test

res_pm <- plot_pairclass_bar(dt = progenitor_wide_pm, metric_name = "GPS-M",
                             diff_col = "diff_pleio10", use_percent = TRUE)
res_pm$plot$layers <- res_pm$plot$layers[
  !sapply(res_pm$plot$layers, function(x) inherits(x$geom, "GeomBracket"))
]
res_pm$plot
res_pm$stat_test

ggsave(plot = res_pn$plot, width = 5, height = 6, device = cairo_pdf,
       filename = sprintf("%s/progenitor_higher_pn.pdf", figure_file))
ggsave(plot = res_pm$plot, width = 5, height = 6, device = cairo_pdf,
       filename = sprintf("%s/progenitor_higher_pm.pdf", figure_file))


#3 compare general gene-level features between the two subsets------------------
##3.1 manage data----------------
load("pleiotropy_maindata.RData")
head(pleiotropy_maindata$pm_ld)
feature_cols_raw <- c("Gene length (bp)", "CDS Length", "Exon Counts",
                      "Number of SNPs (Gene)", "Nonsense/Synonymous ratio",
                      "Transcript count", "GC content", "Missense/Synonymous ratio", "LOEUF", "pLI")
cols_use <- c("gene", feature_cols_raw)
feature_data <- pleiotropy_maindata$pm_ld[, ..cols_use]
head(feature_data)

gene_features <- feature_data[, `:=`(
  progenitor_higher_pro_pm  = gene %in% progenitor_wide_pm[diff_pleio10 > 0, gene.x],
  progenitor_lower_pro_pm   = gene %in% progenitor_wide_pm[diff_pleio10 < 0, gene.x],
  progenitor_higher_copy_pm = gene %in% progenitor_wide_pm[diff_pleio10 > 0, gene.y],
  progenitor_lower_copy_pm  = gene %in% progenitor_wide_pm[diff_pleio10 < 0, gene.y],
  
  progenitor_higher_pro_pn  = gene %in% progenitor_wide_pn[diff_pleio10 > 0, gene.x],
  progenitor_lower_pro_pn   = gene %in% progenitor_wide_pn[diff_pleio10 < 0, gene.x],
  progenitor_higher_copy_pn = gene %in% progenitor_wide_pn[diff_pleio10 > 0, gene.y],
  progenitor_lower_copy_pn  = gene %in% progenitor_wide_pn[diff_pleio10 < 0, gene.y]
)]

head(gene_features)

##3.2 compare gene metrics between two groups-------------
feature_progenitor_pn <- rbindlist(lapply(feature_cols_raw, function(f) {
  x <- gene_features[progenitor_higher_pro_pn == T , get(f)]
  y <- gene_features[progenitor_lower_pro_pn == T , get(f)]
  
  wt <- wilcox.test(x, y)
  
  data.table(feature = f, n_high = sum(!is.na(x)), n_low  = sum(!is.na(y)),
             median_high = median(x, na.rm = TRUE),
             median_low  = median(y, na.rm = TRUE), p = wt$p.value)}), fill = TRUE)

feature_progenitor_pn[, p_adj := p.adjust(p, method = "BH")]
feature_progenitor_pn[order(p_adj)]

feature_progenitor_pm <- rbindlist(lapply(feature_cols_raw, function(f) {
  x <- gene_features[progenitor_higher_pro_pm == T , get(f)]
  y <- gene_features[progenitor_lower_pro_pm == T , get(f)]
  
  wt <- wilcox.test(x, y)
  
  data.table(feature = f, n_high = sum(!is.na(x)), n_low  = sum(!is.na(y)),
             median_high = median(x, na.rm = TRUE),
             median_low  = median(y, na.rm = TRUE), p = wt$p.value)}), fill = TRUE)

feature_progenitor_pm[, p_adj := p.adjust(p, method = "BH")]
feature_progenitor_pm[order(p_adj)]

##compare copy between two groups
feature_copy_pn <- rbindlist(lapply(feature_cols_raw, function(f) {
  x <- gene_features[progenitor_higher_copy_pn == T, get(f)]
  y <- gene_features[progenitor_lower_copy_pn == T , get(f)]
  
  wt <- wilcox.test(x, y)
  
  data.table(feature = f, n_high = sum(!is.na(x)), n_low  = sum(!is.na(y)),
             median_high = median(x, na.rm = TRUE), 
             median_low  = median(y, na.rm = TRUE),p = wt$p.value)}), fill = TRUE)

feature_copy_pn[, p_adj := p.adjust(p, method = "fdr")]
feature_copy_pn[order(p_adj)]


feature_copy_pm <- rbindlist(lapply(feature_cols_raw, function(f) {
  x <- gene_features[progenitor_higher_copy_pm == T, get(f)]
  y <- gene_features[progenitor_lower_copy_pm == T, get(f)]
  
  wt <- wilcox.test(x, y)
  
  data.table(feature = f, n_high = sum(!is.na(x)), n_low  = sum(!is.na(y)),
             median_high = median(x, na.rm = TRUE), 
             median_low  = median(y, na.rm = TRUE),p = wt$p.value)}), fill = TRUE)

feature_copy_pm[, p_adj := p.adjust(p, method = "fdr")]
feature_copy_pm[order(p_adj)]

##compare all between two groups
feature_all_pn <- rbindlist(lapply(feature_cols_raw, function(f) {
  x <- gene_features[progenitor_higher_pro_pn == T | progenitor_higher_copy_pn == T, get(f)]
  y <- gene_features[progenitor_lower_pro_pn == T | progenitor_lower_copy_pn == T , get(f)]
  
  wt <- wilcox.test(x, y)
  
  data.table(feature = f, n_high = sum(!is.na(x)), n_low  = sum(!is.na(y)),
             median_high = median(x, na.rm = TRUE), 
             median_low  = median(y, na.rm = TRUE),p = wt$p.value)}), fill = TRUE)

feature_all_pn[, p_adj := p.adjust(p, method = "fdr")]
feature_all_pn[order(p_adj)]


feature_all_pm <- rbindlist(lapply(feature_cols_raw, function(f) {
  x <- gene_features[progenitor_higher_pro_pm == T | progenitor_higher_copy_pm == T, get(f)]
  y <- gene_features[progenitor_lower_pro_pm == T | progenitor_lower_copy_pm == T, get(f)]
  
  wt <- wilcox.test(x, y)
  
  data.table(feature = f, n_high = sum(!is.na(x)), n_low  = sum(!is.na(y)),
             median_high = median(x, na.rm = TRUE), 
             median_low  = median(y, na.rm = TRUE),p = wt$p.value)}), fill = TRUE)

feature_all_pm[, p_adj := p.adjust(p, method = "fdr")]
feature_all_pm[order(p_adj)]


feature_all <- rbindlist(list(feature_all_pm = feature_all_pm[, class := "feature_all_pm"],
                              feature_all_pn = feature_all_pn[, class := "feature_all_pn"],
                              
                              feature_progenitor_pm = feature_progenitor_pm[, class := "feature_progenitor_pm"],
                              feature_progenitor_pn = feature_progenitor_pn[, class := "feature_progenitor_pn"],
                              
                              feature_copy_pm = feature_copy_pm[, class := "feature_copy_pm"],
                              feature_copy_pn = feature_copy_pn[, class := "feature_copy_pn"]))

head(feature_all)

feature_all_output <- feature_all[, .(feature, n_high, n_low, 
                                      median_high_format = formatC(median_high, format = "f", digits = 2, drop0trailing = TRUE), 
                                      median_low_format = formatC(median_low, format = "f", digits = 2, drop0trailing = TRUE),
                                      P = ifelse(p > 0.001, sprintf("%.3f", p), sprintf("%.2e", p)), 
                                      adjusted_P = ifelse(p_adj > 0.001, sprintf("%.3f", p_adj), sprintf("%.2e", p_adj)), 
                                      comparison = ifelse(startsWith(class, "feature_all"), "All-gene comparison",
                                                          ifelse(startsWith(class, "feature_progenitor"), "Ancestral-only comparison", "Recent-only comparison")),
                                      GPS = ifelse(endsWith(class, "pm"), "GPS-M", "GPS-N"))]
feature_all_output

fwrite(feature_all_output, file = sprintf("%s/feature_progenitor_pleio10.csv", figure_file))


################Duplication (Ensembl Compara)####################--------------------
#1 load data------------------------
paralog_tbl <- fread("Duplication/biomaRt_paralogs.txt")
names(paralog_tbl) <- c("gene_id", "gene_name",
                        "paralog_gene_id", "paralog_gene_name",
                        "lca", "pid_target", "pid_query")
head(paralog_tbl) #3560494
para <- paralog_tbl[!is.na(gene_name) & gene_name != "" &
                      !is.na(paralog_gene_name) & paralog_gene_name != "" &
                      gene_name != paralog_gene_name]
nrow(para) #2945162
head(para)

recent_lca_list <- list(recent_human = c("Homo sapiens"),
                        recent_primate = c( "Homo sapiens", "Homininae", "Hominidae", "Hominoidea", "Catarrhini", "Simiiformes", "Haplorrhini", "Primates"),
                        recent_mammal = c("Homo sapiens", "Homininae", "Hominidae", "Hominoidea", "Catarrhini", "Simiiformes", "Haplorrhini", "Primates",
                                          "Euarchontoglires", "Boreoeutheria", "Eutheria", "Theria", "Mammalia"),
                        ancient_vertebrate_or_older = c("Amniota", "Tetrapoda", "Sarcopterygii", "Euteleostomi",
                                                        "Gnathostomata", "Vertebrata", "Chordata", "Bilateria", "Opisthokonta"))

dupgene_list <- list(ens_dup_human  = para[lca %in% recent_lca_list$recent_human, unique(na.omit(c(gene_name, paralog_gene_name)))],
                     ens_dup_primate = para[lca %in% recent_lca_list$recent_primate, unique(na.omit(c(gene_name, paralog_gene_name)))],
                     ens_dup_mammal = para[lca %in% recent_lca_list$recent_mammal, unique(na.omit(c(gene_name, paralog_gene_name)))],
                     ens_dup_vertebrate_older = para[lca %in% recent_lca_list$ancient_vertebrate_or_older, unique(na.omit(c(gene_name, paralog_gene_name)))],
                     ens_dup_all = para[, unique(na.omit(c(gene_name, paralog_gene_name)))])


##merge with GPS data 
GPS_feature_all <- fread("GPS_feature_all.csv", na.strings = c("", "NA"))
head(GPS_feature_all)

GPS_hsd <- copy(GPS_feature_all)
GPS_hsd[, ens_dup_human := gene %in% dupgene_list$ens_dup_human]
GPS_hsd[, ens_dup_primate := gene %in% dupgene_list$ens_dup_primate]
GPS_hsd[, ens_dup_mammal := gene %in% dupgene_list$ens_dup_mammal]
GPS_hsd[, ens_dup_vertebrate_older := gene %in% dupgene_list$ens_dup_vertebrate_older]
GPS_hsd[, ens_dup_all := gene %in% dupgene_list$ens_dup_all]
table(GPS_hsd$ens_dup_human)
head(GPS_hsd)

#2 compare GPS between singletons and duplicated genes-------------------
library(data.table)
summarize_dup_wilcox <- function(dt, flag_cols = c("ens_dup_human", "ens_dup_primate", "ens_dup_mammal", 
                                                   "ens_dup_vertebrate_older", "ens_dup_all"),
                                 value_cols = c("pm_ld", "pn_ld"),
                                 exact = FALSE) {
  dt <- as.data.table(copy(dt))
  
  # check columns
  stopifnot(all(flag_cols %in% names(dt)))
  stopifnot(all(value_cols %in% names(dt)))
  
  res <- rbindlist(lapply(flag_cols, function(flag_col) {
    flag <- as.logical(dt[[flag_col]])
    
    rbindlist(lapply(value_cols, function(val_col) {
      x <- dt[flag == TRUE,  get(val_col)]
      y <- dt[flag == FALSE, get(val_col)]
      
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      
      wt <- tryCatch(
        wilcox.test(x = x, y = y, exact = exact),
        error = function(e) NULL
      )
      
      data.table(
        dup_def        = flag_col,
        variable       = val_col,
        n_dup          = length(x),
        n_nondup       = length(y),
        mean_dup       = if (length(x) > 0) mean(x) else NA_real_,
        mean_nondup    = if (length(y) > 0) mean(y) else NA_real_,
        median_dup     = if (length(x) > 0) median(x) else NA_real_,
        median_nondup  = if (length(y) > 0) median(y) else NA_real_,
        p_value        = if (!is.null(wt)) wt$p.value else NA_real_
      )
    }), use.names = TRUE)
  }), use.names = TRUE)
  
  res[, p_adj_BH := p.adjust(p_value, method = "BH")]
  setcolorder(res, c(
    "dup_def", "variable",
    "n_dup", "n_nondup",
    "mean_dup", "mean_nondup",
    "median_dup", "median_nondup",
    "p_value", "p_adj_BH"
  ))
  
  return(res[])
}
res_dup <- summarize_dup_wilcox(GPS_hsd)
res_dup


##output table
format_dup_result <- function(res, variable_labels = c(pm_ld = "GPS-M", pn_ld = "GPS-N"),
                              dup_labels = c(ens_dup_human = "Ensembl paralogs: Homo sapiens",
                                             ens_dup_primate = "Ensembl paralogs: Primates or younger",
                                             ens_dup_mammal = "Ensembl paralogs: Mammalia or younger",
                                             ens_dup_vertebrate_older = "Ensembl paralogs: Vertebrata or older",
                                             ens_dup_all = "Ensembl paralogs: all LCA levels")) {
  
  dt <- copy(as.data.table(res))
  
  dt[, score := variable_labels[variable]]
  dt[is.na(score), score := variable]
  
  dt[, duplication_definition := dup_labels[dup_def]]
  dt[is.na(duplication_definition), duplication_definition := dup_def]
  
  # optional ordering
  dt[, score := factor(score, levels = unique(variable_labels))]
  dt[, duplication_definition := factor(duplication_definition, levels = unique(dup_labels))]
  
  setorder(dt, duplication_definition, score)
  
  out <- dt[, .(duplication_definition = as.character(duplication_definition),
                score = as.character(score), n_dup, n_nondup,
                median_dup = sprintf("%.2f", median_dup),
                median_nondup = sprintf("%.2f", median_nondup),
                p_value = sprintf("%.2e", p_value),
                p_adj_BH = sprintf("%.2e", p_adj_BH))]
  
  return(out)
}

res_dup_formatted <- format_dup_result(res_dup)
res_dup_formatted

fwrite(res_dup_formatted, file = sprintf("%s/gpscompare_ensembl_dup_lca.csv", figure_file))

#3 plot data-----------------
plot_dup_violin_box2 <- function(dt,
                                 flag_cols = c("ens_dup_human", "ens_dup_primate", "ens_dup_mammal", "ens_dup_vertebrate_older", "ens_dup_all"),
                                 value_cols = c("pm_ld", "pn_ld"),
                                 dup_labels = c(ens_dup_human = "Homo sapiens",
                                                ens_dup_primate = "Primates or younger",
                                                ens_dup_mammal = "Mammalia or younger",
                                                ens_dup_vertebrate_older = "Vertebrata or older",
                                                ens_dup_all = "all LCA levels"),
                                 value_labels = c(pm_ld = "GPS-M", pn_ld = "GPS-N"),
                                 exact = FALSE) {
  
  dt <- as.data.table(copy(dt))
  
  long_dt <- rbindlist(lapply(flag_cols, function(flag_col) {
    rbindlist(lapply(value_cols, function(val_col) {
      data.table(
        dup_def_raw  = flag_col,
        variable_raw = val_col,
        group        = fifelse(as.logical(dt[[flag_col]]), "Duplicates", "Singletons"),
        value        = dt[[val_col]]
      )
    }))
  }), use.names = TRUE)
  
  long_dt <- long_dt[!is.na(value)]
  
  summary_dt <- long_dt[, .(
    n = .N,
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE)
  ), by = .(dup_def_raw, variable_raw, group)]
  
  test_dt <- rbindlist(lapply(flag_cols, function(flag_col) {
    rbindlist(lapply(value_cols, function(val_col) {
      x <- dt[get(flag_col) == TRUE,  get(val_col)]
      y <- dt[get(flag_col) == FALSE, get(val_col)]
      x <- x[!is.na(x)]
      y <- y[!is.na(y)]
      wt <- wilcox.test(x, y, exact = exact)
      data.table(
        dup_def_raw  = flag_col,
        variable_raw = val_col,
        p_value      = wt$p.value
      )
    }))
  }), use.names = TRUE)
  
  stats_wide <- merge(
    summary_dt[group == "Duplicates"][, .(
      dup_def_raw, variable_raw,
      n_dup = n, mean_dup = mean, median_dup = median
    )],
    summary_dt[group == "Singletons"][, .(
      dup_def_raw, variable_raw,
      n_singleton = n, mean_singleton = mean, median_singleton = median
    )],
    by = c("dup_def_raw", "variable_raw"),
    all = TRUE
  )
  
  stats_wide <- merge(stats_wide, test_dt, by = c("dup_def_raw", "variable_raw"), all = TRUE)
  
  long_dt[, dup_def := factor(dup_def_raw, levels = flag_cols, labels = dup_labels[flag_cols])]
  long_dt[, variable := factor(variable_raw, levels = value_cols, labels = value_labels[value_cols])]
  long_dt[, group := factor(group, levels = c("Singletons", "Duplicates"))]
  
  test_dt[, dup_def := factor(dup_def_raw, levels = flag_cols, labels = dup_labels[flag_cols])]
  test_dt[, variable := factor(variable_raw, levels = value_cols, labels = value_labels[value_cols])]
  test_dt[, p_label := paste0("P = ", formatC(p_value, format = "e", digits = 2))]
  
  plot_dt <- copy(long_dt)
  plot_dt[variable == "GPS-M" & (value < 0 | value > 10), value := NA_real_]
  
  ann_pos <- rbindlist(list(
    plot_dt[variable == "GPS-M", .(y = 9.6), by = .(dup_def, variable)],
    plot_dt[variable == "GPS-N", {
      rng <- range(value, na.rm = TRUE)
      span <- rng[2] - rng[1]
      if (!is.finite(span) || span == 0) span <- max(abs(rng[2]), 1)
      .(y = rng[2] + 0.10 * span)
    }, by = .(dup_def, variable)]
  ), use.names = TRUE)
  
  ann_dt <- merge(
    ann_pos,
    test_dt[, .(dup_def, variable, p_label)],
    by = c("dup_def", "variable"),
    all.x = TRUE
  )
  ann_dt[, x := 1.5]
  
  ## median labels
  median_lab <- summary_dt[, .(
    dup_def_raw, variable_raw, group,
    median,
    median_label = sprintf("%.2f", median)
  )]
  median_lab[, dup_def := factor(dup_def_raw, levels = flag_cols, labels = dup_labels[flag_cols])]
  median_lab[, variable := factor(variable_raw, levels = value_cols, labels = value_labels[value_cols])]
  median_lab[, group := factor(group, levels = c("Singletons", "Duplicates"))]
  median_lab[, x := c(1, 2)[group]]
  
  median_lab <- merge(
    median_lab,
    plot_dt[, {
      rng <- range(value, na.rm = TRUE)
      span <- rng[2] - rng[1]
      if (!is.finite(span) || span == 0) span <- max(abs(rng[2]), 1)
      .(offset = 0.2 * span)
    }, by = .(dup_def, variable)],
    by = c("dup_def", "variable"),
    all.x = TRUE
  )
  
  median_lab[, y := median + offset]
  median_lab[variable == "GPS-M" & y > 9.8, y := 9.8]
  
  p <- ggplot(plot_dt, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, scale = "width", width = 0.9, alpha = 0.7, color = NA, na.rm = TRUE) +
    geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.35, alpha = 0.95, na.rm = TRUE) +
    stat_summary(fun = median, geom = "point", size = 1.2, color = "black", na.rm = TRUE) +
    geom_text(
      data = median_lab,
      aes(x = x, y = y, label = median_label),
      inherit.aes = FALSE,
      size = 3.5, hjust = -0.5, color = "darkred"
    ) +
    geom_text(
      data = ann_dt,
      aes(x = x, y = y, label = p_label),
      inherit.aes = FALSE,
      size = 3
    ) +
    facet_grid(variable ~ dup_def, scales = "free_y") +
    scale_fill_manual(values = c("Singletons" = "#BDBDBD", "Duplicates" = "#4C78A8")) +
    labs(x = NULL, y = NULL, fill = NULL) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "white", colour = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "top"
    )
  
  return(list(
    plot = p,
    long_data = long_dt,
    plot_data = plot_dt,
    summary_stats = stats_wide
  ))
}

res_plot <- plot_dup_violin_box2(GPS_hsd)
res_plot$plot

ggsave(plot = res_plot$plot, width = 10, height = 5, device = cairo_pdf,
       filename = sprintf("%s/plot_dup_box_ensembl.pdf", figure_file))

#4 proportion across groups---------------------
plot_dup_prop_from_raw <- function(
    data,
    group_var,
    facet_var = "threshold",
    dup_var = "duplicated",
    yes_value = "Yes",
    group_order = NULL,
    group_labels = NULL,
    facet_order = NULL,
    compare_levels = NULL,
    title = NULL,
    ylab = "Proportion of duplicated genes (%)",
    fill_values = NULL,
    free_y = TRUE,
    show_p = TRUE,
    p_test = c("fisher", "chisq"),
    percent_digits = 1,
    p_digits = 1,
    show_bar_label = TRUE,
    bar_label_inside = TRUE,
    bracket_pad = 0.12,
    text_pad = 0.20,
    tick_frac = 0.04,
    bar_width = 0.78,
    text_size = 3.4,
    return_data = FALSE
) {
  
  p_test <- match.arg(p_test)
  dt <- as.data.table(copy(data))
  
  required_cols <- c(group_var, facet_var, dup_var)
  miss_cols <- setdiff(required_cols, names(dt))
  if (length(miss_cols) > 0) {
    stop("Missing columns in data: ", paste(miss_cols, collapse = ", "))
  }
  
  # group order
  if (is.null(group_order)) {
    if (is.factor(dt[[group_var]])) {
      group_order <- levels(dt[[group_var]])
    } else {
      group_order <- unique(as.character(dt[[group_var]]))
    }
  }
  
  # facet order
  if (is.null(facet_order)) {
    if (is.factor(dt[[facet_var]])) {
      facet_order <- levels(dt[[facet_var]])
    } else {
      facet_order <- unique(as.character(dt[[facet_var]]))
    }
  }
  
  # comparison levels for bracket
  if (is.null(compare_levels)) {
    compare_levels <- c(group_order[1], group_order[length(group_order)])
  }
  
  if (!all(compare_levels %in% group_order)) {
    stop("compare_levels must be included in group_order.")
  }
  
  # labels
  if (is.null(group_labels)) {
    group_labels <- setNames(group_order, group_order)
  } else if (is.null(names(group_labels))) {
    group_labels <- setNames(group_labels, group_order)
  }
  
  # convert duplicated variable to logical
  if (is.logical(yes_value)) {
    dt[, .dup_yes := get(dup_var) == yes_value]
  } else {
    dt[, .dup_yes := as.character(get(dup_var)) %in% as.character(yes_value)]
  }
  
  dt[is.na(get(dup_var)), .dup_yes := NA]
  
  dt[, .group := factor(as.character(get(group_var)), levels = group_order)]
  dt[, .facet := factor(as.character(get(facet_var)), levels = facet_order)]
  dt <- dt[!is.na(.group) & !is.na(.facet) & !is.na(.dup_yes)]
  
  # calculate proportion
  prop_dt <- dt[, .(
    n_total = .N,
    n_dup = sum(.dup_yes),
    prop_yes = mean(.dup_yes)
  ), by = .(.facet, .group)]
  
  prop_dt[, .x := as.numeric(.group)]
  
  # format functions
  fmt_pct <- function(x) {
    formatC(x * 100, format = "f", digits = percent_digits)
  }
  
  fmt_p <- function(x) {
    out <- rep(NA_character_, length(x))
    out[is.na(x)] <- NA_character_
    out[!is.na(x) & x == 0] <- "< 2.2e-16"
    out[!is.na(x) & x >= 0.001] <- formatC(
      x[!is.na(x) & x >= 0.001],
      format = "f",
      digits = 3
    )
    out[!is.na(x) & x != 0 & x < 0.001] <- formatC(
      x[!is.na(x) & x != 0 & x < 0.001],
      format = "e",
      digits = p_digits
    )
    out
  }
  
  # p value comparing first and last specified groups within each facet
  ann_dt <- NULL
  
  if (show_p) {
    g1 <- compare_levels[1]
    g2 <- compare_levels[2]
    
    ann_dt <- dt[, {
      sub <- .SD[as.character(.group) %in% c(g1, g2)]
      
      n1_yes <- sum(as.character(sub$.group) == g1 & sub$.dup_yes)
      n1_no  <- sum(as.character(sub$.group) == g1 & !sub$.dup_yes)
      n2_yes <- sum(as.character(sub$.group) == g2 & sub$.dup_yes)
      n2_no  <- sum(as.character(sub$.group) == g2 & !sub$.dup_yes)
      
      tab <- matrix(
        c(n1_yes, n1_no, n2_yes, n2_no),
        nrow = 2,
        byrow = TRUE
      )
      
      pval <- if (all(rowSums(tab) > 0) && all(colSums(tab) > 0)) {
        if (p_test == "fisher") {
          fisher.test(tab)$p.value
        } else {
          suppressWarnings(chisq.test(tab)$p.value)
        }
      } else {
        NA_real_
      }
      
      .(p_value = pval, label = fmt_p(pval))
    }, by = .facet]
    
    y_dt <- prop_dt[, .(
      max_y = max(prop_yes, na.rm = TRUE)
    ), by = .facet]
    
    y_dt[!is.finite(max_y) | max_y <= 0, max_y := 0.01]
    
    ann_dt <- merge(ann_dt, y_dt, by = ".facet", all.x = TRUE)
    
    x1 <- match(g1, group_order)
    x2 <- match(g2, group_order)
    
    ann_dt[, `:=`(
      x1 = x1,
      x2 = x2,
      x_mid = mean(c(x1, x2)),
      y = max_y * (1 + bracket_pad),
      y_text = max_y * (1 + text_pad),
      tick = max_y * tick_frac
    )]
  }
  
  # base plot
  p <- ggplot(prop_dt, aes(x = .x, y = prop_yes, fill = .group)) +
    geom_col(width = bar_width, color = "white", linewidth = 0.4)
  
  if (show_bar_label) {
    p <- p +
      geom_text(
        aes(label = fmt_pct(prop_yes)),
        vjust = ifelse(bar_label_inside, 1.4, -0.25),
        size = 3.4,
        fontface = "bold"
      )
  }
  
  p <- p +
    facet_wrap(
      ~ .facet,
      nrow = 5,
      scales = ifelse(free_y, "free_y", "fixed")
    ) +
    scale_x_continuous(
      breaks = seq_along(group_order),
      labels = unname(group_labels[group_order])
    ) +
    scale_y_continuous(
      labels = function(x) fmt_pct(x),
      expand = expansion(mult = c(0, ifelse(show_p, 0.25, 0.08)))
    ) +
    labs(
      x = "",
      y = ylab,
      title = title
    )
  
  if (!is.null(fill_values)) {
    p <- p + scale_fill_manual(values = fill_values, drop = FALSE)
  }
  
  # add p-value bracket
  if (show_p && !is.null(ann_dt)) {
    p <- p +
      geom_segment(
        data = ann_dt,
        aes(x = x1, xend = x2, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.5
      ) +
      geom_segment(
        data = ann_dt,
        aes(x = x1, xend = x1, y = y, yend = y - tick),
        inherit.aes = FALSE,
        linewidth = 0.5
      ) +
      geom_segment(
        data = ann_dt,
        aes(x = x2, xend = x2, y = y, yend = y - tick),
        inherit.aes = FALSE,
        linewidth = 0.5
      ) +
      geom_text(
        data = ann_dt,
        aes(x = x_mid, y = y_text, label = label),
        inherit.aes = FALSE,
        size = text_size
      )
  }
  
  p <- p +
    theme_bw() +
    theme(
      text = element_text(size = 12, color = "black", face = "bold"),
      axis.title = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 12, color = "black", face = "bold"),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "grey50"),
      panel.grid = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      panel.border = element_blank(),
      plot.title = element_text(size = 15, hjust = 0.5),
      strip.background = element_rect(fill = "#F2F2F2", color = "#D9D9D9"),
      strip.text = element_text(
        size = 12,
        face = "bold",
        margin = margin(t = 4, r = 4, b = 4, l = 4)
      )
    )
  
  if (return_data) {
    return(list(plot = p, prop_dt = prop_dt, ann_dt = ann_dt))
  } else {
    return(p)
  }
}


plot_long <- rbindlist(list(
  GPS_hsd[, .(pn_pleio3, pm_pleio3, age_stage4, threshold = "Homo sapiens", duplicated = ens_dup_human)],
  GPS_hsd[, .(pn_pleio3, pm_pleio3, age_stage4, threshold = "Primates or younger", duplicated = ens_dup_primate)],
  GPS_hsd[, .(pn_pleio3, pm_pleio3, age_stage4, threshold = "Mammalia or younger", duplicated = ens_dup_mammal)],
  GPS_hsd[, .(pn_pleio3, pm_pleio3, age_stage4, threshold = "Vertebrata or older", duplicated = ens_dup_vertebrate_older)],
  GPS_hsd[, .(pn_pleio3, pm_pleio3, age_stage4, threshold = "all LCA levels", duplicated = ens_dup_all)]
))

plot_long[, duplicated := ifelse(duplicated, "Yes", "No")]
head(plot_long)

plot_pn <- plot_dup_prop_from_raw(data = plot_long, group_var = "pn_pleio3", facet_var = "threshold", bar_label_inside = F,
                                  group_order = pleio_order <- c("Low pleiotropy", "Intermediate pleiotropy", "High pleiotropy"),
                                  group_labels = c("Low pleiotropy" = "Low", "Intermediate pleiotropy" = "Intermediate", "High pleiotropy" = "High"),
                                  facet_order = c("Homo sapiens", "Primates or younger", "Mammalia or younger", "Vertebrata or older", "all LCA levels"),
                                  compare_levels = c("Low pleiotropy", "High pleiotropy"), title = "GPS-N",  ylab = "Proportion of duplicated genes (%)",
                                  fill_values = c("Low pleiotropy" = "#9EEADD", "Intermediate pleiotropy" = "#A4E3FE", "High pleiotropy" = "#7F8FCF"), free_y = TRUE)

plot_pn

plot_pm <- plot_dup_prop_from_raw(data = plot_long, group_var = "pm_pleio3", facet_var = "threshold", bar_label_inside = F,
                                  group_order = c("Low pleiotropy", "Intermediate pleiotropy", "High pleiotropy"),
                                  group_labels = c("Low pleiotropy" = "Low", "Intermediate pleiotropy" = "Intermediate", "High pleiotropy" = "High"),
                                  facet_order = c("Homo sapiens", "Primates or younger", "Mammalia or younger", "Vertebrata or older", "all LCA levels"),
                                  compare_levels = c("Low pleiotropy", "High pleiotropy"), title = "GPS-N",  ylab = "Proportion of duplicated genes (%)",
                                  fill_values = c("Low pleiotropy" = "#9EEADD", "Intermediate pleiotropy" = "#A4E3FE", "High pleiotropy" = "#7F8FCF"), free_y = TRUE)

plot_pm

plot_age <- plot_dup_prop_from_raw(data = plot_long, group_var = "age_stage4", facet_var = "threshold", bar_label_inside = F,
                                   group_order = c("Euteleostomi", "Tetrapoda", "Amniota" ,"Eutheria"),
                                   group_labels = c("Euteleostomi", "Tetrapoda", "Amniota" ,"Eutheria"),
                                   facet_order = c("Homo sapiens", "Primates or younger", "Mammalia or younger", "Vertebrata or older", "all LCA levels"),
                                   compare_levels = c("Eutheria", "Euteleostomi"), title = "Gene age",  ylab = "Proportion of duplicated genes (%)",
                                   fill_values = c("Euteleostomi" = "#E79DCB", "Tetrapoda" = "#FFAB8C", "Amniota" = "#A6DD87", "Eutheria" = "#FFED77"), free_y = TRUE)

plot_age


ggsave(plot = plot_pn, width = 3.2, height = 16, device = cairo_pdf ,
       filename = sprintf("%s/ensembl_duplication_pn_lca.pdf", figure_file))
ggsave(plot = plot_pm, width = 3.2, height = 16, device = cairo_pdf ,
       filename = sprintf("%s/ensembl_duplication_pm_lca.pdf", figure_file))
ggsave(plot = plot_age, width = 4.8, height = 16, device = cairo_pdf ,
       filename = sprintf("%s/ensembl_duplication_age_lca.pdf", figure_file))


#5 compare major gene features between duplicates and singletons------------------
fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))
  )
}

median_iqr_str <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  med <- median(x)
  q1  <- quantile(x, 0.25, names = FALSE)
  q3  <- quantile(x, 0.75, names = FALSE)
  paste0(fmt_num(med, digits), " (", fmt_num(q1, digits), " - ", fmt_num(q3, digits), ")")
}

n_pct_str <- function(x, level, digits = 1) {
  denom <- sum(!is.na(x))
  n <- sum(x == level, na.rm = TRUE)
  if (denom == 0) return(NA_character_)
  paste0(format(n, big.mark = ","), " (", fmt_num(100 * n / denom, digits), "%)")
}

epsilon2_kruskal <- function(x, g) {
  keep <- complete.cases(x, g)
  x <- x[keep]
  g <- droplevels(factor(g[keep]))
  
  n <- length(x)
  k <- nlevels(g)
  if (n == 0 || k < 2 || n <= k) return(NA_real_)
  
  H <- unname(kruskal.test(x ~ g)$statistic)
  (H - k + 1) / (n - k)
}

low_high_rrb <- function(x, g, low_level, high_level, digits = 2) {
  keep <- !is.na(x) & !is.na(g) & g %in% c(low_level, high_level)
  x <- x[keep]
  g <- droplevels(factor(g[keep], levels = c(low_level, high_level)))
  
  if (length(unique(g)) < 2) {
    return(list(
      p = NA_real_,
      effect = NA_character_
    ))
  }
  
  # Wilcoxon rank-sum P
  wt <- suppressWarnings(wilcox.test(x ~ g, exact = FALSE))
  
  # rank-biserial correlation; order = high vs low so positive means higher in high group
  x_high <- x[g == high_level]
  x_low  <- x[g == low_level]
  rb <- effectsize::rank_biserial(x_high, x_low, ci = 0.95)
  
  est_col <- grep("biserial", names(rb), value = TRUE)[1]
  est <- rb[[est_col]][1]
  cil <- rb$CI_low[1]
  cih <- rb$CI_high[1]
  
  list(
    p = wt$p.value,
    effect = paste0(fmt_num(est, digits), " (", fmt_num(cil, digits), " to ", fmt_num(cih, digits), ")")
  )
}

overall_chisq_cramers_v <- function(x, g, digits = 2) {
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]
  
  tab <- table(g, x)
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(list(
      p = NA_real_,
      effect = NA_character_
    ))
  }
  
  ct <- suppressWarnings(chisq.test(tab, correct = FALSE))
  cv <- suppressWarnings(effectsize::cramers_v(tab))
  est_col <- grep("Cramers", names(cv), value = TRUE)[1]
  v <- cv[[est_col]][1]
  
  list(
    p = ct$p.value,
    effect = fmt_num(v, digits)
  )
}

low_high_or <- function(x, g, level, low_level, high_level, digits = 2) {
  keep <- !is.na(x) & !is.na(g) & g %in% c(low_level, high_level)
  x <- x[keep]
  g <- g[keep]
  
  if (length(unique(g)) < 2) {
    return(list(
      p = NA_real_,
      effect = NA_character_
    ))
  }
  
  a <- sum(g == high_level & x == level, na.rm = TRUE)
  b <- sum(g == high_level & x != level, na.rm = TRUE)
  c <- sum(g == low_level  & x == level, na.rm = TRUE)
  d <- sum(g == low_level  & x != level, na.rm = TRUE)
  
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  # rows = high, low; cols = level present, other
  rownames(tab) <- c("high", "low")
  colnames(tab) <- c("present", "other")
  
  if (sum(tab) == 0 || any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
    return(list(
      p = NA_real_,
      effect = NA_character_
    ))
  }
  
  ft <- fisher.test(tab)
  or  <- unname(ft$estimate)
  cil <- ft$conf.int[1]
  cih <- ft$conf.int[2]
  
  list(
    p = ft$p.value,
    effect = paste0(fmt_num(or, digits), " (", fmt_num(cil, digits), " to ", fmt_num(cih, digits), ")")
  )
}


make_supp_table <- function(data,
                            group_var,
                            group_order,
                            low_level,
                            high_level,
                            cont_vars = NULL,
                            cat_specs = NULL,
                            row_labels = NULL,
                            digits_cont = 2,
                            digits_eff = 2,
                            digits_pct = 1) {
  
  df <- data %>%
    dplyr::mutate(.group = factor(.data[[group_var]], levels = group_order)) %>%
    dplyr::filter(!is.na(.group))
  
  out <- list()
  
  # -------------------------
  # Continuous variables
  # -------------------------
  if (!is.null(cont_vars)) {
    for (v in cont_vars) {
      x <- df[[v]]
      g <- df$.group
      
      overall_p <- suppressWarnings(kruskal.test(x ~ g)$p.value)
      eps2 <- epsilon2_kruskal(x, g)
      low_high <- low_high_rrb(x, g, low_level, high_level, digits = digits_eff)
      
      row <- tibble(
        Characteristic = paste0(v, ", median (IQR)"),
        Total = median_iqr_str(x, digits = digits_cont),
        Overall_P = fmt_p(overall_p),
        Omnibus_effect = paste0("ε² = ", fmt_num(eps2, digits_eff)),
        Low_vs_high_P = fmt_p(low_high$p),
        Low_vs_high_effect = low_high$effect
      )
      
      for (lev in group_order) {
        row[[lev]] <- median_iqr_str(df[df$.group == lev, v, drop = TRUE], digits = digits_cont)
      }
      
      row <- row %>%
        dplyr::select(Characteristic, Total, all_of(group_order), Overall_P, Omnibus_effect,
                      Low_vs_high_P, Low_vs_high_effect)
      
      out[[length(out) + 1]] <- row
    }
  }
  
  # -------------------------
  # Categorical variables
  # -------------------------
  if (!is.null(cat_specs)) {
    for (v in names(cat_specs)) {
      x <- df[[v]]
      g <- df$.group
      
      overall <- overall_chisq_cramers_v(x, g, digits = digits_eff)
      
      for (lev in cat_specs[[v]]) {
        low_high <- low_high_or(x, g, lev, low_level, high_level, digits = digits_eff)
        
        key <- paste(v, lev, sep = "::")
        label <- if (!is.null(row_labels) && key %in% names(row_labels)) {
          row_labels[[key]]
        } else {
          paste0(v, ": ", lev, ", n (%)")
        }
        
        row <- tibble(
          Characteristic = label,
          Total = n_pct_str(x, lev, digits = digits_pct),
          Overall_P = fmt_p(overall$p),
          Omnibus_effect = paste0("V = ", overall$effect),
          Low_vs_high_P = fmt_p(low_high$p),
          Low_vs_high_effect = low_high$effect
        )
        
        for (gr in group_order) {
          row[[gr]] <- n_pct_str(df[df$.group == gr, v, drop = TRUE], lev, digits = digits_pct)
        }
        
        row <- row %>%
          dplyr::select(Characteristic, Total, all_of(group_order), Overall_P, Omnibus_effect,
                        Low_vs_high_P, Low_vs_high_effect)
        
        out[[length(out) + 1]] <- row
      }
    }
  }
  
  dplyr::bind_rows(out)
}

table_dup_ensembl_list <- lapply(list("ens_dup_human", "ens_dup_primate", "ens_dup_mammal", 
                                      "ens_dup_vertebrate_older", "ens_dup_all"), function(x){
                                        make_supp_table(data = as.data.frame(GPS_hsd),
                                                        group_var = x,
                                                        group_order = c("TRUE", "FALSE"),
                                                        low_level = "FALSE",
                                                        high_level = "TRUE",
                                                        cont_vars = c("Gene length (bp)"),
                                                        digits_cont = 0,
                                                        cat_specs = list("LOEUF_class" = c("Intolerant", "Tolerant"),
                                                                         "GES_class" = c("Essentialty", "Non_Essentialty"),
                                                                         "ifspecific" = c("Specific", "Non_specific")))
                                      })


names(table_dup_ensembl_list) <- c("ens_dup_human", "ens_dup_primate", "ens_dup_mammal", "ens_dup_vertebrate_older", 
                                   "ens_dup_all")
table_dup_ensembl_list$ens_dup_human[c(1,2,4,6),]


for (i in names(table_dup_ensembl_list)) {
  fwrite(table_dup_ensembl_list[[i]][c(1,2,4,6),], file = sprintf("%s/%s.csv", figure_file, i))
}

##6 compare duplicates vs singletons using DGD dataset-------------------
cont_vars <- c("Gene length (bp)", "Transcript length (bp)", "CDS Length", "CDS/Transcript Length ratio", "Protein size (aa)",
               "Transcript count", "Exon Counts", "Intron Counts", "GC content",
               "Number of SNPs (Gene)", "Number of SNPs (Transcript)",
               "dN/dS ratio mouse", "dN/dS ratio Chimp", "dN/dS ratio Gorilla",
               "Missense/Synonymous ratio", "Nonsense/Synonymous ratio")

table_dup <- make_supp_table(
  data = as.data.frame(GPS_hsd),
  group_var = "ifhsd2",
  group_order = c("Duplicates", "Singletons"),
  low_level = "Singletons",
  high_level = "Duplicates",
  cont_vars = c(cont_vars, "GES"),
  digits_cont = 1,
  cat_specs = list("LOEUF_class" = c("Intolerant", "Tolerant"),
                   "pLI_class"   = c("Intolerant", "Intermediate", "Tolerant"),
                   "GES_class" = c("Essentialty", "Non_Essentialty"),
                   "ifspecific" = c("Specific", "Non_specific"))
)

table_dup

write.csv(table_dup, file = sprintf("%s/table_dup.csv", figure_file), row.names = FALSE)

#############Duplication (paper 2026)#########################
#1 load data----------------
GPS_hsd_merge <- fread("Duplication/GPS_hsd_merge.csv", na.strings = c("", "NA"))
data_paper2026 <- fread("Duplication/Gene_Age_DataFrame.tsv")
head(data_paper2026)

count_bioprocess <- function(x) {
  sapply(x, function(z) {
    if (is.null(z) || is.na(z)) return(NA_integer_)
    
    z <- as.character(z)
    z <- trimws(z)
    
    if (z %in% c("", "[]", "NA", "NULL")) return(0L)
    
    # remove brackets
    z <- gsub("^\\[|\\]$", "", z)
    
    # split by comma
    items <- unlist(strsplit(z, ","))
    items <- trimws(items)
    items <- gsub("^['\"]|['\"]$", "", items)
    items <- items[items != ""]
    
    length(unique(items))
  })
}

data_paper2026[, n_bioprocess := count_bioprocess(Assc_BioProcesses)]
data_paper2026[, ifhsd := fifelse(
  is.na(Paralogs), NA_character_,
  fifelse(Paralogs == 0, "Singletons",
          fifelse(Paralogs > 0, "Duplicates", NA_character_))
)]

table(data_paper2026$AgeBin)

data_paper2026[, AgeBin_group5 := cut(
  as.numeric(AgeBin),
  breaks = c(0, 1, 3, 5, 7, 9),
  labels = c(1,2,3,4,5),
  include.lowest = TRUE,
  right = TRUE
)]

table(data_paper2026$AgeBin_group5)

data_paper2026_merge <- merge(GPS_hsd_merge, data_paper2026, by.x = "gene", by.y = "Gene_Symbol", all = T)
nrow(data_paper2026_merge) #22136
length(intersect(GPS_hsd_merge$gene, data_paper2026$Gene_Symbol)) #15294

#2 plot data-----------------
plot_feature_by_age_dup <- function(data, feature_var, title_text = "H. sapiens",
                                    age_var = "acrossAgeBin_500",
                                    dup_var = "ifhsd",
                                    age_order = NULL,
                                    dup_levels = c("FALSE", "TRUE"),
                                    dup_labels = c("Singletons", "Duplicates"),
                                    log_transform = TRUE,
                                    ylab = NULL,
                                    add_p = TRUE,
                                    test_method = "wilcox",
                                    p_adjust_method = "BH",
                                    p_label_type = "both") {
  
  library(data.table)
  library(ggplot2)
  
  dt <- copy(as.data.table(data))
  
  
  setnames(dt, old = names(dt), new = trimws(names(dt)))
  age_var <- trimws(age_var)
  dup_var <- trimws(dup_var)
  feature_var <- trimws(feature_var)
  
  dt <- dt[
    !is.na(get(age_var)) &
      !is.na(get(dup_var)) &
      !is.na(get(feature_var))
  ]
  
  if (is.null(age_order)) {
    age_order <- sort(unique(dt[[age_var]]))
  }
  
  dt[, age_group := factor(get(age_var), levels = age_order)]
  
  
  dt[, dup_group := factor(
    as.character(get(dup_var)),
    levels = as.character(dup_levels),
    labels = dup_labels
  )]
  
  dt <- dt[!is.na(age_group) & !is.na(dup_group)]
  
  dt[, y := as.numeric(get(feature_var))]
  
  if (log_transform) {
    dt[, y := log1p(y)]
  }
  
  if (is.null(ylab)) {
    ylab <- if (log_transform) {
      paste0("log(", feature_var, " + 1)")
    } else {
      feature_var
    }
  }
  
  
  stat_dt <- dt[, {
    
    x1 <- y[dup_group == dup_labels[1]]
    x2 <- y[dup_group == dup_labels[2]]
    
    if (length(x1) < 2 | length(x2) < 2) {
      data.table(
        n1 = length(x1),
        n2 = length(x2),
        p = NA_real_
      )
    } else {
      
      pval <- if (test_method == "wilcox") {
        wilcox.test(x1, x2)$p.value
      } else if (test_method == "t.test") {
        t.test(x1, x2)$p.value
      } else {
        stop("test_method should be 'wilcox' or 't.test'")
      }
      
      data.table(
        n1 = length(x1),
        n2 = length(x2),
        p = pval
      )
    }
    
  }, by = age_group]
  
  stat_dt[, p_adj := p.adjust(p, method = p_adjust_method)]
  
  stat_dt[, p_adj_label := ifelse(
    is.na(p_adj), "NA",
    ifelse(p_adj < 0.001, sprintf("%.2e", p_adj), sprintf("%.3f", p_adj))
  )]
  
  stat_dt[, signif_label := fifelse(
    is.na(p_adj), "NA",
    fifelse(p_adj < 0.001, "***",
            fifelse(p_adj < 0.01, "**",
                    fifelse(p_adj < 0.05, "*", "ns")))
  )]
  
  stat_dt[, label := p_adj_label]
  
  
  y_range <- range(dt$y, na.rm = TRUE)
  y_offset <- diff(y_range) * 0.08
  
  ypos_dt <- dt[, .(
    y_pos = max(y, na.rm = TRUE) + y_offset
  ), by = age_group]
  
  stat_dt <- merge(stat_dt, ypos_dt, by = "age_group", all.x = TRUE)
  
  
  p <- ggplot(
    dt,
    aes(x = age_group, y = y, fill = dup_group, group = interaction(age_group, dup_group))
  ) +
    geom_violin(
      position = position_dodge(width = 0.8),
      width = 0.85,
      alpha = 0.45,
      color = NA,
      trim = TRUE
    ) +
    geom_boxplot(
      position = position_dodge(width = 0.8),
      width = 0.18,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.5
    ) +
    stat_summary(
      aes(group = dup_group, color = dup_group),
      fun = median,
      geom = "line",
      position = position_dodge(width = 0.8),
      linewidth = 0.8
    ) +
    stat_summary(
      aes(group = dup_group, color = dup_group),
      fun = median,
      geom = "point",
      position = position_dodge(width = 0.8),
      size = 1.6
    ) +
    scale_fill_manual(
      values = setNames(c("#8EC1E8", "#F2C38B"), dup_labels)
    ) +
    scale_color_manual(
      values = setNames(c("#4C72B0", "#E49B39"), dup_labels)
    ) +
    labs(x = "", y = ylab) +
    scale_x_discrete(
      breaks = c(levels(dt$age_group)[1], levels(dt$age_group)[length(levels(dt$age_group))]),
      labels = c("Oldest", "Youngest")
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.03, 0.18))
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 13),
      legend.title = element_blank(),
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "italic", size = 15)
    ) +
    ggtitle(title_text)
  
  if (add_p) {
    p <- p +
      geom_text(
        data = stat_dt,
        aes(x = age_group, y = y_pos, label = label),
        inherit.aes = FALSE,
        size = 4,
        lineheight = 0.85
      )
  }
  
 
  attr(p, "stat_test") <- stat_dt
  
  return(p)
}

##group5---------------
plot_paper_list5 <- list(
  paper = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "n_bioprocess", title_text = "Martin and Tate 2025",
                                  dup_levels = c("Singletons", "Duplicates"), dup_labels = c("Singletons", "Duplicates"),
                                  age_var = "AgeBin_group5", dup_var = "ifhsd", log_transform = TRUE, ylab = "log(Biological Process Count)"),
  human = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "n_bioprocess", title_text = "human-lineage duplication",
                                  age_var = "AgeBin_group5", dup_var = "ens_dup_human", log_transform = TRUE, ylab = "log(Biological Process Count)"),
  primate = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "n_bioprocess", title_text = "primate-or-younger duplication",
                                    age_var = "AgeBin_group5", dup_var = "ens_dup_primate", log_transform = TRUE, ylab = "log(Biological Process Count)"),
  mammal = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "n_bioprocess", title_text = "mammalian-or-younger duplication",
                                   age_var = "AgeBin_group5", dup_var = "ens_dup_mammal", log_transform = TRUE, ylab = "log(Biological Process Count)"),
  vertebrate = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "n_bioprocess", title_text = "vertebrate-lineage-or-older duplication",
                                       age_var = "AgeBin_group5", dup_var = "ens_dup_vertebrate_older", log_transform = TRUE, ylab = "log(Biological Process Count)"))

plot_paper_list5$paper



plot_pn5 <- list(
  paper = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "pn_ld", title_text = "Martin and Tate 2025",
                                  dup_levels = c("Singletons", "Duplicates"), dup_labels = c("Singletons", "Duplicates"),
                                  age_var = "AgeBin_group5", dup_var = "ifhsd", log_transform = F, ylab = "GPS-N"),
  human = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "pn_ld", title_text = "human-lineage duplication",
                                  age_var = "AgeBin_group5", dup_var = "ens_dup_human", log_transform = F, ylab = "GPS-N"),
  primate = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "pn_ld", title_text = "primate-or-younger duplication",
                                    age_var = "AgeBin_group5", dup_var = "ens_dup_primate", log_transform = F, ylab = "GPS-N"),
  mammal = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "pn_ld", title_text = "mammalian-or-younger duplication",
                                   age_var = "AgeBin_group5", dup_var = "ens_dup_mammal", log_transform = F, ylab = "GPS-N"),
  vertebrate = plot_feature_by_age_dup(data = data_paper2026_merge, feature_var = "pn_ld", title_text = "vertebrate-lineage-or-older duplication",
                                       age_var = "AgeBin_group5", dup_var = "ens_dup_vertebrate_older", log_transform = F, ylab = "GPS-N"))

plot_pn5$mammal


for(i in names(plot_paper_list5)){
  ggsave(plot = plot_paper_list5[[i]], width = 5, height = 4, device = cairo_pdf,
         filename = sprintf("%s/group5_paper2026_%s.pdf", figure_file, i))
}

for(i in names(plot_pn5)){
  ggsave(plot = plot_pn5[[i]], width = 5, height = 4, device = cairo_pdf,
         filename = sprintf("%s/group5_pn_%s.pdf", figure_file, i))
}

