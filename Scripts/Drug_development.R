#1 load packages and set variables-------
sapply(c("data.table", "dplyr", "ggplot2", "corrplot", "ggpubr", "Gmisc", "openxlsx", "readxl", "paletteer", "MuMIn",
         "ggsci", "scales", "RColorBrewer", "gridExtra", "dplyr", "tidyr", "stringr", "colorspace", "cowplot",
         "ggbreak", "ggstatsplot", "viridis"), require, character.only = TRUE)
mycolor <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
pleio_color <- lighten(c("#80CDC1", "#85C6E0", "#223E84"), 0.45)
age_color <- lighten(c("#AB6191", "#DA7235", "#7FB65D", "#F2DD33"), 0.45)
age_color7 <- lighten(mycolor, 0.45)
age7_name <- c("Euteleostomi", "Tetrapoda", "Amniota", "Mammalia", "Theria", "Eutheria", "Primate")
x_text <-  c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria")
age_test <- c("Euteleostomi", "Eutheria")

setwd("/home/liumy/pleiotropy/github/data/")
figure_file <- "/home/liumy/pleiotropy/github/Figure"
#2 load data--------------------------
load("pleiotropy_maindata.RData")
#3 figure7--------------------
##figure7A: Proportion of FDA approved drug targets across groups----------------------
###load data---------------
fda_druggenes <- read.xlsx("fda_2023.xlsx") #drug genes from Nat. Rev. Drug Discov. 22, 864 (2023).
fda_target <- strsplit(fda_druggenes$targetIds, split = ",") %>%
  unlist() %>%  unique() %>%  trimws() %>%  .[!is.na(.) & . != "NA"]  #428 drugs with 504 targets
length(intersect(fda_target, pleiotropy_maindata$pm_ld$ensemblid)) #407

drug_data <- lapply(list("pm_ld", "pn_ld"), function(x, data = pleiotropy_maindata){
  data_merge <- data[[x]][, .(gene, ensemblid, pleio_class3, pleio10, pleio_class5, age_stage4)] %>%
    .[, iffda := ifelse(ensemblid %in% fda_target, 1, 0)] 
  return(data_merge)
})

names(drug_data) <- c("pm_ld", "pn_ld")

drug_percent_plot <- function(data = egenes_prop[1:3,49], 
                              group = factor(c("Low", "Intermediate", "High"), levels = c("Low", "Intermediate", "High")),
                              color_values = pleio_color, x_label = "", x_text = c("Low", "Intermediate", "High"),
                              y_label = "Proportion of eGenes", stat_name = "egene_pleio",
                              stat_test = test, test_list = c("Low", "High")){
  
  plot_data <- data.table(group = group, percent = data*100)
  
  title_wrap <- str_wrap(gsub("_", " ", x_label), width = 25)
  
  plot <- ggplot(data = plot_data, aes(x = group, y = percent, fill = group)) + theme_bw() +
    geom_bar(stat = "identity") +
    geom_signif(comparisons = list(test_list), 
                annotations = stat_test[[stat_name]]) +
    scale_fill_manual(values = color_values) +
    geom_text(aes(label = paste0(sprintf("%.2f",percent), "%")), vjust = 1.5, 
              color = "grey20", fontface = "bold", size = 6) +
    labs(y = y_label, title = title_wrap, x = "") +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_text(size = 15, face = "bold", ), 
          axis.text.x = element_text(size = 12, color = "black", face = "bold"), 
          axis.text.y = element_blank(), 
          axis.ticks = element_blank(), axis.line.x = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", 
          panel.border = element_blank(),
          plot.title = element_text(size=15, hjust=0.5),
          panel.grid.major = element_blank()) 
  
  return(plot)
}  
###plot data-------------------------
library(rstatix)
drug_prop <- list(
  pleio_n_fda = table(drug_data$pn_ld[, .(pleio_class5, iffda)]) %>% 
    prop.table(., margin = 1) %>% as.data.table(.) %>% .[iffda == 1,],
  pleio_m_fda = table(drug_data$pm_ld[, .(pleio_class5, iffda)]) %>% 
    prop.table(., margin = 1) %>% as.data.table(.) %>% .[iffda == 1,],
  age_n_fda = table(drug_data$pn_ld[, .(age_stage4, iffda)]) %>% 
    prop.table(., margin = 1) %>% as.data.table(.) %>% .[iffda == 1,])

test <- list(
  pleio_n = prop_trend_test(as.table(rbind(
    c(table(drug_data$pn_ld[iffda == 1,]$pleio_class5)),
    c(table(drug_data$pn_ld[iffda == 0,]$pleio_class5)))))$p %>% sprintf("%.2e", .),
  pleio_m = prop_trend_test(as.table(rbind(
    c(table(drug_data$pm_ld[iffda == 1,]$pleio_class5)),
    c(table(drug_data$pm_ld[iffda == 0,]$pleio_class5)))))$p %>% sprintf("%.2e", .),
  age = prop_trend_test(as.table(rbind(
    c(table(drug_data$pm_ld[iffda == 1,]$age_stage4)),
    c(table(drug_data$pm_ld[iffda == 0,]$age_stage4)))))$p
)

test

percent_fda <- list(pleio_n = drug_percent_plot(data = drug_prop$pleio_n_fda$N, 
                                                group = factor(c("First", "Second", "Third", "Fourth", "Fifth"), 
                                                               levels = c("First", "Second", "Third", "Fourth", "Fifth")),
                                                color_values = lighten(c("#6baed6","#4292c6","#2171b5","#08519c","#08306b"), 0.6),
                                                stat_test = test, test_list = c("First", "Fifth"), stat_name = "pleio_n",
                                                y_label = "Proportion of genes\ntargeted by FDA-approved drugs (%)"),
                    pleio_m = drug_percent_plot(data = drug_prop$pleio_m_fda$N, 
                                                group = factor(c("First", "Second", "Third", "Fourth", "Fifth"), 
                                                               levels = c("First", "Second", "Third", "Fourth", "Fifth")),
                                                color_values = lighten(c("#6baed6","#4292c6","#2171b5","#08519c","#08306b"), 0.6),
                                                stat_test = test, test_list = c("First", "Fifth"), stat_name = "pleio_m",
                                                y_label = "Proportion of genes\ntargeted by FDA-approved drugs (%)"),
                    age = drug_percent_plot(data = drug_prop$age_n_fda$N, 
                                            group = factor(c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria"), 
                                                           levels = c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria")),
                                            color_values = lighten(c("#FFC3AF", "#FF935F", "#DA7235", "#944F2B"), 0.4),
                                            stat_test = test, test_list = c("Euteleostomi", "Eutheria"), stat_name = "age",
                                            y_label = "Proportion of genes\ntargeted by FDA-approved drugs (%)"))

percent_fda$pleio_n
ggsave(plot = percent_fda$pleio_n, width = 6, height = 5, device = cairo_pdf,
       filename = sprintf("%s/percent_fda_pn.pdf",figure_file))


#figure7BC: Proportion of T–I pairs based on current development status--------------------
###load packages--------------------
options(stringsAsFactors=F)
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(binom))
suppressMessages(library(glue))
suppressMessages(library(lawstat))
suppressMessages(library(weights))
suppressMessages(library(epitools))
suppressMessages(library(DescTools))
suppressMessages(library(openxlsx))
suppressMessages(library(optparse))
suppressMessages(library(MASS)); summarize=dplyr::summarize; select=dplyr::select;
library(data.table)

###load functions-----------------
pipeline_best = function(merged_table,
                         basis='ti',
                         phase='combined',
                         require_insight=TRUE,
                         share_mode='L2G', # other option is V2G
                         min_share=0.5,
                         max_share=1,
                         worst_rank=Inf, # set to Inf if you want to include all
                         min_h4 = 0.9,
                         include_missing=FALSE,
                         associations=c('OMIM','GWAS'),
                         otg_subcat=c(''),
                         genebass_subcat=NULL,
                         mendelian_mechanism='',
                         min_year=2005,
                         max_year=2022,
                         firstyear=F,
                         minusomim=F,
                         lacking=NULL, # association sources required to be lacked by the T-I
                         andalso=NULL, # association sources required to *also* endorse the T-I
                         minusothersubcat=F,
                         mingenecount=0,
                         maxgenecount=Inf,
                         mapping_basis='all',
                         min_beta=0,
                         max_beta=Inf,
                         min_or=1,
                         max_or=Inf,
                         min_maf=0,
                         max_maf=1,
                         threshold=0.8,
                         network_list=NA,
                         pleio_name = "pleio10",
                         pleio_cat = c("1", "10"),
                         verbose=T) {
  
  start_time = Sys.time()
  
  mtable = merged_table
  mtable$pleio_class <- mtable[[pleio_name]]
  
  if (verbose) {
    cat(file=stderr(),'Starting row count: ',nrow(mtable),'\n')
    flush.console()
  }
  
  # add & select unique ID
  if (basis %in% c('di_mesh','drug-indication')) {
    mtable$di_uid = paste0(mtable$drugid,'-',mtable$indication_mesh_id)
    mtable$uid = mtable$di_uid
  } else if (basis %in% c('ti','target-indication')) {
    mtable$ti_uid = paste0(mtable$gene,'-',mtable$indication_mesh_id)
    mtable$uid = mtable$ti_uid
  } else if (basis=='drug') {
    mtable$uid = mtable$drugid
  }
  
  # assign highest level of advancement depending on phase specified
  if (phase == 'active') {
    meta = meta_acat
    mtable$cat = mtable$acat
  } else if (phase == 'historical') {
    meta = meta_hcat
    mtable$cat = mtable$hcat
  } else if (phase == 'combined') {
    meta = meta_ccat
    mtable$cat = mtable$ccat
  }
  
  mtable$catnum = meta$num[match(mtable$cat, meta$cat)]
  # remove "Other"
  mtable = mtable[mtable$cat != 'Other' ,]
  # map L2G share
  mtable$assoc_share = mtable$l2g_share
  mtable$assoc_rank = mtable$l2g_rank
  
  if (verbose) {
    cat(file=stderr(),'Selecting user-specified filters...')
    flush.console()
  }
  
  # by default, require non-missing target & indication
  if (!include_missing) {
    mtable = mtable[mtable$gene != '' & mtable$indication_mesh_id != '' & !is.na(mtable$gene) & !is.na(mtable$indication_mesh_id),]
  }
  # genetic insight requirement
  if (require_insight) {
    mtable = mtable[mtable$indication_mesh_id %in% indic$indication_mesh_id[indic$genetic_insight != 'none'],]
  } else {
    # otherwise simply require the indication be present in the indic table
    mtable = mtable[mtable$indication_mesh_id %in% indic$indication_mesh_id,]
  }
  
  # remove omim-supported associations if desired. only works in T-I mode
  # note that order of operations is important - this must come before associations filter
  if (minusomim) {
    # look for first year in which a target-*indication* pair was genetically supported
    mtable %>%
      filter(comb_norm >= threshold) %>%
      filter(assoc_source %in% c('OMIM')) %>%
      select(ti_uid) -> omim_supported_ti
    # retain the null rows (i.e. no association) or those where hte T-I is not in OMIM
    # what gets removed? e.g. OTG associations that were already established by OMIM
    mtable %>%
      filter(is.na(assoc_source) | !(ti_uid %in% omim_supported_ti$ti_uid)) -> mtable
  }
  
  # remove any association sources required to be lacked
  if (!is.null(lacking)) {
    mtable %>%
      filter(comb_norm >= threshold) %>%
      filter(assoc_source %in% lacking) %>%
      filter(!(assoc_source %in% 'OTG' & l2g_share >= min_share)) %>%
      select(ti_uid) -> lackable_supported_ti
    mtable %>%
      filter(is.na(assoc_source) | !(ti_uid %in% lackable_supported_ti$ti_uid)) -> mtable
  }
  
  if (!is.null(andalso)) {
    mtable %>%
      filter(comb_norm >= threshold) %>%
      filter(assoc_source %in% andalso) %>%
      filter(!(assoc_source %in% 'OTG' & l2g_share < min_share)) %>%
      select(ti_uid) -> andalso_supported_ti
    mtable %>%
      filter(is.na(assoc_source) | (ti_uid %in% andalso_supported_ti$ti_uid)) -> mtable
  }
  
  # user-specified sources of genetic associations
  # allow user to specify either grouping terms like "GWAS", or specific sources
  source_map = tibble(source=c("OTG", "PICCOLO", "Genebass", "OMIM", "intOGen"),
                      source_name=c('GWAS','GWAS','GWAS','OMIM','Somatic'))
  if (!(identical(associations, c('OMIM','GWAS','Somatic'))) ) {
    mtable %>%
      left_join(source_map, by=c('assoc_source'='source')) %>%
      filter(is.na(source_name) | source_name %in% associations | assoc_source %in% associations) -> mtable
    
  }
  
  # further filter of subtype of OTG association
  if (otg_subcat != '') {
    # first subset to just OTG
    mtable %>%
      filter(is.na(assoc_source) | assoc_source %in% 'OTG') -> mtable
    # now pick the types
    mtable %>%
      mutate(gwas_source = case_when(grepl('GCST', original_link) ~ 'GWAS Catalog',
                                     grepl('FINNGEN', original_link) ~ 'FinnGen',
                                     grepl('NEALE',original_link) ~ 'Neale UKBB',
                                     TRUE ~ 'Other')) %>%
      filter(is.na(assoc_source) | gwas_source %in% otg_subcat) -> mtable
  }
  
  # further filter of annotation & test in Genebass
  if (!is.null(genebass_subcat)) {
    grepstring = paste(genebass_subcat, collapse='|')
    # mtable %>%
    #   filter(!is.na(assoc_source) & assoc_source %in% 'Genebass') %>%
    #   filter(grepl(grepstring, assoc_info)) -> genebass_hits
    mtable %>%
      filter(is.na(assoc_source) | !(assoc_source %in% 'Genebass') | grepl(grepstring, assoc_info)) -> mtable
  }
  
  # apply user-specified filter of OMIM disease mechanism
  if (mendelian_mechanism != '') {
    mtable %>%
      filter(!(mtable$assoc_source %in% 'OMIM') | is.na(mtable$assoc_info) | grepl(mendelian_mechanism,mtable$assoc_info)) -> mtable
  }
  
  # apply user-specified OTG gene mapping share & rank minimum/maximum
  if (share_mode == 'V2G') {
    mtable$assoc_share = mtable$v2g_share
    mtable$assoc_rank = mtable$v2g_rank
    assoc$assoc_share = assoc$v2g_share # needed in assoc table too for genecount section below
  } else if (share_mode == 'L2G') {
    mtable$assoc_share = mtable$l2g_share
    mtable$assoc_rank = mtable$l2g_rank
    assoc$assoc_share = assoc$l2g_share
  }
  
  # worst rank
  if (worst_rank < Inf) {
    mtable %>%
      filter(!(assoc_source %in% 'OTG') | mtable$assoc_rank <= worst_rank) -> mtable
  }
  
  # note that among OTG associations, throw out any with NA share, as these would be zeroes (does not occur in Dec 2021 dataset anyway)
  # and note that with L2G a significant number of associations actually have 100% share, so only delete > max_share and not >= max_share
  mtable %>%
    filter(!(assoc_source %in% 'OTG') | (!is.na(assoc_share) & assoc_share >= min_share & assoc_share <= max_share)) -> mtable
  
  # apply user-specified H4 minimum / maximum
  if (min_h4 > .9) {
    mtable %>%
      filter(!(assoc_source %in% 'PICCOLO') | (!is.na(mtable$pic_h4) & mtable$pic_h4 >= min_h4)) -> mtable
  }
  
  # apply user-specified genecount minimum/maximum
  if (mingenecount > 0 | maxgenecount < Inf) {
    assoc %>%
      filter(source %in% associations) %>%
      filter(source!='OTG' | (assoc_share >= min_share & assoc_share <= max_share)) %>%
      group_by(mesh_id) %>%
      summarize(.groups='keep', n_genes=length(unique(gene))) -> gene_counts
    mtable$gene_count = gene_counts$n_genes[match(mtable$assoc_mesh_id, gene_counts$mesh_id)]
    mtable %>%
      filter(is.na(gene_count) | gene_count >= mingenecount & gene_count <= maxgenecount) -> mtable
  }
  
  
  # apply "first year" criteria if applicable
  if (firstyear) {
    # look for first year in which a target-*indication* pair was genetically supported
    # only works in T-I mode
    mtable %>%
      filter(comb_norm >= threshold) %>%
      group_by(ti_uid) %>%
      summarize(.groups='keep', min_assoc_year=min(assoc_year)) -> ti_first_sup
    mtable$min_assoc_year = ti_first_sup$min_assoc_year[match(mtable$ti_uid, ti_first_sup$ti_uid)]
    # keep entries with no assoc year (the null rows), or where assoc year = the min assoc year, i.e. this is
    # the first report of this genetic association (or tied for first)
    mtable %>%
      filter(is.na(assoc_year) | assoc_year == min_assoc_year) -> mtable
  }
  
  # avoid -Inf values in comparisons by hard coding in case of all missing values:
  if (sum(!is.na(mtable$assoc_year)) == 0) {
    mtable_intrinsic_max_year = 2021
    mtable_intrinsic_min_year = 2000
  } else {
    mtable_intrinsic_max_year = max(mtable$assoc_year, na.rm=T)
    mtable_intrinsic_min_year = min(mtable$assoc_year, na.rm=T)
  }
  
  # apply user-specified filter of association years - for OTG only
  if (min_year > mtable_intrinsic_min_year | max_year < mtable_intrinsic_max_year) {
    mtable %>%
      filter(is.na(assoc_source) | assoc_source != 'OTG' | assoc_year >= min_year & assoc_year <= max_year) -> mtable
  }
  
  # join back in beta
  mtable$abs_beta = abs(assoc$beta[match(mtable$arow, assoc$arow)])
  if (min_beta > 0 | max_beta < Inf) {
    stopifnot(associations=='OTG') # only supported for OTG-only mode
    mtable %>%
      filter(is.na(assoc_source) | (!is.na(mtable$abs_beta) & mtable$abs_beta >= min_beta & mtable$abs_beta <= max_beta)) -> mtable
  }
  
  # same as beta but for OR
  mtable$abs_or = abs_or(assoc$odds_ratio[match(mtable$arow, assoc$arow)])  
  if (min_or > 1 | max_or < Inf) {
    stopifnot(associations=='OTG') # only supported for OTG-only mode
    # leave null rows (with no association source) but delete all those mapped to OTG that do not have or, or have or outside range
    mtable %>%
      filter(is.na(assoc_source) | (!is.na(mtable$abs_or) & mtable$abs_or >= min_or & mtable$abs_or < max_or)) -> mtable
  }
  
  
  # lead SNP maf
  if (min_maf > 0 | max_maf < 1) {
    mtable$lead_maf = pmin(mtable$af_gnomad_nfe, 1-mtable$af_gnomad_nfe)
    mtable$lead_maf[!is.na(mtable$lead_maf) & mtable$lead_maf < 0] = NA
    # lead_maf >= min_maf & lead_maf < max_maf
    # >= and < gets you "[, )" logic
    # also remove those that are GWAS where lead_maf is NA - likely in non-European populations so af_gnomad_nfe is not relevant
    mtable %>% 
      filter(is.na(assoc_source) | !(assoc_source %in% c('OTG','PICCOLO','Genebass')) | (!is.na(lead_maf) & (lead_maf >= min_maf & lead_maf < max_maf))) -> mtable 
  }
  
  
  
  
  if (verbose) {
    cat(file=stderr(),'Selecting highest phase reached and best genetic similarity...')
    flush.console()
  }
  
  suppressWarnings(mtable %>% group_by(uid) %>% summarize(maxsim = max(comb_norm, na.rm=T), maxcat=max(catnum, na.rm=T)) -> step1)
  
  if (verbose) {
    cat(file=stderr(),nrow(step1),'rows remain.\n')
    flush.console()
  }
  
  if (verbose) {
    cat(file=stderr(),'Joining back in program details...')
    flush.console()
  }
  # add a filter first - only slightly reduces row count
  mtable %>%
    filter(uid %in% step1$uid & comb_norm %in% unique(step1$maxsim) & catnum %in% unique(step1$maxcat)) -> mtable
  # use tidy to left join
  step1 %>%
    left_join(mtable, by = c("uid" = "uid", "maxsim" = "comb_norm", "maxcat" = "catnum")) %>%
    rename(similarity=maxsim, catnum=maxcat) -> step2
  
  if (verbose) {
    cat(file=stderr(),nrow(step2),'rows found after join.\n')
    flush.console()
  }
  if (verbose) {
    cat(file=stderr(),'Removing duplicates resulting from ties...')
    flush.console()
  }
  
  # prioritize rows with known outcome at most advanced phase, then de-dup
  step2 %>%
    mutate(highest_phase_with_known_outcome = case_when(!is.na(succ_3_a) ~ 3,
                                                        !is.na(succ_2_3) ~ 2,
                                                        !is.na(succ_1_2) ~ 1,
                                                        !is.na(succ_p_1) ~ 0)) %>%
    arrange(uid, desc(highest_phase_with_known_outcome)) %>%
    group_by(uid) %>%
    slice(1) %>%
    ungroup() -> step2 # row count should drop back to ~ that of step1
  
  if (verbose) {
    cat(file=stderr(),nrow(step2),'rows remain.\n')
    flush.console()
  }
  
  # annotate in additional info:
  step2$areas = indic$areas[match(step2$indication_mesh_id, indic$indication_mesh_id)]
  step2$genetic_insight = replace_na(indic$genetic_insight[match(step2$indication_mesh_id, indic$indication_mesh_id)],'none')
  step2$target_status = ''
  step2$target_status[step2$similarity >= threshold] = 'genetically supported target'
  step2$target_status[step2$similarity <  threshold] = 'unsupported target'
  step2$target_status[require_insight & step2$genetic_insight == 'none'] = 'indication lacks genetic insight'
  step2$target_status[is.na(step2$gene) | step2$gene == ''] = 'no target annotated'
  step2$target_status[is.na(step2$indication_mesh_id) | step2$indication_mesh_id == ''] = 'no indication annotated'
  
  step2$pleio <- ""
  step2$pleio[step2$pleio_class == pleio_cat[2]] = "high"
  step2$pleio[step2$pleio_class == pleio_cat[1]] = "low"
  
  if (verbose) {
    cat(file=stderr(),paste0('Using sim threshold ',threshold,', "genetically supported target" rows: ',sum(step2$target_status=='genetically supported target'),'....\n'))
    time_elapsed = (Sys.time() - start_time)
    cat(file=stderr(),'pipeline_best completed in',round(time_elapsed,1),units(time_elapsed),'.\n')
    flush.console()
  }
  
  return (step2)
}


advancement_forest = function(best_table, phase='combined', pleio_cat = "1", pleio_name = "pleio10") {
  if (phase == 'active') {
    meta = meta_acat
  } else if (phase == 'historical') {
    meta = meta_hcat
  } else if (phase == 'combined') {
    meta = meta_ccat
  }
  
  best_table$pleio_class <- best_table[[pleio_name]]
  
  meta %>%
    left_join(best_table, by=c('num'='catnum', 'cat'='cat')) %>%
    filter(cat != 'Other') %>%
    filter(!(target_status %in% c('no indication annotated','no target annotated'))) %>%
    group_by(catnum=num, cat) %>%
    summarize(.groups='keep',
              n_total = sum(!is.na(pleio_class)),
              n_gensup = sum(pleio_class %in% pleio_cat, na.rm=T)) -> forest_data
  bconf_obj = binom.confint(x=forest_data$n_gensup, n=forest_data$n_total, method='wilson', conf.level = .95)
  forest_data$proportion = bconf_obj$mean
  forest_data$l95 = bconf_obj$lower
  forest_data$u95 = bconf_obj$upper
  forest_data$y = max(forest_data$catnum) - forest_data$catnum + 1
  colnames(forest_data) = c('num','label','denominator','numerator','mean','l95','u95','y')
  return (forest_data)
}

abs_or = function(odds_ratio) {
  abs_odds_ratio = odds_ratio
  flip_indices = odds_ratio < 1 & !is.na(odds_ratio)
  abs_odds_ratio[flip_indices] = 1/odds_ratio[flip_indices]
  return (abs_odds_ratio)
}
###load data----------------
merge2_raw = read_tsv('genetic_support-main/merge2.tsv.gz', col_types=cols())
assoc = read_tsv('genetic_support-main/assoc.tsv.gz', col_types=cols())
indic = read_tsv("genetic_support-main/indic.tsv", col_types=cols())
meta_ccat = read_tsv('genetic_support-main/meta_ccat.tsv', col_types=cols())

drug_merge_list <- list(
  merge2_pn = merge(merge2_raw, pleiotropy_maindata$pn_ld[
    , .(gene, pleio_class3, pleio_class5, pleio10, pleio100)], by = "gene", all.x = T) %>% as_tibble(.),
  merge2_pm = merge(merge2_raw, pleiotropy_maindata$pm_ld[
    , .(gene, pleio_class3, pleio_class5, pleio10, pleio100)], by = "gene", all.x = T) %>% as_tibble(.)
)

drug_ti_combined <- list(
  combined_pleio5_pn = pipeline_best(drug_merge_list$merge2_pn,  phase='combined', basis='ti', verbose=F, 
                                     pleio_name = "pleio_class5", pleio_cat = c("1st", "5th")),
  combined_pleio5_pm = pipeline_best(drug_merge_list$merge2_pm,  phase='combined', basis='ti', verbose=F, 
                                     pleio_name = "pleio_class5", pleio_cat = c("1st", "5th"))
)

###Figure for L-GPS and H-GPS------------------
forest_data_list <- list(
  high_pn = advancement_forest(drug_ti_combined$combined_pleio5_pn, phase='combined', pleio_cat = "5th", pleio_name = "pleio_class5"),
  low_pn = advancement_forest(drug_ti_combined$combined_pleio5_pn,phase='combined', pleio_cat = "1st", pleio_name = "pleio_class5"),
  high_pm = advancement_forest(drug_ti_combined$combined_pleio5_pm, phase='combined', pleio_cat = "5th", pleio_name = "pleio_class5"),
  low_pm = advancement_forest(drug_ti_combined$combined_pleio5_pm,phase='combined', pleio_cat = "1st", pleio_name = "pleio_class5")
)

### plot
library(colorspace); library(cowplot)
forest_func <- function(data = hist_ti_forest, legend_label = c("L-GPS-N", "H-GPS-N")){
  
  data[[1]]$pleio <- "high"; data[[2]]$pleio <- "low" 
  data_long <- do.call(rbind, data)
  data_long[,5:7] <- data_long[,5:7]*100
  
  
  plot <- ggplot(data_long, aes(x = mean, y = y-as.numeric(as.factor(pleio))*0.5+0.7, 
                                color = factor(pleio, levels = c("low", "high"), labels = c("low", "high")))) + theme_bw() +
    geom_point(size = 3.5) +
    geom_errorbarh(aes(xmin = l95, xmax = u95), height = 0.1, size = 1.2) +
    labs(x = "Proportion of T–I pairs (%)", y = NULL) +
    scale_color_manual(values = lighten(c("#9E6531", "#346981"), 0.5),
                       labels = legend_label) +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"), 
          axis.title.x = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), 
          legend.title = element_blank(), legend.text = element_text(size = 12, face = "bold"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.position = "top",
          panel.border = element_blank(), axis.line.x.bottom = element_line(color = 'grey30')) +
    geom_hline(yintercept = c(1.5,2.5,3.5,4.5,5.5), color = "gray", size = 0.5) +
    scale_y_continuous(breaks = c(1:5), labels = c("Launched", "Phase III", "Phase II", "Phase I", "Preclinical"))
  
  plot_text <- ggplot(data_long, aes(x = 1, y = y-as.numeric(as.factor(pleio))*0.5+0.7, 
                                     color = factor(pleio, levels = c("low", "high"), labels = c("low", "high")))) + theme_bw() +
    geom_text(aes(label = paste0(formatC(numerator,big.mark=','),'/',formatC(denominator,big.mark=','))), 
              size = 4.5, fontface = "bold", show.legend = F) +
    scale_color_manual(values = lighten(c("#9E6531", "#346981"), 0.5),
                       labels = c("L-GPS-N", "H-GPS-N")) +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_blank(), axis.text = element_blank(), 
          axis.ticks = element_blank(), legend.title = element_blank(),
          #legend.position = "none",
          legend.text = element_blank(),
          legend.key = element_blank(),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.position = "top",
          panel.border = element_blank(), axis.line.x.bottom = element_line(color = 'grey30')) +
    geom_hline(yintercept = c(1.5,2.5,3.5,4.5,5.5), color = "gray", size = 0.5) +
    scale_x_continuous(limits = c(0,2), breaks = 1)
  
  plot_merge <- plot_grid(plot, plot_text, align = "h", rel_widths = c(4, 1))
  
  return(plot_merge)
}

forest_plot_list <- list(
  pn = forest_func(data = list(forest_data_list$high_pn, forest_data_list$low_pn), legend_label = c("L-GPS-N", "H-GPS-N")),
  pm = forest_func(data = list(forest_data_list$high_pm, forest_data_list$low_pm), legend_label = c("L-GPS-M", "H-GPS-M"))
)

forest_plot_list$pn


ggsave(plot = forest_plot_list$pn, width = 7, height = 4.5, device = cairo_pdf,
       filename = sprintf("%s/forest_combined_pn.pdf",figure_file))


prop.test(x = c(2208, 1415),n = c(10051, 10051))$p.value %>% sprintf("%.2e", .) #"7.52e-48"
prop.test(x = c(1074, 575),n = c(4127, 4127))$p.value %>% sprintf("%.2e", .) #"8.93e-43"
prop.test(x = c(1662, 927),n = c(6888, 6888))$p.value %>% sprintf("%.2e", .) #"1.13e-57"
prop.test(x = c(328, 185),n = c(1464, 1464))$p.value %>% sprintf("%.2e", .) #"5.08e-12"
prop.test(x = c(323, 181),n = c(1435, 1435))$p.value %>% sprintf("%.2e", .) #"4.60e-12"



##figure7B--------------
plotline_func <- function(data = forest_pleio5_data, x_label = "Quintiles of GPS-N"){
  
  data <- data[, x := rep(1:5, each = 5)]
  data$label <- factor(data$label, levels = c("Preclinical", "Phase I", "Phase II", "Phase III", "Launched"),
                       labels = c("Preclinical", "Phase I", "Phase II", "Phase III", "Launched"))
  
  plot <- ggplot(data, aes(x = x, y = mean, 
                           group = label, color = label)) + theme_bw() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = l95, ymax = u95), width = 0.1) +
    geom_line() +
    labs(y = "Proportion of T–I pairs (%)", x = x_label) +
    scale_color_manual(values = mycolor) +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"), 
          axis.title.x = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), 
          legend.title = element_blank(), legend.text = element_text(size = 12, face = "bold"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.position.inside = c(5,10),
          panel.border = element_blank(), axis.line = element_line(color = 'grey30')) +
    scale_x_continuous(breaks = c(1:5), labels = c("First", "Second", "Third", "Fourth", "Fifth"))
  
  return(plot) 
}

plotline_data_list <- lapply(list(1,2), function(score_name){
  plotline_data <- lapply(list("1st", "2nd", "3rd", "4th", "5th"), function(x){
    data <- advancement_forest(best_table = drug_ti_combined[[score_name]], phase='combined', 
                               pleio_name = "pleio_class5", pleio_cat = x) 
    data$group <- x
    return(data)}) %>%
    do.call(rbind, .) %>% as.data.table(.) 
  plotline_data[, 5:7] <- plotline_data[, 5:7]*100
  return(plotline_data)
})
names(plotline_data_list) <- c("pn", "pm")

plotline_plot_list <- list(
  pn = plotline_func(plotline_data_list$pn, x_label = "Quintiles of GPS-N"),
  pm = plotline_func(plotline_data_list$pm, x_label = "Quintiles of GPS-M")
)

plotline_plot_list$pn
plotline_plot_list$pm

ggsave(plot = plotline_plot_list$pn, width = 6, height = 4, device = cairo_pdf,
       filename = sprintf("%s/plotline_pn.pdf",figure_file))


## p for trend 
library(rstatix)

pn_pfortrend <- lapply(as.list(1:5), function(x, data = plotline_data_list$pn){
  
  use_data <- as.data.frame(data)
  xtab <- as.table(rbind(
    c(use_data[use_data$num == x, "numerator"]),
    c(use_data[use_data$num == x, "denominator"] - use_data[use_data$num == x, "numerator"])
  ))
  
  test_stat <- sprintf("%.2e",prop_trend_test(xtab)$p) 
  
  return(test_stat)
}) %>% unlist()
pm_pfortrend <- lapply(as.list(1:5), function(x, data = plotline_data_list$pm){
  
  use_data <- as.data.frame(data)
  xtab <- as.table(rbind(
    c(use_data[use_data$num == x, "numerator"]),
    c(use_data[use_data$num == x, "denominator"] - use_data[use_data$num == x, "numerator"])
  ))
  
  test_stat <- sprintf("%.2e",prop_trend_test(xtab)$p) 
  
  return(test_stat)
}) %>% unlist()

pn_pfortrend
pm_pfortrend
