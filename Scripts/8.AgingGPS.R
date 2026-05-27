#1 load packages and set variables-------
setwd("/path/to/Gene_pleiotropy-main/Data")
figure_file <- "/path/to/Gene_pleiotropy-main/Output/AgingGPS"
code_file <- "/path/to/Gene_pleiotropy-main/Scripts"
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

#1 compare with pan-trait GPS and aging-GPS--------------------
##1.1 load data------------
agingGPS <- fread("AgingGPS/agingGPS.csv", na.strings = c("", "NA"))
head(agingGPS)

load("pleiotropy_maindata.RData")
head(pleiotropy_maindata$pm_ld)
GPS <- merge(pleiotropy_maindata$pm_ld,
             pleiotropy_maindata$pn_ld[, .(gene, pn_ld = use_score, pn_pleio3 = pleio_class3, pn_pleio5 = pleio_class5, pn_pleio10 = pleio10)],
             by = "gene")
head(GPS)

setnames(GPS, old = c("use_score", "pleio_class3", "pleio_class5", "pleio10"),
         new = c("pm_ld", "pm_pleio3", "pm_pleio5", "pm_pleio10"))
head(GPS)

merge_GPS_agingGPS <- merge(agingGPS, GPS, by = "gene")
head(merge_GPS_agingGPS)

cor(merge_GPS_agingGPS[, .(pm_ld, pm_ld_aging, pn_ld, pn_ld_aging)], method = "sp", use = "com")
cor.mtest(merge_GPS_agingGPS[, .(pm_ld, pm_ld_aging, pn_ld, pn_ld_aging)], method = "sp", use = "com")$p

##1.2 Jaccard as the effect size; Fisher (hypergeometric) p for significance (use)---------
overlap_jaccard_test <- function(U, A, B, alternative = c("greater", "less", "two.sided")) {
  
  alternative <- match.arg(alternative)
  U <- unique(U)
  A <- unique(intersect(A, U))
  B <- unique(intersect(B, U))
  
  m <- length(A)
  n <- length(B)
  N <- length(U)
  k <- length(intersect(A, B))
  
  denom <- m + n - k
  J <- if (denom > 0) k / denom else NA_real_
  
  #           In B   Not in B
  # In A        k     m - k
  # Not in A  n - k  N - m - n + k
  a <- k
  b <- m - k
  c <- n - k
  d <- N - m - n + k
  cont <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                 dimnames = list(A = c("inA", "notInA"), B = c("inB", "notInB")))
  
  
  fisher_p <- tryCatch(
    fisher.test(cont, alternative = alternative)$p.value,
    error = function(e) NA_real_
  )
  
  fisher_or <- tryCatch(
    fisher.test(cont, alternative = alternative)$estimate,
    error = function(e) NA_real_
  )
  
  
  tab <- data.table(Jaccard = J, k = k, m = m, n = n, N = N, fisher_p = fisher_p, fisher_or = fisher_or)
  return(tab)
}

pn_overlap_jaccard <- lapply(as.list(1:10), function(pleio10){
  
  U = merge_GPS_agingGPS$gene
  A = merge_GPS_agingGPS[pn_pleio10_aging == pleio10,]$gene
  B = merge_GPS_agingGPS[pn_pleio10 == pleio10,]$gene
  
  test <- overlap_jaccard_test(U, A, B, alternative = "greater")
  test[, pleio_group := pleio10]
}) %>% rbindlist()

pn_overlap_jaccard

pm_overlap_jaccard <- lapply(as.list(1:10), function(pleio10){
  
  U = merge_GPS_agingGPS$gene
  A = merge_GPS_agingGPS[pm_pleio10_aging == pleio10,]$gene
  B = merge_GPS_agingGPS[pm_pleio10 == pleio10,]$gene
  
  test <- overlap_jaccard_test(U, A, B, alternative = "greater")
  test[, pleio_group := pleio10]
  
}) %>% rbindlist()

pm_overlap_jaccard

##visualization
library(forcats)
pn_overlap_jaccard_plotdata <- pn_overlap_jaccard %>%
  mutate(minus_log10p = -log10(fisher_p),
         or = fisher_or,
         pleio_group = fct_rev(factor(pleio_group, levels = sort(unique(pleio_group)))))

pn_overlap_jaccard_plot <- ggplot(pn_overlap_jaccard_plotdata, aes(x = Jaccard, y = pleio_group)) +  theme_bw() +
  geom_point(aes(size = minus_log10p, color = or)) +
  scale_size_area(max_size = 12, name = expression(-log[10](P))) +
  scale_color_gradient2(low = "steelblue", mid = "grey80", high = "firebrick",
                        midpoint = 0, name = "Enrichment") +
  labs(x = "Jaccard index", y = "Pleiotropy decile (GPS-N)") +
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        legend.position = "right",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid.minor = element_blank())

pn_overlap_jaccard_plot

ggsave(plot = pn_overlap_jaccard_plot, width = 5, height = 8, device = cairo_pdf,
       filename = sprintf("%s/pn_overlap_jaccard_plot.pdf", figure_file))

pm_overlap_jaccard_plotdata <- pm_overlap_jaccard %>%
  mutate(minus_log10p = -log10(fisher_p),
         or = fisher_or,
         pleio_group = fct_rev(factor(pleio_group, levels = sort(unique(pleio_group)))))

pm_overlap_jaccard_plot <- ggplot(pm_overlap_jaccard_plotdata, aes(x = Jaccard, y = pleio_group)) +  theme_bw() +
  geom_point(aes(size = minus_log10p, color = or)) +
  scale_size_area(max_size = 12, name = expression(-log[10](P))) +
  scale_color_gradient2(low = "steelblue", mid = "grey80", high = "firebrick",
                        midpoint = 0, name = "Enrichment") +
  labs(x = "Jaccard index", y = "Pleiotropy decile (GPS-M)") +
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        legend.position = "right",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid.minor = element_blank())

pm_overlap_jaccard_plot

ggsave(plot = pm_overlap_jaccard_plot, width = 5, height = 8, device = cairo_pdf,
       filename = sprintf("%s/pm_overlap_jaccard_plot.pdf", figure_file))

#2. GO gsea analysis---------------
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(org.Hs.eg.db)

agingGPS <- fread("AgingGPS/agingGPS.csv", na.strings = c("", "NA"))

##GPS-N------------
pn_geneList <- agingGPS$pn_ld_aging
names(pn_geneList) <- agingGPS$gene
head(pn_geneList)
pn_geneList <- sort(pn_geneList, decreasing = TRUE)

set.seed(1234)
pn_gsego <- gseGO(geneList     = pn_geneList, 
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  keyType      = "SYMBOL",
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  scoreType = "pos",
                  eps = 0,
                  seed = TRUE,
                  verbose      = FALSE)


#gene-centered network
library(igraph)
library(ggraph)
library(ggplot2)
library(stringr)

source(sprintf("%s/gene_centered_network_function.R", code_file))
net_res_pn <- make_gene_centered_gsea_network(
  gsea_obj = pn_gsego,
  gene_score = pn_geneList,
  top_n_pathways = 30,
  high_cutoff = min(agingGPS[pn_pleio5_aging == "5th"]$pn_ld_aging),
  top_n_genes = 20,
  min_pathways_per_gene = 2,
  pathway_label_width = 24,
  
  gene_radius = 0.3,
  pathway_radius = 0.60,
  pathway_label_offset = 0.07,
  
  pathway_text_size = 3.5,
  pathway_text_color = "grey20",
  pathway_node_fill = "#4A76A8",
  pathway_node_color = "white",
  
  gene_text_size = 4,
  edge_color = "grey10",
  edge_alpha = 0.22,
  edge_width_range = c(0.12, 0.65),
  
  gps_low = "#E9D8CC",
  gps_high = "#B8572A",
  gene_node_color = "white",
  gene_node_stroke = 0.55,
  gene_label_offset = 0.05,
  
  plot_margin = ggplot2::margin(10, 10, 10, 100),
  
  save_prefix = "Temp/pn_GSEA_high_pleiotropy"
)

net_res_pn$plot

ggsave(net_res_pn$plot, filename = sprintf("%s/pn_GSEA_gene_network.pdf",figure_file), 
       device = cairo_pdf, width = 10, height = 6.5, dpi = 600, bg = "white")




##GPS-M------------
pm_geneList <- agingGPS$pm_ld_aging
names(pm_geneList) <- agingGPS$gene
head(pm_geneList)
pm_geneList <- sort(pm_geneList, decreasing = TRUE)

set.seed(1234)
pm_gsego <- gseGO(geneList     = pm_geneList, 
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  keyType      = "SYMBOL",
                  minGSSize    = 10,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  scoreType = "pos",
                  eps = 0,
                  seed = T,
                  verbose      = FALSE)

##gene-centered network plot
net_res_pm <- make_gene_centered_gsea_network(
  gsea_obj = pm_gsego,
  gene_score = pm_geneList,
  top_n_pathways = 30,
  high_cutoff = min(agingGPS[pm_pleio5_aging == "5th"]$pm_ld_aging),
  top_n_genes = 20,
  min_pathways_per_gene = 2,
  pathway_label_width = 24,
  
  gene_radius = 0.3,
  pathway_radius = 0.60,
  pathway_label_offset = 0.07,
  
  pathway_text_size = 3.5,
  pathway_text_color = "grey20",
  pathway_node_fill = "#4A76A8",
  pathway_node_color = "white",
  
  gene_text_size = 4,
  edge_color = "grey10",
  edge_alpha = 0.22,
  edge_width_range = c(0.12, 0.65),
  
  gps_low = "#E9D8CC",
  gps_high = "#B8572A",
  gene_node_color = "white",
  gene_node_stroke = 0.55,
  gene_label_offset = 0.05,
  gps_legend_title = "GPS-M",
  
  plot_margin = ggplot2::margin(10, 10, 10, 100),
  
  save_prefix = "Temp/pm_GSEA_high_pleiotropy"
)

net_res_pm$plot


##output table--------------
add_high_core_prop <- function(gsea_obj, gene_score, high_cutoff = 20) {
  
  res <- as.data.table(gsea_obj@result)
  
  if (!"core_enrichment" %in% colnames(res)) {
    stop("No core_enrichment column found. This should be a GSEA result object.")
  }
  
  if (is.null(names(gene_score))) {
    stop("gene_score must be a named vector, and names(gene_score) should match genes in core_enrichment.")
  }
  
  high_genes <- names(gene_score)[
    !is.na(gene_score) & gene_score >= high_cutoff
  ]
  
  res[, row_id := .I]
  
  res[, c(
    "Number of core genes",
    "Number of core genes with high agingGPS",
    "Proportion of high agingGPS genes"
  ) := {
    
    genes <- strsplit(core_enrichment, "/", fixed = TRUE)[[1]]
    genes <- unique(genes)
    
    genes_with_score <- genes[genes %in% names(gene_score)]
    high_core <- genes_with_score[genes_with_score %in% high_genes]
    
    list(
      length(genes),
      length(high_core),
      ifelse(
        length(genes_with_score) > 0,
        length(high_core) / length(genes_with_score),
        NA_real_
      )
    )
    
  }, by = row_id]
  
  res[, row_id := NULL]
  
  result <- as.data.table(res)
  
  return(result)
}

pn_gsego_output <- add_high_core_prop(
  gsea_obj = pn_gsego,
  gene_score = pn_geneList,
  high_cutoff = min(agingGPS[pn_pleio5_aging == "5th"]$pn_ld_aging)
)
head(pn_gsego_output)
nrow(pn_gsego_output)

pm_gsego_output <- add_high_core_prop(
  gsea_obj = pm_gsego,
  gene_score = pm_geneList,
  high_cutoff = min(agingGPS[pm_pleio5_aging == "5th"]$pm_ld_aging)
)
head(pm_gsego_output)
nrow(pm_gsego_output)

gsego_outtable <- rbind(pn_gsego_output[, GPS := "agingGPS-N"],
                        pm_gsego_output[, GPS := "agingGPS-M"])

nrow(gsego_outtable)
head(gsego_outtable)

#3 disease enrichment-------
library(DOSE)
##GPS-N--------------
agingGPS <- fread("AgingGPS/agingGPS.csv", na.strings = c("", "NA"))
pn_geneList <- agingGPS$pn_ld_aging
names(pn_geneList) <- agingGPS$geneid
head(pn_geneList)
pn_geneList <- sort(pn_geneList, decreasing = TRUE)

gsea_do_pn <- gseDO(pn_geneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = 123,
                    verbose = FALSE, eps = 0, scoreType = "pos")


##new plot
source(sprintf("%s/disease_gsea_plot.R", code_file))
disease_plot_data_pn <- prepare_disease_gsea_contribution_data(
  gsea_obj = gsea_do_pn,
  gene_score = pn_geneList,
  high_cutoff = min(agingGPS[pn_pleio5_aging == "5th"]$pn_ld_aging),
  top_n_terms = 20,
  top_n_shared_genes = 20
)

p_disease_pn <- plot_disease_gsea_contribution(
  plot_data = disease_plot_data_pn
)

p_disease_pn

ggsave(p_disease_pn, filename = sprintf("%s/p_disease_pn.pdf",figure_file), 
       device = cairo_pdf, width = 11, height = 6.5, dpi = 600, bg = "white")


##GPS-M------------
pm_geneList <- agingGPS$pm_ld_aging
names(pm_geneList) <- agingGPS$geneid
head(pm_geneList)
pm_geneList <- sort(pm_geneList, decreasing = TRUE)

gsea_do_pm <- gseDO(pm_geneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = 123,
                    verbose = FALSE, eps = 0, scoreType = "pos")

##new plot
disease_plot_data_pm <- prepare_disease_gsea_contribution_data(
  gsea_obj = gsea_do_pm,
  gene_score = pm_geneList,
  high_cutoff = min(agingGPS[pm_pleio5_aging == "5th"]$pm_ld_aging),
  top_n_terms = 20,
  top_n_shared_genes = 20
)

p_disease_pm <- plot_disease_gsea_contribution(
  plot_data = disease_plot_data_pm
)

p_disease_pm


##output table--------------
pn_do_ouput <- add_high_core_prop(
  gsea_obj = gsea_do_pn,
  gene_score = pn_geneList,
  high_cutoff = min(agingGPS[pn_pleio5_aging == "5th"]$pn_ld_aging)
)

head(pn_do_ouput)
nrow(pn_do_ouput)

pm_do_ouput <- add_high_core_prop(
  gsea_obj = gsea_do_pm,
  gene_score = pm_geneList,
  high_cutoff = min(agingGPS[pm_pleio5_aging == "5th"]$pm_ld_aging)
)

head(pm_do_ouput)
nrow(pm_do_ouput)

do_outtable <- rbind(pn_do_ouput[, GPS := "agingGPS-N"],
                     pm_do_ouput[, GPS := "agingGPS-M"])

nrow(do_outtable)
head(do_outtable)


#4 gerogenes and gerosuppressor genes----------------------
##4.1 enrichment analysis------------------
gerogenes <- c("APOB", "APOE", "CDKN2A", "CDKN2B", "DBI", "ERVs", "GHR", "IGF1", "IL11", "LMNA", "LINE-1", "miR-29", "PCSK9")
gerosuppressor <- c("APOE", "ATM", "BNP", "BRCA1", "BRCA2", "DNMT3A", "ERCC6", "ERCC8",
                    "FOXO3A", "KL", "SH2B3", "SIRT6", "TET2", "TERT", "WRN", "ZMPSTE24")

aging_drug <- fread("AgingGPS/EFO_0022597-known-drugs-aging.tsv")
head(aging_drug)
drug_targets <- unique(aging_drug$symbol) #97

length(unique(c(gerogenes, gerosuppressor, drug_targets))) #124

##aginggene_new
aginggene_new <- fread("AgingGPS/gene-aging-mechanisms.tsv", sep = "\t", quote = "")
nrow(aginggene_new)
names(aginggene_new) <- c("gene", "hallmarks_of_aging")

known_aging_genes <- unique(c(trimws(gsub('"', '', aginggene_new$gene)),
                              gerogenes, gerosuppressor, drug_targets))
length(known_aging_genes) #2482


##agingGPS
agingGPS <- fread("/home/liumy/pleiotropy/aging_revision/results/AgingGPS/agingGPS.csv", na.strings = c("", "NA"))
agingGPS[, `:=`(if_gerogenes = ifelse(gene %in% gerogenes, 1, 0),
                if_gerosuppressor = ifelse(gene %in% gerosuppressor, 1, 0),
                if_gero = ifelse(gene %in% c(gerogenes, gerosuppressor), 1, 0),
                ifdrug_OTG = ifelse(gene %in% unique(aging_drug$symbol), 1, 0),
                ifknow_aging = ifelse(gene %in% known_aging_genes, 1, 0))]


fisher_test_func <- function(geneset1, geneset2, total_gene) {
  if(!all(c(geneset1, geneset2) %in% total_gene)) {
    stop("Some genes in the input lists are not present in the total_gene list!")
  }
  
  intersection_genes <- intersect(geneset1, geneset2)
  union_genes <- union(geneset1, geneset2)
  
  
  table_matrix <- matrix(c(
    length(intersection_genes),
    length(setdiff(geneset1, geneset2)),
    length(setdiff(geneset2, geneset1)),
    length(total_gene) - length(union_genes)
  ), nrow = 2, byrow = TRUE)
  
  
  print(table_matrix)
  fisher_test_result <- fisher.test(table_matrix)
  return(fisher_test_result)
}


total_gene <- agingGPS$gene #17631
geneset1 <- agingGPS[ifknow_aging == 1,]$gene

aging_enrichment_list_pn <- lapply(list(100,200,500,1000), function(x){
  fisher_test_func(geneset1 = geneset1,
                   geneset2 = agingGPS[order(-pn_ld_aging)]$gene[1:x],
                   total_gene = total_gene)
})

names(aging_enrichment_list_pn) <- paste0("pn_aging_top", c(100,200,500,1000))
aging_enrichment_list_pn$pn_aging_top500

aging_enrichment_list_pm <- lapply(list(100,200,500,1000), function(x){
  fisher_test_func(geneset1 = geneset1,
                   geneset2 = agingGPS[order(-pm_ld_aging)]$gene[1:x],
                   total_gene = total_gene)
})

names(aging_enrichment_list_pm) <- paste0("pm_aging_top", c(100,200,500,1000))
aging_enrichment_list_pm$pm_aging_top500

##output table
format_p <- function(p) {ifelse( p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p))}
extract_fisher <- function(x, name) {
  
  score <- ifelse(grepl("^pn_", name), "agingGPS-N", "agingGPS-M")
  top <- name
  
  data.table(
    top = top,
    score = score,
    OR = unname(x$estimate),
    OR_l95 = x$conf.int[1],
    OR_u95 = x$conf.int[2],
    p = x$p.value
  )
}

aging_enrichment_table_pn <- rbindlist(
  lapply(names(aging_enrichment_list_pn), function(nm) {
    extract_fisher(aging_enrichment_list_pn[[nm]], nm)
  })
)
head(aging_enrichment_table_pn)

aging_enrichment_table_pm <- rbindlist(
  lapply(names(aging_enrichment_list_pm), function(nm) {
    extract_fisher(aging_enrichment_list_pm[[nm]], nm)
  })
)
head(aging_enrichment_table_pm)

aging_enrichment_table <- rbind(aging_enrichment_table_pn, aging_enrichment_table_pm)
head(aging_enrichment_table)

aging_enrichment_table_format <- copy(aging_enrichment_table)
aging_enrichment_table_format[, OR := sprintf("%.2f", OR)]
aging_enrichment_table_format[, OR_l95 := sprintf("%.2f", OR_l95)]
aging_enrichment_table_format[, OR_u95 := sprintf("%.2f", OR_u95)]
aging_enrichment_table_format[, p := format_p(p)]
aging_enrichment_table_format

##4.2 module analysis (pn_ld)-------------
###4.2.1 load data------------
library(STRINGdb); library(igraph); library(clusterProfiler); library(org.Hs.eg.db);library(ReactomePA)
genes <- c(agingGPS[if_gero == 1 & pn_pleio3_aging %in% c("High pleiotropy"),]$gene,
           agingGPS[ifdrug_OTG == 1 & pn_pleio3_aging %in% c("High pleiotropy")]$gene)
length(genes) #25

string_db <- STRINGdb(version="12", species=9606, score_threshold=700, input_directory="AgingGPS/string_db") #high confidence (>0.7) interactions
map <- string_db$map(data.frame(gene=genes), "gene", removeUnmappedRows = TRUE)
head(map)
nrow(map) #29
hits <- string_db$get_interactions(map$STRING_id)
head(hits)
nrow(hits) #92
ids_set <- unique(map$STRING_id) #29
hits <- subset(hits, from %in% ids_set & to %in% ids_set)
head(hits)
nrow(hits) #92

###4.2.2 Build graph & detect modules-------------------
library(igraph); library(ggraph); library(tidygraph); library(scales); library(scatterpie)

id2sym <- setNames(map$gene, map$STRING_id)
head(id2sym)

edge_df <- data.frame(
  a = unname(id2sym[hits$from]),
  b = unname(id2sym[hits$to]),
  score = if ("combined_score" %in% names(hits)) hits$combined_score/1000 else 1
)
edge_df <- subset(edge_df, !is.na(a) & !is.na(b) & a != b)
edge_df <- unique(edge_df)
head(edge_df)

edge_df2 <- edge_df
names(edge_df2)[names(edge_df2)=="score"] <- "weight"
edge_df2$a <- as.character(edge_df2$a)
edge_df2$b <- as.character(edge_df2$b)
head(edge_df2)

g <- igraph::graph_from_data_frame(edge_df2, directed = FALSE)
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                      edge.attr.comb = list(weight = "max"))

print(g, e = TRUE, v = TRUE)
cl <- cluster_louvain(g)

V(g)$module <- as.integer(membership(cl))
V(g)$degree <- 1
V(g)$module_lab <- paste0("Module ", V(g)$module)

sym <- V(g)$name
is_gero <- as.integer(sym %in% gerogenes)
is_gs <- as.integer(sym %in% gerosuppressor)
is_target <- as.integer(sym %in% drug_targets)

xy <- igraph::layout_with_fr(g, weights = E(g)$weight, niter = 2000, grid = "nogrid")
xy = xy*2

node_df <- data.frame(name = sym, module = V(g)$module, module_lab = V(g)$module_lab, degree = V(g)$degree,
                      x = xy[,1], y = xy[,2], gero = is_gero, gs = is_gs, target = is_target, stringsAsFactors = FALSE)
node_df <- node_df[match(V(g)$name, node_df$name), ]
head(node_df)

ea <- as_ids(head_of(g, E(g)))
eb <- as_ids(tail_of(g, E(g)))
mod_map <- setNames(node_df$module_lab, node_df$name)
E(g)$edge_col <- ifelse(mod_map[ea] == mod_map[eb], mod_map[ea], "cross")
E(g)$weight   <- ifelse(is.na(E(g)$weight), 0.5, E(g)$weight)

g_tbl <- as_tbl_graph(g)
mod_levels <- sort(unique(node_df$module_lab))
mod_pal <- setNames(RColorBrewer::brewer.pal(max(3, length(mod_levels)), "Set1")[seq_along(mod_levels)], mod_levels)
mod_pal <- c(mod_pal, cross = "white")
pie_cols <- c(gero="#67C1ECFF", gs="#7BBC53FF", target="#DE6736FF")

set.seed(123)
node_df$degree_use <- 5
plot_pn <- ggraph(g_tbl, layout = "manual", x = node_df$x, y = node_df$y) + theme_void() +
  geom_edge_link(aes(colour = edge_col, width = weight),
                 alpha = 0.5, lineend = "round", show.legend = TRUE) +
  scale_edge_width(range = c(0.2, 1.2), guide = "none") +
  scale_edge_colour_manual(values = rep("grey50", 6), name = "") +
  geom_scatterpie(data = node_df,
                  aes(x = x, y = y, r = scales::rescale(degree_use, to = c(0.1, 0.3))),
                  cols = c("gero","gs","target"),
                  color = "white", lwd = 0.3, alpha = 0.95) +
  scale_fill_manual(values = pie_cols, name = "",
                    labels = c(gero = "Gerogene", gs = "Gerosuppressor", target = "Drug target")) +
  
  geom_text(data = node_df, aes(x = x, y = y, label = name), size = 4, vjust = -1.5) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 4)))

plot_pn
ggsave(plot = plot_pn, width = 10, height = 10, device = cairo_pdf,
       filename = sprintf("%s/PPI_plot_pn.pdf", figure_file))



##4.3 module analysis (pm_ld)-------------
###4.3.1 load data------------
genes <- c(agingGPS[if_gero == 1 & pm_pleio3_aging %in% c("High pleiotropy"),]$gene,
           agingGPS[ifdrug_OTG == 1 & pm_pleio3_aging %in% c("High pleiotropy")]$gene)
length(genes) #25

string_db <- STRINGdb(version="12", species=9606, score_threshold=700, input_directory="AgingGPS/string_db") #high confidence (>0.7) interactions
map <- string_db$map(data.frame(gene=genes), "gene", removeUnmappedRows = TRUE)
head(map)
nrow(map) #25
hits <- string_db$get_interactions(map$STRING_id)
head(hits)
nrow(hits) #64
ids_set <- unique(map$STRING_id) #31
hits <- subset(hits, from %in% ids_set & to %in% ids_set)
head(hits)

###4.3.2 Build graph & detect modules------------------
id2sym <- setNames(map$gene, map$STRING_id)
head(id2sym)

edge_df <- data.frame(
  a = unname(id2sym[hits$from]),
  b = unname(id2sym[hits$to]),
  score = if ("combined_score" %in% names(hits)) hits$combined_score/1000 else 1
)
edge_df <- subset(edge_df, !is.na(a) & !is.na(b) & a != b)
edge_df <- unique(edge_df)
head(edge_df)

edge_df2 <- edge_df
names(edge_df2)[names(edge_df2)=="score"] <- "weight"
edge_df2$a <- as.character(edge_df2$a)
edge_df2$b <- as.character(edge_df2$b)
head(edge_df2)

g <- igraph::graph_from_data_frame(edge_df2, directed = FALSE)
g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                      edge.attr.comb = list(weight = "max"))

print(g, e = TRUE, v = TRUE)
cl <- cluster_louvain(g)

V(g)$module <- as.integer(membership(cl))
V(g)$degree <- degree(g)
V(g)$module_lab <- paste0("Module ", V(g)$module)

sym <- V(g)$name
is_gero <- as.integer(sym %in% gerogenes)
is_gs <- as.integer(sym %in% gerosuppressor)
is_target <- as.integer(sym %in% drug_targets)

xy <- igraph::layout_with_fr(g, weights = E(g)$weight, niter = 2000, grid = "nogrid")
xy <- xy * 2

node_df <- data.frame(name = sym, module = V(g)$module, module_lab = V(g)$module_lab, degree = V(g)$degree,
                      x = xy[,1], y = xy[,2], gero = is_gero, gs = is_gs, target = is_target, stringsAsFactors = FALSE)
node_df <- node_df[match(V(g)$name, node_df$name), ]
head(node_df)

ea <- as_ids(head_of(g, E(g)))
eb <- as_ids(tail_of(g, E(g)))
mod_map <- setNames(node_df$module_lab, node_df$name)
E(g)$edge_col <- ifelse(mod_map[ea] == mod_map[eb], mod_map[ea], "cross")
E(g)$weight   <- ifelse(is.na(E(g)$weight), 0.5, E(g)$weight)

g_tbl <- as_tbl_graph(g)

mod_levels <- sort(unique(node_df$module_lab))
mod_pal <- setNames(RColorBrewer::brewer.pal(max(3, length(mod_levels)), "Set1")[seq_along(mod_levels)], mod_levels)
mod_pal <- c(mod_pal, cross = "white")
pie_cols <- c(gero="#67C1ECFF", gs="#7BBC53FF", target="#DE6736FF")


set.seed(123)
node_df$degree_use <- 5
plot_pm <- ggraph(g_tbl, layout = "manual", x = node_df$x, y = node_df$y) + theme_void() +
  
  geom_edge_link(aes(colour = edge_col, width = weight),
                 alpha = 0.5, lineend = "round", show.legend = TRUE) +
  scale_edge_width(range = c(0.2, 1.2), guide = "none") +
  scale_edge_colour_manual(values = rep("grey50",7), name = "") +
  geom_scatterpie(data = node_df,
                  aes(x = x, y = y, r = scales::rescale(degree_use, to = c(0.1, 0.3))),
                  cols = c("gero","gs","target"),
                  color = "white", lwd = 0.3, alpha = 0.95) +
  scale_fill_manual(values = pie_cols, name = "",
                    labels = c(gero = "Gerogene", gs = "Gerosuppressor", target = "Drug target")) +
  
  geom_text(data = node_df, aes(x = x, y = y, label = name), size = 5, vjust = -1.5) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 5)))

plot_pm
ggsave(plot = plot_pm, width = 10, height = 10, device = cairo_pdf,
       filename = sprintf("%s/PPI_plot_pm.pdf", figure_file))


