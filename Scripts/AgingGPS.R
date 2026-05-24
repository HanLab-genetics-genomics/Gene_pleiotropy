#1 load packages and set variables-------
setwd("/path/to/Gene_pleiotropy-main/Data")
figure_file <- "/path/to/Gene_pleiotropy-main/Output/Aging"
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
load("AgingGPS/hopsgene_aging12new_cut.RData")
head(hopsgene_ageburden_cut$pm_ld)
agingGPS <- merge(hopsgene_ageburden_cut$pm_ld[, .(gene, pm_ld_aging = use_score, pm_pleio3_aging = pleio_class3,
                                                   pm_pleio5_aging = pleio_class5, pm_pleio10_aging = pleio10)],
                  hopsgene_ageburden_cut$pn_ld[, .(gene, pn_ld_aging = use_score, pn_pleio3_aging = pleio_class3,
                                                   pn_pleio5_aging = pleio_class5, pn_pleio10_aging = pleio10)],
                  by = "gene")
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
# Note: tied gene-level scores may still cause minor fluctuations in GSEA p-values,
# adjusted p-values and leading-edge genes. We set a random seed to improve
# reproducibility, but the emapplot layout may still vary slightly across environments.
library(clusterProfiler)
library(enrichplot)
library(GOSemSim)
library(ggplot2)

##GPS-N------------
set.seed(1234)
pn_geneList <- hopsgene_ageburden_cut$pn_ld$use_score
names(pn_geneList) <- hopsgene_ageburden_cut$pn_ld$gene
head(pn_geneList)
pn_geneList <- sort(pn_geneList, decreasing = TRUE)

pn_gsego <- gseGO(geneList     = pn_geneList, 
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  keyType      = "SYMBOL",
                  scoreType = "pos",
                  minGSSize    = 100,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  eps = 0,
                  verbose      = FALSE)

colnames(pn_gsego@result)


pn_gsego_plot <- pn_gsego
pn_gsego_plot@result <- pn_gsego_plot@result[!is.na(pn_gsego_plot@result$p.adjust) & pn_gsego_plot@result$p.adjust < 0.05, ]
pn_gsego_plot@result <- pn_gsego_plot@result[order(pn_gsego_plot@result$p.adjust), ]
bp_sem <- godata("org.Hs.eg.db", ont = "BP")

pn_gsego_plot2 <- pairwise_termsim(pn_gsego_plot, method = "Wang", semData = bp_sem, showCategory = 30)
set.seed(1234)
pn_gsego_emapplot <- emapplot(pn_gsego_plot2, showCategory = 30, color = "p.adjust")
pn_gsego_emapplot

ggsave(plot = pn_gsego_emapplot, width = 10, height = 10, device = cairo_pdf,
       filename = sprintf("%s/pn_gsego_emapplot.pdf", figure_file))



##GPS-M------------
pm_geneList <- hopsgene_ageburden_cut$pm_ld$use_score
names(pm_geneList) <- hopsgene_ageburden_cut$pm_ld$gene
head(pm_geneList)
pm_geneList <- sort(pm_geneList, decreasing = TRUE)

pm_gsego <- gseGO(geneList     = pm_geneList, 
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  keyType      = "SYMBOL",
                  minGSSize    = 100,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  scoreType = "pos",
                  eps = 0,
                  verbose      = FALSE)

colnames(pm_gsego@result)

pm_gsego_plot <- pm_gsego
pm_gsego_plot@result <- pm_gsego_plot@result[!is.na(pm_gsego_plot@result$p.adjust) & pm_gsego_plot@result$p.adjust < 0.05, ]
pm_gsego_plot@result <- pm_gsego_plot@result[order(pm_gsego_plot@result$p.adjust), ]
bp_sem <- godata("org.Hs.eg.db", ont = "BP")

pm_gsego_plot2 <- pairwise_termsim(pm_gsego_plot, method = "Wang", semData = bp_sem, showCategory = 30)
pm_gsego_emapplot <- emapplot(pm_gsego_plot2, showCategory = 30, color = "p.adjust")
pm_gsego_emapplot

ggsave(plot = pm_gsego_emapplot, width = 10, height = 10, device = cairo_pdf,
       filename = sprintf("%s/pm_gsego_emapplot.pdf", figure_file)) 

#3 disease enrichment-------
##3.1 DOSE analysis--------------
# Note: tied gene-level scores may still cause minor fluctuations in GSEA p-values,
# adjusted p-values and leading-edge genes. We set a random seed to improve
# reproducibility, but the emapplot layout may still vary slightly across environments.
library(DOSE)
###pn_ld--------------
usedata_pn <- hopsgene_ageburden_cut$pn_ld
genelist_gsea_pn <- usedata_pn[order(-use_score)]$use_score
names(genelist_gsea_pn) <- usedata_pn[order(-use_score)]$geneid

gsea_do_pn <- gseDO(genelist_gsea_pn, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = 123,
                    verbose = FALSE, eps = 0, scoreType = "pos")
#head(gsea_do_pn)
gsea_do_result_pn <- setDT(gsea_do_pn@result)
nrow(gsea_do_result_pn[qvalue < 0.05,]) #147 
#gsea_do_result[qvalue < 0.05,] %>% head()

gsea_do_result_pn$GPS <- "agingGPS-N"
head(gsea_do_result_pn)


####dotplot by ggplot (use)---------------------
library(stringr)

df_top_pn <- gsea_do_result_pn[order(p.adjust)][1:20] %>%
  .[order(-NES)]
head(df_top_pn)
df_top_pn$Description <- factor(df_top_pn$Description, levels = rev(df_top_pn$Description))

pn_do_dotplot <- ggplot(df_top_pn, aes(x = NES, y = Description, color = p.adjust)) +
  geom_point(size = 6) +
  labs(x = "Normalized enrichment scores", y = NULL, color = "p.adjust") +
  theme_bw() +
  scale_color_gradient(low='#DA7235', high='#7FB65D') +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25)) +
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        legend.position = "right",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

pn_do_dotplot

ggsave(plot = pn_do_dotplot, width = 6, height = 8, device = cairo_pdf,
       filename = sprintf("%s/pn_do_dotplot_defultgseDO_ggplot.pdf", figure_file))


###pm_ld------------
usedata_pm <- hopsgene_ageburden_cut$pm_ld
genelist_gsea_pm <- usedata_pm[order(-use_score)]$use_score
names(genelist_gsea_pm) <- usedata_pm[order(-use_score)]$geneid

gsea_do_pm <- gseDO(genelist_gsea_pm, pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = 123,
                    verbose = FALSE, eps = 0, scoreType = "pos")
#head(gsea_do)
gsea_do_result_pm <- setDT(gsea_do_pm@result)
nrow(gsea_do_result_pm[qvalue < 0.05,]) #92 --> 102
#gsea_do_result[qvalue < 0.05,] %>% head()

gsea_do_result_pm$GPS <- "agingGPS-M"

gsea_do_result_pmpn <- rbind(gsea_do_result_pn, gsea_do_result_pm)
fwrite(gsea_do_result_pmpn, file = sprintf("%s/gsea_do_result_defultgseDO.csv", figure_file))


####dotplot by ggplot (use)---------------------
df_top_pm <- gsea_do_result_pm[order(p.adjust)][1:20] %>%
  .[order(-NES)]
head(df_top_pm)
df_top_pm$Description <- factor(df_top_pm$Description, levels = rev(df_top_pm$Description))

pm_do_dotplot <- ggplot(df_top_pm, aes(x = NES, y = Description, color = p.adjust)) +
  geom_point(size = 6) +
  labs(x = "Normalized enrichment scores", y = NULL, color = "p.adjust") +
  theme_bw() +
  scale_color_gradient(low='#DA7235', high='#7FB65D') +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25)) +
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        legend.position = "right",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

pm_do_dotplot

ggsave(plot = pm_do_dotplot, width = 6, height = 8, device = cairo_pdf,
       filename = sprintf("%s/pm_do_dotplot_defultgseDO_ggplot.pdf", figure_file))


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
aginggene_new <- fread("AgingGPS/gene-aging-mechanisms.tsv", 
                       sep = "\t", quote = "")
nrow(aginggene_new)
names(aginggene_new) <- c("gene", "hallmarks_of_aging")

known_aging_genes <- unique(c(trimws(gsub('"', '', aginggene_new$gene)),
                              gerogenes, gerosuppressor, drug_targets))
length(known_aging_genes) #2482

for(i in c("pm_ld", "pn_ld")) {
  hopsgene_ageburden_cut[[i]][, `:=`(
    if_gerogenes = ifelse(gene %in% gerogenes, 1, 0),
    if_gerosuppressor = ifelse(gene %in% gerosuppressor, 1, 0),
    if_gero = ifelse(gene %in% c(gerogenes, gerosuppressor), 1, 0)
  )]}

for(i in c("pm_ld", "pn_ld")) {
  hopsgene_ageburden_cut[[i]][, `:=`(ifdrug_OTG = ifelse(gene %in% unique(aging_drug$symbol), 1, 0))]}

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

#total_gene <- unique(c(known_aging_genes, hopsgene_ageburden_cut$pn_ld$gene)) #17928
total_gene <- hopsgene_ageburden_cut$pn_ld$gene #17631

aging_enrichment <- list(
  pn_aging_top500 = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                     geneset2 = hopsgene_ageburden_cut$pn_ld[order(-use_score)]$gene[1:500],
                                     total_gene = total_gene),
  pm_aging_top500 = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                     geneset2 = hopsgene_ageburden_cut$pm_ld[order(-use_score)]$gene[1:500],
                                     total_gene = total_gene),
  
  pn_aging_top200 = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                     geneset2 = hopsgene_ageburden_cut$pn_ld[order(-use_score)]$gene[1:200],
                                     total_gene = total_gene),
  pm_aging_top200 = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                     geneset2 = hopsgene_ageburden_cut$pm_ld[order(-use_score)]$gene[1:200],
                                     total_gene = total_gene),
  
  pn_aging_top100 = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                     geneset2 = hopsgene_ageburden_cut$pn_ld[order(-use_score)]$gene[1:100],
                                     total_gene = total_gene),
  pm_aging_top100 = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                     geneset2 = hopsgene_ageburden_cut$pm_ld[order(-use_score)]$gene[1:100],
                                     total_gene = total_gene),
  
  pn_aging_top1000 = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                      geneset2 = hopsgene_ageburden_cut$pn_ld[order(-use_score)]$gene[1:1000],
                                      total_gene = total_gene),
  pm_aging_top1000 = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                      geneset2 = hopsgene_ageburden_cut$pm_ld[order(-use_score)]$gene[1:1000],
                                      total_gene = total_gene),
  
  pn_aging_high = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                   geneset2 = hopsgene_ageburden_cut$pn_ld[pleio_class3 == "High pleiotropy",]$gene,
                                   total_gene = total_gene),
  pm_aging_high = fisher_test_func(geneset1 = intersect(known_aging_genes, hopsgene_ageburden_cut$pm_ld$gene),
                                   geneset2 = hopsgene_ageburden_cut$pm_ld[pleio_class3 == "High pleiotropy",]$gene,
                                   total_gene = total_gene))


aging_enrichment$pn_aging_top1000$p.value
aging_enrichment$pm_aging_top1000$p.value


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

aging_enrichment_table <- rbindlist(
  lapply(names(aging_enrichment), function(nm) {
    extract_fisher(aging_enrichment[[nm]], nm)
  })
)

aging_enrichment_table_format <- copy(aging_enrichment_table)
aging_enrichment_table_format[, OR := sprintf("%.2f", OR)]
aging_enrichment_table_format[, OR_l95 := sprintf("%.2f", OR_l95)]
aging_enrichment_table_format[, OR_u95 := sprintf("%.2f", OR_u95)]
aging_enrichment_table_format[, p := format_p(p)]
aging_enrichment_table_format


##4.2 module analysis (pm_ld)-------------
###4.2.1 load data------------
library(STRINGdb); library(igraph); library(clusterProfiler); library(org.Hs.eg.db);library(ReactomePA)
genes <- c(hopsgene_ageburden_cut$pm_ld[if_gero == 1 & pleio_class3 %in% c("High pleiotropy"),]$gene,
           hopsgene_ageburden_cut$pm_ld[ifdrug_OTG == 1 & pleio_class3 %in% c("High pleiotropy")]$gene)
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

###4.2.2 Build graph & detect modules------------------
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


##4.3 module analysis (pn_ld)-------------
###4.3.1 load data------------
genes <- c(hopsgene_ageburden_cut$pn_ld[if_gero == 1 & pleio_class3 %in% c("High pleiotropy"),]$gene,
           hopsgene_ageburden_cut$pn_ld[ifdrug_OTG == 1 & pleio_class3 %in% c("High pleiotropy")]$gene)
length(genes) #29

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

###4.3.2 Build graph & detect modules-------------------
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


