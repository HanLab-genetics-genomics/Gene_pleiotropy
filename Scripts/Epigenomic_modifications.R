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
bed_file <- "/home/liumy/pleiotropy/github/data/bed_data"
lola_out_file <- "/home/liumy/pleiotropy/github/lola"
#2 load data--------------------------
load("pleiotropy_maindata.RData")
#3 figure5--------------------
##enrichment by lola----------------- 
###Extract gene coordinates ± 1000bp (could skip and load the data provided)-------------------
genecoordinate_func <- function(x, data = pleiotropy_maindata$pm_ld, 
                                var_name = "pleio100", out_name = "pm", ifidentical = F,
                                bedfile = bed_file){
  
  if(ifidentical == F & x[1] <= 50){
    genes <- data[, use_var := get(var_name) %>% as.numeric(.)] %>% .[use_var <= x[1], gene]
  } else if(ifidentical == F & x[1] > 50){genes <- data[, use_var := get(var_name) %>% as.numeric(.)] %>% .[use_var >= x[1], gene]
  } else if(ifidentical == T){genes <- data[, use_var := get(var_name) %>% as.numeric(.)] %>% .[use_var %in% x, gene]}
  
  gene_coords <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
                       filters = "hgnc_symbol",
                       values = genes,
                       mart = mart)
  bed_data <- as.data.table(gene_coords) %>% 
    .[, .(chrom = as.numeric(chromosome_name),
          start = start_position - 1 - 1000,  # BED format uses 0-based start
          end = end_position + 1000,
          name = hgnc_symbol)] %>%
    .[, start := ifelse(start < 0, 0, start)] %>%
    na.omit(.) %>%
    .[!duplicated(name)]
  
  
  bed_data$chrom <- paste0("chr", bed_data$chrom, sep = "")
  
  fwrite(bed_data, file = sprintf("%s/%s.bed", bedfile, out_name), 
         quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(bed_data)
}
####GRCh38---------------------
library(biomaRt); library(dbplyr)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

pm_bed <- list()
for(i in c(20)){
  top_bed_name = paste0("pm_top", i, "_hg38")
  pm_bed[[top_bed_name]] <- genecoordinate_func(x = 101-i, var_name =  "pleio100", bedfile = bed_file,
                                                data = pleiotropy_maindata$pm_ld, out_name = top_bed_name)
  bottom_bed_name = paste0("pm_bottom", i, "_hg38")
  pm_bed[[bottom_bed_name]] <- genecoordinate_func(x = i, var_name =  "pleio100", 
                                                   data = pleiotropy_maindata$pm_ld, out_name = bottom_bed_name)
}

pn_bed <- list()
for(i in c(20)){
  top_bed_name = paste0("pn_top", i, "_hg38")
  pn_bed[[top_bed_name]] <- genecoordinate_func(x = 101-i, var_name =  "pleio100", 
                                                data = pleiotropy_maindata$pn_ld, out_name = top_bed_name)
  bottom_bed_name = paste0("pn_bottom", i, "_hg38")
  pn_bed[[bottom_bed_name]] <- genecoordinate_func(x = i, var_name =  "pleio100", 
                                                   data = pleiotropy_maindata$pn_ld, out_name = bottom_bed_name)
}

age_bed <- list()
for(i in c(1:4)){
  bed_name = paste0("age4_", i, "_hg38")
  age_bed[[bed_name]] <- genecoordinate_func(x = i, var_name =  "age_stage4_num", 
                                             data = pleiotropy_maindata$pn_ld, out_name = bed_name, ifidentical = T)
}

####GRCh37---------------------
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
pm_bed <- list()
for(i in c(10, 20, 30, 40 ,50)){
  top_bed_name = paste0("pm_top", i, "_hg19")
  pm_bed[[top_bed_name]] <- genecoordinate_func(x = 101-i, var_name =  "pleio100", 
                                                data = pleiotropy_maindata$pm_ld, out_name = top_bed_name)
  bottom_bed_name = paste0("pm_bottom", i, "_hg19")
  pm_bed[[bottom_bed_name]] <- genecoordinate_func(x = i, var_name =  "pleio100", 
                                                   data = pleiotropy_maindata$pm_ld, out_name = bottom_bed_name)
}

pn_bed <- list()
for(i in c(10, 20, 30, 40 ,50)){
  top_bed_name = paste0("pn_top", i, "_hg19")
  pn_bed[[top_bed_name]] <- genecoordinate_func(x = 101-i, var_name =  "pleio100", 
                                                data = pleiotropy_maindata$pn_ld, out_name = top_bed_name)
  bottom_bed_name = paste0("pn_bottom", i, "_hg19")
  pn_bed[[bottom_bed_name]] <- genecoordinate_func(x = i, var_name =  "pleio100", 
                                                   data = pleiotropy_maindata$pn_ld, out_name = bottom_bed_name)
}

age_bed <- list()
for(i in c(1:4)){
  bed_name = paste0("age4_", i, "_hg19")
  age_bed[[bed_name]] <- genecoordinate_func(x = i, var_name =  "age_stage4_num", 
                                             data = pleiotropy_maindata$pn_ld, out_name = bed_name, ifidentical = T)
}

##lola enrichment for LOLA core dataset------------------
library(GenomicRanges); library(LOLA)
regionDB = loadRegionDB(dbLocation = "LOLACore") 

lola_enrich_func <- function(regiondb_use, ifout = TRUE, 
                             bed1_file = "pn_1_hg19.bed",
                             bed2_file = "pn_100_hg19.bed",
                             out_file = "pn1to100_3d"){
  
  firstregions <- readBed(bed1_file)
  secondregions <- readBed(bed2_file)
  
  userSets = GRangesList(firstregions, secondregions)
  restrictedUniverse = buildRestrictedUniverse(userSets)
  
  locresults = runLOLA(userSets, restrictedUniverse, regionDB = regiondb_use, redefineUserSets=TRUE, cores = 15)
  
  if(ifout == TRUE) {
    writeCombinedEnrichment(locresults, outFolder = out_file , includeSplits=TRUE)
  }
  
  return(locresults)
}

lola_enrichment <- list(
  pm_topbottom = lapply(as.list(c(20)), function(x){
    lola_enrich_func(regiondb_use = regionDB, ifout = TRUE, 
                     bed1_file = sprintf("%s/pm_bottom%s_hg38.bed", bed_file, x),
                     bed2_file = sprintf("%s/pm_top%s_hg38.bed", bed_file, x),
                     out_file = sprintf("%s/pm_topbottom%s_lola", lola_out_file, x))}),
  pn_topbottom = lapply(as.list(c(20)), function(x){
    lola_enrich_func(regiondb_use = regionDB, ifout = TRUE, 
                     bed1_file = sprintf("%s/pn_bottom%s_hg38.bed", bed_file, x),
                     bed2_file = sprintf("%s/pn_top%s_hg38.bed", bed_file, x),
                     out_file = sprintf("%s/pn_topbottom%s_lola", lola_out_file, x))}),
  age4_class = lapply(as.list(c(1,2,3)), function(x){
    lola_enrich_func(regiondb_use = regionDB, ifout = TRUE, 
                     bed1_file = sprintf("%s/age4_4_hg38.bed", bed_file),
                     bed2_file = sprintf("%s/age4_%s_hg38.bed", bed_file, x),
                     out_file = sprintf("%s/age4_%sto4_lola", lola_out_file, x))})
  
)

###supplementary figrue 24: extract enrichment results by lola for UCSC database------------------
plotTopLOLAEnrichments_lmy <- function(data = ucsc_enrich, list_name, 
                                       color_pair = c("#FFAB8C", "#FFD4AE"), label_pair = c("Euteleostomi", "Eutheria")) {
  
  data <- data[!duplicated(description), class := ifelse(userSet == 1, "low","high")] %>%
    .[, .SD[order(-oddsRatio)], by = "class"]
  
  plot <- ggplot(data, aes(x = factor(description, levels = rev(unique(description))), y=oddsRatio, fill = class)) + theme_bw() +
    geom_bar(stat = "identity", position = "identity", color = "grey50") +
    scale_fill_manual(values = c("high" = color_pair[1], "low" = color_pair[2]),
                      labels = c("high" = label_pair[1], "low" = label_pair[2])) +  
    coord_flip() +
    labs(y = "Odds ratio", title = "", x = "") + 
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          plot.title = element_text(size = 12, face = "bold", color = "black"),
          axis.title.y = element_text(size = 15, face = "bold"), 
          axis.title.x = element_text(size = 15, face = "bold"), 
          axis.text.y = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.title = element_blank(),
          legend.position = "top") 
  
  return(plot)
}

ucsc_enrich_list <- list(
  pm = fread(sprintf("%s/pm_topbottom20_lola/col_ucsc_features.tsv", lola_out_file)) %>% .[qValue < 0.05,], 
  pn = fread(sprintf("%s/pn_topbottom20_lola/col_ucsc_features.tsv", lola_out_file)) %>% .[qValue < 0.05,],
  age = fread(sprintf("%s/age4_1to4_lola/col_ucsc_features.tsv", lola_out_file)) %>% .[qValue < 0.05,])

ucsc_enrich_plot <- list(
  pn = plotTopLOLAEnrichments_lmy(data = ucsc_enrich_list$pn, list_name, 
                                  color_pair = c("#7F8FCF", "#A4E3FE"), label_pair = c("H-GPS-N", "L-GPS-N")),
  pm = plotTopLOLAEnrichments_lmy(data = ucsc_enrich_list$pm, list_name, 
                                  color_pair = c("#7F8FCF", "#A4E3FE"), label_pair = c("H-GPS-M", "L-GPS-M")),
  age = plotTopLOLAEnrichments_lmy(data = ucsc_enrich_list$age, list_name, 
                                   color_pair = c("#FFAB8C", "#FFD4AE"), label_pair = c("Euteleostomi", "Eutheria"))
)

ucsc_enrich_plot$pn
ucsc_enrich_plot$age

ggsave(plot = ucsc_enrich_plot$pn, width = 7, height = 5, device = cairo_pdf,
       filename = sprintf("%s/ucsc_enrich_pn20.pdf",figure_file))
ggsave(plot = ucsc_enrich_plot$age, width = 7, height = 5, device = cairo_pdf,
       filename = sprintf("%s/ucsc_enrich_age.pdf",figure_file))

###figrue5A&5B: extract enrichment results by lola for Roadmap epigenomics region using beeswarm plot----------------------
Roadmapepi_enrich_list <- list(
  pm = fread(sprintf("%s/pm_topbottom20_lola/col_roadmap_epigenomics.tsv", lola_out_file)) %>% .[qValue < 0.05,], 
  pn = fread(sprintf("%s/pn_topbottom20_lola/col_roadmap_epigenomics.tsv", lola_out_file)) %>% .[qValue < 0.05,],
  age = fread(sprintf("%s/age4_1to4_lola/col_roadmap_epigenomics.tsv", lola_out_file)) %>% .[qValue < 0.05,]
)

library(ggbeeswarm); library(RColorBrewer)
beeswarm_epifunc <- function(data = Roadmapepi_enrich_list$pm, userset_num = 1){
  
  data <- data[userSet == userset_num, ]
  
  top_10_epi <- data[, `:=`(Significance = round(-log10(qValue), digits = 2),
                            family_1 = ifelse(is.na(antibody) & substr(filename, 6, 10) == "DNase", "DHS", antibody))] %>%
    .[, .N, by = family_1] %>% .[order(-N)] %>% .[1:10, family_1] 
  
  data_cut <- data[, `:=`(Significance = round(-log10(qValue), digits = 2),
                          family_1 = ifelse(is.na(antibody) & substr(filename, 6, 10) == "DNase", "DHS", antibody))] %>%
    .[, family_2 := ifelse(family_1 %in% top_10_epi, family_1, "Other")] %>%
    .[, family_2 := ifelse(qValue < 0.05 & family_2 == "Other", "NS", family_2)] %>%
    .[order(-Significance)]
  
  point1_markname <- data_cut[, .N, by = family_2][N == 1, family_2]
  data_cut <- data_cut[!family_2 %in% point1_markname, ]
  
  xlabel <- unique(data_cut[!family_2 %in% c("NS", "Other"),]$family_2)
  class.colors <- lighten(colorRampPalette(brewer.pal(12, "Paired"), space="Lab")(10),0.3)
  names(class.colors) <- c(xlabel)
  
  top10_enriched_epiplot <- data_cut %>% .[!family_2 %in% c("NS", "Other"),] %>%
    .[, family_2 := factor(family_2, levels = xlabel)] %>%
    ggplot(aes(x = family_2, y = Significance, color = family_2, label = family_2)) +
    geom_beeswarm(size = 3.5, dodge.width = 0.333, priority='density') +
    scale_color_manual(values = class.colors, labels = names(class.colors)) +
    # scale_colour_manual(values = cell.type.colors) +
    theme_classic() +
    xlab("") +
    ylab("-log10(qvalue) for enrichment") +
    theme(legend.position = "none") +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"), 
          axis.title.x = element_text(size = 15, face = "bold"), 
          axis.text.x = element_text(size = 12, color = "black", face = "bold"),
          axis.text.y = element_text(hjust=1, size = 12),
          axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.title = element_blank())
  
  return(top10_enriched_epiplot)
}

roadmapepi_beeplot <- list(
  pm_epi_bottom20 = beeswarm_epifunc(data = Roadmapepi_enrich_list$pm, userset_num = 1),
  pm_epi_top20 = beeswarm_epifunc(data = Roadmapepi_enrich_list$pm, userset_num = 2),
  pn_epi_bottom20 = beeswarm_epifunc(data = Roadmapepi_enrich_list$pn, userset_num = 1),
  pn_epi_top20 = beeswarm_epifunc(data = Roadmapepi_enrich_list$pn, userset_num = 2),
  age_epi_young = beeswarm_epifunc(data = Roadmapepi_enrich_list$age, userset_num = 1),
  age_epi_old = beeswarm_epifunc(data = Roadmapepi_enrich_list$age, userset_num = 2)
)

roadmapepi_beeplot$pn_epi_bottom20
roadmapepi_beeplot$pn_epi_top20

ggsave(plot = roadmapepi_beeplot$pn_epi_top20, width = 10, height = 5, device = cairo_pdf,
       filename = sprintf("%s/roadmap_pn_top20.pdf",figure_file))
ggsave(plot = roadmapepi_beeplot$pn_epi_bottom20, width = 3, height = 5, device = cairo_pdf,
       filename = sprintf("%s/roadmap_pn_bottom20.pdf",figure_file))
##3D genomic features---------------
###lola enrichment for 3d genomic features-------------------
regionDB = loadRegionDB("LOLA3D")

lola_enrichment <- list(
  pm_topbottom = lapply(as.list(c(10,20,30,40,50)), function(x){
    lola_enrich_func(regiondb_use = regionDB, ifout = TRUE, 
                     bed1_file = sprintf("%s/pm_bottom%s_hg19.bed", bed_file, x),
                     bed2_file = sprintf("%s/pm_top%s_hg19.bed", bed_file, x),
                     out_file = sprintf("%s/pm_topbottom%s_3d", lola_out_file, x))}),
  
  pm_topbottom = lapply(as.list(c(10,20,30,40,50)), function(x){
    lola_enrich_func(regiondb_use = regionDB, ifout = TRUE, 
                     bed1_file = sprintf("%s/pn_bottom%s_hg19.bed", bed_file, x),
                     bed2_file = sprintf("%s/pn_top%s_hg19.bed", bed_file, x),
                     out_file = sprintf("%s/pn_topbottom%s_3d", lola_out_file, x))}),
  
  age4_class = lapply(as.list(1:3), function(x){
    lola_enrich_func(regiondb_use = regionDB, ifout = TRUE, 
                     bed1_file = sprintf("%s/age4_4_hg19.bed", bed_file),
                     bed2_file = sprintf("%s/age4_%s_hg19.bed", bed_file, x),
                     out_file = sprintf("%s/age4_%sto4_3d", lola_out_file, x))})
  
)

####figure5C: plot heatmap for subcompartment--------------------------
lola_subcompartment_results <- list(
  pn = fread(sprintf("%s/pn_topbottom20_3d/col_subcompartment.tsv", lola_out_file)) %>%
    .[, Class := ifelse(userSet == 1, "Bottom 20% GPS-N genes", "Top 20% GPS-N genes")] %>% .[, c(24,3:15,18,21:23)],
  pm = fread(sprintf("%s/pm_topbottom20_3d/col_subcompartment.tsv", lola_out_file)) %>%
    .[, Class := ifelse(userSet == 1, "Bottom 20% GPS-M genes", "Top 20% GPS-M genes")] %>% .[, c(24,3:15,18,21:23)],
  age = fread(sprintf("%s/age4_1to4_3d/col_subcompartment.tsv", lola_out_file)) %>%
    .[, Class := ifelse(userSet == 1,  "Eutheria genes", "Euteleostomi genes")] %>% .[, c(24,3:15,18,21:23)]
)

lola_subcompartment_results_all <- do.call(rbind, lola_subcompartment_results)[, oddsRatio := ifelse(oddsRatio == "Inf", support*d/(b+c), oddsRatio)]


##GPS-N and age
library(ggpattern)
plot_data <- lola_subcompartment_results_all[Class %in% c("Top 20% GPS-N genes", "Bottom 20% GPS-N genes", "Euteleostomi genes", "Eutheria genes"),
                                             .(Class, oddsRatio, description, qValue)][order(oddsRatio)] %>%
  .[, `:=`(`Odds ratio` = c(1:24),
           Class = factor(Class, levels = c("Top 20% GPS-N genes", "Bottom 20% GPS-N genes", "Euteleostomi genes", "Eutheria genes"),
                          labels = c("H-GPS-N", "L-GPS-N", "Euteleostomi", "Eutheria")),
           description = factor(description, levels = c(paste0("Subcompartment_", c("A2", "A1", "B4", "B3", "B2", "B1"))),
                                labels = c("A2", "A1", "B4", "B3", "B2", "B1")),
           facet_x = ifelse(Class %in% c("Top 20% GPS-N genes", "Bottom 20% GPS-N genes"), "GPS-N group", "Age category") %>%
             factor(., levels = c("GPS-N group", "Age category")),
           facet_y = ifelse(description %in% c("Subcompartment_A1",  "Subcompartment_A2"), "Compartment A", "Compartment B"),
           marker = ifelse(qValue < 1e-10, "***", 
                           ifelse(qValue < 1e-5, "**",
                                  ifelse(qValue < 0.05, "*",""))))]

plot_compartment <- ggplot(data = plot_data, 
                           aes(x = Class, y = description, fill = `Odds ratio`)) + theme_bw() +
  geom_tile(color = "grey30") +
  scale_pattern_manual(values = c("Significant" = "stripe", "Non-significant" = "none")) +
  scale_pattern_fill_manual(values = c("Significant" = "black", "Non-significant" = "white")) +
  scale_fill_gradient2(low = "#FFFBF5", mid =  "#FFEBCA", high = "#fed789", midpoint = 13,
                       breaks = c(5,10,15,20), labels = c(5,10,15,20)) +
  geom_text(aes(label = marker), color = "black", size = 6) +
  facet_grid(facet_y~facet_x, space = "free", scales = "free") +
  labs(y = "", title = "", x = "") + 
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(size = 12, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold"), 
        axis.title.x = element_text(size = 15, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.ticks = element_blank(), axis.line = element_blank(),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill=NA, linewidth=0.5),
        legend.position = "top",
        strip.text = element_text(size=0), strip.background = element_blank()) 


plot_compartment
ggsave(plot = plot_compartment, width = 5, height = 5, device = cairo_pdf,
       filename = sprintf("%s/enrich_compartment_20%%.pdf",figure_file))



####figure5F&5G: plot heatmap for loop, atac and domain--------------------------
pm_results_series <- lapply(as.list(c(10,20,30,40,50)), function(x){
  file_path <- sprintf("%s/pm_topbottom%s_3d/allEnrichments.tsv", lola_out_file, x)
  data <- fread(file_path) %>%
    .[, class := x]
  return(data)
}) %>% do.call(rbind, .) 

pn_results_series <- lapply(as.list(c(10,20,30,40,50)), function(x){
  file_path <- sprintf("%s/pn_topbottom%s_3d/allEnrichments.tsv", lola_out_file, x)
  data <- fread(file_path) %>%
    .[, class := x]
  return(data)
}) %>% do.call(rbind, .) 

age4_results <- lapply(as.list(1:3), function(x){
  file_path <- sprintf("%s/age4_%sto4_3d/allEnrichments.tsv", lola_out_file, x)
  data <- fread(file_path) %>%
    .[, class := x]
  return(data)
}) %>% do.call(rbind, .) 

heatmap_func <- function(data = pn_results_series, color_value = c("#DEEBF7", "#08519C"), feature = "loop", 
                         str_keep_num = 2, split_type = "_", x_label = "GPS-N", 
                         pleio_num = c(10,20,30,40,50), 
                         xaxis_break = c(10,20,30,40,50), ifage = F,
                         xaxis_label = c("10%", "20%", "30%", "40%", "50%"),
                         title_name = "Loop anchors"){
  data <- data[, cell_type := fifelse(grepl("domain", description, ignore.case = TRUE), tstrsplit(description, "_", keep = 2)[[1]], 
                                      fifelse(grepl("loop", description, ignore.case = TRUE), tstrsplit(description, "_", keep = 2)[[1]], 
                                              fifelse(grepl("ATAC", description, ignore.case = TRUE), tstrsplit(description, " ", keep = 1)[[1]], "None")))] %>%
    .[, cell_type := ifelse(cell_type %in% c("loop","pooled"), "Pooled", cell_type)] 
  
  
  plot_data <- data[grepl(feature, description, ignore.case = TRUE) & userSet == 2 & 
                      class %in% pleio_num & cell_type %in% c("GM12878", "HeLa", "IMR90", "K562"), 
                    .(oddsRatio, cell_type, qValue, class, description)]  %>%
    .[description !="Rao_IMR90_Domain",]
  
  if(ifage == T){
    Eutheria_data <- plot_data[class == 1, ] %>%
      .[, `:=`(class = 4, oddsRatio = 1, qValue = 0.06)]
    plot_data <- rbind(plot_data, Eutheria_data)
  }
  
  plot_data[, marker := ifelse(qValue < 1e-10, "***", 
                               ifelse(qValue < 1e-5, "**",
                                      ifelse(qValue < 0.05, "*","")))]
  
  plot <- ggplot(data = plot_data, 
                 aes(x = class, y = cell_type, fill = oddsRatio)) + theme_bw() + 
    geom_tile() +
    scale_fill_gradient(low = color_value[1], high = color_value[2], n.breaks=3) +
    labs(y = "", title = title_name, x = x_label) + 
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          plot.title = element_text(size = 12, face = "bold", color = "black"),
          axis.title.y = element_text(size = 15, face = "bold"), 
          axis.title.x = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line = element_blank(),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), 
          panel.border = element_rect(colour = "grey50", fill=NA, linewidth=0.5),
          legend.position = "top",
          legend.title = element_blank()) +
    scale_x_continuous(breaks = xaxis_break, labels = xaxis_label) +
    geom_text(aes(label = marker), color = "black", size = 6) 
  
  return(list(data = plot_data, plot = plot))
}

loop_domian_atac_heatmap <- list(
  pn_loop = heatmap_func(data = pn_results_series, color_value = c("#DEEBF7", lighten("#08519C", 0.45)), feature = "loop", 
                         str_keep_num = 2, split_type = "_", 
                         pleio_num = c(10,20,30,40,50), 
                         xaxis_break = c(10,20,30,40,50), ifage = F,
                         xaxis_label = c("10%", "20%", "30%", "40%", "50%"),
                         title_name = "Loop anchor"),
  pn_domian = heatmap_func(data = pn_results_series, color_value = c( "#C7EAE5",lighten("#01665E", 0.45)), feature = "domain", 
                           str_keep_num = 2, split_type = "_", title_name = "TAD",
                           pleio_num = c(10,20,30,40,50), 
                           xaxis_break = c(10,20,30,40,50), ifage = F,
                           xaxis_label = c("10%", "20%", "30%", "40%", "50%")),
  pn_atac = heatmap_func(data = pn_results_series, color_value = c("#B2ABD2",lighten("#542788", 0.45)), feature = "ATAC", 
                         str_keep_num = 1, split_type = " ", title_name = "ATAC peak",
                         pleio_num = c(10,20,30,40,50), 
                         xaxis_break = c(10,20,30,40,50), ifage = F,
                         xaxis_label = c("10%", "20%", "30%", "40%", "50%")),
  
  pm_loop = heatmap_func(data = pm_results_series, color_value = c("#DEEBF7", lighten("#08519C", 0.45)), feature = "loop", 
                         str_keep_num = 2, split_type = "_",
                         pleio_num = c(10,20,30,40,50), 
                         xaxis_break = c(10,20,30,40,50), ifage = F,
                         xaxis_label = c("10%", "20%", "30%", "40%", "50%"),
                         title_name = "Loop anchor"),
  pm_domian = heatmap_func(data = pm_results_series, color_value = c( "#C7EAE5",lighten("#01665E", 0.45)), feature = "domain", 
                           str_keep_num = 2, split_type = "_", title_name = "TAD",
                           pleio_num = c(10,20,30,40,50), 
                           xaxis_break = c(10,20,30,40,50), ifage = F,
                           xaxis_label = c("10%", "20%", "30%", "40%", "50%")),
  pm_atac = heatmap_func(data = pm_results_series, color_value = c("#B2ABD2",lighten("#542788", 0.45)), feature = "ATAC", 
                         str_keep_num = 1, split_type = " ", title_name = "ATAC peak",
                         pleio_num = c(10,20,30,40,50), 
                         xaxis_break = c(10,20,30,40,50), ifage = F,
                         xaxis_label = c("10%","20%", "30%", "40%", "50%")),
  
  age4_loop = heatmap_func(data = age4_results, color_value = c("#FFD2C5", lighten("#DA7235", 0.3)), feature = "loop", 
                           str_keep_num = 2, split_type = "_", xaxis_break = c(1,2,3,4), pleio_num = c(1,2,3), ifage = T,
                           xaxis_label = c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria"), 
                           x_label = "Gene age", title_name = "Loop anchor"),
  age4_domian = heatmap_func(data = age4_results, color_value = c("#FFD8EF", lighten("#AB6191", 0.3)), feature = "domain", 
                             str_keep_num = 2, split_type = "_", title_name = "TAD", x_label = "Gene age", 
                             xaxis_break = c(1,2,3,4), pleio_num = c(1,2,3), ifage = T,
                             xaxis_label = c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria")),
  age4_atac = heatmap_func(data = age4_results, color_value = c("#FFF8D7", lighten("#F2DD33", 0.3)), feature = "ATAC", 
                           str_keep_num = 1, split_type = " ", title_name = "ATAC peak",
                           xaxis_break = c(1,2,3,4), pleio_num = c(1,2,3), ifage = T, x_label = "Gene age",
                           xaxis_label = c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria"))
)

loop_domian_atac_heatmap$pn_loop$plot
loop_domian_atac_heatmap$pn_domian$plot
loop_domian_atac_heatmap$age4_loop$plot


ggsave(plot = loop_domian_atac_heatmap$pn_loop$plot, width = 4, height = 4, device = cairo_pdf,
       filename = sprintf("%s/enrich_loop_pn_heatmap.pdf",figure_file))
ggsave(plot = loop_domian_atac_heatmap$pn_domian$plot, width = 4, height = 4, device = cairo_pdf,
       filename = sprintf("%s/enrich_domain_pn_heatmap.pdf",figure_file))
ggsave(plot = loop_domian_atac_heatmap$pn_atac$plot, width = 4, height = 4, device = cairo_pdf,
       filename = sprintf("%s/enrich_atac_pn_heatmap.pdf",figure_file)) 

ggsave(plot = loop_domian_atac_heatmap$age4_loop$plot, width = 4, height = 4, device = cairo_pdf,
       filename = sprintf("%s/enrich_loop_age_heatmap.pdf",figure_file))
ggsave(plot = loop_domian_atac_heatmap$age4_domian$plot, width = 4, height = 4, device = cairo_pdf,
       filename = sprintf("%s/enrich_domain_age_heatmap.pdf",figure_file))
ggsave(plot = loop_domian_atac_heatmap$age4_atac$plot, width = 4, height = 4, device = cairo_pdf,
       filename = sprintf("%s/enrich_atac_age_heatmap.pdf",figure_file)) 

ggsave(plot = loop_domian_atac_heatmap$pm_loop$plot, width = 4, height = 4, device = cairo_pdf,
       filename = sprintf("%s/enrich_loop_pm_heatmap.pdf",figure_file))
ggsave(plot = loop_domian_atac_heatmap$pm_domian$plot, width = 4, height = 4, device = cairo_pdf,
       filename = sprintf("%s/enrich_domain_pm_heatmap.pdf",figure_file))
ggsave(plot = loop_domian_atac_heatmap$pm_atac$plot, width = 4, height = 4, device = cairo_pdf,
       filename = sprintf("%s/enrich_atac_pm_heatmap.pdf",figure_file)) 
