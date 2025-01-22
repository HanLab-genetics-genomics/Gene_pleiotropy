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
#3 figure6--------------------------
##figure6A-------------------------
library(genekitr); library(igraph); library(ggraph); library(patchwork); library(rrvgo); library(clusterProfiler)
ego_bp <- list(pm10 = genORA(id = pleiotropy_maindata$pm_ld[pleio10 == 10, gene], 
                             geneset = geneset::getGO(org = "human",ont = "bp")),
               pm1 = genORA(id = pleiotropy_maindata$pm_ld[pleio10 == 1, gene], 
                            geneset = geneset::getGO(org = "human",ont = "bp")),
               pn10 = genORA(id = pleiotropy_maindata$pn_ld[pleio10 == 10, gene], 
                             geneset = geneset::getGO(org = "human",ont = "bp")),
               pn1 = genORA(id = pleiotropy_maindata$pn_ld[pleio10 == 1, gene], 
                            geneset = geneset::getGO(org = "human",ont = "bp")),
               old = genORA(id = pleiotropy_maindata$pn_ld[age_stage4_num == 1 & gene_age_num >= 1500 , gene],  #1106
                            geneset = geneset::getGO(org = "human",ont = "bp")),
               young = genORA(id = pleiotropy_maindata$pn_ld[age_stage4_num == 4, gene], 
                              geneset = geneset::getGO(org = "human",ont = "bp")))

gobp_plot <- lapply(list("pn10", "pn1", "pm10", "pm1", "old", "young"), function(x){
  go_analysis <- ego_bp[[x]]
  simMatrix <- calculateSimMatrix(go_analysis$Hs_BP_ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  scores <- setNames(-log10(go_analysis$qvalue), go_analysis$Hs_BP_ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  
  cut_term <- as.data.table(reducedTerms) %>% .[, .N, by = parentTerm] %>% .[N == 1, parentTerm]
  plot_reducedata <- reducedTerms[!reducedTerms$parentTerm %in% cut_term,]
  plot_reducedata$size <- plot_reducedata$size*0.1
  
  return(plot_reducedata)
})
names(gobp_plot) <- c("pn10", "pn1", "pm10", "pm1", "old", "young")

treemapPlot(gobp_plot$pn10, overlap.labels = 0.99, force.print.labels = T, 
            lowerbound.cex.labels = 0.01, border.lwds = c(0.8,0.4))

lapply(list("pn10"), function(x){
  pdf(file = sprintf("%s/bp_%s_enrich.pdf", figure_file, x), width=8, height=8)
  treemapPlot(gobp_plot[[x]], overlap.labels = 0.99, force.print.labels = T, 
              lowerbound.cex.labels = 0.01, border.lwds = c(0.8,0.4))
  dev.off()
})

##figure6B&6C-------------------------
###load data-----------------
chronos_data <- fread("CRISPRGeneEffect.csv")
chronos_essential <- data.table(gene = names(chronos_data)[-1],
                                essential_score = colMeans(chronos_data[,-1], na.rm = TRUE)) %>%
  separate(gene, into = c("gene", "id"), sep = " ") %>% .[,-2] %>%
  mutate(., essential_num = apply(chronos_data[,-1], 2, function(x){sum(x < -1)})) %>%
  mutate(., essential_chronos = ifelse(essential_num >= 1, 1, 2))

GES_data <- lapply(list("pn_ld", "pm_ld"), function(x, data = pleiotropy_maindata){
  data_merge <- merge(data[[x]][, .(gene, pleio_class3, age_stage4)], chronos_essential, by = "gene", all.x = T) %>%
    .[, GES := -essential_score]
  return(data_merge)
})

names(GES_data) <- c("pn_ld", "pm_ld")

###GESs across groups--------------------------
GES_box_func <- function(data = hopsgene_ageburden_cut_genemetrics$pn_ld, x_label = "", 
                         y_label = "Gene length (bp)" , group = "pleio_class3", x_var = "Gene length (bp)",
                         x_text = c("Low", "Intermediate", "High"), color_value = pleio_color,
                         stat_test = test, test_list = c("Low pleiotropy", "High pleiotropy"), stat_name){
  
  plot_data <- data[, .(median = sapply(.SD, function(x) as.numeric(mean(x, na.rm = TRUE)))), by = group, .SDcols = x_var, with = T] 
  
  plot <- ggplot() + theme_bw() +
    stat_boxplot(data = data, aes(x = !!sym(group), y = !!sym(x_var)), geom ='errorbar', width = 0.2) +
    #geom_jitter(data = data, aes(x = !!sym(group), y = !!sym(x_var), fill = !!sym(group),
    #                color = !!sym(group)), width = 0.4, size = 0.2, alpha = 0.2) +
    geom_boxplot(data = data, aes(x = !!sym(group), y = !!sym(x_var), fill = !!sym(group)),
                 color="grey50", lwd = 0.6, outlier.alpha = 0.3) +
    geom_signif(data = data, aes(x = !!sym(group), y = !!sym(x_var)),
                comparisons = list(test_list), color = "grey20",
                map_signif_level = TRUE,
                annotations = stat_test[[stat_name]]) +
    geom_point(data = plot_data, aes(x = !!sym(group), y = median), color = "grey20", size = 2.5) +
    geom_line(data = plot_data, aes(x = !!sym(group), y = median, group = 1), 
              linetype = "dashed", color = "grey20") +
    scale_color_manual(values = color_value) +
    scale_fill_manual(values = color_value) +
    labs(x=x_label, y= y_label) +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          legend.position = "none",
          axis.title.y = element_text(size = 15, face = "bold"), 
          axis.title.x = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          axis.line.y.right = element_blank(),
          axis.text.y.right = element_blank()) +
    scale_x_discrete(labels = x_text) +
    scale_y_break(c(10e-4, 10e-3), scales = 5) +
    scale_y_log10(breaks = c(10e-4, 10e-3,10e-2, 10e-1, 1, 10),
                  labels = trans_format("log10", math_format(10^.x))) 
  
  return(list(plot_data = plot_data, plot = plot))
}

test <- list(
  ges_pleio_n = pairwise.wilcox.test(GES_data$pn_ld$GES, 
                                     GES_data$pn_ld$pleio_class3, 
                                     p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .), 
  ges_pleio_m = pairwise.wilcox.test(GES_data$pm_ld$GES, 
                                     GES_data$pm_ld$pleio_class3, 
                                     p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .), 
  ges_age = pairwise.wilcox.test(GES_data$pn_ld$GES, 
                                 GES_data$pn_ld$age_stage4, 
                                 p.adjust.method = "fdr")$p.value[3,1] %>% sprintf("%.2e", .) 
)

GES_box <- list(pleio_n = GES_box_func(data = GES_data$pn_ld, color_value = pleio_color,
                                       y_label = "Gene essential scores", x_var = "GES",
                                       stat_name = "ges_pleio_n"),
                pleio_m = GES_box_func(data = GES_data$pm_ld, color_value = pleio_color,
                                       y_label = "Gene essential scores", x_var = "GES",
                                       stat_name = "ges_pleio_m"),
                age_n = GES_box_func(data = GES_data$pn_ld[!is.na(age_stage4),], group = "age_stage4",
                                     color_value = age_color, y_label = "Gene essential scores", x_var = "GES", 
                                     x_text = x_text, test_list = age_test, stat_name = "ges_age")) 


GES_box$pleio_n$plot
GES_box$age_n$plot

ggsave(plot = GES_box$pleio_n$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/ges_pn.pdf",figure_file))
ggsave(plot = GES_box$age_n$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/ges_age.pdf",figure_file))



##figure6D-------------------------
library(dorothea); library(OmnipathR); library(decoupleR)
progeny <- progeny::model_human_full %>%
  group_by(pathway) %>%
  slice_min(order_by = p.value, n = 500)

names(progeny) <- c("target", "source", "weight", "p_value")
###mirrorbar----------------------------
mirrorbar_func <- function(data = pleiotropy_maindata$pm_ld, pleio_group = "pleio_class3", 
                           high_low_pair = list("High pleiotropy", "Low pleiotropy"),
                           high_low_name = c("H-GPS-N", "L-GPS-N"),
                           high_low_color = c(lighten("#65AC9C",0.45), lighten("#3C3E7B",0.7))){
  
  progeny_num <- progeny %>% as.data.table() %>%
    .[p_value < 0.05, .(total = .N), by = source]
  
  data$pleio_use <- data[,get(pleio_group)]
  
  percent_high <- progeny %>% as.data.table() %>%
    .[target %in% data[pleio_use %in% high_low_pair[[1]], gene] & p_value < 0.05, .N, by = source] %>%
    merge(., progeny_num, by = "source") %>%
    .[, percent := N/total*100] %>%
    .[order(-percent)] %>%
    .[, group := "pleio_high"]
  
  percent_low <- progeny %>% as.data.table() %>%
    .[target %in% data[pleio_use %in% high_low_pair[[2]], gene] & p_value < 0.05, .N, by = source] %>%
    merge(., progeny_num, by = "source") %>%
    .[, percent := N/total*100] %>%
    .[order(-percent)] %>%
    .[, group := "pleio_low"]
  
  df_percent <- rbind(percent_high, percent_low) %>%
    .[, source := factor(source, levels = percent_high$source)] %>%
    .[order(source)] %>%
    .[, percent := ifelse(group == "pleio_high", percent, -percent)] 
  
  diff_test = lapply(as.list(unique(df_percent$source)), function(x){
    prop.test(c(df_percent[source == x]$N[1], df_percent[source == x]$N[2]),
              c(df_percent[source == x]$total[1], df_percent[source == x]$total[2]))$p.value
  }) %>% unlist(.) 
  ajust_pvalue = p.adjust(diff_test, method = "fdr")
  
  df_percent[, `:=`(p_value = rep(diff_test, each = 2),
                    fdr = rep(ajust_pvalue, each = 2))] %>%
    .[, text := ifelse(group == "pleio_high" & fdr < 0.001, "***", 
                       ifelse(group == "pleio_high" & fdr < 0.01, "**",
                              ifelse(group == "pleio_high" & fdr < 0.05, "*", "")))] %>%
    .[, color := ifelse(source == "Trail", "2", "1")]
  
  
  plot <- ggplot(df_percent, aes(x = source, y = percent, fill = group)) + theme_bw() +
    geom_bar(stat = "identity", position = "identity", color = "grey30", size = 0.2) +
    scale_fill_manual(values = c("pleio_high" = high_low_color[1], "pleio_low" = high_low_color[2]),
                      labels = c("pleio_high" = high_low_name[1], "pleio_low" = high_low_name[2])) +  
    coord_flip() +
    labs(y = "Percent overlap with PROGENy pathway target genes (%)", title = " ", x = "") + 
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          plot.title = element_text(family = "Lobster Two", size = 12, face = "bold", color = "black"),
          axis.title.x = element_text(size = 14, face = "bold"), 
          axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 10, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.title = element_blank(),
          legend.position = "top") +
    geom_text(aes(x = source, y = 30, label = text, color = color),
              fontface = "bold", size = 5, show.legend = FALSE) +
    scale_color_manual(values = c("1" = darken(high_low_color[1], 0.5), "2" = darken(high_low_color[2], 0.5)),
                       labels = c("1" = "", "2" = ""))  
  
  return(list(data = df_percent, plot = plot))
}

progeny_plot_list <- list(
  pn = mirrorbar_func(data = pleiotropy_maindata$pn_ld,
                      pleio_group = "pleio_class3", 
                      high_low_pair = list("High pleiotropy", "Low pleiotropy"),
                      high_low_name = c("H-GPS-N", "L-GPS-N"),
                      high_low_color = c(lighten("#65AC9C",0.45), lighten("#3C3E7B",0.7))),
  pm = mirrorbar_func(data = pleiotropy_maindata$pm_ld,
                      pleio_group = "pleio_class3", 
                      high_low_pair = list("High pleiotropy", "Low pleiotropy"),
                      high_low_name = c("H-GPS-M", "L-GPS-M"),
                      high_low_color = c(lighten("#65AC9C",0.45), lighten("#3C3E7B",0.7)))
)

progeny_plot_list$pn$plot

ggsave(plot = progeny_plot_list$pn$plot, width = 8, height = 7, device = cairo_pdf,
       filename = sprintf("%s/progeny_pn.pdf",figure_file))

##figure6E--------------------------
kegg_func <- function(data = pleiotropy_maindata$pn_ld, pleio_group = 10){
  
  kegg_pleio10 <- genORA(id = data[pleio10 == pleio_group, gene],
                         geneset = geneset::getKEGG(org = "hsa",category = "pathway"))
  kegg_pleio10$pathway <-sub("hsa", "map", kegg_pleio10$ID)
  
  map_ko <- lapply(as.list(intersect(kegg_pleio10$pathway, unlist(luca_map))), function(x){
    ko <- luca[Probable == 1, ][grep(x, luca_map),]$KEGG_ko %>% as.character(.)
    return(ko)}) %>% unlist(.) %>% unique(.)                       
  
  map_pathway <- kegg_pleio10[kegg_pleio10$pathway %in% intersect(kegg_pleio10$pathway, unlist(luca_map)), "Description"]
  
  return(list(kegg = kegg_pleio10, map_pathway = map_pathway, map_ko =  map_ko))
}
luca <- fread("41559_2024_2461_MOESM4_ESM.tsv")
luca_map <- as.list(luca[Probable == 1,]$Map) %>%
  lapply(., function(x){gsub("path:", "", unlist(strsplit(x, ",")))}) 

kegg_list <- list(
  pn10 = kegg_func(data = pleiotropy_maindata$pn_ld, pleio_group = 10),
  pn1 = kegg_func(data = pleiotropy_maindata$pn_ld, pleio_group = 1),
  pm10 = kegg_func(data = pleiotropy_maindata$pm_ld, pleio_group = 10),
  pm1 = kegg_func(data = pleiotropy_maindata$pm_ld, pleio_group = 1)
)


library(enrichplot)
keggluca_plot_func <- function(data = kegg_list$pn10, plot_height = 10, plot_name = "pn10", plot_type = "dot") {
  
  kegg_data <- data$kegg %>% as.data.table(.) %>% .[order(FoldEnrich)] 
  ylabel_color <- ifelse(kegg_data$Description %in% data$map_pathway, "#D18859", "black")
  
  plot <- plotEnrich(kegg_data, plot_type = plot_type, scale_ratio = 0.3, up_color = "#FFB287", down_color = "#7391E0",
                     stats_metric = "qvalue", term_metric = "FoldEnrich", wrap_length = 25) +
    theme(axis.text.y = element_text(color = ylabel_color, face = "bold", size = 12),
          axis.text.x = element_text(color = "black", face = "bold", size = 12),
          text = element_text(family = "Arial",size = 12, color = "black", face = "bold"),
          #legend.position = "bottom",
          legend.text = element_text(family = "Arial", size = 12, face = "bold"),
          legend.title = element_text(family = "Arial", size = 12, face = "bold"),
          axis.ticks = element_line(colour = "grey50"), 
          axis.line = element_line(colour = "grey50"),
          panel.border = element_rect(colour = "grey50"),
          axis.title.y = element_blank(), 
          axis.title.x = element_text(family = "Arial", size = 15, face = "bold")) +
    xlab("Enrichment")
  
  ggsave(plot = plot, width = 10, height = plot_height, device = cairo_pdf,
         filename = sprintf("kegg_luca_%s.pdf", plot_name))
}

keggluca_plot_list <- list(
  pn10 = keggluca_plot_func(data = kegg_list$pn10, plot_name = "pn10", plot_height = 10, plot_type = "dot"),
  pm10 = keggluca_plot_func(data = kegg_list$pm10, plot_name = "pm10", plot_height = 10, plot_type = "dot"),
  pn1 = keggluca_plot_func(data = kegg_list$pn1, plot_name = "pn1", plot_height = 2, plot_type = "bar"),
  pm1 = keggluca_plot_func(data = kegg_list$pm1, plot_name = "pm1", plot_height = 3, plot_type = "bar")
)
