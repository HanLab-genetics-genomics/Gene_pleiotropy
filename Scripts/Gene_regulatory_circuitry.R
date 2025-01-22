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
#3 figure4-----------------------
##figure4A----------------------
###load data and functions--------------------
qtl_color <- lighten(brewer.pal(n = 10, name = 'PuOr'), 0.8)
qtl_color2 <- lighten(brewer.pal(n = 10, name = 'PiYG'), 0.8)
percent_plot <- function(data = egenes_prop[1:3,49], 
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
###cis-egenes and cis-sgenes across groups---------------------
#wget -c https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar
#wget -c https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_sQTL.tar
####load data (download the data or load the revised data provided)-----------------
prop_data <- lapply(list("/home/liumy/pleiotropy/data/qtl/GTEx_Analysis_v8_eQTL",
                         "/home/liumy/pleiotropy/data/qtl/GTEx_Analysis_v8_sQTL"), 
                    function(file_path, score_data = pleiotropy_maindata$pn_ld){
                      
                      file_egenes <- list.files(path = file_path, pattern = "genes\\.txt$", full.names = TRUE)
                      egenes_gtex_list <- lapply(unlist(file_egenes), function(x, data = score_data){
                        egenes <- fread(x) %>% 
                          .[qval < 0.05, ] %>% 
                          .[, .(ensemblid = sub("\\.\\d+$", "", gene_id))] %>% 
                          setnames(., names(.), sub(".*\\/([^.]+)\\.v8.*", "\\1", x))
                        
                        return(egenes)
                      })
                      
                      egenes_gtex <- do.call(cbind, egenes_gtex_list) 
                      
                      egene_prop_list <- lapply(unlist(colnames(egenes_gtex)), function(x, data = pleiotropy_maindata$pn_ld){
                        
                        data_merge <- data[, .(ensemblid, gene, pleio_class3, pleio_class5, age_stage4)][, ifegenes := ifelse(ensemblid %in% egenes_gtex[, get(x)], 1, 2)]
                        sum_table_pleio3 <- table(data_merge[, .(ifegenes, pleio_class3)]) %>% prop.table(., margin = 2) %>% .[1,] %>% as.data.frame(.) %>%
                          setnames(., names(.), x)
                        sum_table_age <- table(data_merge[, .(ifegenes, age_stage4)]) %>% prop.table(., margin = 2) %>% .[1,] %>% as.data.frame(.) %>%
                          setnames(., names(.), x)
                        
                        sum_table <- rbind(sum_table_pleio3, sum_table_age)
                        return(sum_table)
                      })
                      
                      egenes_prop <- do.call(cbind, egene_prop_list) 
                      
                      return(egenes_prop)
                    })
names(prop_data) <- c("cis-eqtls", "cis-sqtls")
####plot data--------------------
load("cis_esqtl.RData")
test_data <- data.table(
  class = c("Low pleiotropy", "Intermediate pleiotropy", "High pleiotropy",
            "Euteleostomi", "Tetrapoda", "Amniota", "Eutheria"),
  total_num = c(3593, 10777, 3593, 10072, 2427, 2578, 1570),
  isegenes_prop = prop_data$`cis-eqtls`[,49],
  issgenes_prop = prop_data$`cis-sqtls`[,49]) %>%
  .[, .(isegenes = total_num*isegenes_prop, 
        isnotegenes = total_num*(1-isegenes_prop),
        issgenes = total_num*issgenes_prop, 
        isnotsgenes = total_num*(1-issgenes_prop))] %>% as.matrix()

test <- list(
  egene_pleio = pairwise.prop.test(test_data[1:3, "isegenes"], c(3593, 10777, 3593),
                                   p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .),
  egene_age = pairwise.prop.test(test_data[4:7, "isegenes"], c(10072, 2427, 2578, 1570),
                                 p.adjust.method = "fdr")$p.value[3,1] %>% sprintf("%.2e", .),
  sgene_pleio = pairwise.prop.test(test_data[1:3, "issgenes"], c(3593, 10777, 3593),
                                   p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .),
  sgene_age = pairwise.prop.test(test_data[4:7, "issgenes"], c(10072, 2427, 2578, 1570),
                                 p.adjust.method = "fdr")$p.value[3,1] %>% sprintf("%.2e", .)
)


percent_egenes <- list(pleio = percent_plot(data = prop_data$`cis-eqtls`[1:3, 49], 
                                            group = factor(c("Low", "Intermediate", "High"), 
                                                           levels = c("Low", "Intermediate", "High")),
                                            color_values = qtl_color[c(8,9,10)], y_label = "",
                                            x_label = "Cis-eGenes", stat_name = "egene_pleio"),
                       age = percent_plot(data = prop_data$`cis-eqtls`[4:7, 49], 
                                          group = factor(c("Euteleostomi", " Tetrapoda", " Amniota", "Eutheria"), 
                                                         levels = c("Euteleostomi", " Tetrapoda", " Amniota", "Eutheria")),
                                          color_values = qtl_color[1:4], y_label = "",
                                          x_label = "Cis-eGenes", stat_name = "egene_age", test_list = age_test))
percent_sgenes <- list(pleio = percent_plot(data = prop_data$`cis-sqtls`[1:3, 49], 
                                            group = factor(c("Low", "Intermediate", "High"), levels = c("Low", "Intermediate", "High")),
                                            color_values = qtl_color[c(8,9,10)], y_label = "",
                                            x_label = "Cis-sGenes", stat_name = "sgene_pleio"),
                       age = percent_plot(data = prop_data$`cis-sqtls`[4:7, 49], 
                                          group = factor(c("Euteleostomi", " Tetrapoda", " Amniota", "Eutheria"), 
                                                         levels = c("Euteleostomi", " Tetrapoda", " Amniota", "Eutheria")),
                                          color_values = qtl_color[1:4], y_label = "",
                                          x_label = "Cis-sGenes", stat_name = "sgene_age", test_list = age_test))

percent_egenes$pleio
percent_egenes$age
percent_sgenes$pleio
percent_sgenes$age

ggsave(plot = percent_egenes$pleio, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/cisegenes_pn.pdf",figure_file))
ggsave(plot = percent_egenes$age, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/cisegenes_age.pdf",figure_file))

ggsave(plot = percent_sgenes$pleio, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/cissgenes_pn.pdf",figure_file))
ggsave(plot = percent_sgenes$age, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/cissgenes_age.pdf",figure_file))

###trans-eGenes----------------------------------
trans_eqtl <- fread("eqtlgen_transeqtl.txt")
transeqtl_prop <- lapply(list("pm_ld", "pn_ld"), function(x, score_data = pleiotropy_maindata){
  
  gene_list <- unique(trans_eqtl$Gene)
  data_merge <- score_data[[x]][, .(ensemblid, gene, pleio_class3, age_stage4)][
    , ifegenes := ifelse(ensemblid %in% gene_list, 1, 2)] 
  sum_table_pleio3 <- table(data_merge[, .(ifegenes, pleio_class3)]) %>% 
    prop.table(., margin = 2) %>% .[1,] %>% as.data.frame(.) %>%
    setnames(., names(.), x)
  sum_table_age <- table(data_merge[, .(ifegenes, age_stage4)]) %>% 
    prop.table(., margin = 2) %>% .[1,] %>% as.data.frame(.) %>%
    setnames(., names(.), x)
  transegenes_prop <- rbind(sum_table_pleio3, sum_table_age)
  
  pleio_count <- table(data_merge[, .(pleio_class3, ifegenes)])
  age_count <- table(data_merge[, .(age_stage4, ifegenes)])
  
  test <- list(
    transegene_pleio = pairwise.prop.test(pleio_count[,"1"], 
                                          c(sum(pleio_count["Low pleiotropy", ]), sum(pleio_count["Intermediate pleiotropy",]), sum(pleio_count["High pleiotropy",])),
                                          p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .),
    transegene_age = pairwise.prop.test(age_count[,"1"], 
                                        c(sum(age_count["Euteleostomi", ]), sum(age_count["Tetrapoda",]), 
                                          sum(age_count["Amniota",]), sum(age_count["Eutheria",])),
                                        p.adjust.method = "fdr")$p.value[3,1] %>% sprintf("%.2e", .)
  )
  
  return(list(transegenes_prop = transegenes_prop, test = test))
})

names(transeqtl_prop) <- c("pm_ld", "pn_ld")

percent_transegenes <- list(pleio_n = percent_plot(data = transeqtl_prop$pn_ld$transegenes_prop[1:3,], 
                                                   group = factor(c("Low", "Intermediate", "High"), levels = c("Low", "Intermediate", "High")),
                                                   color_values = qtl_color[c(8,9,10)], y_label = "",
                                                   x_label = "Trans-eGenes", 
                                                   stat_test = transeqtl_prop$pn_ld$test, stat_name = "transegene_pleio"),
                            pleio_m = percent_plot(data = transeqtl_prop$pm_ld$transegenes_prop[1:3,], 
                                                   group = factor(c("Low", "Intermediate", "High"), levels = c("Low", "Intermediate", "High")),
                                                   color_values = qtl_color[c(8,9,10)], y_label = "",
                                                   x_label = "Trans-eGenes", 
                                                   stat_test = transeqtl_prop$pn_ld$test, stat_name = "transegene_pleio"),
                            age = percent_plot(data = transeqtl_prop$pn_ld$transegenes_prop[4:7, ], 
                                               group = factor(c("Euteleostomi", " Tetrapoda", " Amniota", "Eutheria"), 
                                                              levels = c("Euteleostomi", " Tetrapoda", " Amniota", "Eutheria")),
                                               color_values = qtl_color[1:4], y_label = "",
                                               x_label = "Trans-eGenes", test_list = age_test,
                                               stat_test = transeqtl_prop$pn_ld$test, stat_name = "transegene_age"))

percent_transegenes$pleio_n
percent_transegenes$age

ggsave(plot = percent_transegenes$pleio_n, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/transegenes_pn.pdf",figure_file))
ggsave(plot = percent_transegenes$age, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/transegene_age.pdf",figure_file))


##figure4B&4C----------------------
###load data
tau_score <- fread("/home/liumy/pleiotropy/data/specificity/Tau_gene_V8.csv")
tau_merge <- lapply(list("pm_ld", "pn_ld"), function(x, data = pleiotropy_maindata){
  data_merge <- merge(data[[x]][, .(gene, ensemblid, pleio_class3, age_stage4)], tau_score, by.x = "ensemblid", by.y = "gene_id", all.x = T)
})
names(tau_merge) <- c("pm_ld", "pn_ld")
tau_merge$pn_ld$ifspecific <- ifelse(tau_merge$pn_ld$tau > 0.8, "1", "2")
tau_merge$pm_ld$ifspecific <- ifelse(tau_merge$pm_ld$tau > 0.8, "1", "2")

##ridge plot
library(ggridges)
tau_ridge_func <- function(data = tau_merge$pn_ld, group = "pleio_class3", color_values = pleio_color){
  data$pleio_class3 <- factor(data$pleio_class3, levels = c("High pleiotropy", "Intermediate pleiotropy", "Low pleiotropy"))
  plot <- ggplot(data = data, aes(x = tau, y = !!sym(group), fill = !!sym(group))) +
    geom_density_ridges() +
    theme_ridges() + 
    scale_fill_manual(values = color_values) +
    labs(y = "", title = "", x = "Tissue-specific τ index") +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_text(size = 15, face = "bold", ), 
          axis.text.x = element_text(size = 12, color = "black", face = "bold"), 
          axis.text.y = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line.x = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", 
          panel.border = element_blank(),
          plot.title = element_text(size=15, hjust=0.5),
          panel.grid.major = element_blank()) +
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), labels = c(0,0.2,0.4,0.6,0.8,1))
  return(plot)
}

tau_pleio_n <- tau_ridge_func(data = tau_merge$pn_ld[!is.na(tau),], group = "pleio_class3", color_values = rev(pleio_color))
tau_pleio_m <- tau_ridge_func(data = tau_merge$pm_ld[!is.na(tau),], group = "pleio_class3", color_values = rev(pleio_color))
tau_age <- tau_ridge_func(data = tau_merge$pn_ld[!is.na(age_stage4) & !is.na(tau),], group = "age_stage4", color_values = age_color) 

tau_pleio_n
tau_age

ggsave(plot = tau_pleio_n, width = 5, height = 5, device = cairo_pdf,
       filename = sprintf("%s/tau_ridge_pn.pdf",figure_file))
ggsave(plot = tau_age, width = 5, height = 5, device = cairo_pdf,
       filename = sprintf("%s/tau_ridge_age.pdf",figure_file))

##figure4D----------------------
express_hist <- function(id = "FTO", data = tau_merge$pn_ld, plot_title = "FTO"){
  
  plot_data <- data.table(tissue = names(data)[-c(1:5, 36)],
                          express = as.numeric(data[gene == id, -c(1:5, 36)])) %>%
    .[order(-express)] %>%
    .[, tissue_num := 1:.N]
  plot_data$tissue <- factor(plot_data$tissue_num, labels = plot_data$tissue)
  
  plot_data <- plot_data[express != 0,]
  plot <- ggplot() + theme_bw() +
    geom_bar(data = plot_data, aes(x = tissue, y = express), stat = "identity", width = 0.7, 
             fill = "#FFECD4", color = "#FFECD4") +
    geom_line(aes(x = seq(1, nrow(plot_data), length.out = nrow(plot_data)), 
                  y = spline(plot_data$tissue_num, plot_data$express, 
                             xout = seq(1, nrow(plot_data), length.out = nrow(plot_data)))$y), 
              color = "#FBD9AA", size = 1) +
    theme(axis.title = element_text(size = 15, face = "bold"),
          plot.title = element_text(hjust = 0.5,face = "bold"),
          axis.text.y = element_text(size = 12,face="bold"),
          axis.text.x = element_text(angle = 45,  vjust = 0.5, size = 12, face="bold"),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.line = element_line(colour = "grey50"),
          legend.position = "none", legend.text = element_text(size = 9,face = "bold"),axis.ticks.length=unit(-0.1, "cm"),
          legend.title = element_blank(),legend.background = element_blank(),legend.key.size = unit(0.5,"cm"),
          plot.margin = unit(c(0.2,0.5,0.2,0.2), "cm")) +
    labs(y = "Gene expression",x = "") 
  #scale_x_discrete(guide = guide_axis(n.dodge=3)) +
  #ggtitle(plot_title)
  
  return(plot)
} 

### high pleiotropy & low specificity
plot_hist <- list(
  PRKAR1A= express_hist(id = "PRKAR1A", data = tau_merge$pn_ld, plot_title = "PRKAR1A (11.80)"),
  CST9L = express_hist(id = "CST9L", data = tau_merge$pn_ld, plot_title = "CST9L (4.93)"))

plot_hist$PRKAR1A
plot_hist$CST9L

ggsave(plot = plot_hist[[1]], width = 10, height = 5, device = cairo_pdf,
       filename = sprintf("%s/tissue_express_fto.pdf",figure_file))
ggsave(plot = plot_hist[[2]], width = 1, height = 5, device = cairo_pdf,
       filename = sprintf("%s/tissue_express_slc12a3.pdf",figure_file))

##figure4E&4F----------------------
###load data-----------------
library(dorothea); library(OmnipathR); library(decoupleR)
tf <- fread("TFC2_16102023b.tsv") %>%
  .[substr(Entrez.Taxa.ID,1,4) == 9606, Associated.Gene.Name] #TFCheckpoint
dorothea <-  get_dorothea(organism = "human",  levels = c("A", "B", "C", "D"),   
                          weight_dict = list(A = 1, B = 2, C = 3, D = 4, E = 5))
dorothea_table <- dorothea %>% as.data.table(.) %>% .[, .N, by = .(confidence, target)]
dorothea_table_wide <- dcast(dorothea_table, target ~ confidence, value.var = "N") %>%
  .[, tf_sum := apply(.SD, MARGIN = 1, sum, na.rm = T), .SDcols = 2:5] %>%
  setnames(., names(.), c("target", "tf_A", "tf_B", "tf_C", "tf_D", "tf_sum"))

hopsgene_ageburden_cut_tf <- lapply(list("pm_ld", "pn_ld"), function(score_name){
  
  data <- pleiotropy_maindata[[score_name]]
  data <- data[, if_tf := ifelse(gene %in% tf, 1, 2)]
  data <- merge(data, dorothea_table_wide, by.x = "gene", by.y = "target", all.x = T)
  return(data)
})
names(hopsgene_ageburden_cut_tf) <- c("pm_ld", "pn_ld")
table(hopsgene_ageburden_cut_tf$pm_ld[, is.na(tf_sum)]) #16045 1918 
sum(hopsgene_ageburden_cut_tf$pm_ld$tf_sum, na.rm = T) #255191

###gene and TF interactions from dorothea-----------------------
tf_scatter_func <- function(data = pleiotropy_maindata$pn_ld, x_label = "", ifscale = T,
                            y_label = "Gene length (bp)" , group = "pleio_class3", x_var = "Gene length (bp)",
                            x_text = c("Low", "Intermediate", "High"), color_value = pleio_color,
                            stat_test = test, test_list = c("Low pleiotropy", "High pleiotropy"), stat_name){
  
  plot_data <- data[, .(median = sapply(.SD, function(x) as.numeric(mean(x, na.rm = TRUE)))), by = group, .SDcols = x_var, with = T] 
  
  if(ifscale){
    plot <- ggplot() + theme_bw() +
      stat_boxplot(data = data, aes(x = !!sym(group), y = !!sym(x_var)), geom ='errorbar', width = 0.4) +
      geom_boxplot(data = data, aes(x = !!sym(group), y = !!sym(x_var), fill = !!sym(group)), 
                   notch=TRUE, notchwidth = 0.6, color = "grey50", fatten = 0.6, lwd = 0.6, 
                   color = "grey50", outlier.shape = NA) +
      geom_signif(data = data, aes(x = !!sym(group), y = !!sym(x_var)),
                  comparisons = list(test_list), color = "grey20",
                  map_signif_level = TRUE,
                  annotations = stat_test[[stat_name]]) +
      geom_point(data = plot_data, aes(x = !!sym(group), y = median), color = "grey20", size = 2.5) +
      geom_line(data = plot_data, aes(x = !!sym(group), y = median, group = 1), 
                linetype = "dashed", color = "grey20") +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
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
            panel.grid.major = element_blank()) +
      scale_x_discrete(labels = x_text)
  } else {
    plot <- ggplot() + theme_bw() +
      stat_boxplot(data = data, aes(x = !!sym(group), y = !!sym(x_var)), geom ='errorbar', width = 0.4) +
      geom_boxplot(data = data, aes(x = !!sym(group), y = !!sym(x_var), fill = !!sym(group)), 
                   notch=TRUE, notchwidth = 0.6, color = "grey50", fatten = 0.6, lwd = 0.6, 
                   color = "grey50", outlier.shape = NA) +
      geom_signif(data = data, aes(x = !!sym(group), y = !!sym(x_var)),
                  comparisons = list(test_list), color = "grey20",
                  map_signif_level = TRUE,
                  annotations = stat_test[[stat_name]]) +
      geom_point(data = plot_data, aes(x = !!sym(group), y = median), color = "grey20", size = 2.5) +
      geom_line(data = plot_data, aes(x = !!sym(group), y = median, group = 1), 
                linetype = "dashed", color = "grey20") +
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
            panel.grid.major = element_blank()) +
      scale_x_discrete(labels = x_text)
  }
  
  return(list(plot_data = plot_data, plot = plot))
}

test <- list(
  tf_pleio_pn = pairwise.wilcox.test(hopsgene_ageburden_cut_tf$pn_ld$tf_sum, 
                                     hopsgene_ageburden_cut_tf$pn_ld$pleio_class3, 
                                     p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .), 
  tf_pleio_pm = pairwise.wilcox.test(hopsgene_ageburden_cut_tf$pm_ld$tf_sum, 
                                     hopsgene_ageburden_cut_tf$pm_ld$pleio_class3, 
                                     p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .), 
  tf_age = pairwise.wilcox.test(hopsgene_ageburden_cut_tf$pn_ld$tf_sum, 
                                hopsgene_ageburden_cut_tf$pn_ld$age_stage4, 
                                p.adjust.method = "fdr")$p.value[3,1] %>% sprintf("%.2e", .)  
)

tf_sum_box <- list(pleio_n = tf_scatter_func(data = hopsgene_ageburden_cut_tf$pn_ld, 
                                             y_label = "Number of TFs per gene", x_var = "tf_sum", ifscale = T,
                                             stat_name = "tf_pleio_pn"),
                   pleio_m = tf_scatter_func(data = hopsgene_ageburden_cut_tf$pm_ld, 
                                             y_label = "Number of TFs per gene", x_var = "tf_sum", ifscale = T,
                                             stat_name = "tf_pleio_pm"),
                   age = tf_scatter_func(data = hopsgene_ageburden_cut_tf$pn_ld[!is.na(age_stage4),], 
                                         group = "age_stage4", color_value = age_color,
                                         y_label = "Number of TFs per gene", x_var = "tf_sum", ifscale = T, 
                                         x_text = x_text, test_list = age_test, stat_name = "tf_age")) 

tf_sum_box$pleio_n$plot
tf_sum_box$age$plot


ggsave(plot = tf_sum_box$pleio_n$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/tf_sum_pn.pdf",figure_file))
ggsave(plot = tf_sum_box$age$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/tf_sum_age.pdf",figure_file))

##figure4G&4H----------------------
###load data----------------------
enhancer_annota <- fread("genehancer_annotation.csv") #244737
length(unique(enhancer_annota$name))
gene_enhancer_all <- fread("genehancer_all.csv") %>%
  merge(., enhancer_annota[, .(name, elementType)], by.x = "geneHancerIdentifier", by.y = "name", all.x = T) #757871 cre-gene pairs
length(unique(gene_enhancer_all$geneName)) #45749 incuding genes encoding RNA
intersect(unique(gene_enhancer_all$geneName), pleiotropy_maindata$pm_ld$gene) %>% length() #17178
gene_enhancer_table <- gene_enhancer_all[elementType == "Enhancer", .(enhancer_n = .N), by = geneName] 

cre_merge_all <- lapply(list("pm_ld", "pn_ld"), function(x){
  data <- pleiotropy_maindata[[x]]
  data_merge <- merge(data[, .(gene, use_score, age_stage4, age_stage7, age_stage4_num, age_stage7_num,
                               pleio_class3, pleio_class5, ensemblid)], 
                      gene_enhancer_table, by.x = "gene", by.y = "geneName", all.x = T) 
  data_merge$ifenhancer <- ifelse(!is.na(data_merge$enhancer_n), 1, 0)
  return(data_merge)
})

names(cre_merge_all) <- c("pm_ld", "pn_ld")
nrow(cre_merge_all$pn_ld[!is.na(enhancer_n)]) #17035
sum(cre_merge_all$pn_ld$enhancer_n, na.rm = T) #368039
unique(gene_enhancer_all[geneName %in% cre_merge_all$pn_ld[!is.na(enhancer_n), gene], geneHancerIdentifier]) %>% length() #199241
###number of enhancers per gene across pleiotropy groups and age_stages----------------------
enhancer_box_func <- function(data = pleiotropy_maindata$pn_ld, x_label = "", scale_y = F,
                              y_label = "Gene length (bp)" , group = "pleio_class3", x_var = "Gene length (bp)",
                              x_text = c("Low", "Intermediate", "High"), color_value = pleio_color,
                              stat_test = test, test_list = c("Low pleiotropy", "High pleiotropy"), stat_name){
  
  plot_data <- data[, .(median = sapply(.SD, function(x) as.numeric(mean(x, na.rm = TRUE)))), by = group, .SDcols = x_var, with = T] 
  
  if(scale_y == F){
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
      #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
      #              labels = trans_format("log10", math_format(10^.x))) +
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
            panel.grid.major = element_blank()) +
      scale_x_discrete(labels = x_text)
  } else {
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
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x))) +
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
            panel.grid.major = element_blank()) +
      scale_x_discrete(labels = x_text)
    
  }
  
  return(list(plot_data = plot_data, plot = plot))
}


test <- list(
  pleio_pn = pairwise.wilcox.test(cre_merge_all$pn_ld$enhancer_n, 
                                  cre_merge_all$pn_ld$pleio_class3,
                                  p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .), 
  pleio_pm = pairwise.wilcox.test(cre_merge_all$pm_ld$enhancer_n, 
                                  cre_merge_all$pm_ld$pleio_class3,
                                  p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .),
  age = pairwise.wilcox.test(cre_merge_all$pm_ld$enhancer_n, 
                             cre_merge_all$pm_ld$age_stage4,
                             p.adjust.method = "fdr")$p.value[3,1] %>% sprintf("%.2e", .)
)

cre_box <- list(pleio_pn = enhancer_box_func(data = cre_merge_all$pn_ld, color_value = pleio_color, scale_y = T,
                                             y_label = "Number of enhancers per gene", x_var = "enhancer_n", group = "pleio_class3",
                                             x_text = c("Low", "Intermediate", "High"), 
                                             stat_name = "pleio_pn", test_list = c("Low pleiotropy", "High pleiotropy")),
                pleio_pm = enhancer_box_func(data = cre_merge_all$pm_ld, color_value = pleio_color, scale_y = T,
                                             y_label = "Number of enhancers per gene", x_var = "enhancer_n", group = "pleio_class3",
                                             x_text = c("Low", "Intermediate", "High"), 
                                             stat_name = "pleio_pm", test_list = c("Low pleiotropy", "High pleiotropy")),
                age = enhancer_box_func(data = cre_merge_all$pm_ld[!is.na(age_stage7),], color_value = age_color, scale_y = T,
                                        y_label = "Number of enhancers per gene", x_var = "enhancer_n", group = "age_stage4",
                                        x_text = c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria"), 
                                        stat_name = "age", test_list = c("Euteleostomi", "Eutheria"))) 


cre_box$pleio_pn$plot
cre_box$age$plot

ggsave(plot = cre_box$pleio_pn$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/enhancer_box_pn.pdf",figure_file))
ggsave(plot = cre_box$age$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/enhancer_box_age.pdf",figure_file))


##figure4I&4J----------------------
###load data------------------
ogee_pip <- fread("connectivity.txt") %>%
  .[taxon_id == 9606, .(gene, ogeepip_score = score, ogee_connectivity = connectivity)]

pip_data <- lapply(list("pn_ld", "pm_ld"), function(x, data = pleiotropy_maindata){
  data_merge <- merge(data[[x]][, .(gene, pleio_class3, age_stage4)], ogee_pip, by = "gene", all.x = T)
  return(data_merge)
})

names(pip_data) <- c("pn_ld", "pm_ld")

###gene connectivities across the L-GPS, I-GPS and H-GPS groups--------------
library(viridis)
pip_violin_func <- function(data = pleiotropy_maindata$pn_ld, x_label = "", 
                            y_label = "Gene length (bp)" , group = "pleio_class3", x_var = "Gene length (bp)",
                            x_text = c("Low", "Intermediate", "High"), color_value = pleio_color,
                            stat_test = test, test_list = c("Low pleiotropy", "High pleiotropy"), stat_name){
  
  plot_data <- data[, .(median = sapply(.SD, function(x) as.numeric(mean(x, na.rm = TRUE)))), by = group, .SDcols = x_var, with = T] 
  
  plot <- ggplot(data = data, aes(x = !!sym(group), y = !!sym(x_var), fill = !!sym(group))) + theme_bw() +
    geom_violin(width=1) +
    geom_boxplot(width=0.4, color="grey50", outlier.shape = NA, lwd = 0.5) +
    geom_signif(data = data, aes(x = !!sym(group), y = !!sym(x_var)),
                comparisons = list(test_list), color = "grey20",
                map_signif_level = TRUE,
                annotations = stat_test[[stat_name]]) +
    scale_fill_viridis(discrete = TRUE) +
    geom_point(data = plot_data, aes(x = !!sym(group), y = median), color = "grey20", size = 2.5) +
    geom_line(data = plot_data, aes(x = !!sym(group), y = median, group = 1), 
              linetype = "dashed", color = "grey20") +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
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
          panel.grid.major = element_blank()) +
    scale_x_discrete(labels = x_text)
  
  return(list(plot_data = plot_data, plot = plot))
}

test <- list(
  pip_pleio_pn = pairwise.wilcox.test(pip_data$pn_ld$ogee_connectivity, 
                                      pip_data$pn_ld$pleio_class3, 
                                      p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .),
  pip_pleio_pm = pairwise.wilcox.test(pip_data$pm_ld$ogee_connectivity, 
                                      pip_data$pm_ld$pleio_class3, 
                                      p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .),
  pip_age = pairwise.wilcox.test(pip_data$pn_ld$ogee_connectivity, 
                                 pip_data$pn_ld$age_stage4, 
                                 p.adjust.method = "fdr")$p.value[3,1] %>% sprintf("%.2e", .) 
)

pip_box <- list(pleio_n = pip_violin_func(data = pip_data$pn_ld, stat_name = "pip_pleio_pn",
                                          y_label = "Number of PPIs per gene", x_var = "ogee_connectivity"),
                pleio_m = pip_violin_func(data = pip_data$pm_ld, stat_name = "pip_pleio_pm",
                                          y_label = "Number of PPIs per gene", x_var = "ogee_connectivity"),
                age_n = pip_violin_func(data = pip_data$pn_ld[!is.na(age_stage4),], group = "age_stage4",
                                        y_label = "Number of PPIs per gene", x_var = "ogee_connectivity", 
                                        color_value = age_color, x_text = x_text,
                                        test_list = age_test, stat_name = "pip_age")) 

pip_box$pleio_n$plot
pip_box$age_n$plot

ggsave(plot = pip_box$pleio_n$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/pip_pn.pdf",figure_file))
ggsave(plot = pip_box$age_n$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/pip_age.pdf",figure_file))
