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
#3 figure3----------------------
##function-------------------
group_scatter_func <- function(data = pleiotropy_maindata$pn_ld, x_label = "", ifscale = T,
                               y_label = "Gene length (bp)" , group = "pleio_class3", x_var = "Gene length (bp)",
                               x_text = c("Low", "Intermediate", "High"), color_value = age_color, 
                               stat_test = test, test_list = c("Low pleiotropy", "High pleiotropy")){
  
  plot_data <- data[, .(median = sapply(.SD, function(x) as.numeric(mean(x, na.rm = TRUE)))), 
                    by = group, .SDcols = x_var, with = T] 
  
  if(ifscale){
    plot <- ggplot() + theme_bw() +
      geom_boxplot(data = data, aes(x = !!sym(group), y = !!sym(x_var), fill = !!sym(group)), 
                   fatten = 0.6, lwd = 0.6, color = "grey50", outlier.alpha = 0.3) +
      geom_signif(data = data, aes(x = !!sym(group), y = !!sym(x_var)),
                  comparisons = list(test_list), color = "grey20",
                  map_signif_level = TRUE, 
                  annotations = stat_test[[x_var]]) +
      geom_point(data = plot_data, aes(x = !!sym(group), y = median), color = "grey20", size = 2.5) +
      geom_line(data = plot_data, aes(x = !!sym(group), y = median, group = 1), linetype = "dashed", color = "grey20") +
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
  } else {
    plot <- ggplot() + theme_bw() +
      geom_boxplot(data = data, aes(x = !!sym(group), y = !!sym(x_var), fill = !!sym(group)), 
                   lwd = 0.6, fatten = 0.6, color = "grey50", outlier.alpha = 0.3) +
      geom_signif(data = data, aes(x = !!sym(group), y = !!sym(x_var)), 
                  comparisons = list(test_list), color = "grey20",
                  map_signif_level = TRUE, 
                  annotations = stat_test[[x_var]]) +
      geom_point(data = plot_data, aes(x = !!sym(group), y = median), color = "grey20", size = 2.5) +
      geom_line(data = plot_data, aes(x = !!sym(group), y = median, group = 1), linetype = "dashed", color = "grey20", linewidth = 1) +
      labs(x=x_label, y= y_label) +
      scale_color_manual(values = color_value) +
      scale_fill_manual(values = color_value) +
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

constraint_plot_func <- function(data = pleiotropy_maindata$pn_ld, x_label = "", ifscale = T,
                                 y_label = "Gene length (bp)" , group = "pleio_class3", x_var = "LOEUF",
                                 x_text = c("Low", "Intermediate", "High"), color_value = age_color,
                                 stat_test = test, test_list = c("Low pleiotropy", "High pleiotropy")){
  
  if(ifscale){
    plot <- ggplot(data = data, aes(x = !!sym(group), y = !!sym(x_var))) + theme_bw() +
      geom_point(aes(color = !!sym(group)), 
                 position = position_jitter(width = 0.35, height = 0), size = 0.5, alpha = 0.15) +
      geom_boxplot(aes(fill = !!sym(group)), outlier.shape = NA, color = "grey40", 
                   fatten = 0.8, linewidth = 0.8, coef = 0, width = 0.3) +
      geom_signif(comparisons = list(test_list), color = "grey20",
                  map_signif_level = TRUE,
                  annotations = stat_test[[x_var]]) +
      stat_summary(fun = mean, geom = "line", group = 1, size = 1, color = "grey20", linetype = "dashed") +
      stat_summary(fun = mean, geom = "point", size = 2.5, color = "grey20") +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    #labels = function(x) sprintf("%.1f", x),
                    labels = trans_format("log10", math_format(10^.x))) +
      scale_color_manual(values = color_value) +
      scale_fill_manual(values = color_value) +
      labs(x=x_label, y= y_label) +
      theme(text = element_text(size = 12, color = "black", face = "bold"),
            legend.position = "none",
            axis.title.y = element_text(size = 15, face = "bold"), 
            axis.title.x = element_blank(), 
            axis.text = element_text(size = 12, color = "black", face = "bold"), 
            axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
            panel.grid = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            panel.grid.major = element_blank()) +
      scale_x_discrete(labels = x_text) 
  } else {
    plot <- ggplot(data = data, aes(x = !!sym(group), y = !!sym(x_var), color = !!sym(group))) + theme_bw() +
      geom_point(position = position_jitter(width = 0.35, height = 0), size = 0.5, alpha = 0.15) +
      geom_boxplot(aes(fill = !!sym(group)), outlier.shape = NA, color = "grey40", 
                   fatten = 0.8, linewidth = 0.8, coef = 0, width = 0.3) +
      geom_signif(comparisons = list(test_list), color = "grey20",
                  map_signif_level = TRUE,
                  annotations = stat_test[[x_var]]) +
      stat_summary(fun = mean, geom = "line", group = 1, color = "grey20", linetype = "dashed") +
      stat_summary(fun = mean, geom = "point", size = 2.5, color = "grey20") +
      scale_color_manual(values = color_value) +
      scale_fill_manual(values = color_value) +
      theme(text = element_text(size = 12, color = "black", face = "bold"),
            legend.position = "none",
            axis.title.y = element_text(size = 15, face = "bold"), 
            axis.title.x = element_blank(), 
            axis.text = element_text(size = 12, color = "black", face = "bold"), 
            axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
            panel.grid = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            panel.grid.major = element_blank()) +
      scale_x_discrete(labels = x_text)
  }
  
  return(plot)
}

group_conpercent_func <- function(data = pleiotropy_maindata$pn_ld, x_label = "", ifscale = T,
                                  y_label = "Percents of intolerant genes classified by LOEUF (%)", 
                                  group = "pleio_class3", x_var = "LOEUF_class",
                                  x_text = c("Low", "Intermediate", "High"),
                                  color_value, 
                                  stat_test = test, test_list = c("Low pleiotropy", "High pleiotropy")){
  
  freq_data <- table(data[,get(x_var)], data[, get(group)]) %>% 
    prop.table(., margin = 2) %>% 
    as.data.table(.) %>% 
    setnames(., names(.), c("feature", "pleio", "freq")) %>%
    .[feature == "Intolerant", list(pleio, freq = freq*100,freq_name = sprintf("%.2f", freq*100))]
  
  freq_data$pleio <- factor(freq_data$pleio, levels = freq_data$pleio) 
  
  plot <- ggplot(data = freq_data) + theme_bw() +
    geom_bar(aes(x = pleio, y = freq, fill = pleio), stat = "identity", width = 0.9) +
    scale_fill_manual(values = color_value) +
    geom_signif(aes(x = pleio, y = freq),
                comparisons = list(test_list), color = "grey20",
                #map_signif_level = TRUE,
                annotations = stat_test[[x_var]]) +
    geom_text(aes(x = pleio, y = freq, label = freq_name), vjust = 1.5, color = "grey30", fontface = "bold", size = 6) +
    labs(y = y_label, title = " ", x = x_label) +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none",
          panel.border = element_blank(), 
          panel.grid.major = element_blank())+
    scale_x_discrete(labels = x_text)
  
  return(list(plot_data = freq_data, plot = plot))
}

dnds_func <- function(data = pleiotropy_maindata$pn_ld, x_label = "GPS-N bins", 
                      y_label = "Percents of intolerant genes classified by LOEUF (%)", 
                      group = "pleio_class3", x_var = "dN/dS ratio mouse", 
                      color_value = pleio_color, x_text = c("Low", "Intermediate", "High"),
                      stat_test = test, test_list = c("Low pleiotropy", "High pleiotropy")){
  
  plot <- ggplot(data = data, aes(x = !!sym(group), y = !!sym(x_var), fill = !!sym(group))) + theme_bw() +
    geom_jitter(aes(color = !!sym(group)), size=0.1, alpha=0.15) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(outlier.shape = NA, color = "grey40", 
                 fatten = 0.7, linewidth = 0.7) +
    geom_signif(comparisons = list(test_list), color = "grey20",
                map_signif_level = TRUE,
                annotations = stat_test[[x_var]]) +
    stat_summary(fun = mean, geom = "line", group = 1, color = "grey20", linetype = "dashed") +
    stat_summary(fun = mean, geom = "point", size = 2.5, color = "grey20") +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values = color_value) +
    scale_fill_manual(values = color_value) +
    labs(y = y_label, title = " ", x = x_label) +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          axis.title = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "none", panel.border = element_blank(), 
          panel.grid.major = element_blank()) +
    scale_x_discrete(labels = x_text)
  
  
  return(plot)
}
##GPS-N (figure 3A, 3C, 3E, 3G)-------------------------
usedata <-  pleiotropy_maindata$pn_ld
###length related metrics + GC content--------------------
length_var <- c("Gene length (bp)", "Transcript length (bp)", "CDS Length",  "Protein size (aa)", "CDS/Transcript Length ratio",
                "Transcript count", "Exon Counts", "Intron Counts", "Number of SNPs (Transcript)", "Number of SNPs (Gene)", 
                "GC content")
test <- lapply(as.list(length_var), function(x, data = usedata){
  pairwise.wilcox.test(data[,get(x)], data$pleio_class3, 
                       p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .)   
})
names(test) <- length_var
test
pn_feature_list <- list(gene_length = group_scatter_func(data = usedata, color_value = pleio_color,  
                                                         y_label = "Gene length (base pair)", x_var = "Gene length (bp)"), 
                        transcript_length = group_scatter_func(data = usedata, color_value = pleio_color,
                                                               y_label = "Transcript length (bp)", x_var = "Transcript length (bp)"),
                        CDS_length = group_scatter_func(data = usedata, color_value = pleio_color,
                                                        y_label = "CDS Length (bp)", x_var = "CDS Length"),
                        protein_length = group_scatter_func(data = usedata, color_value = pleio_color,
                                                            y_label = "Protein size (aa)", x_var = "Protein size (aa)"),
                        CDS_transcript_ratio = group_scatter_func(data = usedata, ifscale = F, color_value = pleio_color,
                                                                  y_label = "CDS/Transcript Length ratio", x_var = "CDS/Transcript Length ratio"),
                        transcript_count = group_scatter_func(data = usedata, color_value = pleio_color,
                                                              y_label = "Transcript count", x_var = "Transcript count"),
                        exon_count = group_scatter_func(data = usedata, color_value = pleio_color,
                                                        y_label = "Exon Counts", x_var = "Exon Counts"),
                        intron_count = group_scatter_func(data = usedata, color_value = pleio_color,
                                                          y_label = "Intron Counts", x_var = "Intron Counts"),
                        number_snp_gene = group_scatter_func(data = usedata, color_value = pleio_color,
                                                             y_label = "Number of SNPs (Gene)", x_var = "Number of SNPs (Gene)"),
                        number_snp_transcript = group_scatter_func(data = usedata, color_value = pleio_color,
                                                                   y_label = "Number of SNPs (Transcript)", x_var = "Number of SNPs (Transcript)"),
                        GC_content = group_scatter_func(data = usedata, ifscale = F, color_value = pleio_color,
                                                        y_label = "GC content (%)", x_var = "GC content"))

pn_feature_list$gene_length$plot

ggsave(plot = pn_feature_list$gene_length$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/gene_length_pn.pdf",figure_file))

###gene constraint metrics--------------------
####intolerance to mutaion-------------------
test <- sapply(list("LOEUF", "Missense/Synonymous ratio", "Nonsense/Synonymous ratio",
                    "dN/dS ratio mouse", "dN/dS ratio Chimp", "dN/dS ratio Gorilla"), 
               function(x, data = usedata){
                 pairwise.wilcox.test(data[,get(x)], data$pleio_class3, p.adjust.method = "fdr", 
                                      exact = TRUE)$p.value[2,1] %>% sprintf("%.2e", .)         
               })
names(test) <- c("LOEUF", "Missense/Synonymous ratio", "Nonsense/Synonymous ratio",
                 "dN/dS ratio mouse", "dN/dS ratio Chimp", "dN/dS ratio Gorilla")

pn_mutation_list <- list(
  lof = constraint_plot_func(data =usedata[!is.na(LOEUF),], 
                             ifscale = F, color_value = pleio_color,
                             y_label = "LOEUF score", x_var = "LOEUF", x_label = ""),
  missense = constraint_plot_func(data =usedata[!is.na(`Missense/Synonymous ratio`),], 
                                  ifscale = T, color_value = pleio_color,
                                  y_label = "Missense/synonymous ratio", x_var = "Missense/Synonymous ratio", x_label = ""),
  nonsense = constraint_plot_func(data =usedata[!is.na(`Nonsense/Synonymous ratio`),], 
                                  ifscale = T, color_value = pleio_color,
                                  y_label = "Nonsense/synonymous ratio", x_var = "Nonsense/Synonymous ratio", x_label = "")
)

####dn/ds ratios-------------------------------
pn_dnds_list <- list(
  mouse = dnds_func(data =usedata, x_label = "", 
                    y_label = "dN/dS ratio for human-mouse pair", 
                    x_text = c("Low", "Intermediate", "High"), 
                    group = "pleio_class3", x_var = "dN/dS ratio mouse", color_value = pleio_color),
  chimp = dnds_func(data =usedata, x_label = "", 
                    y_label = "dN/dS ratio for human-chimpanzee pair", 
                    x_text = c("Low", "Intermediate", "High"), 
                    group = "pleio_class3", x_var = "dN/dS ratio Chimp", color_value = pleio_color),
  gorilla = dnds_func(data =usedata, x_label = "", 
                      y_label = "dN/dS ratio for human-gorilla pair", 
                      x_text = c("Low", "Intermediate", "High"), 
                      group = "pleio_class3", x_var = "dN/dS ratio Gorilla", color_value = pleio_color)
)


pn_dnds_list$mouse
ggsave(plot = pn_dnds_list$mouse, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/dnds_mouse_pn.pdf",figure_file))

####intolerant gene percent----------------------
loeuf_count <- table(usedata[, .(pleio_class3, LOEUF_class)])
test <- list(
  LOEUF_class = pairwise.prop.test(loeuf_count[,"Intolerant"], 
                                   c(sum(loeuf_count["Low pleiotropy", ]), sum(loeuf_count["Intermediate pleiotropy",]), sum(loeuf_count["High pleiotropy",])),
                                   p.adjust.method = "fdr")$p.value[2,1] %>% sprintf("%.2e", .))

pn_constrait_percent_list <- list(lof_tolerant = group_conpercent_func(data =usedata,
                                                                       y_label = "Proportion of intolerant genes\ndefined by LOEUF (%)", 
                                                                       x_var = "LOEUF_class", color_value = pleio_color))

pn_constrait_percent_list$lof_tolerant$plot

ggsave(plot = pn_constrait_percent_list$lof_tolerant$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/lof_tolerant_pn.pdf",figure_file))

##AGE (figure 3B, 3D, 3F, 3H)-------------------------
usedata <- pleiotropy_maindata$pn_ld[!is.na(age_stage4),]
###length related metrics + GC content--------------------
length_var <- c("Gene length (bp)", "Transcript length (bp)", "CDS Length",  "Protein size (aa)", "CDS/Transcript Length ratio",
                "Transcript count", "Exon Counts", "Intron Counts", "Number of SNPs (Transcript)", "Number of SNPs (Gene)", 
                "GC content")
test <- lapply(as.list(length_var), function(x, data = usedata){
  pairwise.wilcox.test(data[,get(x)], data$age_stage4, 
                       p.adjust.method = "fdr", exact = TRUE)$p.value[3,1] %>% sprintf("%.2e", .)   
})
names(test) <- length_var

age_feature_list <- list(gene_length = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color,
                                                          y_label = "Gene length (base pair)", x_var = "Gene length (bp)", 
                                                          x_text = x_text, test_list = age_test),
                         transcript_length = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color,
                                                                y_label = "Transcript length (bp)", x_var = "Transcript length (bp)", 
                                                                x_text = x_text, test_list = age_test),
                         CDS_length = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color,
                                                         y_label = "CDS Length (bp)", x_var = "CDS Length", 
                                                         x_text = x_text, test_list = age_test),
                         protein_length = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color,
                                                             y_label = "Protein size (aa)", x_var = "Protein size (aa)", 
                                                             x_text = x_text, test_list = age_test),
                         CDS_transcript_ratio = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color, ifscale = F, 
                                                                   y_label = "CDS/Transcript Length ratio", x_var = "CDS/Transcript Length ratio",
                                                                   x_text = x_text, test_list = age_test),
                         transcript_count = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color,
                                                               y_label = "Transcript count", x_var = "Transcript count",
                                                               x_text = x_text, test_list = age_test),
                         exon_count = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color,
                                                         y_label = "Exon Counts", x_var = "Exon Counts",
                                                         x_text = x_text, test_list = age_test),
                         intron_count = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color,
                                                           y_label = "Intron Counts", x_var = "Intron Counts",
                                                           x_text = x_text, test_list = age_test),
                         number_snp_gene = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color,
                                                              y_label = "Number of SNPs (Gene)", x_var = "Number of SNPs (Gene)",
                                                              x_text = x_text, test_list = age_test),
                         number_snp_transcript = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color,
                                                                    y_label = "Number of SNPs (Transcript)", x_var = "Number of SNPs (Transcript)",
                                                                    x_text = x_text, test_list = age_test),
                         GC_content = group_scatter_func(data = usedata, group = "age_stage4", color_value = age_color, ifscale = F, 
                                                         y_label = "GC content (%)", x_var = "GC content", 
                                                         x_text = x_text, test_list = age_test))

age_feature_list$gene_length$plot

ggsave(plot = age_feature_list$gene_length$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/gene_length_age.pdf",figure_file))


###gene constraint metrics--------------------
####intolerance to mutaion-------------------
test <- sapply(list("LOEUF", "Missense/Synonymous ratio", "Nonsense/Synonymous ratio",
                    "dN/dS ratio mouse", "dN/dS ratio Chimp", "dN/dS ratio Gorilla"), 
               function(x, data = usedata){
                 pairwise.wilcox.test(data[,get(x)], data$age_stage4, p.adjust.method = "fdr")$p.value[3,1] %>% sprintf("%.2e", .)         
               })
names(test) <- c("LOEUF", "Missense/Synonymous ratio", "Nonsense/Synonymous ratio",
                 "dN/dS ratio mouse", "dN/dS ratio Chimp", "dN/dS ratio Gorilla")

age_mutation_list <- list(
  lof = constraint_plot_func(data =usedata[!is.na(LOEUF),], group = "age_stage4", 
                             ifscale = F, color_value = age_color,
                             y_label = "LOEUF score", x_var = "LOEUF", x_label = "",
                             x_text = x_text, test_list = age_test), 
  
  missense = constraint_plot_func(data =usedata[!is.na(`Missense/Synonymous ratio`),], group = "age_stage4",
                                  ifscale = T, color_value = age_color,
                                  y_label = "Missense/synonymous ratio", x_var = "Missense/Synonymous ratio", x_label = "",
                                  x_text = x_text, test_list = age_test),
  
  nonsense = constraint_plot_func(data =usedata[!is.na(`Nonsense/Synonymous ratio`),], group = "age_stage4",
                                  ifscale = T, color_value = age_color,
                                  y_label = "Nonsense/synonymous ratio", x_var = "Nonsense/Synonymous ratio", x_label = "",
                                  x_text = x_text, test_list = age_test)
)

####dn/ds ratios-------------------------------
age_dnds_list <- list(
  mouse = dnds_func(data =usedata, group = "age_stage4", color_value = age_color, 
                    y_label = "dN/dS ratio for human-mouse pair", x_label = "", x_var = "dN/dS ratio mouse",
                    x_text = x_text, test_list = age_test),
  
  chimp = dnds_func(data =usedata, group = "age_stage4", color_value = age_color,
                    y_label = "dN/dS ratio for human-chimpanzee pair", x_label = "", x_var = "dN/dS ratio Chimp", 
                    x_text = x_text, test_list = age_test),
  
  gorilla = dnds_func(data =usedata, group = "age_stage4",  color_value = age_color,
                      y_label = "dN/dS ratio for human-gorilla pair", x_label = "", x_var = "dN/dS ratio Gorilla", 
                      x_text = x_text, test_list = age_test)
)


age_dnds_list$mouse

ggsave(plot = age_dnds_list$mouse, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/dnds_mouse_age.pdf",figure_file))


####intolerant gene percent----------------------
loeuf_count <- table(usedata[, .(age_stage4, LOEUF_class)])
test <- list(
  LOEUF_class = pairwise.prop.test(loeuf_count[,"Intolerant"], 
                                   c(sum(loeuf_count["Euteleostomi", ]), sum(loeuf_count["Tetrapoda",]), 
                                     sum(loeuf_count["Amniota",]), sum(loeuf_count["Eutheria",])),
                                   p.adjust.method = "fdr")$p.value[3,1] %>% sprintf("%.2e", .)
)

age_constrait_percent_list <- list(lof_tolerant = group_conpercent_func(data =usedata, group = "age_stage4", color_value = age_color,
                                                                        y_label = "Proportion of intolerant genes\ndefined by LOEUF (%)", 
                                                                        x_var = "LOEUF_class", 
                                                                        x_text = x_text, test_list = age_test))

age_constrait_percent_list$lof_tolerant$plot

ggsave(plot = age_constrait_percent_list$lof_tolerant$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = sprintf("%s/lof_tolerant_age.pdf",figure_file))

##figure3I-----------------------
hsd_box_func <- function(data = pleiotropy_maindata$pn_ld, x_label = "", scale_y = F,
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
      scale_y_break(c(10,20), scales = 0.3) +
      scale_y_continuous(limits = c(0,38), breaks = c(0,5,10,20,30), labels = c(0,5,10,20,30)) +
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
      scale_x_discrete(labels = x_text)}
  
  return(list(plot_data = plot_data, plot = plot))
}


test <- list(
  pleio_pn = pairwise.wilcox.test(pleiotropy_maindata$pn_ld$use_score, 
                                  pleiotropy_maindata$pn_ld$ifhsd2,
                                  p.adjust.method = "fdr")$p.value %>% sprintf("%.2e", .), 
  pleio_pm = pairwise.wilcox.test(pleiotropy_maindata$pm_ld$use_score, 
                                  pleiotropy_maindata$pm_ld$ifhsd2,
                                  p.adjust.method = "fdr")$p.value %>% sprintf("%.2e", .) 
)

hsd_box <- list(pleio_pn = hsd_box_func(data = pleiotropy_maindata$pn_ld[, ifhsd2 := factor(ifhsd2, levels = c("Singletons", "Duplicates"))], 
                                        color_value = mycolor,
                                        y_label = "Gene pleiotropic score for number", x_var = "use_score", group = "ifhsd2",
                                        x_text = c("Singletons", "Duplicates"),
                                        stat_name = "pleio_pn", test_list = c("Singletons", "Duplicates")),
                pleio_pm = hsd_box_func(data = pleiotropy_maindata$pm_ld[, ifhsd2 := factor(ifhsd2, levels = c("Singletons", "Duplicates"))], 
                                        color_value = mycolor, scale_y = T,
                                        y_label = "Gene pleiotropic scores for magnitude", x_var = "use_score", group = "ifhsd2",
                                        x_text = c("Singletons", "Duplicates"),
                                        stat_name = "pleio_pm", test_list = c("Singletons", "Duplicates"))) 


hsd_box$pleio_pn$plot

ggsave(plot = hsd_box$pleio_pn$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/hsd_pleio_pn.pdf",figure_file))

##figure3H-----------------------
###manage data---------------------
dup_pairs <- pleiotropy_maindata$pn_ld[!is.na(group_id), 
                                       .(N = .N, gene, age_stage7, age_stage7_mya, age_stage4,
                                         age_genetree_num, gene_age_num), by = group_id][N>1,] #2439 genes with 655 groups

selected_genes <- dup_pairs %>%
  group_by(group_id) %>%
  summarize(
    selected_gene = {
      if (!all(is.na(age_stage7_mya))) {
        max_age_stage7_mya <- max(age_stage7_mya, na.rm = TRUE)
        genes_agestage7 <- gene[age_stage7_mya == max_age_stage7_mya]
        
        if (length(genes_agestage7) == 1) {
          genes_agestage7
        } else {
          # Resolve ties using age_genetree_num
          if (!all(is.na(age_genetree_num[age_stage7_mya == max_age_stage7_mya]))) {
            min_age_genetree_num <- min(age_genetree_num[age_stage7_mya == max_age_stage7_mya], na.rm = TRUE)
            genes_genetree <- na.omit(gene[age_stage7_mya == max_age_stage7_mya & age_genetree_num == min_age_genetree_num])
            
            if (length(genes_genetree) == 1) {
              genes_genetree
            } else {
              # Resolve ties further using gene_age_num
              if (!all(is.na(gene_age_num[age_stage7_mya == max_age_stage7_mya & age_genetree_num == min_age_genetree_num]))) {
                max_gene_age_num <- max(gene_age_num[age_stage7_mya == max_age_stage7_mya & age_genetree_num == min_age_genetree_num], na.rm = TRUE)
                genes_origin <- na.omit(gene[age_stage7_mya == max_age_stage7_mya & age_genetree_num == min_age_genetree_num & gene_age_num == max_gene_age_num])
                
                if (length(genes_origin) == 1) {
                  genes_origin
                } else {
                  NA
                }
              } else {
                NA
              }
            }
          } else {
            NA
          }
        }
      } else if (!all(is.na(age_genetree_num))) {
        # Step 2: Select by smallest age_genetree_num
        min_age_genetree_num <- min(age_genetree_num, na.rm = TRUE)
        genes_genetree <- gene[age_genetree_num == min_age_genetree_num]
        
        if (length(genes_genetree) == 1) {
          genes_genetree
        } else {
          # Step 3: Select by largest gene_age_num
          if (!all(is.na(gene_age_num[age_genetree_num == min_age_genetree_num]))) {
            max_gene_age_num <- max(gene_age_num[age_genetree_num == min_age_genetree_num], na.rm = TRUE)
            genes_origin <- na.omit(gene[age_genetree_num == min_age_genetree_num & gene_age_num == max_gene_age_num])
            
            if (length(genes_origin) == 1) {
              genes_origin
            } else {
              NA
            }
          } else {
            NA
          }
        }
      } else if (!all(is.na(gene_age_num))) {
        # Step 3: Select by largest gene_age_num (if all above fail)
        max_gene_age_num <- max(gene_age_num, na.rm = TRUE)
        genes_origin <- gene[gene_age_num == max_gene_age_num]
        
        if (length(genes_origin) == 1) {
          genes_origin
        } else {
          NA
        }
      } else {
        # All columns are NA
        NA
      }
    },
    .groups = "drop"
  )

dup_pairs_cut <- dup_pairs[group_id %in% selected_genes[!is.na(selected_genes$selected_gene),]$group_id,]
length(unique(dup_pairs_cut$group_id))

progenitor_data <- lapply(list("pm_ld", "pn_ld"), function(score_name){
  data_merge <- merge(dup_pairs_cut[, .(gene)], pleiotropy_maindata[[score_name]], by = "gene") %>%
    .[, progenitor := ifelse(gene %in% selected_genes$selected_gene, 1, 2)] %>%
    .[, progenitor := factor(progenitor, levels = c(1,2), labels = c("Progenitor", "Copy"))]
  return(data_merge)
})
names(progenitor_data) <- c("pm_ld", "pn_ld")

####plot data---------------
progenitor_box_func <- function(data = pleiotropy_maindata$pn_ld, x_label = "", scale_y = F,
                                y_label = "Gene length (bp)" , group = "pleio_class3", x_var = "Gene length (bp)",
                                x_text = c("Low", "Intermediate", "High"), color_value = pleio_color,
                                stat_test = test, test_list = c("Low pleiotropy", "High pleiotropy"), stat_name){
  
  plot_data <- data[, .(median = sapply(.SD, function(x) as.numeric(mean(x, na.rm = TRUE)))), by = group, .SDcols = x_var, with = T] 
  
  if(scale_y == F){
    plot <- ggplot() + theme_bw() +
      geom_point(data = data, aes(x = !!sym(group), y = !!sym(x_var), group = group_id), color = "grey90") +
      geom_line(data = data, aes(x = !!sym(group), y = !!sym(x_var), group = group_id), color = "grey90", linewidth = 0.2) +
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
      geom_point(data = data, aes(x = !!sym(group), y = !!sym(x_var), group = group_id), color = "grey90") +
      geom_line(data = data, aes(x = !!sym(group), y = !!sym(x_var), group = group_id), color = "grey90", linewidth = 0.2) +
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
      scale_y_break(c(10,20), scales = 0.3) +
      scale_y_continuous(limits = c(3,38), breaks = c(3,5,10,20,30), labels = c(3,5,10,20,30)) +
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
      scale_x_discrete(labels = x_text)}
  
  return(list(plot_data = plot_data, plot = plot))
}
test <- lapply(list("pm_ld", "pn_ld"), function(score_name){
  data <- progenitor_data[[score_name]]
  progenitor_wide <- merge(data[progenitor == "Progenitor", .(use_score, group_id, gene, progenitor)], 
                           data[progenitor == "Copy", .(use_score, group_id, gene, progenitor)],
                           by = "group_id", all.x = T, all.y = T) %>%
    .[!is.na(use_score.x),]
  
  pvalue <- wilcox.test(x = progenitor_wide$use_score.x,
                        y = progenitor_wide$use_score.y,
                        paired = T)$p.value %>% sprintf("%.2e",.)
  return(pvalue)
})

names(test) <- c("progenitor_pm", "progenitor_pn")
test

hsd_progenitor_box <- list(pleio_pn = progenitor_box_func(data = progenitor_data$pn, color_value = mycolor[c(3,4)],
                                                          y_label = "Gene pleiotropic score for number", x_var = "use_score", group = "progenitor",
                                                          x_text = c("Ancestral\nparalogs", "Recent\nparalogs"),
                                                          stat_name = "progenitor_pn", test_list = c("Progenitor", "Copy")),
                           pleio_pm = progenitor_box_func(data = progenitor_data$pm, color_value = mycolor[c(3,4)], scale_y = T,
                                                          y_label = "Gene pleiotropic scores for magnitude", x_var = "use_score", group = "progenitor",
                                                          x_text = c("Ancestral\nparalogs", "Recent\nparalogs"),
                                                          stat_name = "progenitor_pm", test_list = c("Progenitor", "Copy"))) 
hsd_progenitor_box$pleio_pn$plot

ggsave(plot = hsd_progenitor_box$pleio_pn$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = sprintf("%s/progenitor_pleio_pn.pdf",figure_file))

##figure3K-----------------------------
###run the random forest (long time, could skip and load the results provided)---------------
library(randomForest)
forest_result <- lapply(list("pn_ld", "pm_ld"), function(x, data = pleiotropy_maindata){
  data = data[[x]]
  covar = c("age_stage7_mya","Gene length (bp)", "CDS Length", 
            "CDS/Transcript Length ratio", 
            "Transcript count", "Exon Counts", "GC content", 
            "Number of SNPs (Gene)", 
            "Missense/Synonymous ratio", "Nonsense/Synonymous ratio", "LOEUF", "pLI", "ifhsd2")
  
  data_random <- data[, c("use_score", covar), with = F] %>% .[complete.cases(.)]
  set.seed(123)
  nsg.rf = randomForest(y = data_random$use_score, x = data_random[, ..covar],
                        mtry=5, proximity = TRUE, importance=TRUE, ntree=1000)
  plot_randomforest <- varImpPlot(nsg.rf, sort=T, n.var= 13, pch=16)
  
  return(list(nsg.rf = nsg.rf, plot_randomforest = plot_randomforest))
})

names(forest_result) <- c("pn_ld", "pm_ld")

###plot results of random forest--------------
load("forest_13var_seed123_250114.RData")
data1 <- forest_result$pm_ld$plot_randomforest %>% as.data.frame(.) %>%
  mutate(., name = rownames(.)) %>%
  mutate(., name = ifelse(name == "age_stage7_mya", "Gene age", 
                          ifelse(name == "ifhsd2", "Gene duplication", name))) %>%
  arrange(`%IncMSE`) %>%   
  mutate(name=factor(name, levels=name))

forest_plot_list <- list(
  forest_pm = forest_result$pm_ld$plot_randomforest %>% as.data.frame(.) %>%
    mutate(., name = rownames(.)) %>%
    mutate(., name = ifelse(name == "age_stage7_mya", "Gene age",
                            ifelse(name == "ifhsd2", "Gene duplication", name))) %>%
    arrange(`%IncMSE`) %>%   
    mutate(name=factor(name, levels=name)) %>%
    ggplot(aes(x=name, y=`%IncMSE`)) +
    geom_segment(aes(xend=name, yend=0)) +
    geom_point(size=4, color="#E5C494") +
    coord_flip() +
    theme_bw() +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          legend.position = "none",
          axis.title.y = element_text(size = 15, face = "bold"), 
          axis.title.x = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +  xlab("") + ylab("GPS-M"),
  forest_pn = forest_result$pn_ld$plot_randomforest %>% as.data.frame(.) %>%
    mutate(., name = rownames(.)) %>%
    mutate(., name = ifelse(name == "age_stage7_mya", "Gene age", 
                            ifelse(name == "ifhsd2", "Gene duplication", name))) %>%
    mutate(name=factor(name, levels=levels(data1$name))) %>%
    ggplot(aes(x=name, y=`%IncMSE`)) +
    geom_segment(aes(xend=name, yend=0)) +
    geom_point(size=4, color="#E5C494") +
    coord_flip() +
    theme_bw() +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          legend.position = "none",
          axis.title.y = element_text(size = 15, face = "bold"), 
          axis.title.x = element_text(size = 15, face = "bold"), 
          axis.text = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +  xlab("") + ylab("GPS-N"))

forest_plot_list$forest_pm
forest_plot_list$forest_pn
ggsave(plot = forest_plot_list$forest_pm, width = 5, height = 5, device = cairo_pdf,
       filename = sprintf("%s/forest_pm.pdf",figure_file))
ggsave(plot = forest_plot_list$forest_pn, width = 5, height = 5, device = cairo_pdf,
       filename = sprintf("%s/forest_pn.pdf",figure_file))

##figure3L & supplementary figure 12A-------------------
###load data--------------
mediation_data <- pleiotropy_maindata$pn_ld[, .(`Gene length (bp)`, `CDS Length`, `CDS/Transcript Length ratio`, `GC content`, 
                                                `Transcript count`, `Exon Counts`, `Number of SNPs (Gene)`,
                                                LOEUF, pLI, `Missense/Synonymous ratio`, `Nonsense/Synonymous ratio`, 
                                                ifhsd2, use_score, age_stage7_mya)] 
mediation_data <- as.data.frame(mediation_data) %>%
  setnames(., names(.), c("gene_length", "cds_length", "cds_ratio", "gc", "transcript", "exon",
                          "snps", "LOEUF", "pli", "missense", "nonsense", "ifhsd2", "use_score", "age_stage7_mya")) %>%
  na.omit()
mediation_data$ifhsd2 <- as.numeric(mediation_data$ifhsd2) 

###singel mediator analysis using mediation (long time, could skip and load the results provided)--------------------
mediation_list <- list()
set.seed(123)
for(var_name in c("gene_length", "cds_length", "transcript", "exon", "cds_ratio", "gc",
                  "snps", "LOEUF", "pli", "missense", "nonsense", "ifhsd2")){
  use_data <- mediation_data 
  use_data$var_med <- use_data[, var_name]
  indirect.model1 <- glm(formula = var_med ~ age_stage7_mya, data = use_data)
  indirect.model2 <- glm(formula = use_score ~ age_stage7_mya + var_med, data = use_data)
  mediation.results <- mediate(indirect.model1, indirect.model2, 
                               treat = 'age_stage7_mya', mediator = "var_med", boot = TRUE, sims = 1000)
  test <- summary(mediation.results)
  mediation_list[[var_name]] <- test
  
}

##forest plot for supplementary figure 12A------------------
load("mediation_list_pn.RData")
library(forestploter)
forest_data <- data.table(
  Variable = c("Gene length", "CDS Length", "Transcript count", "Exon count", "CDS/transcript length ratio",  "GC content", 
               "Number of SNPs per gene", "LOEUF score", "pLI index", "Missense/synonymous ratio", 
               "Nonsense/synonymous ratio", "Gene duplication"),
  estimate =  sprintf("%.2f", sapply(as.list(1:12), function(x){mediation_list[[x]]$n.avg})*100) %>% as.numeric(),
  lower = sprintf("%.2f", sapply(as.list(1:12), function(x){mediation_list[[x]]$n.avg.ci[1]})*100) %>% as.numeric(),
  upper = sprintf("%.2f", sapply(as.list(1:12), function(x){mediation_list[[x]]$n.avg.ci[2]})*100) %>% as.numeric(),
  `P value` = ifelse(sprintf("%.3f", sapply(as.list(1:12), function(x){mediation_list[[x]]$n.avg.p})) == "0.000", "< 0.001", 
                     sprintf("%.3f", sapply(as.list(1:12), function(x){mediation_list[[x]]$n.avg.p})))) %>%
  .[, `Prop.med (95% CI)` := paste0(estimate, " (", lower, " to ", upper, ")")] %>%
  .[, ` ` := paste(rep(" ", 20), collapse = " ")] %>% .[order(-estimate)]

plot_singegmed <- forest(forest_data[, c(1,5,7,6)],
                         est = forest_data$estimate,
                         lower = forest_data$lower, 
                         upper = forest_data$upper,
                         ci_column = 3,
                         ref_line = 0,
                         xlim = c(-15,40),
                         ticks_at = c(-10, 0, 10, 20, 30, 40))
plot_singegmed

###multiple mediator analysis using sem (long time, could skip and load the results provided)---------------
library(sem); library(lavaan);library(bruceR)
model <- '
  use_score ~ b1*ifhsd2 + b2*gene_length + b3*cds_length + b5*transcript + b6*exon + b7*snps + b9*LOEUF + b10*pli + b11*missense + b12*nonsense + c*age_stage7_mya
  ifhsd2 ~ a1*age_stage7_mya
  gene_length ~ a2*age_stage7_mya
  cds_length ~ a3*age_stage7_mya
  transcript ~ a5*age_stage7_mya
  exon ~ a6*age_stage7_mya
  snps ~ a7*age_stage7_mya
  LOEUF ~ a9*age_stage7_mya
  pli ~  a10*age_stage7_mya
  missense ~ a11*age_stage7_mya
  nonsense ~ a12*age_stage7_mya
  indirect := a1*b1 + a2*b2 + a3*b3 + a5*b5 + a6*b6 + a7*b7 + a9*b9 + a10*b10 + a11*b11 + a12*b12
  ind_hsd := a1*b1 
  ind_genelength := a2*b2 
  ind_cds := a3*b3 
  ind_transcript := a5*b5
  ind_exon := a6*b6
  ind_snps := a7*b7
  ind_loeuf := a9*b9
  ind_pli := a10*b10
  ind_missense := a11*b11
  ind_nonsense := a12*b12
  total := a1*b1 + a2*b2 + a3*b3 + a5*b5 + a6*b6 + a7*b7 + a9*b9 + a10*b10 + a11*b11 + a12*b12 + c
'

fit <- lavaan::sem(model, data = mediation_data, se="bootstrap", bootstrap = 1000, iseed = 123) 
###plot the results for figure 3L-------------------------
load("lavaan_multiple_mediation_singlesig_pn.RData")
summary_fit <- lavaan_summary(fit, digits = 5, print = F)
summary_fit$effect$fdr <- p.adjust(summary_fit$effect$pval, method = "fdr")
summary_fit$effect
library(tidyverse); library(patchwork)
forest_data <- summary_fit$effect %>% 
  mutate(., mediator = rownames(.)) %>%
  as.data.table(.) %>% .[fdr < 0.05 & Estimate > 0 & mediator != "  total",] %>% .[order(Estimate)] %>%
  .[, `:=`(prop = Estimate/4.218134e-03*100,
           label = c("Nonsense/Synonymous ratio", "CDS Length", "Missense/Synonymous ratio",
                     "Transcript count","LOEUF score", "Gene duplication", "Total indirect effects"))]
forest_data

med_plot <-  ggplot(forest_data, aes(x = factor(label, levels = label), y = prop)) + theme_bw() +
  geom_bar(stat = "identity", position = "identity", color = "grey50", fill = "#85C6E0", width = 0.5) +
  coord_flip() +
  geom_text(aes(label = paste0(label, " (", sprintf("%.2f", prop),"%)"), 
                y = 0, x = c(1:7+0.5)), fontface = "bold", size = 4, hjust = 0) +
  labs(y = "Mediation effect (%)", title = "", x = "") + 
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(family = "Lobster Two", size = 12, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold"), 
        axis.text.y = element_blank(), axis.text.x = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line.x = element_line(colour = "grey50"),
        axis.line.y = element_blank(),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.border = element_blank()) 
med_plot
ggsave(med_plot, width = 4, height = 5, device = cairo_pdf,
       filename = sprintf("%s/mediation_prop_sigsingel_pn.pdf",figure_file))
