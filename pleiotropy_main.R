#1 load packages and set variables-------
sapply(c("data.table", "dplyr", "ggplot2", "corrplot", "ggpubr", "Gmisc", "openxlsx", "readxl", "paletteer", "MuMIn",
         "ggsci", "scales", "RColorBrewer", "gridExtra", "dplyr", "tidyr", "stringr", "colorspace", "cowplot",
         "ggbreak", "ggstatsplot", "viridis"), require, character.only = TRUE)
setwd("/home/liumy/pleiotropy/github/Figure")
mycolor <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494")
pleio_color <- lighten(c("#80CDC1", "#85C6E0", "#223E84"), 0.45)
age_color <- lighten(c("#AB6191", "#DA7235", "#7FB65D", "#F2DD33"), 0.45)
age_color7 <- lighten(mycolor, 0.45)
age7_name <- c("Euteleostomi", "Tetrapoda", "Amniota", "Mammalia", "Theria", "Eutheria", "Primate")
x_text <-  c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria")
age_test <- c("Euteleostomi", "Eutheria")
#2 load data--------------------------
load("/home/liumy/pleiotropy/github/data/pleiotropy_maindata.RData")
#3 figure 1-----------------------
##figure1B--------------
median(pleiotropy_maindata$pn_ld$use_score) %>% sprintf("%.2f", .) #9.38
median(pleiotropy_maindata$pm_ld$use_score) %>% sprintf("%.2f", .) #6.16

hist_pn <- ggplot(data = pleiotropy_maindata$pn_ld, aes(x = use_score)) +
  geom_histogram(fill = lighten("#148ABB", 0.65, space = "HLS"), binwidth = 1) + theme_cowplot() +
  labs(y = "Gene counts\n",x = "GPS-N") +
  scale_x_break(breaks = c(20,30), scales = 0.2) +
  scale_x_continuous(limits = c(0,30), breaks = c(0,5,10,15,20,30), labels = c(0,5,10,15,20,30)) +
  geom_vline(xintercept=median(pleiotropy_maindata$pn_ld$use_score), 
             color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.8) +
  theme(axis.title = element_text(size = 15, face = "bold"), plot.title = element_text(hjust = 0.5,face = "bold"),
        axis.text = element_text(size=12,face="bold"),
        panel.grid = element_blank(),
        axis.text.x.top = element_blank(), axis.line.x.top = element_blank(), axis.ticks.x.top = element_blank(),
        legend.position = "none", legend.text = element_text(size = 12,face = "bold"),axis.ticks.length=unit(-0.1, "cm"),
        legend.title = element_blank(),legend.background = element_blank(),legend.key.size = unit(0.5,"cm"))  

hist_pm <- ggplot(data = pleiotropy_maindata$pm_ld, aes(x = use_score)) +
  geom_histogram(fill = lighten("#F5C51E", 0.65, space = "HLS"), binwidth = 1) + theme_cowplot() +
  theme(axis.title = element_text(size = 15, face = "bold"), plot.title = element_text(hjust = 0.5,face = "bold"),
        axis.text = element_text(size=12,face="bold"),
        panel.grid = element_blank(),
        axis.text.x.top = element_blank(), axis.line.x.top = element_blank(), axis.ticks.x.top = element_blank(),
        legend.position = "none", legend.text = element_text(size = 12,face = "bold"),axis.ticks.length=unit(-0.1, "cm"),
        legend.title = element_blank(),legend.background = element_blank(),legend.key.size = unit(0.5,"cm")) +
  labs(y = "Gene counts\n",x = "GPS-M") +
  scale_x_break(breaks = c(20,30), scales = 0.2) +
  scale_x_continuous(limits = c(0,30), breaks = c(0,5,10,15,20,30), labels = c(0,5,10,15,20,30)) +
  geom_vline(xintercept=median(pleiotropy_maindata$pm_ld$use_score), color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.8) 


hist_pn
hist_pm

ggsave(plot = hist_pn, width = 5, height = 2.5, device = cairo_pdf ,
       filename = "score_distribution_histograms_pn.pdf") 
ggsave(plot = hist_pm, width = 5, height = 2.5, device = cairo_pdf ,
       filename = "score_distribution_histograms_pm.pdf") 

##figure1C--------------
score_data <- merge(pleiotropy_maindata$pm_ld[, .(gene, use_score, pleio_class3)], 
                    pleiotropy_maindata$pn_ld[, .(gene, use_score, pleio_class3)], by = "gene") %>%
  setnames(., names(.), c("gene", "pm", "pm_class3", "pn", "pn_class3")) %>%
  .[, pleio_class := ifelse(pm_class3 == pn_class3, pm_class3, "4")]

cor(score_data[, .(pm, pn)], method = "sp")[1,2] %>% sprintf("%.2f",.)
cor.test(score_data$pm, score_data$pn, method = "sp")$p.value #p-value < 2.2e-16

sp <- ggplot(data = score_data, aes(x = pm, y = pn)) + theme_bw() +
  geom_point(aes(color = pleio_class), size = 3, alpha = 0.6) +
  scale_color_manual(values = c(mycolor[1:3], "#868686ff"), 
                     labels = c("Low pleiotropy", "Intermediate pleiotropy", "High pleiotropy", "Incosistent classification")) +
  theme(axis.title = element_text(size = 15, face = "bold"), plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(), legend.position = "none",
        legend.text = element_text(size = 12, face = "bold"),
        axis.ticks.length=unit(-0.1, "cm"),
        legend.title = element_blank(),legend.background = element_blank(),legend.key.size = unit(0.5,"cm"),
  ) +
  labs(y = "GPS-N",x = "GPS-M") + 
  geom_smooth(method = "lm", formula = y ~ x, color = "black", alpha = 0.5, linetype = "dashed") +
  ylim(c(0,40)) + xlim(c(0,40))

yplot <- ggboxplot(score_data, x = "pn_class3", y = "pn", alpha = 0.5, fill = "pn_class3", color = "pn_class3", ggtheme = theme_bw())+
  scale_color_manual(values = mycolor[1:3]) +
  scale_fill_manual(values = mycolor[1:3]) +
  theme(panel.grid = element_blank()) +
  ylim(c(0,40)) +
  clean_theme() + rremove("legend") 

xplot <- ggboxplot(score_data, x = "pm_class3", y = "pm", alpha = 0.5, fill = "pm_class3", color = "pm_class3", ggtheme = theme_bw())+
  scale_color_manual(values = mycolor[1:3]) +
  scale_fill_manual(values = mycolor[1:3]) +
  theme(panel.grid = element_blank()) + 
  ylim(c(0,40)) +
  rotate() + clean_theme() +rremove("legend")

sp
xplot
yplot

ggsave(plot = sp, width = 5, height = 5, device = cairo_pdf ,
       filename = "score_distribution_sp.pdf") 
ggsave(plot = yplot, width = 2, height = 5, device = cairo_pdf ,
       filename = "score_distribution_yplot.pdf") 
ggsave(plot = xplot, width = 5, height = 2, device = cairo_pdf ,
       filename = "score_distribution_xplot.pdf") 

##figure1D--------------
mapping <- c("1" = 249, "2" = 242, "3" = 198, "4" = 190, "5" = 182, "6" = 171, "7" = 159, "8" = 145,
             "9" = 138, "10" = 134, "11" = 135, "12" = 133, "13" = 114, "14" = 107, "15" = 102, "16" = 90,
             "17" = 83, "18" = 80, "19" = 59, "20" = 64, "21" = 47, "22" = 51)
chrom_data <- merge(score_data, pleiotropy_maindata$pm_ld[, .(gene, chr)], by = "gene", all.x = T) 
chrom_table <- lapply(list("pm_class3", "pn_class3"), function(x){
  data <- chrom_data[, use_class3 := get(x)]
  data_all <- data[, .N, by = chr]
  data_table <- data[use_class3 == "High pleiotropy", .N, by = chr] %>%
    merge(., data_all, by = "chr") %>%
    .[, chr_length := mapping[chr]] %>%
    setnames(., names(.), c("chr", "high", "total", "chr_length")) %>%
    .[, percent := high/total * 100] %>% .[order(percent)] %>%
    .[, rank := c(1:22)]

  return(data_table)
})

names(chrom_table) <- c("pm_ld", "pn_ld")
head(chrom_table$pm_ld)

chrpercent_func <- function(score_name, x_label = "Proportion of H-GPS-N genes\nper chromosome (%)"){
  data <- chrom_table[[score_name]]
  percent_data <- data.table(rank = as.factor(data$rank), percent = data$percent, 
                             label = sprintf("%.2f", data$percent),
                             chr_name = paste0("Chr", data$chr))
  chr_label <- percent_data[order(rank)]$chr_name
  chr_percent <- ggplot(data = percent_data) + theme_bw() +
    geom_col(aes(percent, rank), fill = lighten("#9E6531",0.6), width = 0.75) +
    scale_x_continuous(limits = c(0, 35), breaks = seq(0, 35, by = 5), expand = c(0, 0), position = "bottom" ) +
    scale_y_discrete(expand = expansion(add = c(0, 0.5)),
                     labels = chr_label) +
    theme(panel.grid = element_blank(), axis.ticks.length = unit(0, "mm"), 
          axis.title.x = element_text(face = "bold", size = 15),
          axis.line = element_line(color = "grey20"), 
          axis.text = element_text(face = "bold", size = 12), 
          panel.border = element_blank()) +
    labs(x = x_label, y = "", title = "") +
    geom_vline(xintercept = 20, color = "grey30", linetype = "dashed", linewidth = 1)
  
  return(list(data = percent_data, plot = chr_percent))
}

chr_percent <- list(
  pn_ld = chrpercent_func(score_name = "pn_ld", x_label = "Proportion of H-GPS-N genes\nper chromosome (%)"),
  pm_ld = chrpercent_func(score_name = "pm_ld", x_label = "Proportion of H-GPS-M genes\nper chromosome (%)")
)

chr_percent$pn_ld$plot
chr_percent$pm_ld$plot

ggsave(plot = chr_percent$pn_ld$plot, width = 5, height = 8, device = cairo_pdf, 
       filename = "chr_percent_pn.pdf")

#4 figure2------------------
##figure2A&2B--------------------
###violin plots of pleiotropic scores among age stages
cor(pleiotropy_maindata$pm_ld[, .(age_stage7_mya, use_score)], method = "sp", use = "pairwise.complete.obs")[2,1] %>% sprintf("%.3f", .)
cor.mtest(pleiotropy_maindata$pm_ld[, .(age_stage7_mya, use_score)], method = "sp")$p[2,1] %>% sprintf("%.2e", .)

cor(pleiotropy_maindata$pn_ld[, .(age_stage7_mya, use_score)], method = "sp", use = "pairwise.complete.obs")[2,1] %>% sprintf("%.3f", .)
cor.mtest(pleiotropy_maindata$pn_ld[, .(age_stage7_mya, use_score)], method = "sp")$p[2,1] %>% sprintf("%.2e", .)

###pm_ld
violin_data <- pleiotropy_maindata$pm_ld
upper_plot <- ggbetweenstats(data = violin_data, centrality.type = "nonparametric", centrality.plotting = F,
                             x = age_stage7, y = use_score, 
                             pairwise.display = "none", bf.message = F, results.subtitle = F) +
  scale_color_manual(values = alpha(mycolor[1:7], 0.01)) +
  labs(y = " ", title = " ") + 
  theme(plot.title = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 15), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line.y = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  coord_cartesian(ylim = c(10, 35)) +
  scale_y_continuous(breaks = c(10, 30))


main_plot <- ggbetweenstats(data = violin_data, centrality.type = "nonparametric",
                            x = age_stage7, y = use_score,
                            pairwise.display = "none", bf.message = F, results.subtitle = F,
                            centrality.point.args = list(size = 4, color = "grey30"),
                            centrality.label.args = list(size = 3, 
                                                         face = "bold",
                                                         nudge_x = 0.1, 
                                                         segment.linetype = 4, 
                                                         min.segment.length = 0)) +
  scale_color_manual(values = alpha(mycolor[1:7], 0.01)) +
  labs(y = "GPS-M") + 
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line.y = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
  coord_cartesian(ylim = c(5, 7.2)) +
  scale_y_continuous(breaks = c(5, 6, 7))


lower_plot <- ggbetweenstats(data = violin_data, centrality.type = "nonparametric", centrality.plotting = F,
                             x = age_stage7, y = use_score, 
                             pairwise.display = "none", bf.message = F, results.subtitle = F) +
  scale_color_manual(values = alpha(mycolor[1:7], 0.01)) +
  labs(x = "", y = "") + 
  theme(axis.text = element_text(size = 12, color = "black", face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + 
  coord_cartesian(ylim = c(1.5, 5)) +
  scale_y_continuous(breaks = c(4)) +
  scale_x_discrete(labels = c("Euteleostomi\n(429 mya)", "Tetrapoda\n(352 mya)", "Amniota\n(319 mya)", 
                              "Mammalia\n(180 mya)", "Theria\n(160 mya)", "Eutheria\n(99 mya)", "Primate\n(66 mya)")) 

combined_plot1 <- plot_grid(upper_plot, main_plot, lower_plot, ncol = 1, align = "v", rel_heights = c(0.27,1))

combined_plot1 

ggsave(plot = combined_plot1, width = 8, height = 7, device = cairo_pdf,
       filename = "pm_ld_violin.pdf")


###pn_ld
violin_data <- pleiotropy_maindata$pn_ld
upper_plot <- ggbetweenstats(data = violin_data, centrality.type = "nonparametric", centrality.plotting = F,
                             x = age_stage7, y = use_score, 
                             pairwise.display = "none", bf.message = F, results.subtitle = F) +
  scale_color_manual(values = alpha(mycolor[1:7], 0.01)) +
  labs(y = " ", title = " ") + 
  theme(plot.title = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"), 
        axis.title.y = element_text(size = 15), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line.y = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  coord_cartesian(ylim = c(15, 30)) +
  scale_y_continuous(breaks = c(20, 30))


main_plot <- ggbetweenstats(data = violin_data, centrality.type = "nonparametric",
                            x = age_stage7, y = use_score, 
                            pairwise.display = "none", bf.message = F, results.subtitle = F,
                            centrality.point.args = list(size = 5, color = "grey30"),
                            centrality.label.args = list(size = 3, 
                                                         face = "bold",
                                                         nudge_x = 0.1, 
                                                         segment.linetype = 4, 
                                                         min.segment.length = 0)) +
  scale_color_manual(values = alpha(mycolor[1:7], 0.01)) +
  labs(y = "GPS-N") + 
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line.y = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(),
        #panel.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4"),
        #plot.background = element_rect(fill = "#fbf9f4", color = "#fbf9f4") 
  ) + coord_cartesian(ylim = c(5, 12)) 


lower_plot <- ggbetweenstats(data = violin_data, centrality.type = "nonparametric", centrality.plotting = F,
                             x = age_stage7, y = use_score, 
                             pairwise.display = "none", bf.message = F, results.subtitle = F) +
  scale_color_manual(values = alpha(mycolor[1:7], 0.01)) +
  labs(x = "", y = "") + 
  theme(axis.text = element_text(size = 12, color = "black", face = "bold"), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + 
  coord_cartesian(ylim = c(0, 3)) +
  scale_y_continuous(breaks = c(2)) +
  scale_x_discrete(labels = c("Euteleostomi\n(429 mya)", "Tetrapoda\n(352 mya)", "Amniota\n(319 mya)", 
                              "Mammalia\n(180 mya)", "Theria\n(160 mya)", "Eutheria\n(99 mya)", "Primate\n(66 mya)"))

combined_plot2 <- plot_grid(upper_plot, main_plot, lower_plot, ncol = 1, align = "v", rel_heights = c(0.3,1))

combined_plot2 

ggsave(plot = combined_plot2, width = 8, height = 7, device = cairo_pdf,
       filename = "pn_ld_violin.pdf")

##figure2C&2D--------------------
pleioage_mirrorbar_func <- function(data = pleiotropy_maindata$pm_ld, 
                                    ylabel = "Percents of pleiotropic genes (%)",
                                    legend_name = c("L-GPS", "H-GPS"),
                                    color_name = c("#FF9874", "#6FCAAD")){
  
  df_percent <- table(data$age_stage7_num, data$pleio_class3) %>% 
    prop.table(., margin = 1) %>% as.data.table(.) %>% 
    setnames(., names(.), c("age_stage", "group", "percent")) %>%
    .[group != "Intermediate pleiotropy",] %>%
    .[, age_stage := as.factor(age_stage)]
  
  df_percent$percent <- ifelse(df_percent$group %in% "High pleiotropy", -df_percent$percent*100, df_percent$percent*100)
  
  plot <- ggplot(df_percent) + theme_bw() +
    geom_bar(aes(x = age_stage, y = percent, fill = group, color = group),
             stat = "identity", position = "identity", width = 0.7) +
    scale_fill_manual(values = c("Low pleiotropy" = color_name[1], "High pleiotropy" = color_name[2]),
                      labels = c("Low pleiotropy" = legend_name[1], "High pleiotropy" = legend_name[2])) + 
    scale_color_manual(values = c("Low pleiotropy" = color_name[1], "High pleiotropy" = color_name[2]),
                       labels = c("Low pleiotropy" = legend_name[1], "High pleiotropy" = legend_name[2])) + 
    coord_flip() +
    geom_text(aes(x = age_stage, y = ifelse(percent<0,-5.5 , 5.5), label = paste0(sprintf("%.2f", abs(percent)), "%")), 
              color = "grey40", size=5, fontface = "bold") +
    labs(y = ylabel, title = " ", x = "") + 
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          axis.title.y = element_text(size = 15, face = "bold"), 
          axis.title.x = element_text(size = 15, face = "bold"), 
          axis.text.y = element_text(size = 12, color = "black", face = "bold"), 
          axis.ticks = element_blank(), axis.line.x = element_line(colour = "grey50"),
          axis.line.y = element_blank(),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), legend.title = element_blank(),
          panel.border = element_blank(),
          legend.position = "top") +
    scale_x_discrete(limits = levels(df_percent$age_stage), 
                     labels = c("Euteleostomi\n(429 mya)", "Tetrapoda\n(352 mya)", "Amniota\n(319 mya)", 
                                "Mammalia\n(180 mya)", "Theria\n(160 mya)", "Eutheria\n(99 mya)", "Primate\n(66 mya)"))
  
  return(list(data = df_percent, plot = plot))
}

pleioage_mirrorbar_pm <- pleioage_mirrorbar_func(data = pleiotropy_maindata$pm_ld,
                                                 ylabel = "Proportion of pleiotropic genes (%)",
                                                 legend_name = c("L-GPS-M", "H-GPS-M"),
                                                 color_name = c("#8ED6BE", "#F7B29A"))
pleioage_mirrorbar_pn <- pleioage_mirrorbar_func(data = pleiotropy_maindata$pn_ld,
                                                 ylabel = "Proportion of pleiotropic genes (%)",
                                                 legend_name = c("L-GPS-N", "H-GPS-N"),
                                                 color_name = c("#B1C0DD", "#F1B5D9"))

pleioage_mirrorbar_pm$plot
pleioage_mirrorbar_pn$plot

ggsave(plot = pleioage_mirrorbar_pm$plot, width = 8, height = 6, device = cairo_pdf,
       filename = "pleioage_mirrorbar_pm.pdf")
ggsave(plot = pleioage_mirrorbar_pn$plot, width = 8, height = 6, device = cairo_pdf,
       filename = "pleioage_mirrorbar_pn.pdf")


###p for trend
library(rstatix)
xtab <- as.table(rbind(
  c(table(pleiotropy_maindata$pn_ld[, .(age_stage7, pleio_class3)])[,"Low pleiotropy"]),
  c(table(pleiotropy_maindata$pn_ld[, .(age_stage7)]))
))
prop_trend_test(xtab)$p %>% sprintf("%.2e", .)

xtab <- as.table(rbind(
  c(table(pleiotropy_maindata$pn_ld[, .(age_stage7, pleio_class3)])[,"High pleiotropy"]),
  c(table(pleiotropy_maindata$pn_ld[, .(age_stage7)]))
))
prop_trend_test(xtab)$p %>% sprintf("%.2e", .)

xtab <- as.table(rbind(
  c(table(pleiotropy_maindata$pm_ld[, .(age_stage7, pleio_class3)])[,"Low pleiotropy"]),
  c(table(pleiotropy_maindata$pm_ld[, .(age_stage7)]))
))
prop_trend_test(xtab)$p %>% sprintf("%.2e", .)

xtab <- as.table(rbind(
  c(table(pleiotropy_maindata$pm_ld[, .(age_stage7, pleio_class3)])[,"High pleiotropy"]),
  c(table(pleiotropy_maindata$pm_ld[, .(age_stage7)]))
))
prop_trend_test(xtab)$p %>% sprintf("%.2e", .)

#5 figure3----------------------
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
                                                         y_label = "Gene length (bp)", x_var = "Gene length (bp)"), 
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
         filename = "gene_length_pn.pdf")

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
         filename = "dnds_mouse_pn.pdf")

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
         filename = "lof_tolerant_pn.pdf")

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
                                                          y_label = "Gene length (bp)", x_var = "Gene length (bp)", 
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
         filename = "gene_length_age.pdf")


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
         filename = "dnds_mouse_age.pdf")


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
         filename = "lof_tolerant_age.pdf")

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
                                        y_label = "GPS-N", x_var = "use_score", group = "ifhsd2",
                                        x_text = c("Singletons", "Duplicates"),
                                        stat_name = "pleio_pn", test_list = c("Singletons", "Duplicates")),
                pleio_pm = hsd_box_func(data = pleiotropy_maindata$pm_ld[, ifhsd2 := factor(ifhsd2, levels = c("Singletons", "Duplicates"))], 
                                        color_value = mycolor, scale_y = T,
                                        y_label = "GPS-M", x_var = "use_score", group = "ifhsd2",
                                        x_text = c("Singletons", "Duplicates"),
                                        stat_name = "pleio_pm", test_list = c("Singletons", "Duplicates"))) 


hsd_box$pleio_pn$plot

ggsave(plot = hsd_box$pleio_pn$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = "hsd_pleio_pn.pdf")

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
                                                          y_label = "GPS-N", x_var = "use_score", group = "progenitor",
                                                          x_text = c("Ancestral\nparalogs", "Recent\nparalogs"),
                                                          stat_name = "progenitor_pn", test_list = c("Progenitor", "Copy")),
                           pleio_pm = progenitor_box_func(data = progenitor_data$pm, color_value = mycolor[c(3,4)], scale_y = T,
                                                          y_label = "GPS-M", x_var = "use_score", group = "progenitor",
                                                          x_text = c("Ancestral\nparalogs", "Recent\nparalogs"),
                                                          stat_name = "progenitor_pm", test_list = c("Progenitor", "Copy"))) 
hsd_progenitor_box$pleio_pn$plot

ggsave(plot = hsd_progenitor_box$pleio_pn$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = "progenitor_pleio_pn.pdf")

##figure3K-----------------------------
###run the random forest
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


###plot results of random forest
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
       filename = "forest_pm.pdf")
ggsave(plot = forest_plot_list$forest_pn, width = 5, height = 5, device = cairo_pdf,
       filename = "forest_pn.pdf")

##figure3L-------------------
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

###singel mediator analysis using mediation--------------------
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

##forest plot
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

ggsave(plot = plot_singegmed, width = 7, height = 5, device = cairo_pdf,
       filename = "singegmed_forest_pn.pdf") 

###multiple mediator analysis using sem---------------
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
save(fit, file = "multiple_mediation_singlesig_pn.RData")

###plot the results
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
       filename = "mediation_prop_sigsingel_pn.pdf")

#6 figure4-----------------------
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
####load data-----------------
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
save(prop_data, file = "cis_esqtl.RData")

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
       filename = "cisegenes_pn.pdf")
ggsave(plot = percent_egenes$age, width = 4.8, height = 5, device = cairo_pdf,
       filename = "cisegenes_age.pdf")

ggsave(plot = percent_sgenes$pleio, width = 3.2, height = 5, device = cairo_pdf,
       filename = "cissgenes_pn.pdf")
ggsave(plot = percent_sgenes$age, width = 4.8, height = 5, device = cairo_pdf,
       filename = "cissgenes_age.pdf")

###trans-eGenes----------------------------------
trans_eqtl <- fread("/home/liumy/pleiotropy/data/qtl/eqtlgen_transeqtl.txt")
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
       filename = "transegenes_pn.pdf")
ggsave(plot = percent_transegenes$age, width = 4.8, height = 5, device = cairo_pdf,
       filename = "transegene_age.pdf")


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
tau_pleio_m
tau_age

ggsave(plot = tau_pleio_n, width = 5, height = 5, device = cairo_pdf,
       filename = "tau_ridge_pn.pdf")
ggsave(plot = tau_age, width = 5, height = 5, device = cairo_pdf,
       filename = "tau_ridge_age.pdf")

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
       filename = "tissue_express_fto.pdf")
ggsave(plot = plot_hist[[2]], width = 1, height = 5, device = cairo_pdf,
       filename = "tissue_express_slc12a3.pdf")

##figure4E&4F----------------------
###load data-----------------
library(dorothea); library(OmnipathR); library(decoupleR)
tf <- fread("/home/liumy/pleiotropy/data/TF/TFC2_16102023b.tsv") %>%
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
tf_sum_box$pleio_m$plot
tf_sum_box$age$plot


ggsave(plot = tf_sum_box$pleio_n$plot, width = 3.2, height = 5, device = cairo_pdf,
       filename = "tf_sum_pn.pdf")
ggsave(plot = tf_sum_box$age$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = "tf_sum_age.pdf")

##figure4G&4H----------------------
###load data----------------------
enhancer_annota <- fread("/home/liumy/pleiotropy/data/CRE/genehancer_annotation.csv") #244737
length(unique(enhancer_annota$name))
gene_enhancer_all <- fread("/home/liumy/pleiotropy/data/CRE/genehancer_all.csv") %>%
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
       filename = "enhancer_box_pn.pdf")
ggsave(plot = cre_box$age$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = "enhancer_box_age.pdf")


##figure4I&4J----------------------
###load data------------------
ogee_pip <- fread("/home/liumy/pleiotropy/data/PIP/connectivity.txt") %>%
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
       filename = "pip_pn.pdf")
ggsave(plot = pip_box$age_n$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = "pip_age.pdf")

#7 figure5--------------------
##enrichment by lola----------------- 
###Extract gene coordinates ± 1000bp-------------------
genecoordinate_func <- function(x, data = pleiotropy_maindata$pm_ld, var_name = "pleio100", out_name = "pm", ifidentical = F){
  
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
  
  fwrite(bed_data, file = sprintf("/home/liumy/pleiotropy/github/bed_file/%s.bed", out_name), 
         quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(bed_data)
}
####GRCh38---------------------
library(biomaRt); library(dbplyr)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

pm_bed <- list()
for(i in c(20)){
  top_bed_name = paste0("pm_top", i, "_hg38")
  pm_bed[[top_bed_name]] <- genecoordinate_func(x = 101-i, var_name =  "pleio100", 
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

bed_file <- "/home/liumy/pleiotropy/github/bed_data"
lola_out_file <- "/home/liumy/pleiotropy/github/lola"

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

###figrue5A&5B: extract enrichment results by lola for UCSC database------------------
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


ggsave(plot = ucsc_enrich_plot$pn, width = 7, height = 5, device = cairo_pdf,
       filename = "ucsc_enrich_pn20.pdf")
ggsave(plot = ucsc_enrich_plot$age, width = 7, height = 5, device = cairo_pdf,
       filename = "ucsc_enrich_age.pdf")

###figrue5C&5D: extract enrichment results by lola for Roadmap epigenomics region using beeswarm plot----------------------
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
    ylab("-log10(qvalue)") +
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
       filename = "roadmap_pn_top20.pdf")
ggsave(plot = roadmapepi_beeplot$pn_epi_bottom20, width = 3, height = 5, device = cairo_pdf,
       filename = "roadmap_pn_bottom20.pdf")
##3D genomic features---------------
###lola enrichment for 3d genomic features-------------------
regionDB = loadRegionDB("/home/liumy/pleiotropy/github/data/lola_3d")

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

####figure5E: plot heatmap for subcompartment--------------------------
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
                       breaks = c(6,12,18), labels = c("25%", "50%", "75%")) +
  geom_text(aes(label = marker), color = "black", size = 6) +
  facet_grid(facet_y~facet_x, space = "free", scales = "free") +
  labs(y = "", title = "", x = "") + 
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(size = 12, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold"), 
        axis.title.x = element_text(size = 15, face = "bold"), 
        axis.text.y = element_text(size = 12, face = "bold", color = c("#AB6191", "#AB6191", "#8DA0CB", "#8DA0CB")),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.ticks = element_blank(), axis.line = element_blank(),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill=NA, linewidth=0.5),
        legend.position = "top",
        strip.text = element_text(size=0), strip.background = element_blank()) 


plot_compartment
ggsave(plot = plot_compartment, width = 5, height = 5, device = cairo_pdf,
       filename = "enrich_compartment_20%%.pdf")



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
       filename = "enrich_loop_pn_heatmap.pdf")
ggsave(plot = loop_domian_atac_heatmap$pn_domian$plot, width = 4, height = 4, device = cairo_pdf,
       filename = "enrich_domain_pn_heatmap.pdf")
ggsave(plot = loop_domian_atac_heatmap$pn_atac$plot, width = 4, height = 4, device = cairo_pdf,
       filename = "enrich_atac_pn_heatmap.pdf") 

ggsave(plot = loop_domian_atac_heatmap$age4_loop$plot, width = 4, height = 4, device = cairo_pdf,
       filename = "enrich_loop_age_heatmap.pdf")
ggsave(plot = loop_domian_atac_heatmap$age4_domian$plot, width = 4, height = 4, device = cairo_pdf,
       filename = "enrich_domain_age_heatmap.pdf")
ggsave(plot = loop_domian_atac_heatmap$age4_atac$plot, width = 4, height = 4, device = cairo_pdf,
       filename = "enrich_atac_age_heatmap.pdf") 

ggsave(plot = loop_domian_atac_heatmap$pm_loop$plot, width = 4, height = 4, device = cairo_pdf,
       filename = "enrich_loop_pm_heatmap.pdf")
ggsave(plot = loop_domian_atac_heatmap$pm_domian$plot, width = 4, height = 4, device = cairo_pdf,
       filename = "enrich_domain_pm_heatmap.pdf")
ggsave(plot = loop_domian_atac_heatmap$pm_atac$plot, width = 4, height = 4, device = cairo_pdf,
       filename = "enrich_atac_pm_heatmap.pdf") 

#8 figure6--------------------------
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
  pdf(file = sprintf("bp_%s_enrich.pdf", x), width=8, height=8)
  treemapPlot(gobp_plot[[x]], overlap.labels = 0.99, force.print.labels = T, 
              lowerbound.cex.labels = 0.01, border.lwds = c(0.8,0.4))
  dev.off()
})

##figure6B&6C-------------------------
###load data-----------------
chronos_data <- fread("/home/liumy/pleiotropy/data/essentiality/CRISPRGeneEffect.csv")
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
       filename = "ges_pn.pdf")
ggsave(plot = GES_box$age_n$plot, width = 4.8, height = 5, device = cairo_pdf,
       filename = "ges_age.pdf")



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
       filename = "progeny_pn.pdf")

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
luca <- fread("/home/liumy/pleiotropy/data/drug/41559_2024_2461_MOESM4_ESM.tsv")
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

#figure7--------------------
##figure7A: Proportion of FDA approved drug targets across groups----------------------
###load data---------------
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
fda_durggenes <- read.xlsx("/home/liumy/pleiotropy/data/drug/fda_2023.xlsx") #drug genes from Nat. Rev. Drug Discov. 22, 864 (2023).
fda_target <- strsplit(fda_durggenes$targetIds, split = ",") %>% unlist(.) %>% unique(.)  #428 drugs with 506 targets

fda_durggenes_long <- fda_durggenes %>%  separate_rows(targetIds, sep = ",\\s*") %>% as.data.table(.) %>%
  .[, drug_gene := paste0(brandDrugName,"_", targetIds)] %>%
  .[!duplicated(drug_gene)] %>%
  .[!is.na(targetIds) & targetIds != "NA",]

drug_data <- lapply(list("pm_ld", "pn_ld"), function(x, data = pleiotropy_maindata){
  target_sum <- fda_durggenes_long[, .(target_pairs = .N), by = targetIds]
  data_merge <- data[[x]][, .(gene, ensemblid, pleio_class3, pleio10, pleio_class5, age_stage4)] %>%
    .[, iffda := ifelse(is.na(ensemblid), NA, ifelse(ensemblid %in% fda_target, 1, 0))] %>%
    merge(., target_sum, by.x = "ensemblid", by.y = "targetIds", all.x = T)
  return(data_merge)
})

names(drug_data) <- c("pm_ld", "pn_ld")

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
                                           y_label = "Proportion of genes\ntargeted by FDA-approved durgs (%)"),
                    pleio_m = drug_percent_plot(data = drug_prop$pleio_m_fda$N, 
                                           group = factor(c("First", "Second", "Third", "Fourth", "Fifth"), 
                                                          levels = c("First", "Second", "Third", "Fourth", "Fifth")),
                                           color_values = lighten(c("#6baed6","#4292c6","#2171b5","#08519c","#08306b"), 0.6),
                                           stat_test = test, test_list = c("First", "Fifth"), stat_name = "pleio_m",
                                           y_label = "Proportion of genes\ntargeted by FDA-approved durgs (%)"),
                    age = drug_percent_plot(data = drug_prop$age_n_fda$N, 
                                       group = factor(c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria"), 
                                                      levels = c("Euteleostomi", "Tetrapoda", "Amniota", "Eutheria")),
                                       color_values = lighten(c("#FFC3AF", "#FF935F", "#DA7235", "#944F2B"), 0.4),
                                       stat_test = test, test_list = c("Euteleostomi", "Eutheria"), stat_name = "age",
                                       y_label = "Proportion of genes\ntargeted by FDA-approved durgs (%)"))

percent_fda$pleio_n
ggsave(plot = percent_fda$pleio_n, width = 6, height = 5, device = cairo_pdf,
       filename = "percent_fda_pn.pdf")


###trend test-------------------------
library(rstatix)
xtab <- as.table(rbind(
  c(table(drug_data$pn_ld[iffda == 1,]$pleio_class5)),
  c(table(drug_data$pn_ld[iffda == 0,]$pleio_class5))
)) 
prop_trend_test(xtab)$p %>% sprintf("%.2e", .)

#figure7BC: Proportion of T–I pairs based on current development status--------------------
#wget -c https://github.com/ericminikel/genetic_support/archive/refs/heads/main.zip
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

setwd("/home/liumy/pleiotropy/genetic_support-main")
output_path = "/home/liumy/pleiotropy/github/Figure"

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
    filter(!(target_status %in% c('indication lacks genetic insight','no indication annotated','no target annotated'))) %>%
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
merge2_raw = read_tsv('data/merge2.tsv.gz', col_types=cols())
assoc = read_tsv('data/assoc.tsv.gz', col_types=cols())
indic = read_tsv("data/indic.tsv", col_types=cols())
meta_hcat = read_tsv('data/meta_hcat.tsv', col_types=cols())
meta_acat = read_tsv('data/meta_acat.tsv', col_types=cols())
meta_ccat = read_tsv('data/meta_ccat.tsv', col_types=cols())


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
       filename = "forest_combined_pn.pdf")


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
       filename = "plotline_pn.pdf")


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