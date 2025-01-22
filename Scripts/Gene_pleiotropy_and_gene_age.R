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
#3 figure2------------------
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
  labs(y = "Gene pleiotropic scores for magnitude") + 
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
       filename = sprintf("%s/pm_ld_violin.pdf",figure_file))


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
  labs(y = "Gene pleiotropic score for number") + 
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
       filename = sprintf("%s/pn_ld_violin.pdf",figure_file))

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
                                                 legend_name = c("Low GPS-M", "High GPS-M"),
                                                 color_name = c("#8ED6BE", "#F7B29A"))
pleioage_mirrorbar_pn <- pleioage_mirrorbar_func(data = pleiotropy_maindata$pn_ld,
                                                 ylabel = "Proportion of pleiotropic genes (%)",
                                                 legend_name = c("Low GPS-N", "High GPS-N"),
                                                 color_name = c("#B1C0DD", "#F1B5D9"))

pleioage_mirrorbar_pm$plot
pleioage_mirrorbar_pn$plot

ggsave(plot = pleioage_mirrorbar_pm$plot, width = 8, height = 6, device = cairo_pdf,
       filename = sprintf("%s/pleioage_mirrorbar_pm.pdf",figure_file))
ggsave(plot = pleioage_mirrorbar_pn$plot, width = 8, height = 6, device = cairo_pdf,
       filename = sprintf("%s/pleioage_mirrorbar_pn.pdf",figure_file))


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
