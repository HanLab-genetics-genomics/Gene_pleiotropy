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
#3 figure 1-----------------------
##figure1B--------------
median(pleiotropy_maindata$pn_ld$use_score) %>% sprintf("%.2f", .) #9.38
median(pleiotropy_maindata$pm_ld$use_score) %>% sprintf("%.2f", .) #6.16

hist_pn <- ggplot(data = pleiotropy_maindata$pn_ld, aes(x = use_score)) +
  geom_histogram(fill = lighten("#148ABB", 0.65, space = "HLS"), binwidth = 1) + theme_cowplot() +
  labs(y = "Gene count\n",x = "Gene pleiotropic score for number") +
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
  labs(y = "Gene count\n",x = "Gene pleiotropic score for magnitude") +
  scale_x_break(breaks = c(20,30), scales = 0.2) +
  scale_x_continuous(limits = c(0,30), breaks = c(0,5,10,15,20,30), labels = c(0,5,10,15,20,30)) +
  geom_vline(xintercept=median(pleiotropy_maindata$pm_ld$use_score), color = "black", linetype = "dashed", linewidth = 0.5, alpha = 0.8) 


hist_pn
hist_pm

ggsave(plot = hist_pn, width = 5, height = 2.5, device = cairo_pdf ,
       filename = sprintf("%s/score_distribution_histograms_pn.pdf",figure_file)) 
ggsave(plot = hist_pm, width = 5, height = 2.5, device = cairo_pdf ,
       filename = sprintf("%s/score_distribution_histograms_pm.pdf",figure_file)) 

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
  labs(y = "Gene pleiotropic score for number",x = "Gene pleiotropic scores for magnitude") + 
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
       filename = sprintf("%s/score_distribution_sp.pdf",figure_file)) 
ggsave(plot = yplot, width = 2, height = 5, device = cairo_pdf ,
       filename = sprintf("%s/score_distribution_yplot.pdf",figure_file)) 
ggsave(plot = xplot, width = 5, height = 2, device = cairo_pdf ,
       filename = sprintf("%s/score_distribution_xplot.pdf",figure_file)) 

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
  pn_ld = chrpercent_func(score_name = "pn_ld", x_label = "Proportion of high GPS-N genes\nper chromosome (%)"),
  pm_ld = chrpercent_func(score_name = "pm_ld", x_label = "Proportion of high GPS-M genes\nper chromosome (%)")
)

chr_percent$pn_ld$plot
chr_percent$pm_ld$plot

ggsave(plot = chr_percent$pn_ld$plot, width = 5, height = 8, device = cairo_pdf, 
       filename = sprintf("%s/chr_percent_pn.pdf",figure_file))
