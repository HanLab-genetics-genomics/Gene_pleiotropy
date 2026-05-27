#1 load packages and set variables-------
setwd("/path/to/Gene_pleiotropy-main/Data")
figure_file <- "/path/to/Gene_pleiotropy-main/Output/Gene_pleiotropy_and_gene_age"
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

#2 load data--------------------------
load("pleiotropy_maindata.RData")
#3 figures------------------
##figure2A--------------------
head(pleiotropy_maindata$pm_ld)

cor(pleiotropy_maindata$pn_ld[, .(gene_age_num, use_score)], method = "sp", use = "pairwise.complete.obs")[2,1] %>% sprintf("%.3f", .)
cor.mtest(pleiotropy_maindata$pn_ld[, .(gene_age_num, use_score)], method = "sp")$p[2,1] %>% sprintf("%.2e", .)

cor(pleiotropy_maindata$pm_ld[, .(gene_age_num, use_score)], method = "sp", use = "pairwise.complete.obs")[2,1] %>% sprintf("%.3f", .)
cor.mtest(pleiotropy_maindata$pm_ld[, .(gene_age_num, use_score)], method = "sp")$p[2,1] %>% sprintf("%.2e", .)


###violin plots for pn_ld (GenOrigin)----------------
violin_data <- pleiotropy_maindata$pn_ld[!is.na(gene_age_num)]
violin_data[, gene_age_group := fcase(
  gene_age_num %in% c(3, 8), "3-8",
  gene_age_num %in% c(12, 18), "12-18",
  gene_age_num %in% c(25, 36), "25-36",
  gene_age_num %in% c(55, 70), "55-70",
  gene_age_num %in% c(78, 86, 93), "78-93",
  gene_age_num %in% c(101, 132), "101-132",
  gene_age_num %in% c(168), "168",
  gene_age_num %in% c(244), "244",
  gene_age_num %in% c(332, 382), "332-382",
  gene_age_num %in% c(424, 454), "424-454",
  gene_age_num %in% c(544), "544",
  gene_age_num %in% c(645), "645",
  gene_age_num %in% c(680), "680",
  gene_age_num %in% c(810, 886), "810-886",
  gene_age_num %in% c(950), "950",
  gene_age_num %in% c(987), "987",
  gene_age_num %in% c(1064, 1292), "1064-1292",
  gene_age_num %in% c(1488), "1488",
  gene_age_num %in% c(1714, 1934), "1714-1934",
  gene_age_num %in% c(4290), ">4290",
  default = NA_character_
)]

group_order <- rev(c("3-8", "12-18", "25-36", "55-70", "78-93", "101-132",
                     "168", "244", "332-382", "424-454", "544", "645", "680",
                     "810-886", "950", "987", "1064-1292", "1488", "1714-1934", ">4290"))

violin_data[, gene_age_group := factor(gene_age_group, levels = group_order)]
table(violin_data$gene_age_group)
violin_data[, gene_age_branch := as.numeric(gene_age_group)]
head(violin_data)

upper_plot <- ggbetweenstats(data = violin_data, centrality.plotting = F, 
                             x = gene_age_branch, y = use_score,
                             pairwise.display = "none", bf.message = F, results.subtitle = F, package = "Polychrome",
                             palette = "dark",
                             #violin.args = list(width = 0, linewidth = 0),
                             box.args = list(color = "grey70", linewidth = 0.6),
                             point.args = list(position = ggplot2::position_jitterdodge(jitter.height = 0.5), alpha = 0.1)) +
  labs(y = " ", title = " ") +
  theme(plot.title = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 15), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line.y = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  coord_cartesian(ylim = c(14, 35)) +
  scale_y_continuous(breaks = c(14, 30))


median <- violin_data[, lapply(.SD, median), by = "gene_age_branch", .SDcols = "use_score"][!is.na(gene_age_branch)]

main_plot <- ggbetweenstats(data = violin_data, centrality.plotting = F, 
                            x = gene_age_branch, y = use_score,
                            pairwise.display = "none", bf.message = F, results.subtitle = F, package = "Polychrome",
                            palette = "dark",
                            #violin.args = list(width = 0, linewidth = 0),
                            box.args = list(color = "grey70", linewidth = 0.6),
                            point.args = list(position = ggplot2::position_jitterdodge(jitter.height = 0.5), alpha = 0.1)) +
  labs(y = "Gene pleiotropic score for number (GPS-N)") +
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line.y = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
  coord_cartesian(ylim = c(6, 12)) +
  scale_y_continuous(breaks = c(6, 8, 10, 12)) +
  geom_smooth(data = median, aes(x = gene_age_branch, y = use_score),
              method = "glm", formula = y ~ x, se = FALSE,
              color = "#002058", linewidth = 0.7, alpha = 0.9, linetype="dashed") 



lower_plot <- ggbetweenstats(data = violin_data, centrality.plotting = F, 
                             x = gene_age_branch, y = use_score,
                             pairwise.display = "none", bf.message = F, results.subtitle = F, package = "Polychrome",
                             palette = "dark",
                             #violin.args = list(width = 0, linewidth = 0),
                             box.args = list(color = "grey70", linewidth = 0.6),
                             point.args = list(position = ggplot2::position_jitterdodge(jitter.height = 0.5), alpha = 0.1)) +
  labs(x = "Gene age (million years ago)", y = "") +
  theme(axis.text.x = element_text(size = 12, color = "black", face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  coord_cartesian(ylim = c(0, 5)) +
  scale_y_continuous(breaks = c(2)) +
  scale_x_discrete(labels = group_order)


combined_plot1 <- plot_grid(upper_plot, main_plot, lower_plot, ncol = 1, align = "v", rel_heights = c(0.27,1))
combined_plot1

ggsave(plot = combined_plot1, width = 10, height = 7, device = cairo_pdf,
       filename = sprintf("%s/figure2A_pn_ld_violin_geneorigin.pdf",figure_file))


#4 Supplementary figures for gene age and gene pleiotropy ----------------------------
##violin plots for pm_ld (GenOrigin)----------------
violin_data <- pleiotropy_maindata$pm_ld[!is.na(gene_age_num)]
violin_data[, gene_age_group := fcase(
  gene_age_num %in% c(3, 8), "3-8",
  gene_age_num %in% c(12, 18), "12-18",
  gene_age_num %in% c(25, 36), "25-36",
  gene_age_num %in% c(55, 70), "55-70",
  gene_age_num %in% c(78, 86, 93), "78-93",
  gene_age_num %in% c(101, 132), "101-132",
  gene_age_num %in% c(168), "168",
  gene_age_num %in% c(244), "244",
  gene_age_num %in% c(332, 382), "332-382",
  gene_age_num %in% c(424, 454), "424-454",
  gene_age_num %in% c(544), "544",
  gene_age_num %in% c(645), "645",
  gene_age_num %in% c(680), "680",
  gene_age_num %in% c(810, 886), "810-886",
  gene_age_num %in% c(950), "950",
  gene_age_num %in% c(987), "987",
  gene_age_num %in% c(1064, 1292), "1064-1292",
  gene_age_num %in% c(1488), "1488",
  gene_age_num %in% c(1714, 1934), "1714-1934",
  gene_age_num %in% c(4290), ">4290",
  default = NA_character_
)]

group_order <- rev(c("3-8", "12-18", "25-36", "55-70", "78-93", "101-132",
                     "168", "244", "332-382", "424-454", "544", "645", "680",
                     "810-886", "950", "987", "1064-1292", "1488", "1714-1934", ">4290"))

violin_data[, gene_age_group := factor(gene_age_group, levels = group_order)]
table(violin_data$gene_age_group)
violin_data[, gene_age_branch := as.numeric(gene_age_group)]
head(violin_data)


upper_plot <- ggbetweenstats(data = violin_data, centrality.plotting = F, 
                             x = gene_age_branch, y = use_score,
                             pairwise.display = "none", bf.message = F, results.subtitle = F, package = "Polychrome",
                             palette = "dark",
                             #violin.args = list(width = 0, linewidth = 0),
                             box.args = list(color = "grey70", linewidth = 0.6),
                             point.args = list(position = ggplot2::position_jitterdodge(jitter.height = 0.5), alpha = 0.1)) +
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


median <- violin_data[, lapply(.SD, median), by = "gene_age_branch", .SDcols = "use_score"][!is.na(gene_age_branch)]
main_plot <- ggbetweenstats(data = violin_data, centrality.plotting = F, 
                            x = gene_age_branch, y = use_score,
                            pairwise.display = "none", bf.message = F, results.subtitle = F, package = "Polychrome",
                            palette = "dark",
                            #violin.args = list(width = 0, linewidth = 0),
                            box.args = list(color = "grey70", linewidth = 0.6),
                            point.args = list(position = ggplot2::position_jitterdodge(jitter.height = 0.5), alpha = 0.1)) +
  labs(y = "Gene pleiotropic score for magnitude (GPS-M)") +
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold"), axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.line.y = element_line(colour = "grey50"),
        panel.grid = element_line(color = "#b4aea9"), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) +
  coord_cartesian(ylim = c(5.5, 7.2)) +
  scale_y_continuous(breaks = c(6, 7)) +
  geom_smooth(data = median, aes(x = gene_age_branch, y = use_score),
              method = "glm", formula = y ~ x, se = FALSE,
              color = "#002058", linewidth = 0.7, alpha = 0.9, linetype="dashed") 



lower_plot <- ggbetweenstats(data = violin_data, centrality.plotting = F, 
                             x = gene_age_branch, y = use_score,
                             pairwise.display = "none", bf.message = F, results.subtitle = F, package = "Polychrome",
                             palette = "dark",
                             #violin.args = list(width = 0, linewidth = 0),
                             box.args = list(color = "grey70", linewidth = 0.6),
                             point.args = list(position = ggplot2::position_jitterdodge(jitter.height = 0.5), alpha = 0.1)) +
  labs(x = "Gene age (million years ago)", y = "") +
  theme(axis.text.x = element_text(size = 12, color = "black", face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, color = "black", face = "bold"),
        axis.title = element_text(size = 15, face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  coord_cartesian(ylim = c(1.5, 5)) +
  scale_y_continuous(breaks = c(2, 4)) +
  scale_x_discrete(labels = group_order)


combined_plot2 <- plot_grid(upper_plot, main_plot, lower_plot, ncol = 1, align = "v", rel_heights = c(0.27,1))
combined_plot2

ggsave(plot = combined_plot2, width = 8, height = 7, device = cairo_pdf,
       filename = sprintf("%s/pm_ld_violin_geneorigin.pdf",figure_file))

rm(combined_plot1, combined_plot2, cor_agewithscore_pm, cor_agewithscore_pn, gene_origin_age, geneorgin_age,
   upper_plot, main_plot, lower_plot, median, p , p2, hsd_merge2, hopsgene_ageburden_cut_genemetrics, violin_data)

##violin plots other gene age datasets----------------
head(pleiotropy_maindata$pm_ld)

violin_func <- function(data = pleiotropy_maindata, age_name = "age_protein_num",
                        score = "pm_ld", ylimit = c(4,10), p_value, cor, x_label){
  violin_data <- data[[score]] %>%
    .[, use_age := get(age_name)] %>% .[!is.na(use_age),]
  
  median <- violin_data[, lapply(.SD, median), by = "use_age", .SDcols = "use_score"]
  
  if(score == "pm_ld"){ylab = "Gene pleiotropic score for magnitude (GPS-M)"} else {ylab = "Gene pleiotropic score for number (GPS-N)"}
  
  plot <- ggbetweenstats(data = violin_data, centrality.plotting = F, outlier.shape = NA,
                         x = use_age, y = use_score,
                         pairwise.display = "none", bf.message = F, results.subtitle = F, package = "Polychrome",
                         palette = "dark",
                         violin.args = list(width = 0, linewidth = 0),
                         box.args = list(color = "grey70", linewidth = 0.6),
                         point.args = list(position = ggplot2::position_jitterdodge(jitter.height = 0.5), alpha = 0.07)) +
    coord_trans(ylim = ylimit) +
    labs(y = ylab, title = " ", x = "Gene age stage (old to young)") +
    theme(text = element_text(size = 12, color = "black", face = "bold"),
          plot.title = element_text(family = "Lobster Two", size = 12, face = "bold", color = "black"),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(size = 14, face = "bold"),
          axis.text.y = element_text(size = 10, color = "black", face = "bold"),
          axis.text.x = element_text(size = 10, color = "black", face = "bold", angle = 45, hjust = 1, vjust = 1),
          axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
          panel.grid = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.margin = margin(t = 10, r = 20, b = 20, l = 20)) +
    scale_x_discrete(labels = x_label) +
    geom_smooth(data = median, aes(x = use_age, y = use_score),
                method = "glm", formula = y ~ x, se = FALSE,
                color = "#002058", linewidth = 0.7, alpha = 0.9, linetype="dashed") +
    geom_text(x = Inf, y = Inf, hjust = 1.1, vjust = 1, size = 5, color = "black",
              label = sprintf("Correlation: %s (P = %s)", cor, p_value))
  
  return(plot)
}

##age_phylos
pm_ld_cor = -cor(pleiotropy_maindata$pm_ld[, .(age_phylos_num, use_score)], method = "sp", use = "pairwise.complete.obs")[2,1] %>% sprintf("%.3f", .)
pm_ld_p = cor.mtest(pleiotropy_maindata$pm_ld[, .(age_phylos_num, use_score)], method = "sp")$p[2,1] %>% sprintf("%.2e", .)

pn_ld_cor = -cor(pleiotropy_maindata$pn_ld[, .(age_phylos_num, use_score)], method = "sp", use = "pairwise.complete.obs")[2,1] %>% sprintf("%.3f", .)
pn_ld_p = cor.mtest(pleiotropy_maindata$pn_ld[, .(age_phylos_num, use_score)], method = "sp")$p[2,1] %>% sprintf("%.2e", .)

pm_ld_cor; pm_ld_p; pn_ld_cor; pn_ld_p


phylostrata <- c("Cellular organisms", "Eukaryota", "Opisthokonta", "Holozoa", "Metazoa", "Eumetazoa", "Bilateria",
                 "Deuterostomia", "Chordata", "Olfactores", "Craniata", "Euteleostomi", "Tetrapoda", "Amniota",
                 "Mammalia", "Eutheria", "Boreoeutheria", "Euarchontoglires", "Primates")
age_phylos <- list(pm_ld = violin_func(data = pleiotropy_maindata,
                                       age_name = "age_phylos_num", score = "pm_ld", ylimit = c(4,10),
                                       p_value = pm_ld_p, cor = pm_ld_cor, x_label = phylostrata),
                   pn_ld = violin_func(data = pleiotropy_maindata,
                                       age_name = "age_phylos_num", score = "pn_ld", ylimit = c(0,20),
                                       p_value = pn_ld_p, cor = pn_ld_cor, x_label = phylostrata))

age_phylos_pnpm <- grid.arrange(age_phylos$pn_ld, age_phylos$pm_ld, ncol = 2)
age_phylos_pnpm

##age_protein
pm_ld_cor = -cor(pleiotropy_maindata$pm_ld[, .(age_protein_num, use_score)], method = "sp", use = "pairwise.complete.obs")[2,1] %>% sprintf("%.3f", .)
pm_ld_p = cor.mtest(pleiotropy_maindata$pm_ld[, .(age_protein_num, use_score)], method = "sp")$p[2,1] %>% sprintf("%.2e", .)

pn_ld_cor = -cor(pleiotropy_maindata$pn_ld[, .(age_protein_num, use_score)], method = "sp", use = "pairwise.complete.obs")[2,1] %>% sprintf("%.3f", .)
pn_ld_p = cor.mtest(pleiotropy_maindata$pn_ld[, .(age_protein_num, use_score)], method = "sp")$p[2,1] %>% sprintf("%.2e", .)

pm_ld_cor; pm_ld_p; pn_ld_cor; pn_ld_p

phylostrata <- c("Cellular organisms", "Eukaryota", "Opisthokonta", "Bilateria", "Deuterostomia", "Chordata", "Euteleostomi", "Tetrapoda",
                 "Amniota", "Mammalia", "Theria", "Eutheria", "Euarchontoglires", "Catarrhini", "Homininae", "Human")
age_protein <- list(pm_ld = violin_func(data = pleiotropy_maindata,
                                        age_name = "age_protein_num", score = "pm_ld", ylimit = c(4,10),
                                        p_value = pm_ld_p, cor = pm_ld_cor, x_label = phylostrata),
                    pn_ld = violin_func(data = pleiotropy_maindata,
                                        age_name = "age_protein_num", score = "pn_ld", ylimit = c(0,20),
                                        p_value = pn_ld_p, cor = pn_ld_cor, x_label = phylostrata))

age_protein_pnpm <- grid.arrange(age_protein$pn_ld, age_protein$pm_ld, ncol = 2)
age_protein_pnpm

ggsave(plot = age_phylos_pnpm, width = 16, height = 8, device = cairo_pdf ,
       filename = sprintf("%s/age_phylos_pnpm.pdf",figure_file))
ggsave(plot = age_protein_pnpm, width = 16, height = 8, device = cairo_pdf ,
       filename = sprintf("%s/age_protein_pnpm.pdf",figure_file))


##violin plots of pleiotropic scores among age stages (age stage 7)------------------------
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
       filename = sprintf("%s/agestage7_pm_ld_violin.pdf",figure_file))


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
       filename = sprintf("%s/agestage7_pn_ld_violin.pdf",figure_file))

##proportion plot--------------------
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
