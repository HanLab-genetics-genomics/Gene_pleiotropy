#1 load packages and set variables-------
setwd("/path/to/Gene_pleiotropy-main/Data")
figure_file <- "/path/to/Gene_pleiotropy-main/Output/Supple"
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

################enrichment of high GPS genes in pleiotropy-enriched genomic hotspots####################--------------------
load("pleiotropy_maindata.RData")
head(pleiotropy_maindata$pm_ld)

library(biomaRt)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
regions <- c("6:28477797:33448354",   # broad MHC/HLA region
             "16:29603941:30198600",  # 16p11.2 BP4-BP5
             "22:18600000:21500000",  # 22q11.2 proximal A-D
             "22:35100000:51304566",  # 22q13 
             "17:43600000:44500000"   # 17q21.31 inversion region
)

genes <- getBM(attributes = c("chromosome_name", "start_position", "end_position",
                              "strand", "hgnc_symbol", "ensembl_gene_id"),
               filters = "chromosomal_region",
               values = regions,  mart = mart)
head(genes)

# Clean a bit
genes <- genes[genes$chromosome_name %in% c("6", "16", "17", "22"), ]
genes <- unique(genes)
head(genes)
setDT(genes)

gene_list <- list(gene_MHC = unique(na.omit(genes[chromosome_name == 6 & hgnc_symbol != "",]$hgnc_symbol)),
                  gene_16 = unique(na.omit(genes[chromosome_name == 16 & hgnc_symbol != "",]$hgnc_symbol)),
                  gene_17 = unique(na.omit(genes[chromosome_name == 17 & hgnc_symbol != "",]$hgnc_symbol)),
                  gene_22 = unique(na.omit(genes[chromosome_name == 22 & hgnc_symbol != "",]$hgnc_symbol)))


fisher_test <- function(U, A, B) {
  
  U <- unique(U)
  A <- unique(intersect(A, U))
  B <- unique(intersect(B, U))
  
  m <- length(A)
  n <- length(B)
  N <- length(U)
  k <- length(intersect(A, B))
  
  a <- k
  b <- m - k
  c <- n - k
  d <- N - m - n + k
  cont <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                 dimnames = list(A = c("inA", "notInA"), B = c("inB", "notInB")))
  
  
  fisher_p <- tryCatch(
    fisher.test(cont)$p.value,
    error = function(e) NA_real_
  )
  
  fisher_or <- tryCatch(
    fisher.test(cont)$estimate,
    error = function(e) NA_real_
  )
  
  
  tab <- data.table(a = a, b = b, c = c, d = d, fisher_p = sprintf("%.2e", fisher_p), fisher_or = sprintf("%.2f", fisher_or))
  return(tab)
}

pm_high <- pleiotropy_maindata$pm_ld[pleio_class3 == "High pleiotropy"]$gene
fisher_test_list_pm <- list(gene_MHC = fisher_test(U = pleiotropy_maindata$pm_ld$gene, A = pm_high, B = gene_list$gene_MHC),
                            gene_16 = fisher_test(U = pleiotropy_maindata$pm_ld$gene, A = pm_high, B = gene_list$gene_16),
                            gene_17 = fisher_test(U = pleiotropy_maindata$pm_ld$gene, A = pm_high, B = gene_list$gene_17),
                            gene_22 = fisher_test(U = pleiotropy_maindata$pm_ld$gene, A = pm_high, B = gene_list$gene_22)) %>% rbindlist

fisher_test_list_pm

pn_high <- pleiotropy_maindata$pn_ld[pleio_class3 == "High pleiotropy"]$gene
fisher_test_list_pn <- list(gene_MHC = fisher_test(U = pleiotropy_maindata$pn_ld$gene, A = pn_high, B = gene_list$gene_MHC),
                            gene_16 = fisher_test(U = pleiotropy_maindata$pn_ld$gene, A = pn_high, B = gene_list$gene_16),
                            gene_17 = fisher_test(U = pleiotropy_maindata$pn_ld$gene, A = pn_high, B = gene_list$gene_17),
                            gene_22 = fisher_test(U = pleiotropy_maindata$pn_ld$gene, A = pn_high, B = gene_list$gene_22)) %>% rbindlist

fisher_test_list_pn


####################PCA analysis#####################--------------------
#1 load data-----------------
GPS_feature_all <- fread("GPS_feature_all.csv", na.strings = c("", "NA"))
head(GPS_feature_all, 2)

#"ogee_connectivity" #6467 NAs
pca_vars <- c("gene_age_num",
              "Gene length (bp)", "CDS Length", "Transcript count", "GC content",  "Number of SNPs (Gene)", "Exon Counts",
              "LOEUF", "Missense/Synonymous ratio","Nonsense/Synonymous ratio",
              "tau", "tf_sum", "enhancer_n", "po_pergene")

dt_pca <- copy(GPS_feature_all[, c("gene", "pm_ld", "pn_ld", pca_vars), with = FALSE])
dt_pca[, lapply(.SD, function(x) sum(is.na(x)))]

dt_pca <- na.omit(dt_pca) #14277
head(dt_pca)
nrow(dt_pca)

###transformations for skewed variables
library(e1071)
skew_stat <- rbindlist(lapply(pca_vars, function(v) {
  x <- dt_pca[[v]]
  x <- x[is.finite(x)]
  
  data.table(variable = v, skewness = skewness(x, na.rm = TRUE, type = 2))}))

skew_stat[order(-abs(skewness))] ##remove skewness > 1

log_vars <- intersect(c("Gene length (bp)", "CDS Length", "Transcript count", "Number of SNPs (Gene)", "Exon Counts", 
                        "Nonsense/Synonymous ratio", "gene_age_num",
                        "tf_sum", "enhancer_n", "po_pergene"), names(dt_pca))

for (v in log_vars) {dt_pca[[v]] <- log1p(dt_pca[[v]])}
colSums(is.na(dt_pca[, ..pca_vars]))
x <- as.matrix(dt_pca[, ..pca_vars])
head(x)

#2 PCA analysis-----------------
pca_fit <- prcomp(x, center = TRUE, scale. = TRUE)

# variance explained
pve <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
pve_dt <- data.table(
  PC = paste0("PC", seq_along(pve)),
  variance_explained = pve,
  cumulative_variance = cumsum(pve)
)

pve_dt

##plot the variance_explained
pve_dt[, PC := factor(PC, levels = paste0('PC', 1:14))]

variance_explained_plot <- ggplot(pve_dt, aes(x = PC, y = variance_explained)) +
  geom_col(width = 0.7, fill = "#8DA0CB", color = "#8DA0CB", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f", variance_explained * 100)),
            vjust = 0.4, hjust = -0.2, size = 4.5) +
  scale_y_continuous(labels = function(x) x * 100,
                     limits = c(0, max(pve_dt$variance_explained) * 1.15),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "Explained variance (%)", title = "") +
  theme_bw() +
  coord_flip() +
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        legend.position = "none",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank())

variance_explained_plot

ggsave(plot = variance_explained_plot, width = 4, height = 5.5, device = cairo_pdf,
       filename = sprintf("%s/variance_explained_plot.pdf", figure_file))

#3 get loadings---------------------
loadings_dt <- as.data.table(pca_fit$rotation, keep.rownames = "feature")
loadings_dt

n_pc <- which(pve_dt$cumulative_variance >= 0.80)[1]
keep_pcs <- pve_dt$PC[1:n_pc]

load_long <- melt(loadings_dt, id.vars = "feature",  variable.name = "PC",  value.name = "loading")
load_long <- load_long[PC %in% keep_pcs]
load_long[, abs_loading := abs(loading)]

head(load_long)

feature_labels <- c("gene_age_num" = "Gene age",
                    "Gene Length (bp)" = "Gene length",
                    "CDS Length" = "CDS length",
                    "Transcript count" = "Transcript count",
                    "GC content" = "GC content",
                    "Number of SNPs (Gene)" = "Number of SNPs per gene",
                    "Exon Counts" = "Exon count",
                    "LOEUF" = "LOEUF",
                    "Missense/Synonymous ratio" = "Missense/synonymous ratio",
                    "Nonsense/Synonymous ratio" = "Nonsense/synonymous ratio",
                    "tau" = "Expression specificity",
                    "tf_sum" = "TF number",
                    "enhancer_n" = "Enhancer number",
                    "po_pergene" = "P-O interactions")

load_long[, feature_label := feature]
for (nm in names(feature_labels)) {
  load_long[feature == nm, feature_label := feature_labels[[nm]]]
}
head(load_long)

pc1_order <- load_long[PC == "PC1"][order(-abs_loading)]$feature_label
load_long[, feature_label := factor(feature_label, levels = rev(pc1_order))]
load_long[, PC := factor(PC, levels = keep_pcs)]
head(load_long)

##plot
PC_loadings_plot <- ggplot(load_long, aes(PC, feature_label)) +
  geom_tile(fill = "white", color = "grey50", linewidth = 0.3) +
  geom_point(aes(size = abs_loading, color = abs_loading), shape = 16) +
  scale_size_area(
    max_size = 13,
    breaks = c(0.2, 0.4, 0.6),
    name = "Absolute\nloading"
  ) +
  scale_color_gradient(
    low = "#F5F5F5",
    high = "#2166AC",
    breaks = c(0.1, 0.2, 0.3, 0.4, 0.5),
    name = "Absolute\nloading"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    text = element_text(size = 12, color = "black", face = "bold"),
    legend.position = "right",
    axis.title.y = element_text(size = 15, face = "bold"),
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank()
  )

PC_loadings_plot


ggsave(plot = PC_loadings_plot, width = 7, height = 5.5, device = cairo_pdf,
       filename = sprintf("%s/PC_loadings.pdf", figure_file))



##3.1 plot PC1 and PC2 loadings----------------------
library(ggrepel)

# PCA loadings
loadings_dt <- as.data.table(pca_fit$rotation, keep.rownames = "feature")
loadings_dt[, feature_label := feature_labels[feature]]
loadings_dt[is.na(feature_label), feature_label := feature]

loadings_dt[, feature_group := fifelse(feature %in% c("Gene length (bp)", "CDS Length", "Transcript count", "GC content",
                                                      "Number of SNPs (Gene)", "Exon Counts"), "Structure",
                                       fifelse(feature %in% c("LOEUF", "Missense/Synonymous ratio", "Nonsense/Synonymous ratio"), "Constraint",
                                               fifelse(feature %in% c("gene_age_num"), "Gene age", "Regulation")))]

loadings_dt[, feature_group := factor(feature_group, levels = c("Structure", "Constraint", "Gene age", "Regulation"))]
loadings_dt

##plot
pc1_var <- round(pve_dt[PC == "PC1", variance_explained] * 100, 1)
pc2_var <- round(pve_dt[PC == "PC2", variance_explained] * 100, 1)

p_loading_map_PC12 <- ggplot(loadings_dt, aes(x = PC1, y = PC2)) +
  ggforce::geom_mark_ellipse(
    aes(fill = feature_group, group = feature_group),
    alpha = 0.10,
    color = NA,
    expand = unit(2.5, "mm"),
    radius = unit(2.5, "mm"),
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.45) +
  geom_point(aes(color = feature_group), size = 3.8, alpha = 0.95) +
  geom_text_repel(
    aes(label = feature_label, color = feature_group),
    size = 4.5,
    box.padding = 0.4,
    point.padding = 0.25,
    segment.color = "grey70",
    segment.linewidth = 0.35,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Structure" = "#4C78A8",
    "Constraint" = "#F58518",
    "Regulation" = "#54A24B",
    "Gene age" = "#B279A2"
  )) +
  scale_fill_manual(values = c(
    "Structure" = "#4C78A8",
    "Constraint" = "#F58518",
    "Regulation" = "#54A24B",
    "Gene age" = "#B279A2"
  )) +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = NULL
  ) +
  coord_equal(xlim = c(-0.5, 0.6), ylim = c(-0.6, 0.75)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank()) 

p_loading_map_PC12

ggsave(plot = p_loading_map_PC12, width = 6, height = 6, device = cairo_pdf,
       filename = sprintf("%s/p_loading_map_PC12.pdf", figure_file))


##3.2 plot PC3 and PC4 loadings----------------------
pc3_var <- round(pve_dt[PC == "PC3", variance_explained] * 100, 1)
pc4_var <- round(pve_dt[PC == "PC4", variance_explained] * 100, 1)

p_loading_map_PC34 <- ggplot(loadings_dt, aes(x = PC3, y = PC4)) +
  ggforce::geom_mark_ellipse(
    aes(fill = feature_group, group = feature_group),
    alpha = 0.10,
    color = NA,
    expand = unit(2.5, "mm"),
    radius = unit(2.5, "mm"),
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.45) +
  geom_point(aes(color = feature_group), size = 3.8, alpha = 0.95) +
  geom_text_repel(
    aes(label = feature_label, color = feature_group),
    size = 4.5,
    box.padding = 0.4,
    point.padding = 0.25,
    segment.color = "grey70",
    segment.linewidth = 0.35,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Structure" = "#4C78A8",
    "Constraint" = "#F58518",
    "Regulation" = "#54A24B",
    "Gene age" = "#B279A2"
  )) +
  scale_fill_manual(values = c(
    "Structure" = "#4C78A8",
    "Constraint" = "#F58518",
    "Regulation" = "#54A24B",
    "Gene age" = "#B279A2"
  )) +
  labs(
    x = paste0("PC3 (", pc3_var, "%)"),
    y = paste0("PC4 (", pc4_var, "%)"),
    color = NULL
  ) +
  coord_equal(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.75)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank()) 

p_loading_map_PC34

ggsave(plot = p_loading_map_PC34, width = 6, height = 6, device = cairo_pdf,
       filename = sprintf("%s/p_loading_map_PC34.pdf", figure_file))


##3.3 plot PC5 and PC6 loadings----------------------
pc5_var <- round(pve_dt[PC == "PC5", variance_explained] * 100, 1)
pc6_var <- round(pve_dt[PC == "PC6", variance_explained] * 100, 1)

p_loading_map_PC56 <- ggplot(loadings_dt, aes(x = PC5, y = PC6)) +
  ggforce::geom_mark_ellipse(
    aes(fill = feature_group, group = feature_group),
    alpha = 0.10,
    color = NA,
    expand = unit(2.5, "mm"),
    radius = unit(2.5, "mm"),
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.45) +
  geom_point(aes(color = feature_group), size = 3.8, alpha = 0.95) +
  geom_text_repel(
    aes(label = feature_label, color = feature_group),
    size = 4.5,
    box.padding = 0.4,
    point.padding = 0.25,
    segment.color = "grey70",
    segment.linewidth = 0.35,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(
    "Structure" = "#4C78A8",
    "Constraint" = "#F58518",
    "Regulation" = "#54A24B",
    "Gene age" = "#B279A2"
  )) +
  scale_fill_manual(values = c(
    "Structure" = "#4C78A8",
    "Constraint" = "#F58518",
    "Regulation" = "#54A24B",
    "Gene age" = "#B279A2"
  )) +
  labs(
    x = paste0("PC5 (", pc5_var, "%)"),
    y = paste0("PC6 (", pc6_var, "%)"),
    color = NULL
  ) +
  coord_equal(xlim = c(-0.6, 0.6), ylim = c(-0.6, 0.4)) +
  theme_bw() +
  theme(legend.position = "right",
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank()) 

p_loading_map_PC56

ggsave(plot = p_loading_map_PC56, width = 6, height = 6, device = cairo_pdf,
       filename = sprintf("%s/p_loading_map_PC56.pdf", figure_file))

#4 check correlations of GPS-N and PCs------------------
scores_dt <- as.data.table(pca_fit$x)
scores_dt[, gene := dt_pca$gene]
head(scores_dt)

dt_pca2 <- cbind(dt_pca[, .(gene, pm_ld, pn_ld)], scores_dt[, !"gene"])
head(dt_pca2)

## original PCs only
pc_cols <- grep("^PC[0-9]+$", names(dt_pca2), value = TRUE)

cor_pn_raw <- rbindlist(lapply(pc_cols, function(pc) {
  ct <- cor.test(dt_pca2[[pc]], dt_pca2$pn_ld, method = "spearman")
  data.table(
    PC = pc,
    rho = unname(ct$estimate),
    p = ct$p.value
  )
}))

cor_pn_raw[, p_adj := p.adjust(p, method = "BH")]
cor_pn_raw[, flip := rho < 0]

cor_pn_raw

## make oriented copies
dt_pca_ori <- copy(dt_pca2)
loadings_ori <- copy(loadings_dt)

for (i in seq_len(nrow(cor_pn_raw))) {
  pc <- cor_pn_raw$PC[i]
  if (cor_pn_raw$flip[i]) {
    dt_pca_ori[, (pc) := -get(pc)]
    loadings_ori[, (pc) := -get(pc)]
  }
}


cor_pn <- rbindlist(lapply(pc_cols, function(pc) {
  ct <- cor.test(dt_pca_ori[[pc]], dt_pca_ori$pn_ld, method = "spearman")
  data.table(
    PC = pc,
    rho = unname(ct$estimate),
    p = ct$p.value
  )
}))

cor_pn[, p_adj := p.adjust(p, method = "BH")]
cor_pn[, log10_padj := -log10(pmax(p_adj, .Machine$double.xmin))]
head(cor_pn)

## sort by effect size
setorder(cor_pn, -rho)
cor_pn[, PC_plot := factor(PC, levels = rev(PC))]
head(cor_pn)


##Horizontal bar plot
pn_cor_bar <- ggplot(cor_pn[PC_plot %in% paste0("PC", 1:7)], aes(x = rho, y = PC_plot, fill = log10_padj)) +
  geom_col(width = 0.75) +
  scale_fill_gradient(low='#7FB65D', high='#DA7235',
                      name = expression(-log[10]("BH-adjusted P"))) +
  labs(x = "Spearman's |ρ|", y = "Principal component", title = "GPS-N") +
  theme_bw() +
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        legend.position = c(0.7, 0.15),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

pn_cor_bar


#5 check correlations of GPS-M and PCs------------------
cor_pm_raw <- rbindlist(lapply(pc_cols, function(pc) {
  ct <- cor.test(dt_pca2[[pc]], dt_pca2$pm_ld, method = "spearman")
  data.table(
    PC = pc,
    rho = unname(ct$estimate),
    p = ct$p.value
  )
}))

cor_pm_raw[, p_adj := p.adjust(p, method = "BH")]
cor_pm_raw[, flip := rho < 0]
cor_pm_raw

## make oriented copies
dt_pca_ori <- copy(dt_pca2)
loadings_ori <- copy(loadings_dt)

for (i in seq_len(nrow(cor_pm_raw))) {
  pc <- cor_pm_raw$PC[i]
  if (cor_pm_raw$flip[i]) {
    dt_pca_ori[, (pc) := -get(pc)]
    loadings_ori[, (pc) := -get(pc)]
  }
}


cor_pm <- rbindlist(lapply(pc_cols, function(pc) {
  ct <- cor.test(dt_pca_ori[[pc]], dt_pca_ori$pm_ld, method = "spearman")
  data.table(
    PC = pc,
    rho = unname(ct$estimate),
    p = ct$p.value
  )
}))

cor_pm[, p_adj := p.adjust(p, method = "BH")]
cor_pm[, log10_padj := -log10(pmax(p_adj, .Machine$double.xmin))]
head(cor_pm)

## sort by effect size
setorder(cor_pm, -rho)
cor_pm[, PC_plot := factor(PC, levels = rev(PC))]
head(cor_pm)


##Horizontal bar plot
pm_cor_bar <- ggplot(cor_pm[PC_plot %in% paste0("PC", 1:7)], aes(x = rho, y = PC_plot, fill = log10_padj)) +
  geom_col(width = 0.75) +
  scale_fill_gradient(low='#7FB65D', high='#DA7235',
                      name = expression(-log[10]("BH-adjusted P"))) +
  labs(x = "Spearman's |ρ|", y = "Principal component", title = "GPS-M") +
  theme_bw() +
  theme(text = element_text(size = 12, color = "black", face = "bold"),
        legend.position = c(0.7, 0.15),
        axis.title.y = element_text(size = 15, face = "bold"),
        axis.title.x = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 12, color = "black", face = "bold"),
        axis.ticks = element_blank(), axis.line = element_line(colour = "grey50"),
        panel.grid = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

pm_cor_bar

pnpm_cor_bar <- grid.arrange(pn_cor_bar, pm_cor_bar, ncol = 2)
pnpm_cor_bar

ggsave(plot = pnpm_cor_bar, width = 8, height = 7, device = cairo_pdf,
       filename = sprintf("%s/PCA_pnpm_cor_bar_PC7.pdf", figure_file))


#6 test duplication---------------
library(broom)

dt_pca_dup <- merge(
  dt_pca2,
  unique(GPS_feature_all[, .(gene, ifhsd2)]),
  by = "gene"
)

dt_pca_dup[, ifhsd2 := factor(ifhsd2, levels = c("Singletons", "Duplicates"))]
dt_pca_dup[, PC2r := -PC2]
dt_pca_dup[, PC5r := -PC5]

## duplication models
fit_dup_logit <- glm(ifhsd2 ~ PC1 + PC2 + PC3 + PC4 +  PC5 + PC6 + PC7, data = dt_pca_dup, family = binomial())

fit_pc1_dup  <- lm(PC1  ~ ifhsd2, data = dt_pca_dup)
fit_pc2_dup  <- lm(PC2r ~ ifhsd2, data = dt_pca_dup)
fit_pc5_dup  <- lm(PC5r ~ ifhsd2, data = dt_pca_dup)

## duplication incremental contribution
fit_pm_base_dup <- lm(pm_ld ~ PC1 + PC2 + PC5 + PC3 + PC4 + PC6 + PC7, data = dt_pca_dup)
fit_pm_add_dup  <- lm(pm_ld ~ PC1 + PC2 + PC5 + PC3 + PC4 + PC6 + PC7 + ifhsd2, data = dt_pca_dup)

fit_pn_base_dup <- lm(pn_ld ~ PC1 + PC2 + PC5 + PC3 + PC4 + PC6 + PC7, data = dt_pca_dup)
fit_pn_add_dup  <- lm(pn_ld ~ PC1 + PC2 + PC5 + PC3 + PC4 + PC6 + PC7 + ifhsd2, data = dt_pca_dup)

#duplication relative to PCs
tab_dup_logit <- as.data.table(
  tidy(fit_dup_logit, conf.int = TRUE, exponentiate = TRUE)
)[term != "(Intercept)"]

setnames(
  tab_dup_logit,
  c("estimate", "conf.low", "conf.high", "p.value"),
  c("OR", "OR_low", "OR_high", "P")
)

tab_dup_logit[, term := factor(term, levels = paste0("PC", 1:7))]
tab_dup_logit

#Group-difference table
tab_dup_pc <- rbindlist(list(
  as.data.table(tidy(fit_pc1_dup, conf.int = TRUE))[term == "ifhsd2Duplicates"][, PC := "PC1"],
  as.data.table(tidy(fit_pc2_dup, conf.int = TRUE))[term == "ifhsd2Duplicates"][, PC := "PC2r"],
  as.data.table(tidy(fit_pc5_dup, conf.int = TRUE))[term == "ifhsd2Duplicates"][, PC := "PC5r"]
), fill = TRUE)

tab_dup_pc <- tab_dup_pc[, .(
  PC,
  beta = estimate,
  CI_low = conf.low,
  CI_high = conf.high,
  P = p.value
)]

tab_dup_pc


#incremental contribution of duplication beyond PCs
anova_pm_dup <- anova(fit_pm_base_dup, fit_pm_add_dup)
anova_pn_dup <- anova(fit_pn_base_dup, fit_pn_add_dup)

tab_dup_increment <- data.table(
  metric = c("pm_ld", "pn_ld"),
  R2_base = c(summary(fit_pm_base_dup)$r.squared,
              summary(fit_pn_base_dup)$r.squared),
  R2_plus_dup = c(summary(fit_pm_add_dup)$r.squared,
                  summary(fit_pn_add_dup)$r.squared),
  delta_R2 = c(summary(fit_pm_add_dup)$r.squared - summary(fit_pm_base_dup)$r.squared,
               summary(fit_pn_add_dup)$r.squared - summary(fit_pn_base_dup)$r.squared),
  P_increment = c(anova_pm_dup$`Pr(>F)`[2],
                  anova_pn_dup$`Pr(>F)`[2])
)

tab_dup_increment

table_all <- list(Duplication_logistic_OR = tab_dup_logit,
                  Duplication_group_diff = tab_dup_pc,
                  Duplication_increment = tab_dup_increment)

table_all$Duplication_logistic_OR
table_all$Duplication_group_diff
table_all$Duplication_increment

####################robustness of GPS framework#####################--------------------
#1 load data-----------------
newGPS_file <- list.files(path = "GPS", pattern = "^hopsgene_.*_cut\\.RData$", full.names = T)
newGPS_file

load_rdata_obj <- function(f) {
  e <- new.env()
  obj_names <- load(f, envir = e)
  e[[obj_names[1]]]
}

gps_obj_list <- lapply(newGPS_file, load_rdata_obj)
names(gps_obj_list) <- sub("_cut\\.RData$", "", sub("^hopsgene_", "", basename(newGPS_file)))
names(gps_obj_list[[1]])
names(gps_obj_list)

get_metric_dt <- function(obj_list, metric = c("pn_ld", "pm_ld")) {
  metric <- match.arg(metric)
  nms <- names(obj_list)
  
  hit <- grep(metric, nms, ignore.case = TRUE, value = TRUE)
  
  if (length(hit) >= 1) {
    dt <- obj_list[[hit[1]]]
  } else {
    if (metric == "pn_ld") {
      dt <- obj_list[[1]]
    } else {
      dt <- obj_list[[2]]
    }
  }
  as.data.table(dt)
}
rename_metric_dt <- function(dt, suffix, metric = c("pn_ld", "pm_ld"), gene_id_col = "gene") {
  metric <- match.arg(metric)
  dt <- copy(as.data.table(dt))
  
  # keep only needed columns if they exist
  keep_cols <- intersect(c(gene_id_col, "use_score", "pleio_class5"), names(dt))
  dt <- dt[, ..keep_cols]
  
  # rename
  if ("use_score" %in% names(dt)) {
    setnames(dt, "use_score", paste0(metric, "_", suffix))
  }
  if ("pleio_class5" %in% names(dt)) {
    setnames(dt, "pleio_class5", paste0(metric, "_class5_", suffix))
  }
  
  dt
}

gps_dt_list <- list()
for (nm in names(gps_obj_list)) {
  suffix <- nm
  
  pn_dt <- get_metric_dt(gps_obj_list[[nm]], "pn_ld")
  pm_dt <- get_metric_dt(gps_obj_list[[nm]], "pm_ld")
  
  pn_dt <- rename_metric_dt(pn_dt, suffix = suffix, metric = "pn_ld", gene_id_col = "gene")
  pm_dt <- rename_metric_dt(pm_dt, suffix = suffix, metric = "pm_ld", gene_id_col = "gene")
  
  gps_dt_list[[suffix]] <- merge(pn_dt, pm_dt, by = "gene", all = TRUE)
}

names(gps_dt_list)

mergeGPS_sensitivity <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), gps_dt_list)
dim(mergeGPS_sensitivity)
head(mergeGPS_sensitivity)

##reproduction GPS
load("hopsgene_repro_cut.RData")
head(hopsgene_repro_cut$pm_ld)
GPS_repro <- merge(hopsgene_repro_cut$pn_ld[, .(gene, pn_ld_repro = use_score, pn_ld_class5_repro = pleio_class5)],
                   hopsgene_repro_cut$pm_ld[, .(gene, pm_ld_repro = use_score, pm_ld_class5_repro = pleio_class5)], by = "gene")
head(GPS_repro)

mergeGPS_repro <- merge(mergeGPS_sensitivity, GPS_repro, by = "gene", all.x = T)
head(mergeGPS_repro)


##original GPS
load("pleiotropy_maindata.RData")
head(pleiotropy_maindata$pm_ld)
GPS <- merge(pleiotropy_maindata$pn_ld[, .(gene, pn_ld_raw = use_score, pn_ld_class5_raw = pleio_class5)],
             pleiotropy_maindata$pm_ld[, .(gene, pm_ld_raw = use_score, pm_ld_class5_raw = pleio_class5, gene_age_num)], by = "gene")

head(GPS)

mergeGPS_all <- merge(mergeGPS_repro, GPS, by = "gene", all.y = T)
head(mergeGPS_all)

#2 correlations between raw GPS and each sensitivity GPS------------------
get_spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  
  tt <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  list(rho = unname(tt$estimate), p = tt$p.value, n = sum(ok))
}

suffixes <- c("mafbin1", "mafbin2", "mafbin3", "pnadjust", "polygenicity", "samplematch", "repro")
cor_res <- rbindlist(lapply(suffixes, function(suf) {
  pn_res <- get_spearman(mergeGPS_all[["pn_ld_raw"]], mergeGPS_all[[paste0("pn_ld_", suf)]])
  pm_res <- get_spearman(mergeGPS_all[["pm_ld_raw"]], mergeGPS_all[[paste0("pm_ld_", suf)]])
  
  rbind(
    data.table(
      analysis = suf,
      metric = "pn",
      rho = pn_res$rho,
      p = pn_res$p,
      n = pn_res$n
    ),
    data.table(
      analysis = suf,
      metric = "pm",
      rho = pm_res$rho,
      p = pm_res$p,
      n = pm_res$n
    )
  )
}))

cor_res

cor_age <- rbindlist(lapply(c("raw", "mafbin1", "mafbin2", "mafbin3", 
                              "pnadjust", "polygenicity", "samplematch", "repro"), function(suf) {
                                pn_res <- get_spearman(mergeGPS_all$gene_age_num, mergeGPS_all[[paste0("pn_ld_", suf)]])
                                pm_res <- get_spearman(mergeGPS_all$gene_age_num, mergeGPS_all[[paste0("pm_ld_", suf)]])
                                
                                rbind(
                                  data.table(
                                    analysis = suf,
                                    metric = "pn",
                                    rho = sprintf("%.2f", pn_res$rho),
                                    p = sprintf("%.2e", pn_res$p),
                                    n = pn_res$n
                                  ),
                                  data.table(
                                    analysis = suf,
                                    metric = "pm",
                                    rho = sprintf("%.2f", pm_res$rho),
                                    p = sprintf("%.2e", pm_res$p),
                                    n = pm_res$n
                                  )
                                )
                              }))

cor_age

#3 preserved gene percent in extreme groups---------------
extreme_preserve <- function(ref_group, new_group, low_label = "1st", high_label = "5th") {
  ok <- complete.cases(ref_group, new_group)
  ref_group <- as.character(ref_group[ok])
  new_group <- as.character(new_group[ok])
  N <- length(ref_group)
  
  one_group_stats <- function(target_label, group_name) {
    ref_idx <- which(ref_group == target_label)
    new_idx <- which(new_group == target_label)
    
    k <- length(intersect(ref_idx, new_idx))   # overlap
    m <- length(ref_idx)                       # size of reference group
    n <- length(new_idx)                       # size of new group
    
    retention <- if (m > 0) k / m else NA_real_
    jaccard <- if ((m + n - k) > 0) k / (m + n - k) else NA_real_
    
    expected_k <- if (N > 0) m * n / N else NA_real_
    expected_retention <- if (N > 0) n / N else NA_real_
    
    # 2x2 table for Fisher exact test
    #                new=in     new=out
    # ref=in           k       m-k
    # ref=out         n-k   N-m-n+k
    mat <- matrix(
      c(k, m - k,
        n - k, N - m - n + k),
      nrow = 2,
      byrow = TRUE
    )
    
    ft <- fisher.test(mat, alternative = "greater")
    
    data.table(
      group = group_name,
      N = N,
      k = k,
      m_ref = m,
      n_new = n,
      retention = retention,
      jaccard = jaccard,
      expected_k = expected_k,
      expected_retention = expected_retention,
      fisher_or = unname(ft$estimate),
      fisher_p = ft$p.value
    )
  }
  
  rbind(
    one_group_stats(low_label, "low"),
    one_group_stats(high_label, "high")
  )
}

preserve_res <- rbindlist(lapply(suffixes, function(suf) {
  pn_dt <- extreme_preserve(
    mergeGPS_all[["pn_ld_class5_raw"]],
    mergeGPS_all[[paste0("pn_ld_class5_", suf)]]
  )[, `:=`(analysis = suf, metric = "pn")]
  
  pm_dt <- extreme_preserve(
    mergeGPS_all[["pm_ld_class5_raw"]],
    mergeGPS_all[[paste0("pm_ld_class5_", suf)]]
  )[, `:=`(analysis = suf, metric = "pm")]
  
  rbind(pn_dt, pm_dt, fill = TRUE)
}), fill = TRUE)

preserve_res

#4 output table-----------------
preserve_wide <- dcast(preserve_res[group %in% c("low", "high")],
                       analysis + metric ~ group,
                       value.var = c("retention", "jaccard", "k", "m_ref", "n_new", "fisher_or", "fisher_p"),
                       sep = "_")
head(preserve_wide)

summary_result <- merge(cor_res, preserve_wide,  by = c("analysis", "metric"),  all = TRUE)
head(summary_result)
fwrite(summary_result, file = sprintf("%s/newGPS_rawGPS_correlations.csv", figure_file))

analysis_order <- c("pnadjust", "polygenicity", "samplematch", "repro", "mafbin1", "mafbin2", "mafbin3")
summary_result_format <- summary_result[analysis %in% analysis_order
                                        , .(analysis, metric, 
                                            rho = sprintf("%.2f", rho),
                                            high_overlap = sprintf("%.1f", retention_high * 100),
                                            low_overlap  = sprintf("%.1f", retention_low * 100))][ 
                                              , analysis := factor(analysis, levels = analysis_order)][order(analysis, metric)]
head(summary_result_format)
summary_result_format[metric == "pm"]$metric <- "GPS-M"
summary_result_format[metric == "pn"]$metric <- "GPS-N"
summary_result_format

#5 gene characters across sensitivity GPS groups------------------
library(dplyr); library(tibble); library(purrr); library(effectsize)
head(mergeGPS_all)

GPS_feature_all <- fread("GPS_feature_all.csv")
head(GPS_feature_all)

##5.1 function--------------
fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))
  )
}

median_iqr_str <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  med <- median(x)
  q1  <- quantile(x, 0.25, names = FALSE)
  q3  <- quantile(x, 0.75, names = FALSE)
  paste0(fmt_num(med, digits), " (", fmt_num(q1, digits), " - ", fmt_num(q3, digits), ")")
}

n_pct_str <- function(x, level, digits = 1) {
  denom <- sum(!is.na(x))
  n <- sum(x == level, na.rm = TRUE)
  if (denom == 0) return(NA_character_)
  paste0(format(n, big.mark = ","), " (", fmt_num(100 * n / denom, digits), "%)")
}

epsilon2_kruskal <- function(x, g) {
  keep <- complete.cases(x, g)
  x <- x[keep]
  g <- droplevels(factor(g[keep]))
  
  n <- length(x)
  k <- nlevels(g)
  if (n == 0 || k < 2 || n <= k) return(NA_real_)
  
  H <- unname(kruskal.test(x ~ g)$statistic)
  (H - k + 1) / (n - k)
}

low_high_rrb <- function(x, g, low_level, high_level, digits = 2) {
  keep <- !is.na(x) & !is.na(g) & g %in% c(low_level, high_level)
  x <- x[keep]
  g <- droplevels(factor(g[keep], levels = c(low_level, high_level)))
  
  if (length(unique(g)) < 2) {
    return(list(
      p = NA_real_,
      effect = NA_character_
    ))
  }
  
  # Wilcoxon rank-sum P
  wt <- suppressWarnings(wilcox.test(x ~ g, exact = FALSE))
  
  # rank-biserial correlation; order = high vs low so positive means higher in high group
  x_high <- x[g == high_level]
  x_low  <- x[g == low_level]
  rb <- effectsize::rank_biserial(x_high, x_low, ci = 0.95)
  
  est_col <- grep("biserial", names(rb), value = TRUE)[1]
  est <- rb[[est_col]][1]
  cil <- rb$CI_low[1]
  cih <- rb$CI_high[1]
  
  list(
    p = wt$p.value,
    effect = paste0(fmt_num(est, digits), " (", fmt_num(cil, digits), " to ", fmt_num(cih, digits), ")")
  )
}

overall_chisq_cramers_v <- function(x, g, digits = 2) {
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]
  
  tab <- table(g, x)
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(list(
      p = NA_real_,
      effect = NA_character_
    ))
  }
  
  ct <- suppressWarnings(chisq.test(tab, correct = FALSE))
  cv <- suppressWarnings(effectsize::cramers_v(tab))
  est_col <- grep("Cramers", names(cv), value = TRUE)[1]
  v <- cv[[est_col]][1]
  
  list(
    p = ct$p.value,
    effect = fmt_num(v, digits)
  )
}

low_high_or <- function(x, g, level, low_level, high_level, digits = 2) {
  keep <- !is.na(x) & !is.na(g) & g %in% c(low_level, high_level)
  x <- x[keep]
  g <- g[keep]
  
  if (length(unique(g)) < 2) {
    return(list(
      p = NA_real_,
      effect = NA_character_
    ))
  }
  
  a <- sum(g == high_level & x == level, na.rm = TRUE)
  b <- sum(g == high_level & x != level, na.rm = TRUE)
  c <- sum(g == low_level  & x == level, na.rm = TRUE)
  d <- sum(g == low_level  & x != level, na.rm = TRUE)
  
  tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  # rows = high, low; cols = level present, other
  rownames(tab) <- c("high", "low")
  colnames(tab) <- c("present", "other")
  
  if (sum(tab) == 0 || any(rowSums(tab) == 0) || any(colSums(tab) == 0)) {
    return(list(
      p = NA_real_,
      effect = NA_character_
    ))
  }
  
  ft <- fisher.test(tab)
  or  <- unname(ft$estimate)
  cil <- ft$conf.int[1]
  cih <- ft$conf.int[2]
  
  list(
    p = ft$p.value,
    effect = paste0(fmt_num(or, digits), " (", fmt_num(cil, digits), " to ", fmt_num(cih, digits), ")")
  )
}


make_supp_table <- function(data,
                            group_var,
                            group_order,
                            low_level,
                            high_level,
                            cont_vars = NULL,
                            cat_specs = NULL,
                            row_labels = NULL,
                            digits_cont = 0,
                            digits_eff = 2,
                            digits_pct = 1) {
  
  df <- data %>%
    dplyr::mutate(.group = factor(.data[[group_var]], levels = group_order)) %>%
    dplyr::filter(!is.na(.group))
  
  out <- list()
  
  # -------------------------
  # Continuous variables
  # -------------------------
  if (!is.null(cont_vars)) {
    for (v in cont_vars) {
      x <- df[[v]]
      g <- df$.group
      
      overall_p <- suppressWarnings(kruskal.test(x ~ g)$p.value)
      eps2 <- epsilon2_kruskal(x, g)
      low_high <- low_high_rrb(x, g, low_level, high_level, digits = digits_eff)
      
      row <- tibble(
        Characteristic = paste0(v, ", median (IQR)"),
        Total = median_iqr_str(x, digits = digits_cont),
        Overall_P = fmt_p(overall_p),
        Omnibus_effect = paste0("ε² = ", fmt_num(eps2, digits_eff)),
        Low_vs_high_P = fmt_p(low_high$p),
        Low_vs_high_effect = low_high$effect
      )
      
      for (lev in group_order) {
        row[[lev]] <- median_iqr_str(df[df$.group == lev, v, drop = TRUE], digits = digits_cont)
      }
      
      row <- row %>%
        dplyr::select(Characteristic, Total, all_of(group_order), Overall_P, Omnibus_effect,
                      Low_vs_high_P, Low_vs_high_effect)
      
      out[[length(out) + 1]] <- row
    }
  }
  
  # -------------------------
  # Categorical variables
  # -------------------------
  if (!is.null(cat_specs)) {
    for (v in names(cat_specs)) {
      x <- df[[v]]
      g <- df$.group
      
      overall <- overall_chisq_cramers_v(x, g, digits = digits_eff)
      
      for (lev in cat_specs[[v]]) {
        low_high <- low_high_or(x, g, lev, low_level, high_level, digits = digits_eff)
        
        key <- paste(v, lev, sep = "::")
        label <- if (!is.null(row_labels) && key %in% names(row_labels)) {
          row_labels[[key]]
        } else {
          paste0(v, ": ", lev, ", n (%)")
        }
        
        row <- tibble(
          Characteristic = label,
          Total = n_pct_str(x, lev, digits = digits_pct),
          Overall_P = fmt_p(overall$p),
          Omnibus_effect = paste0("V = ", overall$effect),
          Low_vs_high_P = fmt_p(low_high$p),
          Low_vs_high_effect = low_high$effect
        )
        
        for (gr in group_order) {
          row[[gr]] <- n_pct_str(df[df$.group == gr, v, drop = TRUE], lev, digits = digits_pct)
        }
        
        row <- row %>%
          dplyr::select(Characteristic, Total, all_of(group_order), Overall_P, Omnibus_effect,
                        Low_vs_high_P, Low_vs_high_effect)
        
        out[[length(out) + 1]] <- row
      }
    }
  }
  
  dplyr::bind_rows(out)
}

##5.2 output table----------------
mergeGPS_all_feature <- merge(mergeGPS_all, GPS_feature_all[, .(gene, `Gene length (bp)`, LOEUF_class, GES_class, ifspecific, ifhsd2)],
                              by = "gene", all.x = T)

head(mergeGPS_all_feature)

table_list <- lapply(list("pnadjust", "polygenicity", "samplematch", "mafbin1", "mafbin2", "mafbin3", "repro"), function(x){
  table_pn <- make_supp_table(data = as.data.frame(mergeGPS_all_feature),
                              group_var = paste0("pn_ld_class5_", x),
                              group_order = c("1st", "5th"),
                              low_level = "1st",
                              high_level = "5th",
                              cont_vars = c("Gene length (bp)"),
                              digits_cont = 0,
                              cat_specs = list("LOEUF_class" = c("Intolerant", "Tolerant"),
                                               "GES_class" = c("Essentialty", "Non_essentialty"),
                                               "ifspecific" = c("Specific", "Non_specific"),
                                               "ifhsd2" = c("Singletons", "Duplicates")))
  write.csv(table_pn[c(1,2,9,4,6),], file = sprintf("%s/table_pn_%s.csv", figure_file, x), row.names = FALSE)
  return(table_pn)
  
})

names(table_list) <- c("pnadjust", "polygenicity", "samplematch", "mafbin1", "mafbin2", "mafbin3", "repro")

table_list$pnadjust[c(1,2,9,4,6),]


