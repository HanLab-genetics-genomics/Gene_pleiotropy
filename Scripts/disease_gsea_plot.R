library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(scales)
library(AnnotationDbi)
library(org.Hs.eg.db)

prepare_disease_gsea_contribution_data <- function(
    gsea_obj,
    gene_score,
    high_cutoff = 20,
    top_n_terms = 15,
    top_n_shared_genes = 5,
    gene_id_type = "ENTREZID",
    OrgDb = org.Hs.eg.db
) {
  
  res <- as.data.table(gsea_obj@result)
  
  if (!"core_enrichment" %in% colnames(res)) {
    stop("The GSEA result does not contain a core_enrichment column.")
  }
  
  if (is.null(names(gene_score))) {
    stop("gene_score must be a named vector.")
  }
  
  # ----------------------------
  # 1. keep valid GSEA results and rank by p.adjust
  # ----------------------------
  res <- res[
    !is.na(p.adjust) &
      !is.na(NES) &
      !is.na(core_enrichment) &
      core_enrichment != ""
  ]
  
  # first select top disease terms by adjusted P value
  setorder(res, p.adjust)
  res_top <- copy(res[1:min(top_n_terms, .N)])
  
  # then display them ordered by NES
  setorder(res_top, -NES)
  
  # ----------------------------
  # 2. define high-pleiotropy genes
  # ----------------------------
  high_genes <- names(gene_score)[
    !is.na(gene_score) & gene_score >= high_cutoff
  ]
  
  # ----------------------------
  # 3. parse core_enrichment genes
  # ----------------------------
  res_top[, core_gene_vec := lapply(
    strsplit(core_enrichment, "/", fixed = TRUE),
    unique
  )]
  
  res_top[, core_gene_scored := lapply(core_gene_vec, function(g) {
    g[g %in% names(gene_score)]
  })]
  
  res_top[, high_core_gene_vec := lapply(core_gene_scored, function(g) {
    g[g %in% high_genes]
  })]
  
  # ----------------------------
  # 4. compute high-pleiotropy core gene proportion
  # ----------------------------
  res_top[, `:=`(
    n_core = lengths(core_gene_vec),
    n_core_with_score = lengths(core_gene_scored),
    n_high_core = lengths(high_core_gene_vec)
  )]
  
  res_top[, prop_high_core := fifelse(
    n_core_with_score > 0,
    n_high_core / n_core_with_score,
    NA_real_
  )]
  
  res_top[, core_ratio_label := paste0(n_high_core, "/", n_core_with_score)]
  
  # disease order: top enriched term at the top
  res_top[, y := rev(seq_len(.N))]
  res_top[, disease_label := Description]
  res_top[, disease_label_wrap := stringr::str_wrap(Description, width = 28)]
  
  # ----------------------------
  # 5. make disease-gene link table
  # ----------------------------
  link_dt <- rbindlist(
    lapply(seq_len(nrow(res_top)), function(i) {
      g <- res_top$high_core_gene_vec[[i]]
      if (length(g) == 0) return(NULL)
      
      data.table(
        Description = res_top$Description[i],
        disease_label = res_top$disease_label[i],
        disease_label_wrap = res_top$disease_label_wrap[i],
        y = res_top$y[i],
        gene = g,
        p.adjust = res_top$p.adjust[i],
        NES = res_top$NES[i]
      )
    }),
    fill = TRUE
  )
  
  if (is.null(link_dt) || nrow(link_dt) == 0) {
    stop("No high-pleiotropy core genes were found in the selected top disease gene sets.")
  }
  
  # ----------------------------
  # 6. select recurrent shared high-pleiotropy genes
  # ----------------------------
  gene_stats <- link_dt[, .(
    n_diseases = uniqueN(Description),
    contribution_score = sum(-log10(p.adjust), na.rm = TRUE)
  ), by = gene]
  
  gene_stats[, GPS := as.numeric(gene_score[gene])]
  
  setorder(gene_stats, -n_diseases, -contribution_score, -GPS)
  gene_stats <- gene_stats[1:min(top_n_shared_genes, .N)]
  
  # ----------------------------
  # 7. convert gene ID to SYMBOL for display
  # ----------------------------
  gene_map <- AnnotationDbi::select(
    OrgDb,
    keys = unique(gene_stats$gene),
    keytype = gene_id_type,
    columns = c("SYMBOL")
  )
  
  gene_map <- as.data.table(gene_map)
  setnames(gene_map, gene_id_type, "gene")
  gene_map <- gene_map[!is.na(SYMBOL)]
  gene_map <- gene_map[!duplicated(gene)]
  
  gene_stats <- merge(
    gene_stats,
    gene_map[, .(gene, SYMBOL)],
    by = "gene",
    all.x = TRUE
  )
  
  gene_stats[, gene_label := fifelse(
    !is.na(SYMBOL) & SYMBOL != "",
    SYMBOL,
    gene
  )]
  
  # keep only selected shared genes
  link_dt <- link_dt[gene %in% gene_stats$gene]
  
  # assign positions for genes in the right panel
  gene_stats[, y_gene := seq(
    from = max(res_top$y),
    to = min(res_top$y),
    length.out = .N
  )]
  gene_stats[, x_gene := 1]
  
  link_dt <- merge(
    link_dt,
    gene_stats[, .(gene, gene_label, y_gene, x_gene, GPS, n_diseases)],
    by = "gene",
    all.x = TRUE
  )
  
  return(list(
    term_summary = res_top,
    link_dt = link_dt,
    gene_stats = gene_stats
  ))
}

plot_disease_gsea_contribution <- function(
    plot_data,
    nes_point_size = 4.2,
    disease_text_size = 11,
    bar_text_size = 3.6,
    gene_text_size = 3.6,
    line_alpha = 0.65
) {
  
  term_summary <- copy(plot_data$term_summary)
  link_dt <- copy(plot_data$link_dt)
  gene_stats <- copy(plot_data$gene_stats)
  
  n_terms <- nrow(term_summary)
  y_limits <- c(0.5, n_terms + 0.5)
  y_breaks <- term_summary$y
  
  # lollipop baseline
  x_min <- min(term_summary$NES, na.rm = TRUE) - 0.05
  
  # for middle bar plot
  term_summary[, bar_ymin := y - 0.32]
  term_summary[, bar_ymax := y + 0.32]
  
  max_prop <- max(term_summary$prop_high_core, na.rm = TRUE)
  if (!is.finite(max_prop)) max_prop <- 1
  
  # ----------------------------
  # Panel 1: disease GSEA NES + p.adjust
  # ----------------------------
  p_left <- ggplot(term_summary, aes(x = NES, y = y)) +
    geom_segment(
      aes(x = x_min, xend = NES, y = y, yend = y),
      color = "grey85",
      linewidth = 0.8
    ) +
    geom_point(
      aes(color = p.adjust),
      size = nes_point_size
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = term_summary$disease_label_wrap,
      limits = y_limits,
      expand = c(0, 0)
    ) +
    scale_color_gradient(
      low  = "#B8572A",
      high = "#E9D8CC",
      #low = "#D07C45",
      #high = "#8EBC66",
      name = "p.adjust"
    ) +
    labs(x = "Normalized enrichment scores", y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_text(size = disease_text_size, face = "bold"),
      axis.title.x = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(5, 10, 5, 5)
    )
  
  # ----------------------------
  # Panel 2: proportion of high-pleiotropy core genes
  # use geom_rect to ensure horizontal bars are shown
  # ----------------------------
  p_mid <- ggplot(term_summary) +
    geom_rect(
      aes(
        xmin = 0,
        xmax = prop_high_core,
        ymin = bar_ymin,
        ymax = bar_ymax,
        fill = prop_high_core
      ),
      color = NA
    ) +
    geom_text(
      aes(
        x = prop_high_core + max_prop * 0.04,
        y = y,
        label = core_ratio_label
      ),
      hjust = 0,
      size = bar_text_size,
      fontface = "bold"
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = NULL,
      limits = y_limits,
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 1, scale = 100),
      limits = c(0, max_prop * 1.35),
      expand = c(0, 0)
    ) +
    scale_fill_gradient(
      low = "#DCE8F5",
      high = "#4A76A8",
      name = "High-GPS\ncore proportion"
    ) +
    labs(x = "Proportion of high-pleiotropy\ncore genes (%)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(face = "bold"),
      legend.position = "none",
      plot.margin = margin(5, 15, 5, 5)
    )
  
  # ----------------------------
  # Panel 3: shared recurrent high-pleiotropy genes + links
  # ----------------------------
  p_right <- ggplot() +
    geom_segment(
      data = link_dt,
      aes(
        x = 0.05, y = y,
        xend = 0.95, yend = y_gene
      ),
      color = "grey80",
      alpha = line_alpha,
      linewidth = 0.5
    ) +
    geom_point(
      data = term_summary,
      aes(x = 0.05, y = y),
      size = 1.6,
      color = "grey65"
    ) +
    geom_point(
      data = gene_stats,
      aes(x = 1, y = y_gene),
      size = 3.8,
      shape = 21,
      stroke = 0.5,
      fill = "#D7301F",
      color = "grey30"
    ) +
    geom_text(
      data = gene_stats,
      aes(x = 1.08, y = y_gene, label = gene_label),
      hjust = 0,
      size = gene_text_size,
      fontface = "bold"
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = NULL,
      limits = y_limits,
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      limits = c(0, 1.70),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_void() +
    theme(
      plot.margin = margin(5, 55, 5, 10)
    )
  
  # ----------------------------
  # Combine
  # ----------------------------
  p <- p_left + p_mid + p_right +
    patchwork::plot_layout(widths = c(1.45, 1.15, 1.35))
  
  return(p)
}