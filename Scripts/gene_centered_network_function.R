make_gene_centered_gsea_network <- function(
    gsea_obj,
    gene_score,
    top_n_pathways = 15,
    high_cutoff = quantile(gene_score, 0.8, na.rm = TRUE),
    top_n_genes = 15,
    min_pathways_per_gene = 2,
    pathway_label_width = 26,
    save_prefix = NULL,
    
    # layout control
    gene_radius = 0.42,
    pathway_radius = 0.90,
    pathway_label_offset = 0.12,
    start_angle = pi / 2,
    
    # pathway aesthetics
    pathway_text_size = 3.8,
    pathway_text_color = "#222222",
    pathway_text_lineheight = 0.88,
    pathway_node_fill = "#E9EEF6",
    pathway_node_color = "#6F86B6",
    pathway_node_stroke = 0.45,
    
    # gene aesthetics
    gene_text_size = 3.9,
    gene_node_color = "black",
    gene_node_stroke = 0.55,
    gene_label_offset = 0.06,
    
    # edge aesthetics
    edge_color = "grey78",
    edge_alpha = 0.18,
    edge_width_range = c(0.10, 0.60),
    
    # GPS color scale
    gps_low = "#FEE8C8",
    gps_high = "#D7301F",
    gps_legend_title = "GPS-N",
    
    # node size scale
    node_size_range = c(3.2, 9.5),
    size_legend_title = "Number of\nconnections",
    
    # figure control
    fig_width = 9.5,
    fig_height = 8.5,
    plot_margin = ggplot2::margin(25, 120, 25, 120)
) {
  
  # ----------------------------
  # 0. Required packages
  # ----------------------------
  require(data.table)
  require(ggplot2)
  require(ggraph)
  require(igraph)
  require(stringr)
  
  # ----------------------------
  # 1. Prepare GSEA result
  # ----------------------------
  res <- data.table::as.data.table(gsea_obj@result)
  
  if (!"core_enrichment" %in% colnames(res)) {
    stop("The GSEA result does not contain a core_enrichment column.")
  }
  
  if (is.null(names(gene_score))) {
    stop("gene_score must be a named vector. names(gene_score) should match core_enrichment gene IDs.")
  }
  
  res <- res[
    !is.na(p.adjust) &
      !is.na(NES) &
      !is.na(core_enrichment) &
      core_enrichment != ""
  ][order(p.adjust)]
  
  if (nrow(res) == 0) {
    stop("No valid enriched pathways found.")
  }
  
  res_top <- res[1:min(top_n_pathways, .N)]
  
  res_top[, pathway_node := paste0("PW__", ID)]
  res_top[, pathway_label := stringr::str_wrap(Description, width = pathway_label_width)]
  
  # ----------------------------
  # 2. Extract core-enrichment genes
  # ----------------------------
  edge_dt <- res_top[
    ,
    {
      genes <- strsplit(core_enrichment, "/", fixed = TRUE)[[1]]
      genes <- genes[genes %in% names(gene_score)]
      .(gene = genes)
    },
    by = .(
      ID,
      Description,
      pathway_node,
      pathway_label,
      NES,
      p.adjust
    )
  ]
  
  if (nrow(edge_dt) == 0) {
    stop("No overlap between core_enrichment genes and names(gene_score). Please check gene ID type.")
  }
  
  edge_dt[, gene_score := as.numeric(gene_score[gene])]
  
  # ----------------------------
  # 3. Keep high-pleiotropy core genes only
  # ----------------------------
  edge_dt <- edge_dt[!is.na(gene_score) & gene_score >= high_cutoff]
  
  if (nrow(edge_dt) == 0) {
    stop("No high-pleiotropy core genes remained after filtering. Try lowering high_cutoff.")
  }
  
  # avoid Inf edge weights if p.adjust = 0
  edge_dt[, p_adj_for_weight := pmax(p.adjust, .Machine$double.xmin)]
  
  # ----------------------------
  # 4. Select recurrent high-pleiotropy genes
  # ----------------------------
  gene_stats <- edge_dt[
    ,
    .(
      n_pathways = data.table::uniqueN(pathway_node),
      contribution_score = sum(-log10(p_adj_for_weight), na.rm = TRUE),
      mean_NES = mean(NES, na.rm = TRUE),
      GPS = unique(gene_score)[1]
    ),
    by = gene
  ]
  
  gene_stats <- gene_stats[n_pathways >= min_pathways_per_gene]
  
  if (nrow(gene_stats) == 0) {
    stop("No gene is connected to at least min_pathways_per_gene pathways. Try setting min_pathways_per_gene = 1.")
  }
  
  data.table::setorder(gene_stats, -n_pathways, -contribution_score, -GPS)
  
  selected_genes <- gene_stats[1:min(top_n_genes, .N), gene]
  
  edge_dt <- edge_dt[gene %in% selected_genes]
  gene_stats <- gene_stats[gene %in% selected_genes]
  
  # ----------------------------
  # 5. Build edge table
  # ----------------------------
  edges <- edge_dt[
    ,
    .(
      from = paste0("G__", gene),
      to = pathway_node,
      weight = -log10(p_adj_for_weight),
      gene = gene,
      pathway = Description,
      NES = NES,
      p.adjust = p.adjust,
      gene_score = gene_score
    )
  ]
  
  if (nrow(edges) == 0) {
    stop("No edges remained after gene selection.")
  }
  
  # ----------------------------
  # 6. Build node table
  # ----------------------------
  gene_nodes <- gene_stats[
    ,
    .(
      name = paste0("G__", gene),
      label = gene,
      type = "Gene",
      GPS = GPS,
      node_size = n_pathways,
      contribution_score = contribution_score,
      NES = NA_real_,
      p.adjust = NA_real_,
      label_x = NA_real_,
      label_y = NA_real_,
      label_hjust = NA_real_
    )
  ]
  
  pathway_stats <- edge_dt[
    ,
    .(
      n_high_core_genes = data.table::uniqueN(gene),
      NES = unique(NES)[1],
      p.adjust = unique(p.adjust)[1],
      pathway_label = unique(pathway_label)[1],
      Description = unique(Description)[1]
    ),
    by = pathway_node
  ]
  
  pathway_nodes <- pathway_stats[
    ,
    .(
      name = pathway_node,
      label = pathway_label,
      type = "Pathway",
      GPS = NA_real_,
      node_size = n_high_core_genes,
      contribution_score = NA_real_,
      NES = NES,
      p.adjust = p.adjust,
      label_x = NA_real_,
      label_y = NA_real_,
      label_hjust = NA_real_
    )
  ]
  
  nodes <- data.table::rbindlist(list(gene_nodes, pathway_nodes), fill = TRUE)
  
  # ----------------------------
  # 7. Manual gene-centered layout
  #    Genes: inner circle
  #    Pathways: outer circle
  #    Pathway labels: outside pathway nodes
  # ----------------------------
  gene_nodes2 <- nodes[type == "Gene"]
  pathway_nodes2 <- nodes[type == "Pathway"][order(p.adjust)]
  
  n_gene <- nrow(gene_nodes2)
  n_pathway <- nrow(pathway_nodes2)
  
  gene_nodes2[, angle := seq(
    start_angle,
    start_angle + 2 * pi,
    length.out = n_gene + 1
  )[1:n_gene]]
  
  gene_nodes2[, x := gene_radius * cos(angle)]
  gene_nodes2[, y := gene_radius * sin(angle)]
  
  # gene label positions: slightly outside gene nodes
  gene_nodes2[, label_x := (gene_radius + gene_label_offset) * cos(angle)]
  gene_nodes2[, label_y := (gene_radius + gene_label_offset) * sin(angle)]
  
  gene_nodes2[, label_hjust := ifelse(cos(angle) >= 0, 0, 1)]
  gene_nodes2[abs(cos(angle)) < 0.15, label_hjust := 0.5]
  
  pathway_nodes2[, angle := seq(
    start_angle,
    start_angle + 2 * pi,
    length.out = n_pathway + 1
  )[1:n_pathway]]
  
  pathway_nodes2[, x := pathway_radius * cos(angle)]
  pathway_nodes2[, y := pathway_radius * sin(angle)]
  
  # Label positions are pushed outward from pathway nodes
  pathway_nodes2[, label_x := (pathway_radius + pathway_label_offset) * cos(angle)]
  pathway_nodes2[, label_y := (pathway_radius + pathway_label_offset) * sin(angle)]
  
  # Align labels based on left/right side
  pathway_nodes2[, label_hjust := ifelse(cos(angle) >= 0, 0, 1)]
  
  # Top/bottom labels are centered
  pathway_nodes2[abs(cos(angle)) < 0.15, label_hjust := 0.5]
  
  nodes_layout <- data.table::rbindlist(list(gene_nodes2, pathway_nodes2), fill = TRUE)
  
  # ----------------------------
  # 8. Make graph
  # ----------------------------
  g <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes_layout,
    directed = FALSE
  )
  table(igraph::V(g)$type)
  
  layout_tbl <- ggraph::create_layout(
    g,
    layout = "manual",
    x = igraph::V(g)$x,
    y = igraph::V(g)$y
  )
  
  # only use data.frame version for pathway labels
  layout_data <- as.data.frame(layout_tbl)
  
  if (!"type" %in% colnames(layout_data)) {
    stop("The layout object does not contain a 'type' column. Please check vertex attributes in nodes_layout.")
  }
  
  pathway_label_df <- layout_data[layout_data$type == "Pathway", ]
  gene_label_df <- layout_data[layout_data$type == "Gene", ]
  
  # ----------------------------
  # 9. Plot
  # ----------------------------
  p <- ggraph::ggraph(layout_tbl) +
    
    # edges
    ggraph::geom_edge_link(
      aes(width = weight),
      color = edge_color,
      alpha = edge_alpha,
      lineend = "round",
      show.legend = FALSE
    ) +
    ggraph::scale_edge_width(range = edge_width_range) +
    
    # pathway nodes
    ggraph::geom_node_point(
      data = function(x) x[x$type == "Pathway", ],
      aes(size = node_size),
      shape = 22,
      fill = pathway_node_fill,
      color = pathway_node_color,
      stroke = pathway_node_stroke
    ) +
    
    # gene nodes
    ggraph::geom_node_point(
      data = function(x) x[x$type == "Gene", ],
      aes(size = node_size, fill = GPS),
      shape = 21,
      color = gene_node_color,
      stroke = gene_node_stroke
    ) +
    
    # gene labels
    ggplot2::geom_text(
      data = gene_label_df,
      aes(
        x = label_x,
        y = label_y,
        label = label,
        hjust = label_hjust
      ),
      size = gene_text_size,
      fontface = "bold",
      color = "black",
      inherit.aes = FALSE
    ) +
    
    # pathway labels: separate coordinates, no overlap with squares
    ggplot2::geom_text(
      data = pathway_label_df,
      aes(
        x = label_x,
        y = label_y,
        label = label,
        hjust = label_hjust
      ),
      size = pathway_text_size,
      color = pathway_text_color,
      lineheight = pathway_text_lineheight,
      inherit.aes = FALSE,
      fontface = "bold"
    ) +
    
    ggplot2::scale_fill_gradient(
      low = gps_low,
      high = gps_high,
      name = gps_legend_title
    ) +
    
    ggplot2::scale_size_continuous(
      range = node_size_range,
      name = size_legend_title
    ) +
    
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      plot.margin = plot_margin,
      legend.box = "vertical",
      legend.box.just = "left",
      legend.spacing.y = grid::unit(10, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 60),
      legend.box.margin = ggplot2::margin(0, 0, 0, 50),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(face = "bold")
    )
  
  # ----------------------------
  # 10. Optional saving
  # ----------------------------
  if (!is.null(save_prefix)) {
    
    out_dir <- dirname(save_prefix)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    ggplot2::ggsave(
      paste0(save_prefix, "_gene_centered_network.pdf"),
      p,
      width = fig_width,
      height = fig_height,
      device = grDevices::cairo_pdf
    )
    
    data.table::fwrite(
      edges,
      paste0(save_prefix, "_gene_pathway_edges.tsv"),
      sep = "\t"
    )
    
    data.table::fwrite(
      nodes_layout,
      paste0(save_prefix, "_gene_pathway_nodes.tsv"),
      sep = "\t"
    )
    
    data.table::fwrite(
      gene_stats,
      paste0(save_prefix, "_gene_summary.tsv"),
      sep = "\t"
    )
    
    data.table::fwrite(
      pathway_stats,
      paste0(save_prefix, "_pathway_summary.tsv"),
      sep = "\t"
    )
  }
  
  return(
    list(
      plot = p,
      edges = edges,
      nodes = nodes_layout,
      gene_stats = gene_stats,
      pathway_stats = pathway_stats,
      high_cutoff = high_cutoff
    )
  )
}