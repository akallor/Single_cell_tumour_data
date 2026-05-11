#R-version of the visualization script

# =============================================================================
# scRNA-seq Reporting Script — R / ggplot2
# Plots: UMAP (cluster, condition, marker gene), dot plot, composition bar chart
# Dataset: NF1 neurofibromatosis — Sus scrofa (Tumor1, Tumor2, NormalAdjacent)
# Dependencies: Seurat, ggplot2, dplyr, patchwork, viridis, RColorBrewer
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridis)
  library(RColorBrewer)
  library(scales)
})

# -----------------------------------------------------------------------------
# 0. Configuration — edit these paths and parameters
# -----------------------------------------------------------------------------

SEURAT_RDS   <- "integrated_seurat.rds"   # path to the saved Seurat object
OUT_DIR      <- "figures/"                # output directory for saved plots
MARKER_GENE  <- "NF1"                     # gene of interest
CONDITION_COL <- "condition"              # metadata column: "tumour" / "normal"
SAMPLE_COL    <- "sample"                # metadata column: sample name
CLUSTER_COL   <- "seurat_clusters"       # metadata column: cluster labels
N_TOP_MARKERS <- 5                        # top markers per cluster for dot plot
PLOT_WIDTH    <- 10                       # inches
PLOT_HEIGHT   <- 8

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load data
#    NOTE: If a Seurat object is not available, replace this block with: (e.g. Read10X → CreateSeuratObject → integrate).
# -----------------------------------------------------------------------------

cat("Loading Seurat object...\n")
sobj <- readRDS(SEURAT_RDS)

# Verify required metadata columns exist
stopifnot(
  CONDITION_COL %in% colnames(sobj@meta.data),
  CLUSTER_COL   %in% colnames(sobj@meta.data),
  SAMPLE_COL    %in% colnames(sobj@meta.data)
)

# Convenience: ensure condition is a factor with consistent level order
sobj@meta.data[[CONDITION_COL]] <- factor(
  sobj@meta.data[[CONDITION_COL]],
  levels = c("normal", "tumour")
)

n_clusters  <- nlevels(factor(sobj@meta.data[[CLUSTER_COL]]))
cluster_pal <- colorRampPalette(brewer.pal(12, "Paired"))(n_clusters)

# -----------------------------------------------------------------------------
# 2. Helper: shared UMAP theme
# -----------------------------------------------------------------------------

umap_theme <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      axis.title       = element_text(size = 9, colour = "grey50"),
      legend.title     = element_text(size = 10, face = "bold"),
      legend.text      = element_text(size = 9),
      plot.title       = element_text(size = 13, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = 10, colour = "grey40", hjust = 0),
      plot.caption     = element_text(size = 8, colour = "grey60"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "grey97", colour = NA)
    )
}

# -----------------------------------------------------------------------------
# 3. Extract UMAP embedding + metadata into a single data frame
# -----------------------------------------------------------------------------

cat("Extracting UMAP coordinates...\n")
umap_coords <- as.data.frame(Embeddings(sobj, reduction = "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

meta_df <- sobj@meta.data[, c(CLUSTER_COL, CONDITION_COL, SAMPLE_COL)]
plot_df  <- cbind(umap_coords, meta_df)

# Fetch normalised expression for the marker gene
if (MARKER_GENE %in% rownames(sobj)) {
  expr <- FetchData(sobj, vars = MARKER_GENE, layer = "data")
  plot_df[[MARKER_GENE]] <- expr[[MARKER_GENE]]
} else {
  warning(paste0(MARKER_GENE, " not found in feature names. Skipping expression UMAP."))
  plot_df[[MARKER_GENE]] <- NA_real_
}

# Sort so high-expressing cells are plotted on top
plot_df_expr <- plot_df[order(plot_df[[MARKER_GENE]], na.last = FALSE), ]

# -----------------------------------------------------------------------------
# 4a. UMAP — coloured by cluster
# -----------------------------------------------------------------------------

cat("Plotting UMAP by cluster...\n")

p_cluster <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2,
                                  colour = .data[[CLUSTER_COL]])) +
  geom_point(size = 0.4, alpha = 0.7) +
  scale_colour_manual(values = cluster_pal, name = "Cluster") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1),
                                ncol = ceiling(n_clusters / 12))) +
  labs(title    = "UMAP — cell clusters",
       subtitle = paste0(format(nrow(plot_df), big.mark = ","), " cells · Leiden clustering"),
       x = "UMAP 1", y = "UMAP 2") +
  umap_theme()

# -----------------------------------------------------------------------------
# 4b. UMAP — coloured by condition (tumour vs normal)
# -----------------------------------------------------------------------------

cat("Plotting UMAP by condition...\n")

condition_pal <- c("normal" = "#4393C3", "tumour" = "#D6604D")

p_condition <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2,
                                    colour = .data[[CONDITION_COL]])) +
  geom_point(size = 0.4, alpha = 0.6) +
  scale_colour_manual(values = condition_pal, name = "Condition") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title    = "UMAP — condition",
       subtitle = "Blue = NormalAdjacent · Red = Tumor1 + Tumor2",
       x = "UMAP 1", y = "UMAP 2") +
  umap_theme()

# -----------------------------------------------------------------------------
# 4c. UMAP — expression of marker gene (NF1)
# -----------------------------------------------------------------------------

cat(paste0("Plotting UMAP for ", MARKER_GENE, " expression...\n"))

p_gene <- ggplot(plot_df_expr, aes(x = UMAP_1, y = UMAP_2,
                                    colour = .data[[MARKER_GENE]])) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_colour_viridis_c(
    option   = "magma",
    name     = paste0(MARKER_GENE, "\n(log-norm)"),
    na.value = "grey85",
    limits   = c(0, NA)
  ) +
  labs(title    = paste0("UMAP — ", MARKER_GENE, " expression"),
       subtitle = "Normalised log expression · high expressors on top",
       x = "UMAP 1", y = "UMAP 2") +
  umap_theme()

# -----------------------------------------------------------------------------
# 4d. Combined UMAP panel (3 plots side by side)
# -----------------------------------------------------------------------------

p_umap_combined <- p_cluster | p_condition | p_gene
p_umap_combined <- p_umap_combined +
  plot_annotation(
    title   = "Single-cell RNA-seq · Sus scrofa NF1 neurofibromatosis",
    caption = paste0("Samples: NormalAdjacent, Tumor1, Tumor2 · Marker: ", MARKER_GENE),
    theme   = theme(
      plot.title   = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 8, colour = "grey60")
    )
  )

ggsave(file.path(OUT_DIR, "umap_combined.pdf"),
       p_umap_combined, width = PLOT_WIDTH * 3, height = PLOT_HEIGHT,
       units = "in", useDingbats = FALSE)
ggsave(file.path(OUT_DIR, "umap_combined.png"),
       p_umap_combined, width = PLOT_WIDTH * 3, height = PLOT_HEIGHT,
       dpi = 200)
cat("Saved umap_combined.pdf/.png\n")

# Save individual UMAPs too
ggsave(file.path(OUT_DIR, "umap_cluster.png"),   p_cluster,   width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 200)
ggsave(file.path(OUT_DIR, "umap_condition.png"), p_condition, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 200)
ggsave(file.path(OUT_DIR, "umap_NF1.png"),       p_gene,      width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = 200)

# -----------------------------------------------------------------------------
# 5. Find cluster marker genes (skip if already computed)
# -----------------------------------------------------------------------------

cat("Finding cluster markers (this may take a few minutes)...\n")

# PrepSCTFindMarkers if SCT was used; otherwise ensure RNA assay is active
DefaultAssay(sobj) <- "RNA"

all_markers <- FindAllMarkers(
  sobj,
  assay          = "RNA",
  only.pos       = TRUE,        # upregulated markers only
  min.pct        = 0.25,
  logfc.threshold = 0.25,
  test.use       = "wilcox",
  verbose        = FALSE
)

top_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = N_TOP_MARKERS) %>%
  ungroup()

write.csv(all_markers, file.path(OUT_DIR, "all_cluster_markers.csv"), row.names = FALSE)
write.csv(top_markers, file.path(OUT_DIR, "top_cluster_markers.csv"), row.names = FALSE)
cat("Saved marker tables.\n")

# -----------------------------------------------------------------------------
# 6. Dot plot — top markers per cluster
# -----------------------------------------------------------------------------

cat("Generating dot plot...\n")

genes_to_plot <- unique(top_markers$gene)

# Build dot plot data: average expression + percent expressed per cluster
dot_data <- do.call(rbind, lapply(unique(top_markers$cluster), function(cl) {
  cells_in_cl <- WhichCells(sobj, idents = cl)
  mat <- GetAssayData(sobj, assay = "RNA", layer = "data")[genes_to_plot, cells_in_cl, drop = FALSE]
  data.frame(
    gene        = genes_to_plot,
    cluster     = cl,
    avg_expr    = rowMeans(expm1(mat)),           # back-transform then average
    pct_express = rowMeans(mat > 0) * 100
  )
}))

# Scale avg_expr within each gene for cross-cluster comparison
dot_data <- dot_data %>%
  group_by(gene) %>%
  mutate(avg_expr_scaled = scales::rescale(avg_expr, to = c(0, 1))) %>%
  ungroup()

# Order genes by cluster of origin for readability
gene_order  <- top_markers %>% arrange(cluster, desc(avg_log2FC)) %>% pull(gene) %>% unique()
dot_data$gene    <- factor(dot_data$gene,    levels = rev(gene_order))
dot_data$cluster <- factor(dot_data$cluster, levels = sort(unique(as.numeric(as.character(dot_data$cluster)))))

p_dot <- ggplot(dot_data, aes(x = cluster, y = gene)) +
  geom_point(aes(size = pct_express, fill = avg_expr_scaled),
             shape = 21, colour = "grey30", stroke = 0.3) +
  scale_fill_gradient2(
    low      = "#EFF3FF",
    mid      = "#6BAED6",
    high     = "#08306B",
    midpoint = 0.5,
    name     = "Scaled\nexpression"
  ) +
  scale_size_continuous(
    range  = c(0.5, 8),
    name   = "% cells\nexpressing",
    breaks = c(10, 25, 50, 75)
  ) +
  labs(
    title    = paste0("Top ", N_TOP_MARKERS, " marker genes per cluster"),
    subtitle = "Dot size = % cells expressing · Colour = scaled mean expression",
    x        = "Cluster",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y      = element_text(size = 8, face = "italic"),
    axis.text.x      = element_text(size = 9),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(size = 13, face = "bold"),
    plot.subtitle    = element_text(size = 9,  colour = "grey40"),
    plot.background  = element_rect(fill = "white", colour = NA)
  )

dot_h <- max(6, length(genes_to_plot) * 0.28)
ggsave(file.path(OUT_DIR, "dotplot_top_markers.pdf"),
       p_dot, width = max(8, n_clusters * 0.9), height = dot_h,
       units = "in", useDingbats = FALSE)
ggsave(file.path(OUT_DIR, "dotplot_top_markers.png"),
       p_dot, width = max(8, n_clusters * 0.9), height = dot_h,
       dpi = 200)
cat("Saved dotplot_top_markers.pdf/.png\n")

# -----------------------------------------------------------------------------
# 7. Cluster composition bar chart — tumour vs normal proportions
# -----------------------------------------------------------------------------

cat("Generating cluster composition bar chart...\n")

comp_df <- sobj@meta.data %>%
  as.data.frame() %>%
  group_by(.data[[CLUSTER_COL]], .data[[CONDITION_COL]]) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(.data[[CLUSTER_COL]]) %>%
  mutate(
    total      = sum(n_cells),
    proportion = n_cells / total
  ) %>%
  ungroup()

comp_df[[CLUSTER_COL]]   <- factor(comp_df[[CLUSTER_COL]],
                                    levels = sort(unique(as.numeric(as.character(comp_df[[CLUSTER_COL]])))))
comp_df[[CONDITION_COL]] <- factor(comp_df[[CONDITION_COL]], levels = c("normal", "tumour"))

# Add cell count annotations at the top of each bar
total_labels <- comp_df %>%
  distinct(.data[[CLUSTER_COL]], total)

p_comp <- ggplot(comp_df,
                 aes(x    = .data[[CLUSTER_COL]],
                     y    = proportion,
                     fill = .data[[CONDITION_COL]])) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.4) +
  geom_text(
    data = total_labels,
    aes(x = .data[[CLUSTER_COL]], y = 1.02, label = scales::comma(total)),
    inherit.aes = FALSE,
    size = 2.8, colour = "grey40", hjust = 0.5
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  scale_fill_manual(
    values = condition_pal,
    name   = "Condition",
    labels = c("Normal adjacent", "Tumour (T1 + T2)")
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    title    = "Cluster composition — tumour vs normal",
    subtitle = "Stacked proportions per cluster · numbers above bars = total cell count",
    x        = "Cluster",
    y        = "Proportion of cells"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
    axis.ticks.x       = element_blank(),
    legend.position    = "top",
    legend.title       = element_text(size = 10, face = "bold"),
    plot.title         = element_text(size = 13, face = "bold"),
    plot.subtitle      = element_text(size = 9,  colour = "grey40"),
    plot.background    = element_rect(fill = "white", colour = NA)
  )

ggsave(file.path(OUT_DIR, "cluster_composition.pdf"),
       p_comp, width = max(8, n_clusters * 0.65), height = 6,
       units = "in", useDingbats = FALSE)
ggsave(file.path(OUT_DIR, "cluster_composition.png"),
       p_comp, width = max(8, n_clusters * 0.65), height = 6,
       dpi = 200)
cat("Saved cluster_composition.pdf/.png\n")

# -----------------------------------------------------------------------------
# 8. Full report panel — all plots on one page
# -----------------------------------------------------------------------------

cat("Assembling full report panel...\n")

p_report <- (p_cluster | p_condition | p_gene) /
            (p_dot | p_comp) +
  plot_annotation(
    title    = "scRNA-seq report · Sus scrofa NF1 neurofibromatosis",
    subtitle = paste0("Samples: NormalAdjacent, Tumor1, Tumor2 · Marker gene: ", MARKER_GENE,
                      " · Top ", N_TOP_MARKERS, " markers per cluster"),
    caption  = paste0("Generated: ", Sys.time()),
    theme = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey40"),
      plot.caption  = element_text(size = 8,  colour = "grey60")
    )
  ) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(file.path(OUT_DIR, "full_report.pdf"),
       p_report, width = 30, height = 22,
       units = "in", useDingbats = FALSE)
cat("Saved full_report.pdf\n")

cat("\nDone. All figures written to:", OUT_DIR, "\n")
cat("Files produced:\n")
cat("  umap_combined.pdf/.png\n")
cat("  umap_cluster.png\n")
cat("  umap_condition.png\n")
cat(paste0("  umap_", MARKER_GENE, ".png\n"))
cat("  dotplot_top_markers.pdf/.png\n")
cat("  cluster_composition.pdf/.png\n")
cat("  full_report.pdf\n")
cat("  all_cluster_markers.csv\n")
cat("  top_cluster_markers.csv\n")
