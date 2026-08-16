suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(ragg)
})

set.seed(20260813)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(root, "11.单细胞", "GSE214695_PPARG_glycolysis_20260809")
object_file <- file.path(analysis_dir, "objects", "GSE214695_HC_UC_scored_seurat.rds")
score_file <- file.path(analysis_dir, "tables", "cell_level_pparg_glycolysis_scores.csv")
out_dir <- file.path(root, "submission_package", "figures", "standalone")
mirror_dir <- file.path(root, "图片1200dpi")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(mirror_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(object_file)
scores <- fread(score_file, select = c("cell", "cell_type_broad"))
umap <- as.data.table(Embeddings(obj, reduction = "umap"), keep.rownames = "cell")
setnames(umap, names(umap)[2:3], c("UMAP1", "UMAP2"))
plot_dt <- merge(umap, scores, by = "cell", all = FALSE)

broad_levels <- c(
  "Epithelial cells", "T/NK/ILC", "B cells", "Plasma cells",
  "Myeloid cells", "Neutrophils", "Mast cells", "Stromal cells",
  "Endothelial cells", "Enteric glia"
)
plot_dt[, cell_type_broad := factor(cell_type_broad, levels = broad_levels)]

type_palette <- c(
  "Epithelial cells" = "#0072B2", "T/NK/ILC" = "#E69F00",
  "B cells" = "#56B4E9", "Plasma cells" = "#CC79A7",
  "Myeloid cells" = "#D55E00", "Neutrophils" = "#F0C808",
  "Mast cells" = "#009E73", "Stromal cells" = "#7A5195",
  "Endothelial cells" = "#2F4B7C", "Enteric glia" = "#8C564B"
)

centroids <- plot_dt[, .(
  UMAP1 = median(UMAP1),
  UMAP2 = median(UMAP2)
), by = cell_type_broad]

p <- ggplot(plot_dt, aes(UMAP1, UMAP2, color = cell_type_broad)) +
  geom_point(size = 0.20, alpha = 0.76) +
  geom_text_repel(
    data = centroids,
    aes(label = cell_type_broad),
    color = "black",
    family = "Arial",
    fontface = "bold",
    size = 3.8,
    box.padding = 0.65,
    point.padding = 0.35,
    force = 12,
    force_pull = 0.35,
    min.segment.length = 0,
    segment.size = 0.45,
    max.overlaps = Inf,
    max.time = 8,
    seed = 20260813,
    show.legend = FALSE
  ) +
  scale_color_manual(values = type_palette, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.035)) +
  scale_y_continuous(expand = expansion(mult = 0.035)) +
  coord_equal(clip = "off") +
  labs(
    title = "Broad cell-type atlas",
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  ) +
  guides(color = guide_legend(
    ncol = 1,
    byrow = TRUE,
    override.aes = list(size = 3.6, alpha = 1)
  )) +
  theme_classic(base_family = "Arial", base_size = 15) +
  theme(
    plot.title = element_text(size = 19, face = "bold", hjust = 0.5, margin = margin(b = 4)),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13, color = "#243746"),
    axis.line = element_line(linewidth = 0.8, color = "#243746"),
    axis.ticks = element_line(linewidth = 0.7, color = "#243746"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.text = element_text(size = 12),
    legend.key.width = unit(0.50, "cm"),
    legend.key.height = unit(0.43, "cm"),
    legend.spacing.x = unit(0.04, "cm"),
    legend.spacing.y = unit(0, "cm"),
    legend.box.spacing = unit(0.1, "cm"),
    plot.margin = margin(7, 9, 3, 7)
  )

stem <- "Fig4A_broad_celltype_UMAP"
ggsave(file.path(out_dir, paste0(stem, ".pdf")), p,
       width = 7.0, height = 5.6, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(out_dir, paste0(stem, ".png")), p,
       width = 7.0, height = 5.6, units = "in", device = ragg::agg_png,
       dpi = 600, bg = "white")
ggsave(file.path(out_dir, paste0(stem, ".tif")), p,
       width = 7.0, height = 5.6, units = "in", device = ragg::agg_tiff,
       dpi = 1200, compression = "lzw", bg = "white")
ggsave(file.path(mirror_dir, paste0(stem, ".png")), p,
       width = 7.0, height = 5.6, units = "in", device = ragg::agg_png,
       dpi = 1200, bg = "white")

message("Wrote compact Fig4A PNG/PDF/TIFF without geometric distortion.")
