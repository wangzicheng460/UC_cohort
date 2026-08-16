options(stringsAsFactors = FALSE, warn = 1)

ml_figure_palette <- function(n = 256L) {
  grDevices::hcl.colors(n, palette = "Cividis")
}

read_ml_figure_data <- function(out_dir) {
  auc_path <- file.path(out_dir, "workflow_auc_by_outer_cohort.csv")
  summary_path <- file.path(out_dir, "workflow_summary.csv")
  stability_path <- file.path(out_dir, "gene_stability_summary.csv")
  stopifnot(file.exists(auc_path), file.exists(summary_path), file.exists(stability_path))

  auc_long <- read.csv(auc_path, check.names = FALSE)
  workflow_summary <- read.csv(summary_path, check.names = FALSE)
  gene_stability <- read.csv(stability_path, check.names = FALSE)

  required_auc <- c("workflow", "heldout_cohort", "auc")
  required_summary <- c("workflow", "mean_auc")
  required_stability <- c("gene", "selection_frequency", "stable_core")
  stopifnot(all(required_auc %in% names(auc_long)))
  stopifnot(all(required_summary %in% names(workflow_summary)))
  stopifnot(all(required_stability %in% names(gene_stability)))

  auc_long$auc <- as.numeric(auc_long$auc)
  gene_stability$selection_frequency <- as.numeric(gene_stability$selection_frequency)
  gene_stability$stable_core <- as.logical(gene_stability$stable_core)
  if (anyNA(auc_long$auc)) stop("Missing AUC values in workflow_auc_by_outer_cohort.csv")
  if (anyNA(gene_stability$selection_frequency) || anyNA(gene_stability$stable_core)) {
    stop("Missing gene-stability values in gene_stability_summary.csv")
  }

  list(
    auc_long = auc_long,
    workflow_summary = workflow_summary,
    gene_stability = gene_stability
  )
}

build_auc_matrix <- function(auc_long, workflow_summary, n_workflows) {
  workflows <- head(workflow_summary$workflow, n_workflows)
  preferred_order <- c("GSE73661", "GSE75214", "GSE87466", "GSE107499")
  cohorts <- c(
    intersect(preferred_order, unique(auc_long$heldout_cohort)),
    setdiff(unique(auc_long$heldout_cohort), preferred_order)
  )
  auc_mat <- matrix(
    NA_real_, nrow = length(workflows), ncol = length(cohorts),
    dimnames = list(workflows, cohorts)
  )
  for (i in seq_len(nrow(auc_long))) {
    workflow <- auc_long$workflow[i]
    cohort <- auc_long$heldout_cohort[i]
    if (workflow %in% workflows && cohort %in% cohorts) {
      auc_mat[workflow, cohort] <- auc_long$auc[i]
    }
  }
  if (anyNA(auc_mat)) stop("Incomplete workflow-by-cohort AUC matrix")
  auc_mat
}

draw_auc_legend <- function(palette, zlim) {
  graphics::par(mar = c(5.2, 0.8, 5.0, 4.2))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = zlim)
  breaks <- seq(zlim[1], zlim[2], length.out = length(palette) + 1L)
  graphics::rect(
    xleft = 0.18, ybottom = head(breaks, -1L),
    xright = 0.48, ytop = tail(breaks, -1L),
    col = palette, border = NA
  )
  legend_ticks <- seq(zlim[1], zlim[2], by = 0.01)
  graphics::axis(
    side = 4, at = legend_ticks, labels = sprintf("%.2f", legend_ticks),
    pos = 0.48, las = 1, tck = -0.12, cex.axis = 0.72, mgp = c(1.2, 0.45, 0)
  )
  graphics::rect(0.18, zlim[1], 0.48, zlim[2], border = "#3F4850", lwd = 0.8)
  graphics::mtext("Held out AUC", side = 3, line = 1.1, adj = 0.15, font = 2, cex = 0.88)
}

plot_auc_heatmap <- function(auc_mat, output_path, annotate_cells) {
  n_workflows <- nrow(auc_mat)
  zlim <- c(0.93, 1.00)
  if (min(auc_mat) < zlim[1] || max(auc_mat) > zlim[2]) {
    stop("AUC values fall outside the approved heatmap display range")
  }
  palette <- ml_figure_palette()
  plot_height <- if (n_workflows <= 20L) 8.4 else 12.2
  grDevices::pdf(output_path, width = 9.4, height = plot_height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(7.3, 1.1))

  graphics::par(mar = c(5.2, 17.5, 5.0, 0.6), xpd = FALSE)
  display_mat <- auc_mat[nrow(auc_mat):1L, , drop = FALSE]
  graphics::image(
    x = seq_len(ncol(display_mat)), y = seq_len(nrow(display_mat)),
    z = t(display_mat), col = palette, zlim = zlim,
    axes = FALSE, xlab = "", ylab = "", useRaster = TRUE
  )
  graphics::abline(
    v = seq(0.5, ncol(display_mat) + 0.5, by = 1),
    h = seq(0.5, nrow(display_mat) + 0.5, by = 1),
    col = grDevices::adjustcolor("white", alpha.f = 0.68), lwd = 0.55
  )
  graphics::axis(
    1, at = seq_len(ncol(display_mat)), labels = colnames(display_mat),
    las = 2, tick = FALSE, line = -0.5, cex.axis = 0.88
  )
  graphics::axis(
    2, at = seq_len(nrow(display_mat)), labels = rownames(display_mat),
    las = 1, tick = FALSE, line = -0.35,
    cex.axis = if (n_workflows <= 20L) 0.72 else 0.56
  )
  graphics::box(col = "#3F4850", lwd = 0.8)

  if (annotate_cells) {
    for (i in seq_len(nrow(display_mat))) {
      for (j in seq_len(ncol(display_mat))) {
        value <- display_mat[i, j]
        label_col <- if (value < 0.965) "white" else "#17202A"
        graphics::text(j, i, sprintf("%.3f", value), col = label_col, cex = 0.57, font = 2)
      }
    }
  }

  graphics::mtext(
    "Leakage controlled leave one cohort out performance",
    side = 3, line = 2.6, adj = 0, font = 2, cex = 1.18
  )
  graphics::mtext(
    "AUC in the held out discovery cohort",
    side = 3, line = 1.15, adj = 0, col = "#52606D", cex = 0.88
  )
  draw_auc_legend(palette, zlim)
  invisible(output_path)
}

plot_gene_stability <- function(gene_stability, output_path) {
  stable_colour <- "#D55E00"
  other_colour <- "#7A8A99"
  threshold_colour <- "#3F4850"
  gene_stability <- gene_stability[order(-gene_stability$selection_frequency), , drop = FALSE]
  stable_n <- sum(gene_stability$stable_core)

  grDevices::pdf(output_path, width = 9.2, height = 7.4, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(1.0, 6.4))

  graphics::par(mar = c(0.1, 7.4, 0.2, 2.0))
  graphics::plot.new()
  graphics::text(
    0, 0.84, "Consensus gene selection stability",
    adj = c(0, 0.5), font = 2, cex = 1.24, xpd = NA
  )
  graphics::legend(
    x = 0, y = 0.48,
    legend = c(
      sprintf("Stable core (n = %d)", stable_n),
      "Other candidates", "Selection threshold = 0.50"
    ),
    col = c(stable_colour, other_colour, threshold_colour),
    pch = c(15, 15, NA), pt.cex = 1.35,
    lty = c(NA, NA, 2), lwd = c(NA, NA, 1.2),
    ncol = 3, bty = "n", cex = 0.78, x.intersp = 0.7, y.intersp = 0.9
  )

  display <- gene_stability[nrow(gene_stability):1L, , drop = FALSE]
  bar_colours <- ifelse(display$stable_core, stable_colour, other_colour)
  graphics::par(mar = c(5.1, 7.4, 0.4, 2.0))
  midpoints <- graphics::barplot(
    display$selection_frequency,
    names.arg = display$gene,
    horiz = TRUE, las = 1, cex.names = 0.72,
    xlim = c(0, 1.0), xaxt = "n",
    xlab = "Selection frequency across selectors and outer LOCO folds",
    col = bar_colours, border = NA
  )
  graphics::axis(1, at = seq(0, 1, by = 0.1), labels = sprintf("%.1f", seq(0, 1, by = 0.1)),
                 cex.axis = 0.78)
  graphics::abline(v = 0.5, lty = 2, lwd = 1.2, col = threshold_colour)
  graphics::text(
    x = display$selection_frequency + 0.012, y = midpoints,
    labels = sprintf("%.3f", display$selection_frequency),
    adj = 0, cex = 0.62, col = "#25313C"
  )
  invisible(output_path)
}

write_ml_figures <- function(out_dir) {
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  data <- read_ml_figure_data(out_dir)
  stable_genes <- data$gene_stability$gene[data$gene_stability$stable_core]
  expected_stable <- c("PDE6A", "PPARG", "ADH6", "LCN2", "VLDLR", "SLC2A3", "TRPM6", "KDELR3")
  if (!setequal(stable_genes, expected_stable)) {
    stop("Stable-core membership differs from the approved eight-gene set")
  }

  auc_top20 <- build_auc_matrix(data$auc_long, data$workflow_summary, 20L)
  auc_top40 <- build_auc_matrix(data$auc_long, data$workflow_summary, 40L)
  plot_auc_heatmap(
    auc_top20, file.path(out_dir, "workflow_auc_heatmap.pdf"), annotate_cells = TRUE
  )
  plot_auc_heatmap(
    auc_top40, file.path(out_dir, "workflow_auc_heatmap_top40_supplement.pdf"), annotate_cells = FALSE
  )
  plot_gene_stability(
    data$gene_stability, file.path(out_dir, "gene_stability.pdf")
  )

  message("Wrote publication figures to: ", out_dir)
  invisible(list(auc_top20 = auc_top20, auc_top40 = auc_top40, gene_stability = data$gene_stability))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_dir <- if (length(args)) args[[1L]] else getwd()
  write_ml_figures(output_dir)
}
