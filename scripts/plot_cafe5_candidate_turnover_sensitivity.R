#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: plot_cafe5_candidate_turnover_sensitivity.R <turnover.tsv> <output.pdf> <output.png>")
}

input_path <- args[[1]]
output_pdf <- args[[2]]
output_png <- args[[3]]
turnover <- read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "N/A", "NA"))

required <- c(
  "candidate_order", "orthogroup_id", "candidate_label", "availability_state",
  "branch_id", "taxon_label", "change", "branch_pvalue", "model", "lambda",
  "maximum_possible_lambda", "model_diagnostic_summary", "interpretation_status"
)
if (!all(required %in% names(turnover))) stop("Turnover table is missing required columns")
if (length(unique(turnover$orthogroup_id)) != 12) stop("Expected exactly 12 candidate orthogroups")
if (!all(turnover$interpretation_status == "sensitivity_only_model_not_accepted")) {
  stop("Every row must retain the rejected-model sensitivity warning")
}

available <- turnover[turnover$availability_state == "available", , drop = FALSE]
if (!nrow(available)) stop("No available turnover rows to plot")
available$change <- as.numeric(available$change)
available$branch_pvalue <- suppressWarnings(as.numeric(available$branch_pvalue))
if (any(!is.finite(available$change))) stop("Available rows require finite branch changes")

available$branch_label <- ifelse(
  !is.na(available$taxon_label) & nzchar(available$taxon_label),
  available$taxon_label,
  paste0("node ", sub(".*<([0-9]+)>", "\\1", available$branch_id))
)
terminal_order <- c("cula", "sese", "psme", "pabi", "plam", "plon", "ptae", "pden", "pmas", "ptab")
observed_branches <- unique(available$branch_label)
branch_order <- c(
  terminal_order[terminal_order %in% observed_branches],
  sort(setdiff(observed_branches, terminal_order))
)

candidate_meta <- unique(turnover[c("candidate_order", "orthogroup_id", "candidate_label")])
candidate_meta <- candidate_meta[order(as.integer(candidate_meta$candidate_order)), ]
candidate_meta$display_label <- ifelse(
  !is.na(candidate_meta$candidate_label) & nzchar(candidate_meta$candidate_label),
  paste(candidate_meta$orthogroup_id, candidate_meta$candidate_label, sep = " | "),
  candidate_meta$orthogroup_id
)

change_matrix <- matrix(
  NA_real_,
  nrow = nrow(candidate_meta),
  ncol = length(branch_order),
  dimnames = list(candidate_meta$orthogroup_id, branch_order)
)
p_matrix <- change_matrix
for (index in seq_len(nrow(available))) {
  row <- available[index, ]
  change_matrix[row$orthogroup_id, row$branch_label] <- row$change
  p_matrix[row$orthogroup_id, row$branch_label] <- row$branch_pvalue
}

palette <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(17)
fill_for_change <- function(value) {
  if (!is.finite(value)) return("#D9D9D9")
  index <- round((max(-8, min(8, value)) + 8) / 16 * 16) + 1
  palette[[index]]
}

draw_turnover <- function() {
  n_candidates <- nrow(candidate_meta)
  n_branches <- length(branch_order)
  layout(matrix(c(1, 2), nrow = 1), widths = c(5.2, 1.35))
  par(mar = c(8.0, 10.0, 5.6, 0.5), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, n_branches), ylim = c(0, n_candidates), xaxs = "i", yaxs = "i")

  for (candidate_index in seq_len(n_candidates)) {
    y <- n_candidates - candidate_index
    orthogroup <- candidate_meta$orthogroup_id[[candidate_index]]
    for (branch_index in seq_len(n_branches)) {
      x <- branch_index - 1
      change <- change_matrix[orthogroup, branch_index]
      rect(x, y, x + 1, y + 1, col = fill_for_change(change), border = "white", lwd = 1.0)
      pvalue <- p_matrix[orthogroup, branch_index]
      if (is.finite(pvalue) && pvalue > 0) {
        score <- min(4, -log10(pvalue))
        if (score > 0) points(x + 0.5, y + 0.5, pch = 21, bg = "black", col = "white", cex = 0.45 + 0.28 * score, lwd = 0.5)
      }
    }
  }

  axis(1, at = seq_len(n_branches) - 0.5, labels = branch_order, las = 2, tick = FALSE, cex.axis = 0.75)
  axis(2, at = n_candidates - seq_len(n_candidates) + 0.5, labels = candidate_meta$display_label, las = 1, tick = FALSE, cex.axis = 0.72)
  title(
    main = "Candidate gene-family turnover across the provisional conifer tree",
    sub = "Terminal taxa and internal branches",
    line = 3.5,
    cex.main = 1.25,
    font.main = 2
  )
  mtext(
    "Turnover sensitivity analysis; current model not accepted for final lineage calls",
    side = 3,
    line = 1.3,
    cex = 0.88,
    font = 2,
    col = "#9B2C2C"
  )

  par(mar = c(4.0, 0.6, 5.6, 1.0), family = "sans")
  plot.new()
  text(0, 0.96, "Figure key", adj = 0, font = 2, cex = 0.9)
  legend_y <- seq(0.86, 0.56, length.out = length(palette))
  for (index in seq_along(palette)) {
    rect(0, legend_y[[index]] - 0.012, 0.13, legend_y[[index]] + 0.012, col = rev(palette)[[index]], border = NA)
  }
  text(0.17, 0.86, "+8 expansion", adj = 0, cex = 0.72)
  text(0.17, 0.71, "0", adj = 0, cex = 0.72)
  text(0.17, 0.56, "-8 contraction", adj = 0, cex = 0.72)
  rect(0, 0.49, 0.13, 0.515, col = "#D9D9D9", border = NA)
  text(0.17, 0.502, "not available", adj = 0, cex = 0.72)
  text(0, 0.42, "Point size: nominal branch P", adj = 0, font = 2, cex = 0.72)
  point_y <- c(0.36, 0.31, 0.25)
  point_p <- c(0.1, 0.01, 0.001)
  for (index in seq_along(point_y)) {
    score <- -log10(point_p[[index]])
    points(0.06, point_y[[index]], pch = 21, bg = "black", col = "white", cex = 0.45 + 0.28 * score, lwd = 0.5)
    text(0.17, point_y[[index]], format(point_p[[index]]), adj = 0, cex = 0.72)
  }
  model <- unique(turnover[c("model", "lambda", "maximum_possible_lambda", "model_diagnostic_summary")])[1, ]
  note <- sprintf(
    "Base model: lambda=%s\nmaximum=%s\n\n%s\n\nNo symbol denotes an accepted lineage call.",
    model$lambda,
    model$maximum_possible_lambda,
    gsub(";", "\n", model$model_diagnostic_summary)
  )
  text(0, 0.16, note, adj = c(0, 1), cex = 0.66, col = "#444444")
  layout(1)
}

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
grDevices::cairo_pdf(output_pdf, width = 15, height = 9)
draw_turnover()
grDevices::dev.off()
grDevices::png(output_png, width = 4500, height = 2700, res = 300, bg = "white")
draw_turnover()
grDevices::dev.off()
message(sprintf("Wrote CAFE5 turnover sensitivity figure: %s, %s", output_pdf, output_png))
