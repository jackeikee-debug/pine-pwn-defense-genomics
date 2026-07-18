#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: plot_fig5_candidate_coevolution.R <summary.tsv> <top.tsv> <output.pdf>")
}

summary_path <- args[[1]]
top_path <- args[[2]]
output_path <- args[[3]]

summary_df <- read.delim(summary_path, stringsAsFactors = FALSE, check.names = FALSE)
top_df <- read.delim(top_path, stringsAsFactors = FALSE, check.names = FALSE)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

format_label <- function(value) {
  gsub("_", " ", value, fixed = TRUE)
}

bias_colors <- c(
  susceptible_enriched = "#c95f4a",
  tolerant_enriched = "#3d8f6f",
  susceptible_only = "#9e3f32",
  tolerant_only = "#23694f",
  balanced_or_mixed = "#777777"
)

pdf(output_path, width = 10.5, height = 7.2, onefile = TRUE)

par(mar = c(1.5, 1.5, 3.2, 1.5), family = "sans")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
title("Fig5. Candidate nematode effector classes linked to pine defense modules", cex.main = 1.15)

modules <- sort(unique(summary_df$module_id))
effectors <- sort(unique(summary_df$effector_class))
module_y <- setNames(seq(0.84, 0.18, length.out = length(modules)), modules)
effector_y <- setNames(seq(0.84, 0.18, length.out = length(effectors)), effectors)
left_x <- 0.22
right_x <- 0.78

for (i in seq_len(nrow(summary_df))) {
  row <- summary_df[i, ]
  y1 <- module_y[[row$module_id]]
  y2 <- effector_y[[row$effector_class]]
  count <- as.numeric(row$orthogroup_count)
  color <- bias_colors[[row$host_bias_direction]]
  if (is.null(color) || is.na(color)) {
    color <- "#777777"
  }
  offset <- ifelse(row$host_bias_direction == "susceptible_enriched", 0.012, -0.012)
  segments(
    left_x + 0.07,
    y1 + offset,
    right_x - 0.07,
    y2 + offset,
    col = color,
    lwd = 1.2 + sqrt(count) / 1.6
  )
  text(
    (left_x + right_x) / 2,
    (y1 + y2) / 2 + offset,
    labels = count,
    cex = 0.68,
    col = color
  )
}

draw_node <- function(x, y, label, fill) {
  rect(x - 0.14, y - 0.035, x + 0.14, y + 0.035, col = fill, border = "#333333", lwd = 0.8)
  text(x, y, labels = format_label(label), cex = 0.75)
}

for (module in modules) {
  draw_node(left_x, module_y[[module]], module, "#f5efe6")
}
for (effector in effectors) {
  draw_node(right_x, effector_y[[effector]], effector, "#e6eef5")
}

text(left_x, 0.94, "Host defense modules", font = 2, cex = 0.9)
text(right_x, 0.94, "Nematode effector classes", font = 2, cex = 0.9)
legend(
  "bottomleft",
  legend = c("susceptible enriched host orthogroups", "tolerant enriched host orthogroups"),
  col = c(bias_colors[["susceptible_enriched"]], bias_colors[["tolerant_enriched"]]),
  lwd = 4,
  bty = "n",
  cex = 0.75
)
text(
  0.5,
  0.035,
  "Edge labels are candidate host orthogroup counts. Links are functional hypotheses, not direct interaction or coevolution evidence.",
  cex = 0.72,
  col = "#444444"
)

par(mar = c(1, 1, 3, 1), family = "sans")
plot.new()
title("Top candidate linked orthogroups for manual review", cex.main = 1.05)
top_display <- head(top_df, 24)
if (nrow(top_display) == 0) {
  text(0.5, 0.5, "No candidate rows available", cex = 1)
} else {
  header <- sprintf("%-4s %-22s %-18s %-20s %-8s %-12s %-18s",
                    "Rank", "Module", "Effector", "Orthogroup", "Score", "Bias", "Keywords")
  lines <- c(header, paste(rep("-", 112), collapse = ""))
  for (i in seq_len(nrow(top_display))) {
    row <- top_display[i, ]
    keyword <- substr(row$matched_keywords, 1, 18)
    lines <- c(
      lines,
      sprintf(
        "%-4s %-22s %-18s %-20s %-8s %-12s %-18s",
        row$top_rank,
        substr(row$module_id, 1, 22),
        substr(row$effector_class, 1, 18),
        substr(row$orthogroup_id, 1, 20),
        row$host_priority_score,
        substr(row$host_bias_direction, 1, 12),
        keyword
      )
    )
  }
  y <- seq(0.92, 0.08, length.out = length(lines))
  text(0.03, y, labels = lines, adj = 0, family = "mono", cex = 0.62)
}

dev.off()
message("Fig5 candidate coevolution PDF written: ", output_path)
