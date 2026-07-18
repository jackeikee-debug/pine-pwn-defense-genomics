#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste(
    "Usage: plot_integrated_candidate_evidence.R",
    "<bubble.tsv> <member_audit.tsv> <candidate_summary.tsv> <output.pdf> <output.png>"
  ))
}

bubble <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
audit <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
candidates <- read.delim(args[[3]], stringsAsFactors = FALSE, check.names = FALSE)
pdf_path <- args[[4]]
png_path <- args[[5]]

if (nrow(bubble) != 36 || length(unique(bubble$orthogroup_id)) != 12 || nrow(candidates) != 12) {
  stop("Integrated figure requires 36 P. massoniana contrast rows and 12 candidates")
}

module_order <- c(
  immune_signaling = 1,
  ros_detoxification = 2,
  phenylpropanoid_lignin = 3,
  hydraulic_xylem = 4,
  wound_periderm = 5
)
module_labels <- c(
  immune_signaling = "Immune signaling",
  ros_detoxification = "ROS detoxification",
  phenylpropanoid_lignin = "Phenylpropanoid/lignin",
  hydraulic_xylem = "Hydraulic/xylem",
  wound_periderm = "Wound/periderm"
)
module_colors <- c(
  immune_signaling = "#4E79A7",
  ros_detoxification = "#C7665B",
  phenylpropanoid_lignin = "#6F9E55",
  hydraulic_xylem = "#3A7D78",
  wound_periderm = "#7E709A"
)

primary_module <- function(value) {
  modules <- strsplit(value, ";", fixed = TRUE)[[1]]
  modules <- modules[modules %in% names(module_order)]
  if (length(modules) == 0) return("hydraulic_xylem")
  modules[which.min(module_order[modules])]
}

candidates$primary_module <- vapply(candidates$module_ids, primary_module, character(1))
candidates$priority_label <- ifelse(candidates$overall_priority_class == "mechanism_leading", "High", "Intermediate")
candidates$mapped_label <- vapply(seq_len(nrow(candidates)), function(i) {
  symbol <- strsplit(candidates$mapped_symbols[i], ";", fixed = TRUE)[[1]][1]
  if (is.na(symbol) || !nzchar(symbol)) candidates$orthogroup_id[i] else paste(candidates$orthogroup_id[i], symbol, sep = " | ")
}, character(1))
candidates <- candidates[order(module_order[candidates$primary_module], candidates$priority_label, candidates$orthogroup_id), ]

finite_number <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  value[!is.finite(value)] <- NA_real_
  value
}

bubble$log2_fold_change <- finite_number(bubble$log2_fold_change)
bubble$whole_transcriptome_padj <- finite_number(bubble$whole_transcriptome_padj)
bubble$candidate_set_simes_padj <- finite_number(bubble$candidate_set_simes_padj)

audit_subset <- audit[
  audit$orthogroup_id %in% candidates$orthogroup_id &
    audit$evidence_unit == "mapped_gene_set_summary" &
    audit$contrast %in% c("pathogen_associated", "bxyl_vs_water", "two_week_vs_zero", "four_week_vs_zero"),
]
audit_subset$log2_fold_change <- finite_number(audit_subset$log2_fold_change)
audit_subset$whole_transcriptome_padj <- finite_number(audit_subset$whole_transcriptome_padj)

lfc_palette <- colorRampPalette(c("#376AA6", "#F7F7F5", "#C54E52"))(201)
lfc_color <- function(value, limit = 4) {
  if (!is.finite(value)) return("#FFFFFF")
  clipped <- max(-limit, min(limit, value))
  lfc_palette[round((clipped + limit) / (2 * limit) * 200) + 1]
}

draw_symbol <- function(x, y, supported, type) {
  if (!supported) {
    points(x, y, pch = 1, col = "#B5BABD", cex = 0.80, lwd = 0.8)
    return()
  }
  pch <- c(network = 21, structure = 24, sequence = 23)[[type]]
  points(x, y, pch = pch, bg = "#2E5968", col = "#2E5968", cex = 0.88)
}

draw_figure <- function() {
  par(mar = c(5.0, 1.0, 3.7, 1.0), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 14.8), ylim = c(0, 15.3), xaxs = "i", yaxs = "i")
  text(0.05, 15.05, "Integrated expression and independent evidence for 12 host orthogroups", adj = c(0, 1), font = 2, cex = 1.05)
  text(0.05, 14.55, "P. massoniana is a balanced 2 x 2 experiment; independent pine contrasts are displayed separately and are not pooled.", adj = c(0, 1), cex = 0.63, col = "#50585C")

  label_x <- 0.25
  module_x <- 3.72
  pmas_x <- c(R = 5.15, S = 5.90, GxI = 6.65)
  independent_x <- c(Pden_BxBt = 8.00, Pden_BxW = 8.75, Pstro_2w = 9.50, Pstro_4w = 10.25)
  priority_x <- 11.35
  support_x <- c(network = 12.35, structure = 13.05, sequence = 13.75)
  row_y <- setNames(seq(13.0, 2.25, length.out = 12), candidates$orthogroup_id)

  text(label_x, 13.85, "Orthogroup | mapped homolog", adj = 0, font = 2, cex = 0.62)
  text(module_x, 13.85, "Module", font = 2, cex = 0.62)
  text(mean(pmas_x), 14.10, "P. massoniana", font = 2, cex = 0.68)
  text(pmas_x, 13.68, names(pmas_x), font = 2, cex = 0.63)
  text(mean(independent_x), 14.10, "Independent expression contrasts", font = 2, cex = 0.68)
  text(independent_x, 13.68, c("P. den.\nBx/Bt", "P. den.\nBx/water", "P. str.\n2w/0w", "P. str.\n4w/0w"), cex = 0.53)
  text(priority_x, 13.85, "Priority", font = 2, cex = 0.60)
  text(support_x, 13.85, c("Network", "Structure", "Sequence"), font = 2, cex = 0.53)

  previous_module <- ""
  for (i in seq_len(nrow(candidates))) {
    candidate <- candidates[i, ]
    y <- row_y[[candidate$orthogroup_id]]
    if (candidate$primary_module != previous_module) {
      rect(0.05, y + 0.39, 14.25, y + 0.57, col = module_colors[[candidate$primary_module]], border = NA)
      previous_module <- candidate$primary_module
    }
    segments(0.05, y - 0.38, 14.25, y - 0.38, col = "#E4E7E8", lwd = 0.7)
    text(label_x, y, candidate$mapped_label, adj = 0, cex = 0.57)
    text(module_x, y, module_labels[[candidate$primary_module]], cex = 0.49, col = module_colors[[candidate$primary_module]])

    rows <- bubble[bubble$orthogroup_id == candidate$orthogroup_id, ]
    for (contrast in names(pmas_x)) {
      row <- rows[rows$contrast_label == contrast, ]
      padj <- row$whole_transcriptome_padj
      if (length(padj) == 0 || !is.finite(padj)) {
        rect(pmas_x[[contrast]] - 0.16, y - 0.16, pmas_x[[contrast]] + 0.16, y + 0.16,
             col = lfc_color(row$log2_fold_change), border = "#9AA0A3", lwd = 0.8)
        text(pmas_x[[contrast]], y, "x", cex = 0.58, font = 2, col = "#62696D")
      } else {
        score <- min(-log10(max(padj, 1e-8)), 6)
        symbols(pmas_x[[contrast]], y, circles = sqrt(max(score, 0.22)), inches = 0.105,
                add = TRUE, bg = lfc_color(row$log2_fold_change),
                fg = ifelse(padj < 0.05, "#171A1C", "#A5AAAD"),
                lwd = ifelse(padj < 0.05, 1.45, 0.75))
      }
      if (length(row$candidate_set_simes_padj) && is.finite(row$candidate_set_simes_padj) && row$candidate_set_simes_padj < 0.05) {
        points(pmas_x[[contrast]] + 0.22, y + 0.22, pch = 18, col = "#7A5A00", cex = 0.50)
      }
    }

    audit_map <- c(Pden_BxBt = "pathogen_associated", Pden_BxW = "bxyl_vs_water", Pstro_2w = "two_week_vs_zero", Pstro_4w = "four_week_vs_zero")
    for (key in names(independent_x)) {
      row <- audit_subset[audit_subset$orthogroup_id == candidate$orthogroup_id & audit_subset$contrast == audit_map[[key]], ]
      if (nrow(row) == 0 || !is.finite(row$log2_fold_change[1])) {
        rect(independent_x[[key]] - 0.16, y - 0.16, independent_x[[key]] + 0.16, y + 0.16, col = "white", border = "#B5BABD")
        text(independent_x[[key]], y, "x", cex = 0.56, col = "#737A7E")
      } else {
        significant <- is.finite(row$whole_transcriptome_padj[1]) && row$whole_transcriptome_padj[1] < 0.05
        rect(independent_x[[key]] - 0.18, y - 0.18, independent_x[[key]] + 0.18, y + 0.18,
             col = lfc_color(row$log2_fold_change[1]), border = ifelse(significant, "#171A1C", "#A5AAAD"),
             lwd = ifelse(significant, 1.35, 0.70))
      }
    }

    rect(priority_x - 0.34, y - 0.20, priority_x + 0.34, y + 0.20,
         col = ifelse(candidate$priority_label == "High", "#DCE9E4", "#E9E6F0"), border = NA)
    text(priority_x, y, candidate$priority_label, cex = 0.48)
    draw_symbol(support_x[["network"]], y, candidate$experimental_or_database_supported == "true", "network")
    draw_symbol(support_x[["structure"]], y, as.integer(candidate$structurally_supported_effector_count) > 0, "structure")
    draw_symbol(support_x[["sequence"]], y, candidate$apoplastic_cell_wall_sequence_support_class == "tier1_candidate_substrate", "sequence")
  }

  legend_y <- 0.95
  text(0.20, legend_y, "Expression fill: log2FC", adj = 0, cex = 0.55)
  gradient_x <- seq(1.85, 3.65, length.out = 201)
  for (i in seq_len(200)) rect(gradient_x[i], legend_y - 0.12, gradient_x[i + 1], legend_y + 0.12, col = lfc_palette[i], border = NA)
  text(c(1.85, 2.75, 3.65), legend_y - 0.30, c("-4", "0", "+4"), cex = 0.47)
  symbols(4.35, legend_y, circles = sqrt(3), inches = 0.105, add = TRUE, bg = "white", fg = "#171A1C", lwd = 1.4)
  text(4.60, legend_y, "whole-transcriptome padj < 0.05", adj = 0, cex = 0.51)
  rect(7.15, legend_y - 0.14, 7.43, legend_y + 0.14, col = "white", border = "#9AA0A3")
  text(7.29, legend_y, "x", cex = 0.52)
  text(7.55, legend_y, "padj unavailable", adj = 0, cex = 0.51)
  points(9.30, legend_y, pch = 18, col = "#7A5A00", cex = 0.60)
  text(9.50, legend_y, "candidate-set Simes FDR < 0.05 (descriptive)", adj = 0, cex = 0.51)
  text(0.20, 0.30, "Bubble area is clipped -log10(whole-transcriptome padj); independent contrasts use square cells. Missing adjusted values are never coded as zero.", adj = 0, cex = 0.54, col = "#50585C")
}

dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
pdf(pdf_path, width = 14.0, height = 8.6, family = "sans", useDingbats = FALSE); draw_figure(); dev.off()
png(png_path, width = 14.0, height = 8.6, units = "in", res = 300, type = "cairo"); draw_figure(); dev.off()
