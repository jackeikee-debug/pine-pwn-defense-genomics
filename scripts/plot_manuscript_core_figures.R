#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 17) {
  stop(paste(
    "Usage: plot_manuscript_core_figures.R <species.tsv> <defense_species.tsv>",
    "<orthology_summary.tsv> <candidate_summary.tsv> <effectors.tsv> <links.tsv>",
    "<og_members.tsv> <og_expression.tsv> <og_tree.nwk>",
    "<fig1.pdf> <fig1.png> <fig2.pdf> <fig2.png>",
    "<fig3.pdf> <fig3.png> <fig5.pdf> <fig5.png>"
  ))
}

species <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
defense <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
orthology <- read.delim(args[[3]], stringsAsFactors = FALSE, check.names = FALSE)
candidates <- read.delim(args[[4]], stringsAsFactors = FALSE, check.names = FALSE)
effectors <- read.delim(args[[5]], stringsAsFactors = FALSE, check.names = FALSE)
links <- read.delim(args[[6]], stringsAsFactors = FALSE, check.names = FALSE)
members <- read.delim(args[[7]], stringsAsFactors = FALSE, check.names = FALSE)
expression <- read.delim(args[[8]], stringsAsFactors = FALSE, check.names = FALSE)
tree_path <- args[[9]]
outputs <- args[10:17]

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The R package 'ape' is required to render the OG0005853 gene tree")
}
if (nrow(species) != 7 || nrow(candidates) != 12 || nrow(members) != 15) {
  stop("Unexpected final-table dimensions for species, candidates, or OG0005853 members")
}
for (path in outputs) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

palette <- c(
  east = "#C75B4B",
  north = "#287D8E",
  outgroup = "#777777",
  gold = "#D9A63B",
  pale = "#E9E9E9",
  ink = "#272727",
  green = "#477A5B"
)

pretty_species <- c(
  pden = "P. densiflora", pmas = "P. massoniana", ptab = "P. tabuliformis",
  ptae = "P. taeda", plam = "P. lambertiana", plon = "P. longaeva",
  pabi = "P. abies"
)
species_order <- c("pden", "pmas", "ptab", "ptae", "plam", "plon", "pabi")
species_colors <- c(
  pden = palette[["east"]], pmas = palette[["east"]], ptab = palette[["east"]],
  ptae = palette[["north"]], plam = palette[["north"]], plon = palette[["north"]],
  pabi = palette[["outgroup"]]
)

module_labels <- c(
  hormone_signaling = "Hormone signaling",
  hydraulic_xylem = "Hydraulic/xylem",
  immune_signaling = "Immune signaling",
  phenylpropanoid_lignin = "Phenylpropanoid/lignin",
  resin_terpenoid = "Resin/terpenoid",
  ros_detoxification = "ROS detoxification",
  wound_periderm = "Wound/periderm"
)

render_pair <- function(pdf_path, png_path, draw_function, width, height) {
  pdf(pdf_path, width = width, height = height, family = "sans", useDingbats = FALSE)
  draw_function()
  dev.off()
  png(png_path, width = width, height = height, units = "in", res = 300, type = "cairo")
  draw_function()
  dev.off()
}

draw_box <- function(x, y, w, h, title, detail, fill, border = "white") {
  rect(x, y, x + w, y + h, col = fill, border = border, lwd = 1.2)
  wrapped_title <- paste(strwrap(title, width = 24), collapse = "\n")
  wrapped_detail <- paste(strwrap(detail, width = 34), collapse = "\n")
  text(x + w / 2, y + h * 0.68, wrapped_title, font = 2, cex = 0.72, col = "white")
  text(x + w / 2, y + h * 0.25, wrapped_detail, cex = 0.56, col = "white")
}

draw_fig1 <- function() {
  par(mar = c(0.6, 0.7, 2.8, 0.7), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title("Study design, quality gates, and candidate attrition", cex.main = 1.25, font.main = 2)

  step_box <- function(x, y, w, h, count, label, fill) {
    rect(x, y, x + w, y + h, col = fill, border = NA)
    text(x + w / 2, y + h * 0.63, count, font = 2, cex = 0.90, col = "white")
    text(x + w / 2, y + h * 0.27, paste(strwrap(label, 27), collapse = "\n"), cex = 0.55, col = "white")
  }

  text(0.16, 0.91, "Host proteomes", font = 2, cex = 0.90, col = palette[["green"]])
  text(0.84, 0.91, "PWN proteome", font = 2, cex = 0.90, col = palette[["north"]])
  step_box(0.04, 0.72, 0.24, 0.13, "29,637", "stable host orthogroups", palette[["green"]])
  step_box(0.04, 0.53, 0.24, 0.13, "1,508", "defense-annotated orthogroups", "#5B8368")
  arrows(0.16, 0.715, 0.16, 0.67, length = 0.055, lwd = 1.3, col = "#5E666A")

  step_box(0.72, 0.72, 0.24, 0.13, "15,860", "BXYJv5 proteins", "#273238")
  step_box(0.72, 0.53, 0.24, 0.13, "2,741", "DeepSig signal-peptide support", palette[["north"]])
  step_box(0.72, 0.34, 0.24, 0.13, "2,164", "secreted-soluble candidates", palette[["gold"]])
  arrows(0.84, 0.715, 0.84, 0.67, length = 0.055, lwd = 1.3, col = "#5E666A")
  arrows(0.84, 0.525, 0.84, 0.48, length = 0.055, lwd = 1.3, col = "#5E666A")
  text(0.84, 0.29, "InterPro / MEROPS / dbCAN /\npublished-secretome support", cex = 0.57, col = "#4E565A")

  arrows(0.285, 0.595, 0.395, 0.595, length = 0.06, lwd = 1.4, col = "#5E666A")
  arrows(0.715, 0.405, 0.605, 0.515, length = 0.06, lwd = 1.4, col = "#5E666A")
  step_box(0.395, 0.48, 0.21, 0.14, "127", "prior-linked host orthogroups", "#7E709A")
  step_box(0.395, 0.29, 0.21, 0.14, "37", "explicitly reviewed candidates", "#6A7B84")
  step_box(0.395, 0.10, 0.21, 0.14, "12", "core host candidates", palette[["east"]])
  arrows(0.50, 0.475, 0.50, 0.44, length = 0.055, lwd = 1.3, col = "#5E666A")
  arrows(0.50, 0.285, 0.50, 0.25, length = 0.055, lwd = 1.3, col = "#5E666A")

  gate_y <- 0.015
  gate_x <- c(0.02, 0.66, 0.83)
  gate_w <- c(0.31, 0.16, 0.15)
  gate_titles <- c("Comparability gate", "Evidence convergence", "Claim boundary")
  gate_details <- c(
    "Stable membership; proteome-quality sensitivity; transcript evidence excluded from copy counts",
    "Expression, network, structure, and sequence kept as separate channels",
    "Candidate mechanisms; no direct targeting or causal resistance claim"
  )
  for (i in seq_along(gate_x)) {
    rect(gate_x[i], gate_y, gate_x[i] + gate_w[i], gate_y + 0.075, col = "#F2F4F4", border = "#8E9598", lwd = 0.8)
    text(gate_x[i] + gate_w[i] / 2, gate_y + 0.052, gate_titles[i], font = 2, cex = 0.54)
    text(gate_x[i] + gate_w[i] / 2, gate_y + 0.022, paste(strwrap(gate_details[i], width = ifelse(i == 1, 50, 31)), collapse = "\n"), cex = 0.40, col = "#4D5559")
  }
}

draw_fig2 <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.0, 2.0))
  ordered_species <- species[match(species_order, species$species_id), ]

  par(mar = c(7, 4.5, 4, 1), family = "sans")
  proteome_max <- max(ordered_species$protein_count / 1000)
  mids <- barplot(
    ordered_species$protein_count / 1000,
    names.arg = pretty_species[species_order],
    col = species_colors[species_order], border = NA, las = 2,
    ylab = "Proteins retained (thousands)", main = "A  Host proteome panel", cex.names = 0.75,
    ylim = c(0, proteome_max * 1.16)
  )
  text(mids, ordered_species$protein_count / 1000 + proteome_max * 0.03, labels = format(ordered_species$protein_count, big.mark = ","), cex = 0.65, srt = 90)
  legend("topright", legend = c("East Asian pine", "North American pine", "Outgroup"), fill = c(palette[["east"]], palette[["north"]], palette[["outgroup"]]), border = NA, bty = "n", cex = 0.72)

  module_order <- names(module_labels)
  matrix_values <- matrix(NA_real_, nrow = length(module_order), ncol = length(species_order), dimnames = list(module_order, species_order))
  for (i in seq_len(nrow(defense))) {
    module <- defense$module_id[[i]]
    sid <- defense$species_id[[i]]
    if (module %in% module_order && sid %in% species_order) {
      protein_count <- species$protein_count[match(sid, species$species_id)]
      matrix_values[module, sid] <- defense$gene_count[[i]] / protein_count * 10000
    }
  }
  breaks <- seq(min(matrix_values, na.rm = TRUE), max(matrix_values, na.rm = TRUE), length.out = 101)
  heat_colors <- colorRampPalette(c("#F4F4F4", "#D9A63B", "#C75B4B"))(100)
  color_index <- matrix(
    pmax(1, pmin(100, findInterval(matrix_values, breaks, all.inside = TRUE))),
    nrow = nrow(matrix_values),
    ncol = ncol(matrix_values),
    dimnames = dimnames(matrix_values)
  )

  par(mar = c(7, 11, 4, 4.5), family = "sans")
  image(seq_along(species_order), seq_along(module_order), t(color_index[nrow(color_index):1, ]), col = heat_colors, axes = FALSE, xlab = "", ylab = "", main = "B  Defense-module density")
  axis(1, at = seq_along(species_order), labels = pretty_species[species_order], las = 2, cex.axis = 0.75)
  axis(2, at = seq_along(module_order), labels = rev(module_labels[module_order]), las = 2, cex.axis = 0.74)
  for (r in seq_along(module_order)) for (c in seq_along(species_order)) {
    value <- matrix_values[module_order[[r]], species_order[[c]]]
    text(c, length(module_order) - r + 1, sprintf("%.0f", value), cex = 0.62, col = ifelse(value > median(matrix_values, na.rm = TRUE), "white", palette[["ink"]]))
  }
  legend("right", inset = c(-0.16, 0), legend = c("Lower", "Higher"), fill = c(heat_colors[[15]], heat_colors[[85]]), border = NA, bty = "n", xpd = NA, cex = 0.7, title = "Genes per\n10,000 proteins")
  layout(1)
}

format_class <- function(x) gsub("_", " ", x, fixed = TRUE)

draw_fig3 <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.0, 1.7))
  class_counts <- sort(table(effectors$effector_class), decreasing = FALSE)
  class_max <- max(class_counts)
  class_colors <- setNames(colorRampPalette(c(palette[["gold"]], palette[["north"]]))(length(class_counts)), names(class_counts))

  par(mar = c(4.5, 10, 4, 1), family = "sans")
  mids <- barplot(class_counts, horiz = TRUE, names.arg = format_class(names(class_counts)), las = 1, col = class_colors, border = NA, xlab = "Predicted candidate proteins", main = "A  Nematode candidate classes", cex.names = 0.75, xlim = c(0, class_max * 1.16))
  text(as.numeric(class_counts) + max(class_counts) * 0.02, mids, labels = as.numeric(class_counts), adj = 0, cex = 0.68)

  effector_levels <- rev(sort(unique(links$effector_class)))
  module_levels <- names(module_labels)
  par(mar = c(8, 10, 4, 2), family = "sans")
  plot(NA, xlim = c(0.5, length(module_levels) + 0.5), ylim = c(0.5, length(effector_levels) + 0.5), axes = FALSE, xlab = "", ylab = "", main = "B  Effector-class to host-module hypotheses")
  axis(1, at = seq_along(module_levels), labels = module_labels[module_levels], las = 2, cex.axis = 0.72)
  axis(2, at = seq_along(effector_levels), labels = format_class(effector_levels), las = 2, cex.axis = 0.72)
  abline(v = seq_along(module_levels), h = seq_along(effector_levels), col = "#EFEFEF", lwd = 0.8)
  for (i in seq_len(nrow(links))) {
    x <- match(links$host_module[[i]], module_levels)
    y <- match(links$effector_class[[i]], effector_levels)
    if (!is.na(x) && !is.na(y)) {
      size <- 0.7 + sqrt(as.numeric(links$host_candidate_count[[i]])) / 3.2
      points(x, y, pch = 21, bg = class_colors[[links$effector_class[[i]]]], col = "white", cex = size)
      text(x, y, links$host_candidate_count[[i]], cex = 0.54, font = 2, col = "white")
    }
  }
  mtext("Circle labels are linked host orthogroup counts; connections are functional hypotheses.", side = 1, line = 6.7, cex = 0.68, col = "#444444")
  layout(1)
}

short_gene_label <- function(gene_id) {
  parts <- strsplit(gene_id, "|", fixed = TRUE)[[1]]
  if (length(parts) > 1) paste0(parts[[1]], " | ", parts[[2]]) else gene_id
}

draw_fig5 <- function() {
  tree <- ape::read.tree(tree_path)
  original_tips <- tree$tip.label
  stable_ids <- sub("^[^_]+_", "", original_tips)
  member_index <- match(stable_ids, members$gene_id)
  if (any(is.na(member_index))) stop("Could not map all OG0005853 tree tips to member annotations")
  tip_species <- members$species_id[member_index]
  tip_classes <- members$annotation_class[member_index]
  tip_numbers <- ave(seq_along(member_index), interaction(tip_species, tip_classes), FUN = seq_along)
  tree$tip.label <- sprintf("%s | %s-%d", tip_species, tip_classes, tip_numbers)
  tip_regions <- members$region[member_index]
  tip_colors <- c(East_Asia = palette[["east"]], North_America = palette[["north"]], Outgroup = palette[["outgroup"]])[tip_regions]
  class_colors <- c(PR4 = palette[["green"]], WIN1 = palette[["gold"]], WHEATWIN = "#8C6BB1")

  layout(matrix(c(1, 2), nrow = 1), widths = c(1.25, 1.0))
  par(mar = c(3, 1, 4, 1), family = "sans")
  ape::plot.phylo(tree, type = "phylogram", show.tip.label = TRUE, tip.color = tip_colors, cex = 0.62, label.offset = 0.006, edge.color = "#555555", main = "A  OG0005853 PR-4/WIN1 family")
  ape::tiplabels(pch = 21, bg = class_colors[members$annotation_class[member_index]], col = "white", cex = 0.85)
  legend("bottomleft", legend = c("East Asia", "North America", "Outgroup"), text.col = c(palette[["east"]], palette[["north"]], palette[["outgroup"]]), bty = "n", cex = 0.68)
  legend("topleft", legend = names(class_colors), pt.bg = class_colors, pch = 21, pt.cex = 1.0, bty = "n", cex = 0.68)
  mtext("All 15 members have the recorded N-terminal hydrophobic flag; full protein IDs are tabulated separately.", side = 1, line = 1.1, cex = 0.58, col = "#4D5559")

  finite_expression <- expression[is.finite(expression$log2FoldChange), ]
  contrast_labels <- c(
    pathogen_associated = "P. densiflora:\npathogenic vs non-pathogenic",
    bxyl_vs_water = "P. densiflora:\nPWN vs water",
    bthai_vs_water = "P. densiflora:\nnon-pathogenic vs water",
    two_week_vs_zero = "P. strobus:\n2 weeks vs 0",
    four_week_vs_zero = "P. strobus:\n4 weeks vs 0"
  )
  contrast_order <- names(contrast_labels)
  x_base <- match(finite_expression$contrast, contrast_order)
  offsets <- ave(seq_along(x_base), x_base, FUN = function(x) seq(-0.18, 0.18, length.out = length(x)))
  x <- x_base + offsets
  significant <- finite_expression$de_interpretation == "significant_up"
  point_fill <- ifelse(finite_expression$species_id == "pden", palette[["east"]], palette[["north"]])

  par(mar = c(9, 4.5, 4, 1), family = "sans")
  plot(x, finite_expression$log2FoldChange, type = "n", xaxt = "n", xlab = "", ylab = expression(log[2]~fold~change), xlim = c(0.55, 5.45), ylim = range(c(-3.5, finite_expression$log2FoldChange + 1.2)), main = "B  Infection-expression evidence")
  abline(h = 0, lty = 2, col = "#777777")
  axis(1, at = seq_along(contrast_order), labels = contrast_labels[contrast_order], las = 2, cex.axis = 0.66)
  abline(v = 3.5, col = "#B9BEC1", lty = 3)
  text(2.0, par("usr")[4] * 0.98, "P. densiflora controlled contrasts", cex = 0.60, font = 2)
  text(4.5, par("usr")[4] * 0.98, "P. strobus temporal contrasts", cex = 0.60, font = 2)
  points(x, finite_expression$log2FoldChange, pch = ifelse(significant, 21, 1), bg = point_fill, col = point_fill, cex = ifelse(significant, 1.25, 0.9), lwd = 1.3)
  text(x[significant], finite_expression$log2FoldChange[significant] + 0.45, labels = sub("TRINITY_", "", finite_expression$expression_gene_id[significant]), cex = 0.52, srt = 18)
  legend("bottomright", legend = c("P. densiflora", "P. strobus", "FDR < 0.05"), pch = c(16, 16, 21), col = c(palette[["east"]], palette[["north"]], palette[["ink"]]), pt.bg = c(palette[["east"]], palette[["north"]], "white"), bty = "n", cex = 0.7)
  mtext("Fold-change magnitudes are shown within each experiment and are not compared statistically across species.", side = 1, line = 7.8, cex = 0.66, col = "#444444")
  layout(1)
}

render_pair(outputs[[1]], outputs[[2]], draw_fig1, 13.0, 7.2)
render_pair(outputs[[3]], outputs[[4]], draw_fig2, 13.0, 7.5)
render_pair(outputs[[5]], outputs[[6]], draw_fig3, 13.5, 7.5)
render_pair(outputs[[7]], outputs[[8]], draw_fig5, 14.0, 8.0)

message("Wrote core manuscript figures: ", paste(outputs, collapse = ", "))
