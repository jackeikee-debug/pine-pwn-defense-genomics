#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste(
    "Usage: plot_effector_module_candidate_network.R <candidates.tsv>",
    "<effector_context.tsv> <edges.tsv> <sectors.tsv> <output.pdf> <output.png>"
  ))
}

candidates <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
effector_context <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
edges <- read.delim(args[[3]], stringsAsFactors = FALSE, check.names = FALSE)
sectors <- read.delim(args[[4]], stringsAsFactors = FALSE, check.names = FALSE)
output_pdf <- args[[5]]
output_png <- args[[6]]

required_candidates <- c(
  "orthogroup_id", "priority_class", "mapped_symbols",
  "network_support", "structure_support", "sequence_support", "cross_species_state"
)
required_context <- c(
  "effector_class", "functional_evaluated_count", "functional_supported_count",
  "functional_annotation_state"
)
required_edges <- c("source_id", "target_id", "edge_type", "weight", "evidence_boundary")
required_sectors <- c("sector_id", "sector_type", "label", "display_order")
if (!all(required_candidates %in% names(candidates))) stop("Candidate table is missing required columns")
if (!all(required_context %in% names(effector_context))) stop("Effector-context table is missing required columns")
if (!all(required_edges %in% names(edges))) stop("Edge table is missing required columns")
if (!all(required_sectors %in% names(sectors))) stop("Sector table is missing required columns")
if (any(edges$evidence_boundary %in% c("direct_interaction", "validated_binding"))) {
  stop("Direct interaction links are outside scope")
}

functional_edges <- edges[edges$edge_type == "functional_prior", , drop = FALSE]
membership_edges <- edges[edges$edge_type == "candidate_membership", , drop = FALSE]
candidate_ids <- candidates$orthogroup_id
effector_ids <- unique(functional_edges$source_id)
module_ids <- unique(c(functional_edges$target_id, membership_edges$source_id))
if (!setequal(effector_ids, effector_context$effector_class)) {
  stop("Every effector class requires one context row")
}
if (!all(membership_edges$target_id %in% candidate_ids)) {
  stop("Every module-membership edge must point to a plotted candidate")
}

sector_order <- function(ids) {
  matched <- sectors[match(ids, sectors$sector_id), , drop = FALSE]
  if (any(is.na(matched$sector_id))) stop("Every plotted node requires a sector row")
  ids[order(as.numeric(matched$display_order), ids)]
}
effector_ids <- sector_order(effector_ids)
module_ids <- sector_order(module_ids)
candidate_ids <- candidates$orthogroup_id
candidates <- candidates[match(candidate_ids, candidates$orthogroup_id), , drop = FALSE]
effector_context <- effector_context[match(effector_ids, effector_context$effector_class), , drop = FALSE]

short_labels <- c(
  cell_wall_modifying = "Cell-wall\nmodifying",
  detoxification = "Detoxification",
  mimicry = "Mimicry",
  protease = "Protease",
  stress_adaptation = "Stress\nadaptation",
  hydraulic_xylem = "Hydraulic /\nxylem",
  immune_signaling = "Immune\nregulation",
  phenylpropanoid_lignin = "Phenylpropanoid /\nlignin",
  ros_detoxification = "ROS\ncontrol",
  wound_periderm = "Wound /\nperiderm"
)
candidate_numbers <- setNames(sprintf("%02d", seq_along(candidate_ids)), candidate_ids)
candidate_labels <- paste0(candidate_numbers, "\n", candidate_ids)
symbols <- candidates$mapped_symbols
candidate_labels[nzchar(symbols)] <- paste0(
  candidate_labels[nzchar(symbols)], "\n", sub(";.*", "", symbols[nzchar(symbols)])
)

palette <- list(
  effector = c("#B68A31", "#867A3D", "#537D55", "#2F7B68", "#267487"),
  module = c("#426B91", "#27818B", "#39806A", "#72754E", "#8D6448"),
  leading = "#BB4E43", supported = "#D78554",
  ink = "#222222", pale = "#F6F6F6", muted = "#666666"
)
effector_colors <- setNames(rep(palette$effector, length.out = length(effector_ids)), effector_ids)
module_colors <- setNames(rep(palette$module, length.out = length(module_ids)), module_ids)
candidate_colors <- setNames(
  ifelse(candidates$priority_class == "mechanism_leading", palette$leading, palette$supported),
  candidate_ids
)
node_colors <- c(effector_colors, module_colors, candidate_colors)

support_count <- rowSums(candidates[c("network_support", "structure_support", "sequence_support")] == "yes")
candidate_cex <- setNames(1.05 + support_count * 0.13, candidate_ids)
effector_cex <- setNames(rep(1.15, length(effector_ids)), effector_ids)
module_cex <- setNames(rep(1.20, length(module_ids)), module_ids)
node_cex <- c(effector_cex, module_cex, candidate_cex)

node_positions <- rbind(
  data.frame(id = effector_ids, type = "effector", x = 0.14, y = seq(0.84, 0.30, length.out = length(effector_ids))),
  data.frame(id = module_ids, type = "module", x = 0.50, y = seq(0.84, 0.30, length.out = length(module_ids))),
  data.frame(id = candidate_ids, type = "candidate", x = 0.86, y = seq(0.94, 0.08, length.out = length(candidate_ids)))
)
rownames(node_positions) <- node_positions$id

draw_curve <- function(x1, y1, x2, y2, color, lwd) {
  xspline(
    c(x1, (x1 + x2) / 2, x2),
    c(y1, y1, y2),
    shape = c(0, 1, 0), open = TRUE, col = color, lwd = lwd
  )
}

draw_node <- function(id) {
  pos <- node_positions[id, ]
  points(pos$x, pos$y, pch = 21, bg = node_colors[[id]], col = "#333333", lwd = 0.7, cex = node_cex[[id]])
}

draw_figure <- function() {
  par(mar = c(1.2, 1.2, 2.6, 1.0), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = NA)
  title("Predicted effector-prior defense-module network", adj = 0, cex.main = 1.05, font.main = 2)

  text(c(0.14, 0.50, 0.86), 0.985, c("PWN effector classes", "Host defense modules", "Candidate orthogroups"),
       font = 2, cex = 0.72)
  abline(v = c(0.32, 0.68), col = "#E6E6E6", lwd = 0.8)

  for (i in seq_len(nrow(functional_edges))) {
    edge <- functional_edges[i, ]
    from <- node_positions[edge$source_id, ]
    to <- node_positions[edge$target_id, ]
    lwd <- 1.2 + 1.8 * as.numeric(edge$weight) / max(as.numeric(functional_edges$weight))
    draw_curve(from$x, from$y, to$x, to$y, grDevices::adjustcolor(module_colors[[edge$target_id]], 0.35), lwd)
  }
  for (i in seq_len(nrow(membership_edges))) {
    edge <- membership_edges[i, ]
    from <- node_positions[edge$source_id, ]
    to <- node_positions[edge$target_id, ]
    draw_curve(from$x, from$y, to$x, to$y, grDevices::adjustcolor(module_colors[[edge$source_id]], 0.55), 1.0)
  }

  invisible(lapply(node_positions$id, draw_node))

  effector_labels <- short_labels[effector_ids]
  module_labels <- short_labels[module_ids]
  text(node_positions[effector_ids, "x"] - 0.025, node_positions[effector_ids, "y"], effector_labels,
       adj = 1, cex = 0.60)
  text(node_positions[module_ids, "x"], node_positions[module_ids, "y"] - 0.035, module_labels,
       adj = c(0.5, 1), cex = 0.55)
  text(node_positions[candidate_ids, "x"] + 0.025, node_positions[candidate_ids, "y"], candidate_labels,
       adj = 0, cex = 0.50)

  legend_x <- 0.02
  legend_y <- 0.07
  points(legend_x, legend_y + 0.09, pch = 21, bg = palette$leading, col = "#333333", cex = 0.95)
  text(legend_x + 0.025, legend_y + 0.09, "leading candidate", adj = 0, cex = 0.55)
  points(legend_x, legend_y + 0.055, pch = 21, bg = palette$supported, col = "#333333", cex = 0.95)
  text(legend_x + 0.025, legend_y + 0.055, "supported candidate", adj = 0, cex = 0.55)
  segments(legend_x - 0.005, legend_y + 0.02, legend_x + 0.04, legend_y + 0.02,
           col = grDevices::adjustcolor("#557A75", 0.55), lwd = 2.6)
  text(legend_x + 0.055, legend_y + 0.02, "functional prior / module annotation; not direct interaction",
       adj = 0, cex = 0.55, col = palette$muted)
  text(0.02, 0.005, "Edges encode class-level hypotheses and module membership only; they do not establish direct effector targeting.",
       adj = 0, cex = 0.55, col = palette$muted)
}

for (path in c(output_pdf, output_png)) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
pdf(output_pdf, width = 10, height = 7.4, pointsize = 16, family = "sans", useDingbats = FALSE)
draw_figure()
dev.off()
png(output_png, width = 10, height = 7.4, units = "in", res = 600, pointsize = 16, type = "cairo")
draw_figure()
dev.off()
message("Wrote effector-module-candidate network: ", output_pdf, ", ", output_png)
