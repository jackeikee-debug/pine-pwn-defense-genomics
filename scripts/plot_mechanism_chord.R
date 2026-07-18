#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: plot_mechanism_chord.R <edges.tsv> <sectors.tsv> <output.pdf> <output.png>")
}

if (!requireNamespace("circlize", quietly = TRUE)) {
  stop("The R package 'circlize' is required")
}

edges <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
sectors <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
output_pdf <- args[[3]]
output_png <- args[[4]]

required_edge_columns <- c("source_id", "target_id", "edge_type", "weight", "evidence_boundary")
required_sector_columns <- c("sector_id", "sector_type", "label", "display_order")
if (!all(required_edge_columns %in% names(edges))) stop("Chord edge table is missing required columns")
if (!all(required_sector_columns %in% names(sectors))) stop("Chord sector table is missing required columns")
if (any(duplicated(edges[c("source_id", "target_id", "edge_type")]))) stop("Chord edges must be unique")
if (any(edges$evidence_boundary %in% c("direct_interaction", "validated_binding"))) {
  stop("Direct-interaction edge types are outside the figure scope")
}

used_ids <- unique(c(edges$source_id, edges$target_id))
sectors <- sectors[sectors$sector_id %in% used_ids, ]
sectors <- sectors[order(as.integer(sectors$display_order)), ]
if (!setequal(used_ids, sectors$sector_id)) stop("Every used edge endpoint must have a sector record")

palette <- list(
  effector = c("#D9A63B", "#B79B46", "#8B9B57", "#67966A", "#4F8A7A", "#337E88", "#287D8E"),
  module = c("#5D77A8", "#4E8AA5", "#3E91A0", "#478F83", "#617F6A", "#7A755D", "#8D684E"),
  leading = "#C75B4B",
  supported = "#D58B62",
  ink = "#2B2B2B"
)

effector_ids <- sectors$sector_id[sectors$sector_type == "effector_class"]
module_ids <- sectors$sector_id[sectors$sector_type == "host_module"]
candidate_ids <- sectors$sector_id[sectors$sector_type == "candidate"]
sector_colors <- setNames(rep("#999999", nrow(sectors)), sectors$sector_id)
sector_colors[effector_ids] <- rep(palette$effector, length.out = length(effector_ids))
sector_colors[module_ids] <- rep(palette$module, length.out = length(module_ids))
candidate_priority <- sectors$priority_class[match(candidate_ids, sectors$sector_id)]
sector_colors[candidate_ids] <- ifelse(candidate_priority == "mechanism_leading", palette$leading, palette$supported)

edge_colors <- ifelse(
  edges$edge_type == "functional_prior",
  grDevices::adjustcolor(sector_colors[edges$source_id], alpha.f = 0.68),
  grDevices::adjustcolor(sector_colors[edges$target_id], alpha.f = 0.58)
)
plot_weight <- ifelse(
  edges$edge_type == "functional_prior",
  1 + sqrt(pmax(1, as.numeric(edges$weight))) / 4,
  1
)
chord_data <- data.frame(from = edges$source_id, to = edges$target_id, value = plot_weight)

candidate_label <- function(row) {
  codes <- character()
  if (row[["expression_support"]] == "cross_species_concordant") codes <- c(codes, "E2")
  if (row[["expression_support"]] == "single_species_concordant") codes <- c(codes, "E1")
  if (row[["network_support"]] == "experimental_or_database_supported") codes <- c(codes, "N")
  if (row[["structure_support"]] == "structurally_supported_anchor") codes <- c(codes, "S")
  if (row[["sequence_support"]] == "tier1_candidate_substrate") codes <- c(codes, "Q")
  suffix <- if (length(codes)) paste0(" [", paste(codes, collapse = ","), "]") else ""
  paste0(row[["sector_id"]], suffix)
}

display_labels <- setNames(sectors$label, sectors$sector_id)
candidate_numbers <- sprintf("%02d", seq_along(candidate_ids))
candidate_keys <- vapply(candidate_ids, function(candidate_id) {
  row <- sectors[sectors$sector_id == candidate_id, , drop = FALSE]
  candidate_label(row[1, ])
}, character(1))
candidate_key_lines <- paste0(candidate_numbers, "  ", candidate_keys)
display_labels[candidate_ids] <- candidate_numbers

group_end <- c(length(effector_ids), length(effector_ids) + length(module_ids), nrow(sectors))
gaps <- rep(2.0, nrow(sectors))
gaps[group_end] <- c(9, 9, 9)

draw_chord <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(4.4, 1.25))
  par(mar = c(2.2, 2.2, 4.0, 1.0), family = "sans", xpd = NA)
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.after = gaps, track.margin = c(0.005, 0.005))
  circlize::chordDiagram(
    x = chord_data,
    order = sectors$sector_id,
    grid.col = sector_colors,
    col = edge_colors,
    transparency = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.18),
    directional = 1,
    direction.type = c("arrows", "diffHeight"),
    link.arr.type = "big.arrow",
    diffHeight = circlize::mm_h(2),
    link.sort = TRUE,
    link.largest.ontop = TRUE
  )
  circlize::circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector_id <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      sector_type <- sectors$sector_type[match(sector_id, sectors$sector_id)]
      label_cex <- if (sector_type == "candidate") 0.62 else 0.56
      circlize::circos.text(
        mean(xlim),
        mean(ylim),
        display_labels[[sector_id]],
        facing = "bending.inside",
        niceFacing = TRUE,
        adj = c(0.5, 0.5),
        cex = label_cex,
        col = palette$ink
      )
    },
    bg.border = NA
  )
  title("Predicted nematode classes, pine defense modules, and prioritized orthogroups", cex.main = 1.15, font.main = 2)
  mtext(
    "Links encode curated functional hypotheses and module membership; they are not direct protein-protein interactions.",
    side = 1,
    line = 0.4,
    cex = 0.70,
    col = "#444444"
  )
  circlize::circos.clear()

  par(mar = c(4, 0.2, 4, 0.5), family = "sans")
  plot.new()
  legend(
    "topleft",
    legend = c(
      "Nematode candidate class", "Host defense module",
      "Leading host candidate", "Supported host candidate"
    ),
    fill = c(palette$effector[[1]], palette$module[[1]], palette$leading, palette$supported),
    border = NA,
    bty = "n",
    cex = 0.76
  )
  text(0, 0.70, "Candidate key", adj = 0, font = 2, cex = 0.78)
  key_y <- seq(0.66, 0.25, length.out = length(candidate_key_lines))
  text(0, key_y, candidate_key_lines, adj = 0, cex = 0.57, col = palette$ink)
  text(
    0,
    0.19,
    "Evidence codes\nE2/E1: concordant expression in 2/1 species\nN: curated or experimental network context\nS: structurally supported nematode anchor\nQ: sequence-resolved host candidate",
    adj = 0,
    cex = 0.60,
    col = "#444444"
  )
  text(0, 0.05, "Link widths aid readability;\nthey do not encode binding strength.", adj = 0, cex = 0.62, col = "#444444")
  layout(1)
}

for (path in c(output_pdf, output_png)) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
pdf(output_pdf, width = 13.5, height = 10.5, family = "sans", useDingbats = FALSE)
draw_chord()
dev.off()
png(output_png, width = 13.5, height = 10.5, units = "in", res = 300, type = "cairo")
draw_chord()
dev.off()

message("Wrote mechanism chord figure: ", output_pdf, ", ", output_png)
