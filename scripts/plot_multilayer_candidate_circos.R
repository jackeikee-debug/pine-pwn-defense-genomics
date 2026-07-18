#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste(
    "Usage: plot_multilayer_candidate_circos.R <candidates.tsv>",
    "<effector_context.tsv> <edges.tsv> <sectors.tsv> <output.pdf> <output.png>"
  ))
}
if (!requireNamespace("circlize", quietly = TRUE)) stop("The R package 'circlize' is required")

candidates <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
effector_context <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
edges <- read.delim(args[[3]], stringsAsFactors = FALSE, check.names = FALSE)
sectors <- read.delim(args[[4]], stringsAsFactors = FALSE, check.names = FALSE)
output_pdf <- args[[5]]
output_png <- args[[6]]

contrast_ids <- c("resistant_pwn_vs_water", "susceptible_pwn_vs_water", "genotype_by_inoculum")
species_ids <- c("pden", "pmas", "ptab", "plam", "plon", "ptae")
contrast_short <- c("R", "S", "GxI")
track_labels <- LETTERS[1:8]
COPY_TRANSFORM_TOLERANCE <- 1e-8
DEVICE_WIDTH_IN <- 10.5
DEVICE_HEIGHT_IN <- 8.0
FINAL_WIDTH_MM <- 180
DEVICE_POINTSIZE <- 16
MIN_CEX <- 7 / (DEVICE_POINTSIZE * FINAL_WIDTH_MM / (DEVICE_WIDTH_IN * 25.4))
figure_title <- "Integrated effector-prior and host-defense evidence landscape"
copy_claim_text <- paste(
  strwrap(
    "Species values are within-orthogroup centered descriptive copy context; not comparable across families; not expansion/contraction calls.",
    width = 55
  ),
  collapse = "\n"
)
functional_annotation_label <- paste(
  strwrap("Effector functional annotation (supported/evaluated)", width = 48),
  collapse = "\n"
)
east_asia_label <- paste(strwrap("East Asia: pden, pmas, ptab", width = 22), collapse = "\n")
north_america_label <- paste(strwrap("North America: plam, plon, ptae", width = 22), collapse = "\n")
copy_scale_label <- "Centered log2 copy"
expression_scale_label <- "Expression log2FC"
significance_scale_label <- "-log10(candidate-set FDR)"
contrast_profile_note <- "D: R, S, and GxI form a categorical contrast profile; not a time series."
required_candidates <- c(
  "orthogroup_id", "priority_class", "regional_log2_ratio", "cross_species_state",
  "network_support", "structure_support", "sequence_support", "mapped_symbols",
  "copy_context_claim_boundary",
  unlist(lapply(species_ids, function(species_id) c(
    paste0("copy_", species_id), paste0("centered_log2_copy_", species_id)
  ))),
  unlist(lapply(contrast_ids, function(id) c(
    paste0(id, "_representative_log2_fold_change"), paste0(id, "_simes_padj"), paste0(id, "_state")
  )))
)
required_context <- c(
  "effector_class", "standard_candidate_count", "functional_evaluated_count",
  "functional_supported_count", "functional_supported_fraction",
  "functional_annotation_state", "claim_boundary"
)
required_edges <- c("source_id", "target_id", "edge_type", "weight", "evidence_boundary")
required_sectors <- c("sector_id", "sector_type", "label", "display_order", "priority_class")
if (!all(required_candidates %in% names(candidates))) stop("Candidate table is missing required columns")
if (!all(required_context %in% names(effector_context))) stop("Effector-context table is missing required columns")
if (!all(required_edges %in% names(edges))) stop("Edge table is missing required columns")
if (!all(required_sectors %in% names(sectors))) stop("Sector table is missing required columns")
if (nrow(candidates) != 12 || anyDuplicated(candidates$orthogroup_id)) stop("Expected exactly 12 unique candidates")
if (anyDuplicated(effector_context$effector_class)) stop("Effector-context classes must be unique")
if (anyDuplicated(edges[c("source_id", "target_id", "edge_type")])) stop("Chord edges must be unique")
allowed_edge_types <- c("functional_prior", "candidate_membership")
if (any(!edges$edge_type %in% allowed_edge_types)) {
  stop("Unsupported chord edge type; only functional_prior and candidate_membership are allowed")
}
if (any(edges$evidence_boundary %in% c("direct_interaction", "validated_binding"))) {
  stop("Direct interaction links are outside scope")
}

copy_columns <- paste0("copy_", species_ids)
centered_copy_columns <- paste0("centered_log2_copy_", species_ids)
raw_copy_values <- as.matrix(data.frame(lapply(candidates[copy_columns], as.numeric), check.names = FALSE))
centered_copy_values <- as.matrix(data.frame(lapply(candidates[centered_copy_columns], as.numeric), check.names = FALSE))
if (any(!is.finite(raw_copy_values)) || any(raw_copy_values < 0)) {
  stop("Raw species copy counts must be finite and nonnegative")
}
if (any(!is.finite(centered_copy_values))) {
  stop("Centered species copy values must be finite")
}
if (any(abs(rowMeans(centered_copy_values)) > 1e-8)) {
  stop("Each candidate's six centered copy values must average to zero")
}
# Allows builder TSV rounding while requiring every centered value to reproduce the documented transform.
transformed_copy_values <- log2(raw_copy_values + 0.5)
expected_centered_copy_values <- transformed_copy_values - rowMeans(transformed_copy_values)
if (any(abs(centered_copy_values - expected_centered_copy_values) > COPY_TRANSFORM_TOLERANCE)) {
  stop("Centered species copy values disagree with log2(raw + 0.5) within-orthogroup centering")
}
if (any(!nzchar(trimws(as.character(candidates$copy_context_claim_boundary))))) {
  stop("Each candidate requires a copy-context claim boundary")
}
# Scope boundary: No chromosome coordinates, Mb labels, synteny, direct targeting, or expansion/contraction inference.

functional_edges <- edges[edges$edge_type == "functional_prior", , drop = FALSE]
membership_edges <- edges[edges$edge_type == "candidate_membership", , drop = FALSE]
effector_ids <- unique(functional_edges$source_id)
module_ids <- unique(c(functional_edges$target_id, membership_edges$source_id))
candidate_ids <- candidates$orthogroup_id
if (!setequal(effector_ids, effector_context$effector_class)) {
  stop("Every linked effector class requires one context row")
}
effector_context <- effector_context[match(effector_ids, effector_context$effector_class), , drop = FALSE]
context_standard <- suppressWarnings(as.numeric(effector_context$standard_candidate_count))
context_evaluated <- suppressWarnings(as.numeric(effector_context$functional_evaluated_count))
context_supported <- suppressWarnings(as.numeric(effector_context$functional_supported_count))
context_fraction <- suppressWarnings(as.numeric(effector_context$functional_supported_fraction))
context_counts <- cbind(context_standard, context_evaluated, context_supported)
if (any(!is.finite(context_counts)) || any(context_counts < 0) ||
    any(context_counts != round(context_counts))) {
  stop("Effector functional-annotation counts must be finite nonnegative integers")
}
if (any(context_supported > context_evaluated | context_evaluated > context_standard)) {
  stop("Effector functional-annotation counts require supported <= evaluated <= standard")
}
expected_context_state <- ifelse(
  context_supported > 0, "supported",
  ifelse(context_evaluated > 0, "evaluated_none", "unavailable")
)
if (any(effector_context$functional_annotation_state != expected_context_state)) {
  stop("Effector functional-annotation state disagrees with its counts")
}
evaluated_context <- context_evaluated > 0
if (any(evaluated_context & (
  !is.finite(context_fraction) |
    abs(context_fraction - context_supported / context_evaluated) > 1e-8
))) {
  stop("Effector functional-annotation fractions must equal supported/evaluated")
}
if (any(!evaluated_context & !is.na(context_fraction))) {
  stop("Effector functional-annotation fractions must be unavailable when no proteins were evaluated")
}

all_sectors <- sectors
used_ids <- candidate_ids
if (anyDuplicated(used_ids)) stop("Candidate sectors must not overlap")
sectors <- sectors[match(used_ids, sectors$sector_id), , drop = FALSE]
if (any(is.na(sectors$sector_id))) stop("Every plotted candidate identifier must have a sector")
if (!all(sectors$sector_type == "candidate")) stop("Candidate-centered Circos requires candidate sectors only")
candidates <- candidates[match(candidate_ids, candidates$orthogroup_id), , drop = FALSE]

palette <- list(
  effector = c("#B68A31", "#867A3D", "#537D55", "#2F7B68", "#267487"),
  module = c("#426B91", "#27818B", "#39806A", "#72754E", "#8D6448"),
  leading = "#BB4E43", supported = "#D78554", ink = "#222222",
  missing = "#FFFFFF", nonsignificant = "#E2E2E2", pale = "#F7F7F7"
)
expression_color <- circlize::colorRamp2(c(-3, 0, 3), c("#2166AC", "#F7F7F7", "#B2182B"))
regional_color <- circlize::colorRamp2(c(-2, 0, 2), c("#2C7FB8", "#F7F7F7", "#D99A32"))
copy_clip <- 2
copy_color <- circlize::colorRamp2(c(-copy_clip, 0, copy_clip), c("#2166AC", "#F7F7F7", "#B2182B"))
cross_colors <- c(
  concordant = "#2F8F68", heterogeneous = "#D89A32", discordant = "#B8504B",
  nonsignificant = "#D7D7D7", unavailable = "#F2F2F2"
)
functional_state_colors <- c(
  supported = "#EEEEEA", evaluated_none = palette$nonsignificant,
  unavailable = palette$missing
)

sector_colors <- setNames(rep("#BBBBBB", length(used_ids)), used_ids)
sector_colors[candidate_ids] <- ifelse(
  candidates$priority_class == "mechanism_leading", palette$leading, palette$supported
)

gaps <- rep(1.5, length(used_ids))
gaps[length(used_ids)] <- 7
candidate_numbers <- setNames(sprintf("%02d", seq_along(candidate_ids)), candidate_ids)
track_label_sector <- tail(candidate_ids, 1)
effector_numbers <- setNames(paste0("E", seq_along(effector_ids)), effector_ids)
module_numbers <- setNames(paste0("M", seq_along(module_ids)), module_ids)
module_member_counts <- table(factor(membership_edges$source_id, levels = module_ids))
candidate_row <- function(sector_id) candidates[match(sector_id, candidates$orthogroup_id), , drop = FALSE]
context_row <- function(sector_id) effector_context[match(sector_id, effector_context$effector_class), , drop = FALSE]

short_labels <- c(
  cell_wall_modifying = "Cell-wall modifying",
  detoxification = "Detoxification",
  mimicry = "Mimicry",
  protease = "Protease",
  stress_adaptation = "Stress adaptation",
  hydraulic_xylem = "Hydraulic/xylem",
  immune_signaling = "Immune regulation",
  phenylpropanoid_lignin = "Phenylpropanoid/lignin",
  ros_detoxification = "ROS control",
  wound_periderm = "Wound/periderm"
)

draw_cell <- function(fill, border = "#C8C8C8", lwd = 0.35) {
  xlim <- circlize::get.cell.meta.data("xlim")
  ylim <- circlize::get.cell.meta.data("ylim")
  circlize::circos.rect(xlim[1], ylim[1], xlim[2], ylim[2], col = fill, border = border, lwd = lwd)
}

# Place one A-H index in the enlarged gap after the final sector so track labels do not mask data.
draw_track_index_label <- function(letter, y) {
  sector_id <- circlize::get.cell.meta.data("sector.index")
  if (sector_id != track_label_sector) return(invisible(NULL))
  xlim <- circlize::get.cell.meta.data("xlim")
  label_x <- xlim[2] + diff(xlim) * 0.15
  circlize::circos.points(
    label_x, y, pch = 21, cex = 0.86,
    col = "#555555", bg = "#FFFFFF"
  )
  circlize::circos.text(
    label_x, y, letter, cex = 0.40, font = 2, col = palette$ink,
    facing = "inside", niceFacing = TRUE
  )
}

contrast_colors <- c("#2C7FB8", "#D95F4C", "#7A5AA6")
names(contrast_colors) <- contrast_ids
module_colors <- setNames(palette$module[seq_along(module_ids)], module_ids)
regional_bar_colors <- c(North_America = "#2C7FB8", East_Asia = "#D99A32")
expression_clip <- 6
significance_clip <- 8
significance_reference <- -log10(0.05)

draw_label_track <- function() {
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.090,
    panel.fun = function(x, y) {
      sector_id <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      circlize::circos.text(
        mean(xlim), 0.48, candidate_numbers[[sector_id]],
        facing = "bending.inside", niceFacing = TRUE, cex = 0.68, col = palette$ink
      )
    },
    bg.col = "#FAFAF8", bg.border = "#B9B9B5"
  )
}

draw_evidence_tier_track <- function() {
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.030,
    panel.fun = function(x, y) {
      sector_id <- circlize::get.cell.meta.data("sector.index")
      draw_cell(sector_colors[[sector_id]], border = "#FFFFFF", lwd = 0.5)
      draw_track_index_label(track_labels[[1]], 0.5)
    },
    bg.border = NA
  )
}

draw_copy_heatmap_tracks <- function() {
  circlize::circos.trackPlotRegion(
    ylim = c(0, length(species_ids)), track.height = 0.150,
    panel.fun = function(x, y) {
      sector_id <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      row <- candidate_row(sector_id)
      for (species_index in seq_along(species_ids)) {
        species_id <- species_ids[[species_index]]
        value <- as.numeric(row[[paste0("centered_log2_copy_", species_id)]])
        y0 <- length(species_ids) - species_index
        circlize::circos.rect(
          xlim[1], y0, xlim[2], y0 + 1,
          col = copy_color(max(-copy_clip, min(copy_clip, value))),
          border = "#F4F4F1", lwd = 0.28
        )
      }
      circlize::circos.lines(xlim, c(3, 3), col = "#575757", lwd = 0.75)
      draw_track_index_label(track_labels[[2]], 3)
    },
    bg.border = "#A8A8A4"
  )
}

draw_regional_bar_track <- function() {
  circlize::circos.trackPlotRegion(
    ylim = c(-2, 2), track.height = 0.065,
    panel.fun = function(x, y) {
      sector_id <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      value <- as.numeric(candidate_row(sector_id)$regional_log2_ratio)
      circlize::circos.lines(xlim, c(0, 0), col = "#8A8A8A", lwd = 0.45)
      if (is.finite(value)) {
        clipped <- max(-2, min(2, value))
        fill <- if (clipped >= 0) regional_bar_colors[["East_Asia"]] else regional_bar_colors[["North_America"]]
        inset <- diff(xlim) * 0.23
        circlize::circos.rect(
          xlim[1] + inset, min(0, clipped), xlim[2] - inset, max(0, clipped),
          col = fill, border = NA
        )
      }
      draw_track_index_label(track_labels[[3]], 0)
    },
    bg.col = "#F8F8F5", bg.border = "#B9B9B5"
  )
}

candidate_contrast_values <- function(row, suffix) {
  vapply(contrast_ids, function(contrast_id) {
    suppressWarnings(as.numeric(row[[paste0(contrast_id, suffix)]]))
  }, numeric(1))
}

draw_expression_profile_track <- function() {
  circlize::circos.trackPlotRegion(
    ylim = c(-expression_clip, expression_clip), track.height = 0.135,
    panel.fun = function(x, y) {
      sector_id <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      row <- candidate_row(sector_id)
      values <- candidate_contrast_values(row, "_representative_log2_fold_change")
      states <- vapply(contrast_ids, function(contrast_id) {
        as.character(row[[paste0(contrast_id, "_state")]])
      }, character(1))
      finite <- is.finite(values) & states != "unavailable"
      xs <- seq(xlim[1] + diff(xlim) * 0.18, xlim[2] - diff(xlim) * 0.18, length.out = 3)
      clipped <- pmax(-expression_clip, pmin(expression_clip, values))
      circlize::circos.lines(xlim, c(0, 0), col = "#9A9A96", lwd = 0.45)
      if (sum(finite) >= 2) {
        circlize::circos.lines(xs[finite], clipped[finite], col = "#595959", lwd = 0.75)
      }
      if (any(finite)) {
        circlize::circos.points(
          xs[finite], clipped[finite], pch = 16, cex = 0.48,
          col = contrast_colors[contrast_ids[finite]]
        )
      }
      if (any(!finite)) {
        circlize::circos.points(xs[!finite], rep(0, sum(!finite)), pch = 4, cex = 0.42, col = "#8F8F8F")
      }
      draw_track_index_label(track_labels[[4]], 0)
    },
    bg.col = "#FBFBF8", bg.border = "#B9B9B5"
  )
}

draw_significance_track <- function() {
  circlize::circos.trackPlotRegion(
    ylim = c(0, significance_clip), track.height = 0.080,
    panel.fun = function(x, y) {
      sector_id <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      row <- candidate_row(sector_id)
      padj <- candidate_contrast_values(row, "_simes_padj")
      finite <- is.finite(padj) & padj > 0
      xs <- seq(xlim[1] + diff(xlim) * 0.18, xlim[2] - diff(xlim) * 0.18, length.out = 3)
      scores <- pmin(significance_clip, -log10(pmax(padj, 10^-significance_clip)))
      circlize::circos.lines(xlim, rep(significance_reference, 2), col = "#777777", lty = 3, lwd = 0.45)
      for (i in which(finite)) {
        circlize::circos.lines(c(xs[[i]], xs[[i]]), c(0, scores[[i]]), col = contrast_colors[[i]], lwd = 0.7)
        circlize::circos.points(
          xs[[i]], scores[[i]], pch = if (padj[[i]] < 0.05) 16 else 1,
          cex = 0.38, col = contrast_colors[[i]]
        )
      }
      draw_track_index_label(track_labels[[5]], significance_clip / 2)
    },
    bg.col = "#F8F8F5", bg.border = "#B9B9B5"
  )
}

draw_cross_track <- function() {
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.033,
    panel.fun = function(x, y) {
      state <- candidate_row(circlize::get.cell.meta.data("sector.index"))$cross_species_state
      draw_cell(cross_colors[[state]], border = "#FFFFFF", lwd = 0.35)
      draw_track_index_label(track_labels[[6]], 0.5)
    },
    bg.border = NA
  )
}

draw_evidence_glyph_track <- function() {
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.045,
    panel.fun = function(x, y) {
      sector_id <- circlize::get.cell.meta.data("sector.index")
      row <- candidate_row(sector_id)
      xlim <- circlize::get.cell.meta.data("xlim")
      xs <- seq(xlim[1] + diff(xlim) * 0.22, xlim[2] - diff(xlim) * 0.22, length.out = 3)
      support <- c(row$network_support, row$structure_support, row$sequence_support) == "yes"
      cols <- ifelse(support, c("#4D7EA8", "#8D6AA8", "#2F8F68"), "#D8D8D8")
      draw_cell("#FFFFFF", border = "#D0D0CC")
      circlize::circos.points(xs, rep(0.5, 3), pch = c(16, 17, 15), cex = 0.42, col = cols)
      draw_track_index_label(track_labels[[7]], 0.5)
    },
    bg.border = NA
  )
}

draw_module_membership_track <- function() {
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.050,
    panel.fun = function(x, y) {
      sector_id <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      memberships <- membership_edges$source_id[membership_edges$target_id == sector_id]
      breaks <- seq(xlim[1], xlim[2], length.out = length(module_ids) + 1)
      for (module_index in seq_along(module_ids)) {
        module_id <- module_ids[[module_index]]
        fill <- if (module_id %in% memberships) module_colors[[module_id]] else "#ECECE8"
        circlize::circos.rect(
          breaks[[module_index]], 0, breaks[[module_index + 1]], 1,
          col = fill, border = "#FFFFFF", lwd = 0.28
        )
      }
      draw_track_index_label(track_labels[[8]], 0.5)
    },
    bg.border = "#B9B9B5"
  )
}

build_shared_module_links <- function() {
  links <- lapply(module_ids, function(module_id) {
    members <- membership_edges$target_id[membership_edges$source_id == module_id]
    members <- candidate_ids[candidate_ids %in% members]
    if (length(members) < 2) return(NULL)
    data.frame(
      module_id = module_id,
      from = members[-length(members)],
      to = members[-1],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, links)
}

draw_shared_module_links <- function() {
  shared_links <- build_shared_module_links()
  if (is.null(shared_links) || nrow(shared_links) == 0) return(invisible(NULL))
  for (i in seq_len(nrow(shared_links))) {
    module_id <- shared_links$module_id[[i]]
    circlize::circos.link(
      shared_links$from[[i]], 0.5, shared_links$to[[i]], 0.5,
      col = grDevices::adjustcolor(module_colors[[module_id]], alpha.f = 0.42),
      border = NA, lwd = 1.5
    )
  }
  invisible(shared_links)
}

draw_effector_context_panel <- function() {
  par(mar = c(0.7, 0.1, 1.8, 0.3), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0, 0.98, "PWN effector-class context", adj = c(0, 1), font = 2, cex = 0.84)
  text(0, 0.905, "Audited class-level counts", adj = c(0, 1), cex = MIN_CEX, col = "#555555")
  text(0.995, 0.905, "standard/evaluated/supported", adj = c(1, 1), cex = MIN_CEX, col = "#555555")
  max_count <- max(1, context_standard)
  bar_left <- 0.38
  bar_width <- 0.50
  ys <- seq(0.78, 0.28, length.out = length(effector_ids))
  for (i in seq_along(effector_ids)) {
    effector_id <- effector_ids[[i]]
    text(0, ys[[i]], paste0(effector_numbers[[effector_id]], "  ", short_labels[[effector_id]]), adj = 0, cex = MIN_CEX)
    rect(bar_left, ys[[i]] - 0.034, bar_left + bar_width, ys[[i]] + 0.034, col = "#F0F0EC", border = NA)
    standard_x <- bar_left + bar_width * context_standard[[i]] / max_count
    evaluated_x <- bar_left + bar_width * context_evaluated[[i]] / max_count
    supported_x <- bar_left + bar_width * context_supported[[i]] / max_count
    rect(
      bar_left, ys[[i]] - 0.034, standard_x, ys[[i]] + 0.034,
      col = grDevices::adjustcolor(palette$effector[[i]], alpha.f = 0.55), border = NA
    )
    segments(bar_left, ys[[i]], evaluated_x, ys[[i]], col = palette$effector[[i]], lwd = 2.2)
    points(supported_x, ys[[i]], pch = 16, cex = 0.48, col = palette$ink)
    text(0.995, ys[[i]], paste(context_standard[[i]], context_evaluated[[i]], context_supported[[i]], sep = "/"), adj = 1, cex = MIN_CEX)
  }
  rect(0.00, 0.13, 0.04, 0.165, col = grDevices::adjustcolor("#537D55", alpha.f = 0.55), border = NA)
  text(0.055, 0.147, "Standard candidates", adj = 0, cex = MIN_CEX)
  segments(0.39, 0.147, 0.45, 0.147, col = "#537D55", lwd = 2.2)
  text(0.465, 0.147, "Functionally evaluated", adj = 0, cex = MIN_CEX)
  points(0.39, 0.065, pch = 16, cex = 0.48)
  text(0.415, 0.065, "Functionally supported", adj = 0, cex = MIN_CEX)
}

draw_legend <- function() {
  par(mar = c(0.6, 0.1, 0.5, 0.3), family = "sans")
  plot.new()
  text(0, 0.98, "Integrated evidence tracks", adj = c(0, 1), font = 2, cex = 0.80)
  track_lines <- c(
    paste0(track_labels[1], "  Candidate priority"),
    paste0(track_labels[2], "  Six-species copy heatmap"),
    paste0(track_labels[3], "  Regional copy-ratio bar"),
    paste0(track_labels[4], "  Expression profile (R / S / GxI)"),
    paste0(track_labels[5], "  ", significance_scale_label),
    paste0(track_labels[6], "  Cross-species expression state"),
    paste0(track_labels[7], "  Network / structure / sequence"),
    paste0(track_labels[8], "  Host-module membership")
  )
  text(0, seq(0.91, 0.69, length.out = length(track_lines)), track_lines, adj = 0, cex = MIN_CEX)
  text(0, 0.650, contrast_profile_note, adj = 0, cex = MIN_CEX, col = "#555555")

  text(0, 0.610, "Species heatmap order", adj = 0, font = 2, cex = 0.66)
  text(0, 0.570, paste(species_ids[1:3], collapse = "  "), adj = 0, cex = MIN_CEX)
  text(0.52, 0.570, paste(species_ids[4:6], collapse = "  "), adj = 0, cex = MIN_CEX)
  text(0, 0.530, copy_scale_label, adj = 0, font = 2, cex = MIN_CEX)
  copy_breaks <- seq(0, 0.39, length.out = 10)
  rect(copy_breaks[1:9], 0.480, copy_breaks[2:10], 0.497,
       col = copy_color(seq(-copy_clip, copy_clip, length.out = 9)), border = NA)
  text(c(0, 0.195, 0.39), 0.453, c("-2", "0", "+2"), cex = MIN_CEX)

  text(0.52, 0.530, expression_scale_label, adj = 0, font = 2, cex = MIN_CEX)
  points(c(0.54, 0.68, 0.82), rep(0.490, 3), pch = 16, cex = 0.55, col = contrast_colors)
  text(c(0.56, 0.70, 0.84), rep(0.490, 3), contrast_short, adj = 0, cex = MIN_CEX)
  segments(0.54, 0.450, 0.60, 0.450, col = "#777777", lty = 3)
  text(0.62, 0.450, "FDR 0.05 reference", adj = 0, cex = MIN_CEX)
  points(c(0.54, 0.72), rep(0.415, 2), pch = c(16, 1), cex = 0.50, col = "#555555")
  text(c(0.56, 0.74), rep(0.415, 2), c("FDR < 0.05", "tested nonsignificant"), adj = 0, cex = MIN_CEX)

  text(0, 0.375, "Host-module colors", adj = 0, font = 2, cex = 0.66)
  module_key <- paste0(module_numbers[module_ids], " ", short_labels[module_ids])
  module_y <- seq(0.330, 0.190, length.out = length(module_ids))
  rect(rep(0.00, length(module_ids)), module_y - 0.012, rep(0.035, length(module_ids)), module_y + 0.012,
       col = module_colors[module_ids], border = NA)
  text(rep(0.05, length(module_ids)), module_y, module_key, adj = 0, cex = MIN_CEX)
  text(0.52, 0.375, "Cross-species state", adj = 0, font = 2, cex = 0.66)
  cross_ids <- c("concordant", "heterogeneous", "discordant", "unavailable")
  cross_y <- seq(0.330, 0.225, length.out = length(cross_ids))
  rect(rep(0.52, length(cross_ids)), cross_y - 0.012, rep(0.555, length(cross_ids)), cross_y + 0.012,
       col = cross_colors[cross_ids], border = "#B0B0B0")
  text(rep(0.57, length(cross_ids)), cross_y, cross_ids, adj = 0, cex = MIN_CEX)

  text(0, 0.145, "Scope boundary", adj = 0, font = 2, cex = 0.66)
  text(0, 0.108, "Center links: shared module membership; not interaction.", adj = 0, cex = MIN_CEX, col = "#555555")
  text(0, 0.068, "No direct effector-target inference or causal resistance claim.", adj = 0, cex = MIN_CEX, col = "#555555")
  text(0, 0.028, copy_claim_text, adj = c(0, 1), cex = MIN_CEX, col = "#555555")
}

draw_candidate_key <- function() {
  par(mar = c(0.7, 0.1, 0.2, 0.3), family = "sans")
  plot.new()
  text(0, 0.98, "Candidate orthogroups", adj = c(0, 1), font = 2, cex = 0.78)
  labels <- paste0(candidate_numbers, "  ", candidate_ids)
  symbols <- candidates$mapped_symbols
  labels[nzchar(symbols)] <- paste0(labels[nzchar(symbols)], " | ", sub(";.*", "", symbols[nzchar(symbols)]))
  text(rep(c(0, 0.52), each = 6), rep(seq(0.86, 0.14, length.out = 6), 2), labels, adj = 0, cex = MIN_CEX)
  rect(0.00, 0.035, 0.035, 0.075, col = palette$leading, border = NA)
  text(0.05, 0.055, "leading", adj = 0, cex = MIN_CEX)
  rect(0.27, 0.035, 0.305, 0.075, col = palette$supported, border = NA)
  text(0.32, 0.055, "supported", adj = 0, cex = MIN_CEX)
}

draw_figure <- function() {
  layout(
    matrix(c(1, 2, 1, 3, 1, 4), nrow = 3, byrow = TRUE),
    widths = c(5.2, 2.8), heights = c(2.15, 3.35, 2.5)
  )
  par(mar = c(1.7, 0.8, 3.4, 0.4), family = "sans", xpd = NA)
  circlize::circos.clear()
  circlize::circos.par(
    start.degree = 92, gap.after = gaps, track.margin = c(0.002, 0.002),
    cell.padding = c(0, 0, 0, 0), points.overflow.warning = FALSE
  )
  circlize::circos.initialize(
    factors = factor(candidate_ids, levels = candidate_ids),
    xlim = c(0, 1)
  )
  draw_label_track()
  draw_evidence_tier_track()
  draw_copy_heatmap_tracks()
  draw_regional_bar_track()
  draw_expression_profile_track()
  draw_significance_track()
  draw_cross_track()
  draw_evidence_glyph_track()
  draw_module_membership_track()
  draw_shared_module_links()
  mtext(figure_title, side = 3, line = 1.35, adj = 0, cex = 0.90, font = 2)
  mtext(
    "Orthogroup-level synthesis; quantitative tracks share the frozen Circos v1 evidence table.",
    side = 1, line = 0.50, cex = MIN_CEX, col = "#444444"
  )
  circlize::circos.clear()
  draw_effector_context_panel()
  draw_legend()
  draw_candidate_key()
  layout(1)
}

for (path in c(output_pdf, output_png)) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
pdf(
  output_pdf, width = DEVICE_WIDTH_IN, height = DEVICE_HEIGHT_IN,
  pointsize = DEVICE_POINTSIZE, family = "sans", useDingbats = FALSE
)
draw_figure()
dev.off()
png(
  output_png, width = DEVICE_WIDTH_IN, height = DEVICE_HEIGHT_IN, units = "in",
  res = 600, pointsize = DEVICE_POINTSIZE, type = "cairo"
)
draw_figure()
dev.off()
message("Wrote refined multilayer candidate Circos: ", output_pdf, ", ", output_png)
