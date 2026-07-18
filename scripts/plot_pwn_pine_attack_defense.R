#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(
    paste(
      "Usage: plot_pwn_pine_attack_defense.R",
      "<bubble.tsv> <nodes.tsv> <edges.tsv> <output.pdf> <output.png>"
    )
  )
}

bubble_path <- args[[1]]
nodes_path <- args[[2]]
edges_path <- args[[3]]
pdf_path <- args[[4]]
png_path <- args[[5]]

read_table <- function(path) {
  if (!file.exists(path)) {
    stop("Input table does not exist: ", path)
  }
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

require_columns <- function(data, required, table_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(table_name, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

bubble <- read_table(bubble_path)
nodes <- read_table(nodes_path)
edges <- read_table(edges_path)

bubble_columns <- c(
  "display_order", "orthogroup_id", "candidate_label", "contrast_id",
  "contrast_label", "log2_fold_change", "whole_transcriptome_padj", "whole_transcriptome_state",
  "candidate_set_simes_padj",
  "priority_class", "cross_species_state", "network_support",
  "structure_support", "sequence_support"
)
node_columns <- c(
  "node_id", "node_label", "node_layer", "display_order", "color_group",
  "evidence_boundary"
)
edge_columns <- c(
  "source_id", "target_id", "source_layer", "target_layer", "weight",
  "edge_class", "color_group", "statistical_state", "evidence_boundary"
)
require_columns(bubble, bubble_columns, "Bubble table")
require_columns(nodes, node_columns, "Node table")
require_columns(edges, edge_columns, "Edge table")

validate_inputs <- function(bubble, nodes, edges) {
  bubble$display_order <- suppressWarnings(as.integer(bubble$display_order))
  bubble$log2_fold_change <- suppressWarnings(as.numeric(bubble$log2_fold_change))
  bubble$whole_transcriptome_padj <- suppressWarnings(as.numeric(bubble$whole_transcriptome_padj))
  nodes$display_order <- suppressWarnings(as.integer(nodes$display_order))
  edges$weight <- suppressWarnings(as.numeric(edges$weight))

  candidate_ids <- unique(bubble$orthogroup_id)
  if (nrow(bubble) != 36 || length(candidate_ids) != 12) {
    stop("Bubble table must contain 36 rows for exactly 12 candidates")
  }
  if (any(is.na(bubble$display_order)) || length(unique(bubble$display_order)) != 12) {
    stop("Bubble candidate display_order values must define 12 ordered candidates")
  }
  if (!setequal(unique(bubble$contrast_label), c("R", "S", "GxI"))) {
    stop("Bubble table must contain only R, S, and GxI contrasts")
  }
  contrast_counts <- table(bubble$orthogroup_id, bubble$contrast_label)
  if (!all(contrast_counts == 1)) {
    stop("Every candidate must have exactly one row for each contrast")
  }
  if (!all(bubble$whole_transcriptome_state %in% c("significant", "not_significant", "unavailable"))) {
    stop("Bubble table contains an unsupported significance state")
  }
  available <- bubble$whole_transcriptome_state != "unavailable"
  if (any(!is.finite(bubble$log2_fold_change[available]))) {
    stop("Available bubble rows require finite log2 fold changes")
  }
  if (any(!is.finite(bubble$whole_transcriptome_padj[available])) ||
      any(bubble$whole_transcriptome_padj[available] < 0 | bubble$whole_transcriptome_padj[available] > 1)) {
    stop("Available bubble rows require whole-transcriptome adjusted P values in [0, 1]")
  }

  if (any(is.na(nodes$display_order)) || any(nodes$node_id == "") ||
      anyDuplicated(nodes$node_id)) {
    stop("Alluvial nodes require unique non-empty IDs and finite display orders")
  }
  expected_layers <- c("effector_class", "host_module", "candidate", "response")
  if (!setequal(unique(nodes$node_layer), expected_layers)) {
    stop("Alluvial nodes must define the four expected layers")
  }
  if (any(!is.finite(edges$weight)) || any(edges$weight <= 0)) {
    stop("Alluvial edge weights must be positive and finite")
  }
  response_state <- edges$statistical_state[edges$edge_class == "response_summary"]
  if (any(!response_state %in% c("significant", "not_significant", "unavailable"))) {
    stop("Response edges contain an unsupported whole-transcriptome FDR status")
  }
  nonresponse_state <- edges$statistical_state[edges$edge_class != "response_summary"]
  if (any(nonresponse_state != "not_applicable")) {
    stop("Non-response edges must use statistical_state=not_applicable")
  }
  if (any(!edges$source_id %in% nodes$node_id) || any(!edges$target_id %in% nodes$node_id)) {
    stop("Alluvial edges reference unknown nodes")
  }
  candidate_nodes <- nodes$node_id[nodes$node_layer == "candidate"]
  if (!setequal(candidate_ids, candidate_nodes)) {
    stop("Panel A and Panel B candidate orthogroup sets must be identical")
  }
  node_layers <- setNames(nodes$node_layer, nodes$node_id)
  if (any(edges$source_layer != unname(node_layers[edges$source_id]))) {
    stop("Edge source_layer does not match source node layer")
  }
  if (any(edges$target_layer != unname(node_layers[edges$target_id]))) {
    stop("Edge target_layer does not match target node layer")
  }
  edge_routes <- paste(edges$source_layer, edges$target_layer, sep = "->")
  allowed_routes <- c(
    "effector_class->host_module",
    "host_module->candidate",
    "candidate->response"
  )
  if (any(!edge_routes %in% allowed_routes)) {
    stop("Unsupported alluvial edge transition")
  }
  prohibited <- grepl("direct_interaction|validated_binding", edges$evidence_boundary)
  if (any(prohibited)) {
    stop("Alluvial edge table contains a prohibited direct-interaction claim")
  }

  incoming <- tapply(edges$weight, edges$target_id, sum)
  outgoing <- tapply(edges$weight, edges$source_id, sum)
  flow <- function(values, ids) {
    result <- values[ids]
    result[is.na(result)] <- 0
    unname(result)
  }
  module_nodes <- nodes$node_id[nodes$node_layer == "host_module"]
  if (length(candidate_nodes) != 12 ||
      any(abs(flow(incoming, candidate_nodes) - 1) > 1e-8) ||
      any(abs(flow(outgoing, candidate_nodes) - 1) > 1e-8)) {
    stop("Each candidate must conserve one unit of normalized flow")
  }
  if (any(abs(flow(incoming, module_nodes) - flow(outgoing, module_nodes)) > 1e-8)) {
    stop("Host-module flow is not conserved")
  }

  list(bubble = bubble, nodes = nodes, edges = edges)
}

validated <- validate_inputs(bubble, nodes, edges)
bubble <- validated$bubble
nodes <- validated$nodes
edges <- validated$edges

module_colors <- c(
  hydraulic_xylem = "#3A7D78",
  immune_signaling = "#4E79A7",
  phenylpropanoid_lignin = "#6F9E55",
  ros_detoxification = "#C7665B",
  wound_periderm = "#7E709A"
)
response_colors <- c(
  r_response_greater = "#3B6FA1",
  s_response_greater = "#C6534F",
  interaction_unresolved = "#A9ADB1"
)
response_labels <- c(
  r_response_greater = "R-response greater",
  s_response_greater = "S-response greater",
  interaction_unresolved = "unresolved / padj unavailable"
)
neutral_ribbon <- "#B8BEC4"

response_nodes <- nodes[nodes$node_layer == "response", c("node_id", "node_label")]
if (!identical(
  setNames(response_nodes$node_label, response_nodes$node_id)[names(response_labels)],
  response_labels
)) {
  stop("Response nodes must use the declared GxI direction labels")
}

lookup_color <- function(values, palette, fallback = "#8B8F93") {
  colors <- unname(palette[values])
  colors[is.na(colors)] <- fallback
  colors
}

lfc_palette <- colorRampPalette(c("#376AA6", "#F7F7F5", "#C54E52"))(201)
lfc_color <- function(values, limit = 4) {
  clipped <- pmax(-limit, pmin(limit, values))
  indices <- round((clipped + limit) / (2 * limit) * 200) + 1
  colors <- lfc_palette[indices]
  colors[!is.finite(values)] <- "#FFFFFF"
  colors
}

draw_key_box <- function(x, y, width, height, fill, label, border = "#FFFFFF",
                         text_color = "#222222", cex = 0.48) {
  rect(x - width / 2, y - height / 2, x + width / 2, y + height / 2,
       col = fill, border = border, lwd = 0.65)
  text(x, y, labels = label, cex = cex, col = text_color)
}

draw_bubble_matrix <- function(bubble) {
  par(mar = c(0.15, 0.2, 0.15, 0.15), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 10.5), ylim = c(0, 14.5), xaxs = "i", yaxs = "i")

  contrast_order <- c("R", "S", "GxI")
  contrast_x <- setNames(c(3.62, 4.30, 4.98), contrast_order)
  candidate_meta <- bubble[order(bubble$display_order), ]
  candidate_meta <- candidate_meta[!duplicated(candidate_meta$orthogroup_id), ]
  row_y <- setNames(seq(12.25, 2.20, length.out = nrow(candidate_meta)), candidate_meta$orthogroup_id)

  text(0.08, 14.05, "A", adj = c(0, 0.5), font = 2, cex = 1.25)
  text(0.48, 14.05, "Host response bubble matrix", adj = c(0, 0.5),
       font = 2, cex = 0.95)
  text(2.70, 13.25, "Prioritized host orthogroup", adj = 1, font = 2, cex = 0.63)
  text(contrast_x, 13.25, labels = contrast_order, font = 2, cex = 0.69)
  text(4.30, 13.78, "P. massoniana contrasts", font = 2, cex = 0.62)
  text(6.18, 13.25, "Priority", font = 2, cex = 0.58)
  text(7.18, 13.25, "Cross-species", font = 2, cex = 0.58)
  text(8.18, 13.25, "Network", font = 2, cex = 0.52)
  text(8.80, 13.25, "Structure", font = 2, cex = 0.52)
  text(9.42, 13.25, "Sequence", font = 2, cex = 0.52)
  text(10.12, 13.25, "log2FC", font = 2, cex = 0.55)

  for (y in row_y) {
    segments(2.86, y, 9.72, y, col = "#ECEDEE", lwd = 0.55)
  }
  text(2.70, row_y[candidate_meta$orthogroup_id], labels = candidate_meta$candidate_label,
       adj = 1, cex = 0.59, col = "#202326")

  ordered <- bubble[order(bubble$display_order, match(bubble$contrast_label, contrast_order)), ]
  point_x <- unname(contrast_x[ordered$contrast_label])
  point_y <- unname(row_y[ordered$orthogroup_id])
  significance_score <- rep(0.35, nrow(ordered))
  finite_fdr <- is.finite(ordered$whole_transcriptome_padj)
  significance_score[finite_fdr] <- pmin(-log10(pmax(ordered$whole_transcriptome_padj[finite_fdr], 1e-8)), 6)
  radii <- sqrt(pmax(significance_score, 0.20))
  point_fill <- lfc_color(ordered$log2_fold_change)
  point_border <- ifelse(ordered$whole_transcriptome_state == "significant", "#202326", "#AEB3B7")
  point_width <- ifelse(ordered$whole_transcriptome_state == "significant", 1.45, 0.70)
  symbols(
    point_x, point_y, circles = radii, inches = 0.115, add = TRUE,
    bg = point_fill, fg = point_border, lwd = point_width
  )
  unavailable <- ordered$whole_transcriptome_state == "unavailable"
  if (any(unavailable)) {
    text(point_x[unavailable], point_y[unavailable], "x", cex = 0.58,
         font = 2, col = "#72777B")
  }

  priority_fill <- c(mechanism_leading = "#DCE9E4", mechanism_supported = "#E9E6F0")
  cross_fill <- c(concordant = "#DCE9E4", heterogeneous = "#F2E1D8")
  for (i in seq_len(nrow(candidate_meta))) {
    row <- candidate_meta[i, ]
    y <- row_y[[row$orthogroup_id]]
    priority_label <- if (row$priority_class == "mechanism_leading") "Lead" else "Support"
    cross_label <- if (row$cross_species_state == "concordant") "Concord." else "Heterog."
    draw_key_box(6.18, y, 0.78, 0.50, priority_fill[[row$priority_class]], priority_label)
    draw_key_box(7.18, y, 0.78, 0.50, cross_fill[[row$cross_species_state]], cross_label)
    support_values <- c(row$network_support, row$structure_support, row$sequence_support)
    support_x <- c(8.18, 8.80, 9.42)
    for (j in seq_along(support_values)) {
      supported <- identical(support_values[[j]], "yes")
      draw_key_box(
        support_x[[j]], y, 0.42, 0.50,
        if (supported) "#4E79A7" else "#F0F1F2",
        if (supported) "+" else "-",
        text_color = if (supported) "#FFFFFF" else "#777B7F",
        cex = 0.56
      )
    }
  }

  gradient_y <- seq(2.20, 12.25, length.out = 201)
  for (i in seq_len(200)) {
    rect(10.01, gradient_y[i], 10.20, gradient_y[i + 1],
         col = lfc_palette[i], border = NA)
  }
  rect(10.01, 2.20, 10.20, 12.25, border = "#8C9093", lwd = 0.6)
  text(10.27, c(2.20, 7.23, 12.25), labels = c("-4", "0", "+4"),
       adj = 0, cex = 0.48, col = "#55595D")

  legend_y <- 1.10
  text(0.48, legend_y, "Bubble area = clipped -log10(whole-transcriptome padj):", adj = 0,
       cex = 0.51, col = "#404448")
  legend_scores <- c(1, 3, 6)
  legend_x <- c(3.36, 3.82, 4.32)
  symbols(legend_x, rep(legend_y, 3), circles = sqrt(legend_scores), inches = 0.115,
          add = TRUE, bg = "#E4E6E7", fg = "#777B7F")
  text(legend_x + 0.18, rep(legend_y, 3), labels = legend_scores, adj = 0,
       cex = 0.48, col = "#55595D")
  symbols(5.20, legend_y, circles = 1, inches = 0.055, add = TRUE,
          bg = "#FFFFFF", fg = "#202326", lwd = 1.4)
  text(5.38, legend_y, "FDR < 0.05", adj = 0, cex = 0.49)
  symbols(6.28, legend_y, circles = 1, inches = 0.055, add = TRUE,
          bg = "#FFFFFF", fg = "#AEB3B7", lwd = 0.7)
  text(6.46, legend_y, "tested nonsignificant", adj = 0, cex = 0.49)
  text(8.18, legend_y, "x", font = 2, cex = 0.55, col = "#72777B")
  text(8.35, legend_y, "unavailable", adj = 0, cex = 0.49)
  text(5.25, 0.28, "R, S, and GxI are categorical contrasts; not time points.",
       cex = 0.55, col = "#4D5155")
}

node_mass <- function(node_ids, edges) {
  incoming <- tapply(edges$weight, edges$target_id, sum)
  outgoing <- tapply(edges$weight, edges$source_id, sum)
  get_mass <- function(values) {
    mass <- values[node_ids]
    mass[is.na(mass)] <- 0
    unname(mass)
  }
  pmax(get_mass(incoming), get_mass(outgoing))
}

calculate_global_mass_scale <- function(nodes, edges, layers, bottom, top, gaps,
                                        zero_mass_height = 0.024) {
  capacities <- vapply(layers, function(layer) {
    layer_nodes <- nodes[nodes$node_layer == layer, ]
    masses <- node_mass(layer_nodes$node_id, edges)
    mass_height <- top - bottom - gaps[[layer]] * (nrow(layer_nodes) - 1) -
      sum(masses == 0) * zero_mass_height
    if (mass_height <= 0 || sum(masses) <= 0) {
      stop("Layer geometry has no space for positive-mass nodes")
    }
    mass_height / sum(masses)
  }, numeric(1))
  min(capacities)
}

calculate_node_spans <- function(layer_nodes, edges, bottom, top, gap,
                                 mass_scale, zero_mass_height = 0.024) {
  layer_nodes <- layer_nodes[order(layer_nodes$display_order), ]
  masses <- node_mass(layer_nodes$node_id, edges)
  zero_mass <- masses == 0
  used_height <- sum(masses) * mass_scale + sum(zero_mass) * zero_mass_height +
    gap * (nrow(layer_nodes) - 1)
  if (!is.finite(mass_scale) || mass_scale <= 0 ||
      used_height > top - bottom + 1e-12) {
    stop("Layer geometry has no space for positive-mass nodes")
  }
  ymin <- ymax <- numeric(nrow(layer_nodes))
  cursor <- (top + bottom + used_height) / 2
  for (i in seq_len(nrow(layer_nodes))) {
    ymax[i] <- cursor
    node_height <- if (zero_mass[i]) zero_mass_height else masses[i] * mass_scale
    ymin[i] <- cursor - node_height
    cursor <- ymin[i] - gap
  }
  data.frame(
    node_id = layer_nodes$node_id,
    node_label = layer_nodes$node_label,
    node_layer = layer_nodes$node_layer,
    display_order = layer_nodes$display_order,
    color_group = layer_nodes$color_group,
    mass = masses,
    ymin = ymin,
    ymax = ymax,
    stringsAsFactors = FALSE
  )
}

span_centers <- function(spans) {
  setNames((spans$ymin + spans$ymax) / 2, spans$node_id)
}

allocate_edge_spans <- function(pair_edges, source_spans, target_spans) {
  pair_edges$source_ymin <- pair_edges$source_ymax <- NA_real_
  pair_edges$target_ymin <- pair_edges$target_ymax <- NA_real_
  source_centers <- span_centers(source_spans)
  target_centers <- span_centers(target_spans)

  for (source_id in source_spans$node_id) {
    indices <- which(pair_edges$source_id == source_id)
    if (length(indices) == 0) next
    indices <- indices[order(target_centers[pair_edges$target_id[indices]])]
    node <- source_spans[source_spans$node_id == source_id, ]
    scale <- (node$ymax - node$ymin) / sum(pair_edges$weight[indices])
    cursor <- node$ymin
    for (index in indices) {
      pair_edges$source_ymin[index] <- cursor
      pair_edges$source_ymax[index] <- cursor + pair_edges$weight[index] * scale
      cursor <- pair_edges$source_ymax[index]
    }
  }

  for (target_id in target_spans$node_id) {
    indices <- which(pair_edges$target_id == target_id)
    if (length(indices) == 0) next
    indices <- indices[order(source_centers[pair_edges$source_id[indices]])]
    node <- target_spans[target_spans$node_id == target_id, ]
    scale <- (node$ymax - node$ymin) / sum(pair_edges$weight[indices])
    cursor <- node$ymin
    for (index in indices) {
      pair_edges$target_ymin[index] <- cursor
      pair_edges$target_ymax[index] <- cursor + pair_edges$weight[index] * scale
      cursor <- pair_edges$target_ymax[index]
    }
  }
  pair_edges
}

draw_ribbon <- function(x0, x1, source_span, target_span, fill, border = NA, lty = 1) {
  t <- seq(0, 1, length.out = 40)
  smooth <- 3 * t^2 - 2 * t^3
  top <- cbind(
    x0 + (x1 - x0) * t,
    source_span[2] + (target_span[2] - source_span[2]) * smooth
  )
  bottom <- cbind(
    x0 + (x1 - x0) * rev(t),
    source_span[1] + (target_span[1] - source_span[1]) * rev(smooth)
  )
  polygon(
    c(top[, 1], bottom[, 1]), c(top[, 2], bottom[, 2]),
    col = fill, border = border, lty = lty
  )
}

edge_color <- function(edge) {
  if (edge$edge_class == "functional_prior") {
    return(adjustcolor(neutral_ribbon, alpha.f = 0.52))
  }
  if (edge$edge_class == "candidate_membership") {
    return(adjustcolor(lookup_color(edge$color_group, module_colors), alpha.f = 0.58))
  }
  alpha <- if (edge$statistical_state == "significant") {
    0.72
  } else if (edge$statistical_state == "not_significant") {
    0.25
  } else {
    0.18
  }
  adjustcolor(lookup_color(edge$color_group, response_colors), alpha.f = alpha)
}

draw_layer_nodes <- function(spans, x, width, palette, fill_alpha = 0.20,
                             fallback = "#8B8F93", cex = 0.50,
                             show_mass = FALSE) {
  for (i in seq_len(nrow(spans))) {
    row <- spans[i, ]
    border <- lookup_color(row$color_group, palette, fallback)
    fill <- adjustcolor(border, alpha.f = fill_alpha)
    rect(
      x - width / 2, row$ymin, x + width / 2, row$ymax,
      col = fill, border = border, lwd = 0.85
    )
    label <- if (show_mass) {
      paste0(row$node_label, "\n", "n = ", format(round(row$mass, 2), trim = TRUE))
    } else {
      row$node_label
    }
    node_height <- row$ymax - row$ymin
    if (show_mass && node_height < 0.035) {
      text(
        x, row$ymax + 0.008,
        labels = paste0(row$node_label, " (n = ", format(round(row$mass, 2), trim = TRUE), ")"),
        cex = 0.34, col = "#202326", adj = c(0.5, 0)
      )
    } else {
      text(x, (row$ymin + row$ymax) / 2, labels = label,
           cex = cex, col = "#202326")
    }
  }
}

draw_alluvial_panel <- function(nodes, edges) {
  par(mar = c(0.15, 0.2, 0.15, 0.15), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  layers <- c("effector_class", "host_module", "candidate", "response")
  layer_x <- c(effector_class = 0.085, host_module = 0.315, candidate = 0.635, response = 0.915)
  layer_width <- c(effector_class = 0.145, host_module = 0.160, candidate = 0.180, response = 0.155)
  layer_gap <- c(effector_class = 0.016, host_module = 0.016, candidate = 0.008, response = 0.020)
  layer_titles <- c(
    effector_class = "Predicted PWN\nfunctional class",
    host_module = "Pine defense\nmodule",
    candidate = "Prioritized host\northogroup",
    response = "Observed GxI\ndirection"
  )

  spans <- setNames(vector("list", length(layers)), layers)
  mass_scale <- calculate_global_mass_scale(
    nodes, edges, layers,
    bottom = 0.205, top = 0.865, gaps = layer_gap
  )
  for (layer in layers) {
    spans[[layer]] <- calculate_node_spans(
      nodes[nodes$node_layer == layer, ], edges,
      bottom = 0.205, top = 0.865, gap = layer_gap[[layer]],
      mass_scale = mass_scale
    )
  }

  text(0.01, 0.982, "Effector-informed functional-prior and host-response map",
       adj = c(0, 0.5), font = 2, cex = 1.02)
  text(layer_x, 0.925, labels = layer_titles, font = 2, cex = 0.58)

  pairs <- list(
    c("effector_class", "host_module"),
    c("host_module", "candidate"),
    c("candidate", "response")
  )
  allocated_pairs <- list()
  for (i in seq_along(pairs)) {
    source_layer <- pairs[[i]][1]
    target_layer <- pairs[[i]][2]
    pair_edges <- edges[
      edges$source_layer == source_layer & edges$target_layer == target_layer,
    ]
    allocated_pairs[[i]] <- allocate_edge_spans(
      pair_edges, spans[[source_layer]], spans[[target_layer]]
    )
  }
  if (sum(vapply(allocated_pairs, nrow, integer(1))) != nrow(edges)) {
    stop("Validated alluvial edge was not assigned to a rendered transition")
  }

  for (i in seq_along(pairs)) {
    source_layer <- pairs[[i]][1]
    target_layer <- pairs[[i]][2]
    pair_edges <- allocated_pairs[[i]]
    draw_order <- order(pair_edges$weight, decreasing = TRUE)
    x0 <- layer_x[[source_layer]] + layer_width[[source_layer]] / 2
    x1 <- layer_x[[target_layer]] - layer_width[[target_layer]] / 2
    for (index in draw_order) {
      edge <- pair_edges[index, ]
      draw_ribbon(
        x0, x1,
        c(edge$source_ymin, edge$source_ymax),
        c(edge$target_ymin, edge$target_ymax),
        fill = edge_color(edge),
        border = if (edge$statistical_state == "unavailable") "#7B8084" else NA,
        lty = if (edge$statistical_state == "unavailable") 2 else 1
      )
    }
  }

  effector_palette <- c(functional_prior = "#7A8791")
  candidate_palette <- c(candidate = "#646A70")
  draw_layer_nodes(
    spans$effector_class, layer_x[["effector_class"]], layer_width[["effector_class"]],
    effector_palette, fill_alpha = 0.16, cex = 0.46, show_mass = TRUE
  )
  draw_layer_nodes(
    spans$host_module, layer_x[["host_module"]], layer_width[["host_module"]],
    module_colors, fill_alpha = 0.19, cex = 0.41, show_mass = TRUE
  )
  draw_layer_nodes(
    spans$candidate, layer_x[["candidate"]], layer_width[["candidate"]],
    candidate_palette, fill_alpha = 0.08, cex = 0.48
  )
  draw_layer_nodes(
    spans$response, layer_x[["response"]], layer_width[["response"]],
    response_colors, fill_alpha = 0.19, cex = 0.49, show_mass = TRUE
  )

  text(0.02, 0.151, "Ribbon width = normalized number of candidate orthogroups, not interaction strength.",
       adj = 0, cex = 0.49, col = "#414549")
  segments(0.035, 0.108, 0.105, 0.108, col = adjustcolor(neutral_ribbon, alpha.f = 0.65), lwd = 7)
  text(0.115, 0.108, "functional prior: hypothesis; not direct interaction", adj = 0,
       cex = 0.46, col = "#414549")
  segments(0.405, 0.108, 0.475, 0.108, col = adjustcolor(module_colors[["hydraulic_xylem"]], alpha.f = 0.75), lwd = 7)
  text(0.485, 0.108, "module membership", adj = 0, cex = 0.46, col = "#414549")
  segments(0.655, 0.108, 0.700, 0.108, col = response_colors[["r_response_greater"]], lwd = 7)
  segments(0.700, 0.108, 0.745, 0.108, col = response_colors[["s_response_greater"]], lwd = 7)
  segments(0.745, 0.108, 0.790, 0.108, col = response_colors[["interaction_unresolved"]], lwd = 7)
  text(0.800, 0.108, "GxI summary", adj = 0, cex = 0.46, col = "#414549")
  text(0.02, 0.058, "Functional priors and module annotations; not direct effector-target interactions.",
       adj = 0, cex = 0.54, font = 2, col = "#303438")
  text(0.98, 0.058, "Response paths encode whole-transcriptome FDR status: deep = significant; pale = nonsignificant; dashed = unavailable.",
       adj = 1, cex = 0.46, col = "#55595D")
}

draw_figure <- function() {
  draw_alluvial_panel(nodes, edges)
}

dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)

pdf(pdf_path, width = 11.5, height = 7.4, pointsize = 14, useDingbats = FALSE)
draw_figure()
invisible(dev.off())

png(
  png_path, width = 11.5, height = 7.4, units = "in", res = 320,
  pointsize = 14, type = if (capabilities("cairo")) "cairo" else "windows"
)
draw_figure()
invisible(dev.off())

message("Functional-prior response map written: ", pdf_path, " and ", png_path)
