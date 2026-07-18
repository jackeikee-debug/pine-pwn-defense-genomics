#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(circlize))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(6, 7)) {
  stop("Usage: plot_evidence_sector_circos.R sectors.tsv tracks.tsv links.tsv gene_sample_mosaic.tsv output.pdf output.png [main|full]")
}

sector_path <- args[[1]]
track_path <- args[[2]]
link_path <- args[[3]]
mosaic_path <- args[[4]]
pdf_path <- args[[5]]
png_path <- args[[6]]
plot_mode <- if (length(args) == 7) args[[7]] else "full"
if (!plot_mode %in% c("main", "full")) stop("Plot mode must be main or full")

read_table <- function(path) {
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "")
  )
}

sectors <- read_table(sector_path)
tracks <- read_table(track_path)
links <- read_table(link_path)
mosaic <- read_table(mosaic_path)

sectors <- sectors[order(sectors$sector_order), , drop = FALSE]
sector_ids <- sectors$sector_id
sector_counts <- pmax(as.numeric(sectors$item_count), 1)
names(sector_counts) <- sector_ids

group_colors <- c(
  "East Asia" = "#D95F59",
  "North America" = "#2A7F9E",
  "P. massoniana expression" = "#2F8F6B",
  "East Asia expression" = "#D49A27",
  "North America expression" = "#7C5AA6",
  "Pine wood nematode" = "#A83B3B"
)
sector_colors <- unname(group_colors[sectors$sector_group])
sector_colors[is.na(sector_colors)] <- "#68737D"
names(sector_colors) <- sector_ids

module_colors <- c(
  "hydraulic_xylem" = "#168A91",
  "immune_signaling" = "#7756A8",
  "phenylpropanoid_lignin" = "#D8842F",
  "ros_detoxification" = "#C84C4C",
  "wound_periderm" = "#5C8E45"
)

alpha_color <- function(color, alpha) {
  grDevices::adjustcolor(color, alpha.f = alpha)
}

signed_color <- function(value, limit) {
  if (!is.finite(value)) return("#ECEDEF")
  limit <- max(limit, 1e-9)
  x <- max(-limit, min(limit, value))
  if (x < 0) {
    grDevices::colorRampPalette(c("#2455A4", "#F7F7F7"))(101)[round((x + limit) / limit * 100) + 1]
  } else {
    grDevices::colorRampPalette(c("#F7F7F7", "#C7333D"))(101)[round(x / limit * 100) + 1]
  }
}

sequential_color <- function(value) {
  if (!is.finite(value)) return("#ECEDEF")
  palette <- grDevices::colorRampPalette(c("#EFF7F1", "#72B7A1", "#165F54"))(101)
  palette[round(max(0, min(1, value)) * 100) + 1]
}

sector_track <- function(sector_id, track_id) {
  tracks[tracks$sector_id == sector_id & tracks$track_id == track_id, , drop = FALSE]
}

mosaic_sector <- function(sector_id, record_type) {
  mosaic[mosaic$sector_id == sector_id & mosaic$record_type == record_type, , drop = FALSE]
}

item_mid <- function(sector_id, item_id) {
  hit <- tracks[tracks$sector_id == sector_id & tracks$item_id == item_id, , drop = FALSE]
  if (nrow(hit) == 0) return(NA_real_)
  as.numeric(hit$item_order[[1]]) - 0.5
}

gap_after <- rep(1.4, nrow(sectors))
named_gaps <- c(
  ptab_proteome = 5.5,
  ptae_proteome = 5.5,
  pmas_gxi = 5.5,
  pthun_expr = 4.0,
  prigtae_expr = 5.5,
  pwn_secretome = 7.0
)
matched_gaps <- match(sectors$sector_id, names(named_gaps))
gap_after[!is.na(matched_gaps)] <- unname(named_gaps[matched_gaps[!is.na(matched_gaps)]])

heat_values <- tracks$value[tracks$track_id %in% c("copy_centered", "expression_log2fc")]
heat_limit <- max(abs(heat_values[is.finite(heat_values)]), 1)
sample_z_limit <- 2.5
de_values <- as.numeric(mosaic$log2fc[mosaic$record_type == "contrast_de"])
de_limit <- max(stats::quantile(abs(de_values[is.finite(de_values)]), 0.95, names = FALSE), 1)
raw_p_significance <- -log10(pmax(as.numeric(tracks$pvalue), 1e-300))
raw_p_limit <- max(raw_p_significance[is.finite(raw_p_significance)], 1)

draw_plot <- function() {
  par(mar = c(0.3, 0.3, 1.8, 0.3), xpd = NA, family = "sans")

  circos.clear()
  circos.par(
    start.degree = 89,
    gap.after = gap_after,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.0035, 0.0035),
    points.overflow.warning = FALSE,
    canvas.xlim = c(-1.15, 1.15),
    canvas.ylim = c(-1.13, 1.13)
  )
  circos.initialize(
    factors = sector_ids,
    xlim = cbind(rep(0, length(sector_ids)), sector_counts)
  )

  # A: evidence-unit sectors.
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.060, bg.border = NA,
    panel.fun = function(x, y) {
      sector_id <- CELL_META$sector.index
      color <- sector_colors[[sector_id]]
      circos.rect(
        CELL_META$xlim[[1]], 0, CELL_META$xlim[[2]], 1,
        col = alpha_color(color, 0.92), border = "white", lwd = 0.8
      )
      sector_row <- sectors[sectors$sector_id == sector_id, , drop = FALSE]
      circos.text(
        CELL_META$xcenter, 1.55, sector_row$sector_label,
        facing = "bending.outside", niceFacing = TRUE,
        cex = 0.47, font = 2, col = "#20262B"
      )
      if (sector_id == sector_ids[[1]]) {
        circos.text(CELL_META$xlim[[1]] + 0.12, 0.5, "A", cex = 0.56, font = 2, col = "white")
      }
    }
  )

  # G: gene x sample expression mosaic. Each narrow column is an exact mapped gene.
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.145, bg.col = "#F2F4F4", bg.border = "white",
    panel.fun = function(x, y) {
      sector_id <- CELL_META$sector.index
      dat <- mosaic_sector(sector_id, "sample_expression")
      if (nrow(dat) == 0) {
        context <- sector_track(sector_id, if (sector_id == "pwn_secretome") "effector_supported_fraction" else "copy_centered")
        for (i in seq_len(nrow(context))) {
          left <- as.numeric(context$item_order[[i]]) - 1
          right <- as.numeric(context$item_order[[i]])
          value <- as.numeric(context$value[[i]])
          color <- if (sector_id == "pwn_secretome") sequential_color(value) else signed_color(value, heat_limit)
          circos.rect(left, 0, right, 1, col = alpha_color(color, 0.34), border = "white", lwd = 0.18)
        }
      } else {
        item_total <- max(as.numeric(sectors$item_count[sectors$sector_id == sector_id]), 1)
        mapped_ogs <- unique(as.numeric(dat$og_order))
        for (j in seq_len(item_total)) {
          if (!(j %in% mapped_ogs)) {
            circos.rect(j - 1, 0, j, 1, col = "#DDE2E4", border = "white", lwd = 0.22)
            circos.lines(c(j - 0.82, j - 0.18), c(0.12, 0.88), col = "#A9B2B6", lwd = 0.28)
          }
        }
        gene_count <- as.numeric(dat$gene_count)
        gene_order <- as.numeric(dat$gene_order)
        sample_count <- as.numeric(dat$sample_count)
        sample_order <- as.numeric(dat$sample_order)
        og_left <- as.numeric(dat$og_order) - 1
        left <- og_left + (gene_order - 1) / gene_count
        right <- og_left + gene_order / gene_count
        bottom <- (sample_order - 1) / sample_count
        top <- sample_order / sample_count
        colors <- vapply(as.numeric(dat$expression_zscore), signed_color, character(1), limit = sample_z_limit)
        circos.rect(left, bottom, right, top, col = colors,
                    border = alpha_color("white", 0.42), lwd = 0.06)
        de_dat <- mosaic_sector(sector_id, "contrast_de")
        representatives <- unique(de_dat[de_dat$is_representative == "yes", c("og_order", "gene_order", "gene_count"), drop = FALSE])
        for (i in seq_len(nrow(representatives))) {
          rep <- representatives[i, ]
          left <- as.numeric(rep$og_order) - 1 + (as.numeric(rep$gene_order) - 1) / as.numeric(rep$gene_count)
          right <- as.numeric(rep$og_order) - 1 + as.numeric(rep$gene_order) / as.numeric(rep$gene_count)
          circos.rect(left, 0, right, 1, col = NA, border = "#111111", lwd = 0.42)
        }
        for (boundary in seq_len(item_total - 1)) {
          circos.lines(c(boundary, boundary), c(0, 1), col = "#697277", lwd = 0.30)
        }
      }
      if (sector_id == "pmas_r") {
        circos.text(CELL_META$xlim[[1]] + 0.16, 0.5, "G", cex = 0.55, font = 2)
      }
    }
  )

  # H: gene x contrast log2FC mosaic; dots mark whole-transcriptome FDR < 0.05.
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.125, bg.col = "#F5F6F6", bg.border = "white",
    panel.fun = function(x, y) {
      sector_id <- CELL_META$sector.index
      dat <- mosaic_sector(sector_id, "contrast_de")
      if (nrow(dat) == 0) {
        context_id <- if (sector_id == "pwn_secretome") "effector_candidate_count" else "copy_count"
        context <- sector_track(sector_id, context_id)
        max_value <- max(as.numeric(context$value), 1, na.rm = TRUE)
        for (i in seq_len(nrow(context))) {
          left <- as.numeric(context$item_order[[i]]) - 1
          right <- as.numeric(context$item_order[[i]])
          level <- as.numeric(context$value[[i]]) / max_value
          circos.rect(left, 0, right, 1, col = alpha_color(sector_colors[[sector_id]], 0.08 + 0.25 * level), border = "white", lwd = 0.18)
        }
      } else {
        item_total <- max(as.numeric(sectors$item_count[sectors$sector_id == sector_id]), 1)
        mapped_ogs <- unique(as.numeric(dat$og_order))
        for (j in seq_len(item_total)) {
          if (!(j %in% mapped_ogs)) {
            circos.rect(j - 1, 0, j, 1, col = "#DDE2E4", border = "white", lwd = 0.22)
            circos.lines(c(j - 0.82, j - 0.18), c(0.12, 0.88), col = "#A9B2B6", lwd = 0.28)
          }
        }
        gene_count <- as.numeric(dat$gene_count)
        gene_order <- as.numeric(dat$gene_order)
        contrast_count <- as.numeric(dat$contrast_count)
        contrast_order <- as.numeric(dat$contrast_order)
        og_left <- as.numeric(dat$og_order) - 1
        left <- og_left + (gene_order - 1) / gene_count
        right <- og_left + gene_order / gene_count
        bottom <- (contrast_order - 1) / contrast_count
        top <- contrast_order / contrast_count
        values <- as.numeric(dat$log2fc)
        colors <- vapply(values, signed_color, character(1), limit = de_limit)
        colors[!is.finite(values)] <- "#E4E7E8"
        circos.rect(left, bottom, right, top, col = colors,
                    border = alpha_color("white", 0.48), lwd = 0.07)
        representative <- dat$is_representative == "yes"
        if (any(representative)) {
          circos.rect(left[representative], bottom[representative], right[representative], top[representative],
                      col = NA, border = "#111111", lwd = 0.48)
        }
        significant <- dat$significant == "yes"
        if (any(significant)) {
          circos.points((left[significant] + right[significant]) / 2,
                        (bottom[significant] + top[significant]) / 2,
                        pch = 16, cex = 0.11, col = "#F4C542")
        }
        for (boundary in seq_len(item_total - 1)) {
          circos.lines(c(boundary, boundary), c(0, 1), col = "#697277", lwd = 0.30)
        }
      }
      if (sector_id == "pmas_r") {
        circos.text(CELL_META$xlim[[1]] + 0.16, 0.72, "H", cex = 0.55, font = 2)
      }
    }
  )

  # B: signed heatmap for copy context or expression response.
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.075, bg.col = "#F7F8F8", bg.border = "white",
    panel.fun = function(x, y) {
      sector_id <- CELL_META$sector.index
      sector_type <- sectors$evidence_type[sectors$sector_id == sector_id]
      track_id <- if (sector_type == "proteome_copy_context") {
        "copy_centered"
      } else if (sector_type == "secretome_functional_context") {
        "effector_supported_fraction"
      } else {
        "expression_log2fc"
      }
      dat <- sector_track(sector_id, track_id)
      for (i in seq_len(nrow(dat))) {
        left <- as.numeric(dat$item_order[[i]]) - 1
        right <- as.numeric(dat$item_order[[i]])
        color <- if (track_id == "effector_supported_fraction") {
          sequential_color(as.numeric(dat$value[[i]]))
        } else {
          signed_color(as.numeric(dat$value[[i]]), heat_limit)
        }
        border <- if (dat$highlight_tier[[i]] == "lead") "#20262B" else "white"
        circos.rect(left, 0, right, 1, col = color, border = border, lwd = 0.45)
      }
      if (sector_id == sector_ids[[1]]) {
        circos.text(CELL_META$xlim[[1]] + 0.12, 0.5, "B", cex = 0.52, font = 2)
      }
    }
  )

  # C: radial bars for copy count, response magnitude, or effector count.
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.075, bg.col = "#F5F6F6", bg.border = "white",
    panel.fun = function(x, y) {
      sector_id <- CELL_META$sector.index
      sector_type <- sectors$evidence_type[sectors$sector_id == sector_id]
      track_id <- if (sector_type == "proteome_copy_context") {
        "copy_count"
      } else if (sector_type == "secretome_functional_context") {
        "effector_candidate_count"
      } else {
        "expression_abs_log2fc"
      }
      dat <- sector_track(sector_id, track_id)
      for (i in seq_len(nrow(dat))) {
        left <- as.numeric(dat$item_order[[i]]) - 0.84
        right <- as.numeric(dat$item_order[[i]]) - 0.16
        height <- as.numeric(dat$plot_value[[i]])
        color <- if (sector_type == "controlled_expression" || sector_type == "cross_species_expression") {
          signed_color(as.numeric(sector_track(sector_id, "expression_log2fc")$value[[i]]), heat_limit)
        } else {
          alpha_color(sector_colors[[sector_id]], 0.88)
        }
        if (is.finite(height)) {
          circos.rect(left, 0, right, height, col = color, border = NA)
        }
      }
      circos.lines(CELL_META$xlim, c(0, 0), col = "#AEB5B9", lwd = 0.35)
      if (sector_id == sector_ids[[1]]) {
        circos.text(CELL_META$xlim[[1]] + 0.12, 0.82, "C", cex = 0.52, font = 2)
      }
    }
  )

  # D: FDR points for expression and evidence-support points elsewhere.
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.060, bg.col = "#FAFAFA", bg.border = "white",
    panel.fun = function(x, y) {
      sector_id <- CELL_META$sector.index
      sector_type <- sectors$evidence_type[sectors$sector_id == sector_id]
      if (sector_type %in% c("controlled_expression", "cross_species_expression")) {
        dat <- sector_track(sector_id, "expression_fdr")
        for (i in seq_len(nrow(dat))) {
          point_y <- as.numeric(dat$plot_value[[i]])
          if (!is.finite(point_y)) next
          significant <- is.finite(as.numeric(dat$fdr[[i]])) && as.numeric(dat$fdr[[i]]) < 0.05
          circos.points(
            as.numeric(dat$item_order[[i]]) - 0.5, point_y,
            pch = if (significant) 16 else 1,
            cex = if (dat$highlight_tier[[i]] == "lead") 0.48 else 0.34,
            col = if (significant) "#B7222E" else "#8B9499"
          )
          raw_p <- as.numeric(dat$pvalue[[i]])
          if (is.finite(raw_p) && raw_p > 0) {
            raw_p_y <- min(1, -log10(max(raw_p, 1e-300)) / raw_p_limit)
            circos.points(
              as.numeric(dat$item_order[[i]]) - 0.5, raw_p_y,
              pch = 5, cex = 0.30, col = "#315D8A"
            )
          }
        }
      } else if (sector_type == "proteome_copy_context") {
        support_ids <- c("network_support", "structure_support", "sequence_support")
        for (support_index in seq_along(support_ids)) {
          dat <- sector_track(sector_id, support_ids[[support_index]])
          dat <- dat[is.finite(dat$value) & dat$value > 0, , drop = FALSE]
          if (nrow(dat) > 0) {
            circos.points(
              as.numeric(dat$item_order) - 0.5,
              rep(0.18 + 0.30 * (support_index - 1), nrow(dat)),
              pch = c(16, 17, 15)[[support_index]], cex = 0.34,
              col = c("#315D8A", "#D17A22", "#4A8C5A")[[support_index]]
            )
          }
        }
      } else {
        dat <- sector_track(sector_id, "functional_support")
        dat <- dat[is.finite(dat$value) & dat$value > 0, , drop = FALSE]
        if (nrow(dat) > 0) {
          circos.points(as.numeric(dat$item_order) - 0.5, rep(0.55, nrow(dat)), pch = 18, cex = 0.56, col = "#1B6E5E")
        }
      }
      if (sector_id == sector_ids[[1]]) {
        circos.text(CELL_META$xlim[[1]] + 0.12, 0.5, "D", cex = 0.52, font = 2)
      }
    }
  )

  # E: within-sector profile line; scaling is already frozen in the input table.
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.070, bg.col = "#F4F5F5", bg.border = "white",
    panel.fun = function(x, y) {
      sector_id <- CELL_META$sector.index
      sector_type <- sectors$evidence_type[sectors$sector_id == sector_id]
      track_id <- if (sector_type == "proteome_copy_context") {
        "copy_centered"
      } else if (sector_type == "secretome_functional_context") {
        "effector_supported_fraction"
      } else {
        "expression_profile"
      }
      dat <- sector_track(sector_id, track_id)
      dat <- dat[order(dat$item_order), , drop = FALSE]
      valid <- is.finite(as.numeric(dat$plot_value))
      if (sum(valid) >= 2) {
        circos.lines(
          as.numeric(dat$item_order[valid]) - 0.5,
          as.numeric(dat$plot_value[valid]),
          col = alpha_color(sector_colors[[sector_id]], 0.92), lwd = 1.15
        )
      }
      if (sum(valid) > 0) {
        circos.points(
          as.numeric(dat$item_order[valid]) - 0.5,
          as.numeric(dat$plot_value[valid]),
          pch = 16, cex = 0.20, col = sector_colors[[sector_id]]
        )
      }
      circos.lines(CELL_META$xlim, c(0.5, 0.5), col = "#BFC5C8", lwd = 0.35, lty = 3)
      if (sector_id == sector_ids[[1]]) {
        circos.text(CELL_META$xlim[[1]] + 0.12, 0.78, "E", cex = 0.52, font = 2)
      }
    }
  )

  # F: mechanism-leading highlights and compact labels at anchor sectors.
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.050, bg.col = "#FFFFFF", bg.border = "white",
    panel.fun = function(x, y) {
      sector_id <- CELL_META$sector.index
      dat <- sector_track(sector_id, "priority_highlight")
      for (i in seq_len(nrow(dat))) {
        if (!is.finite(as.numeric(dat$value[[i]])) || as.numeric(dat$value[[i]]) <= 0) next
        midpoint <- as.numeric(dat$item_order[[i]]) - 0.5
        circos.rect(midpoint - 0.32, 0.25, midpoint + 0.32, 0.75,
                    col = "#F2C84B", border = "#5A4B12", lwd = 0.45)
      }
      label_sector <- sector_id %in% c("pmas_proteome", "pwn_secretome")
      if (label_sector) {
        label_dat <- tracks[tracks$sector_id == sector_id & tracks$track_id %in% c("copy_count", "effector_candidate_count"), , drop = FALSE]
        for (i in seq_len(nrow(label_dat))) {
          is_lead <- label_dat$highlight_tier[[i]] == "lead"
          if (!is_lead && sector_id != "pwn_secretome") next
          circos.text(
            as.numeric(label_dat$item_order[[i]]) - 0.5, 0.05,
            label_dat$item_label[[i]],
            facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
            cex = if (sector_id == "pwn_secretome") 0.30 else 0.27,
            col = "#252A2D"
          )
        }
      }
      if (sector_id == sector_ids[[1]]) {
        circos.text(CELL_META$xlim[[1]] + 0.12, 0.5, "F", cex = 0.52, font = 2)
      }
    }
  )

  # Pale orthogroup bridges first, then functional hypotheses on top.
  evidence_links <- links[links$link_type == "orthogroup_evidence", , drop = FALSE]
  if (nrow(evidence_links) > 0) {
    for (i in seq_len(nrow(evidence_links))) {
      row <- evidence_links[i, ]
      source_mid <- item_mid(row$source_sector, row$source_item)
      target_mid <- item_mid(row$target_sector, row$target_item)
      if (!is.finite(source_mid) || !is.finite(target_mid)) next
      width <- 0.035 + 0.045 * max(0, min(1, as.numeric(row$weight)))
      color <- if (row$highlight_tier == "lead") "#5078A5" else "#8A949C"
      circos.link(
        row$source_sector, c(source_mid - width, source_mid + width),
        row$target_sector, c(target_mid - width, target_mid + width),
        col = alpha_color(color, if (row$highlight_tier == "lead") 0.16 else 0.07),
        border = NA
      )
    }
  }

  hypothesis_links <- links[links$link_type == "functional_hypothesis", , drop = FALSE]
  if (nrow(hypothesis_links) > 0) {
    hypothesis_links <- hypothesis_links[order(hypothesis_links$highlight_tier != "lead"), , drop = FALSE]
    for (i in seq_len(nrow(hypothesis_links))) {
      row <- hypothesis_links[i, ]
      source_mid <- item_mid(row$source_sector, row$source_item)
      target_mid <- item_mid(row$target_sector, row$target_item)
      if (!is.finite(source_mid) || !is.finite(target_mid)) next
      weight <- max(0, min(3, as.numeric(row$weight)))
      width <- 0.045 + 0.035 * weight
      color <- module_colors[[as.character(row$color_group)]]
      if (is.null(color) || is.na(color)) color <- "#A64B4B"
      circos.link(
        row$source_sector, c(source_mid - width, source_mid + width),
        row$target_sector, c(target_mid - width, target_mid + width),
        col = alpha_color(color, if (row$highlight_tier == "lead") 0.42 else 0.20),
        border = alpha_color(color, 0.20), lwd = 0.35
      )
    }
  }

  figure_title <- if (plot_mode == "main") {
    "Integrated cross-species evidence landscape for prioritized pine orthogroups"
  } else {
    "Full cross-species evidence landscape and mapping provenance"
  }
  title(
    main = figure_title,
    cex.main = 0.88, font.main = 2, col.main = "#1F262B", line = 0.45
  )

  circos.clear()
}

draw_gradient <- function(y, limit, label) {
  palette <- vapply(seq(-limit, limit, length.out = 80), signed_color, character(1), limit = limit)
  x <- seq(0.08, 0.88, length.out = length(palette) + 1)
  for (i in seq_along(palette)) {
    rect(x[[i]], y, x[[i + 1]], y + 0.018, col = palette[[i]], border = NA)
  }
  text(0.08, y - 0.012, sprintf("-%s", format(limit, digits = 2)), adj = c(0, 1), cex = 0.57, col = "#596267")
  text(0.48, y - 0.012, "0", adj = c(0.5, 1), cex = 0.57, col = "#596267")
  text(0.88, y - 0.012, sprintf("+%s", format(limit, digits = 2)), adj = c(1, 1), cex = 0.57, col = "#596267")
  text(0.08, y + 0.031, label, adj = c(0, 0), cex = 0.68, font = 2, col = "#263036")
}

draw_side_legend <- function() {
  par(mar = c(1.0, 0.2, 2.2, 1.0), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.04, 0.965, "TRACK GUIDE", adj = c(0, 1), cex = 0.88, font = 2, col = "#1F292E")
  text(0.04, 0.925, "Each angular unit is one candidate OG; each narrow column in G/H is one mapped gene.",
       adj = c(0, 1), cex = 0.62, col = "#4D585E")

  ring_text <- c(
    "A  evidence sectors",
    "G  gene x sample expression",
    "H  gene x contrast log2FC + FDR",
    "B  OG copy/response heatmap",
    "C  copy/response magnitude",
    "D  OG-level significance/support",
    "E  within-sector profile",
    "F  mechanism-leading candidates"
  )
  y <- 0.855
  for (label in ring_text) {
    text(0.06, y, label, adj = c(0, 0.5), cex = 0.64, col = "#313A3F")
    y <- y - 0.032
  }

  draw_gradient(0.555, sample_z_limit, "G: within-gene expression z-score")
  draw_gradient(0.465, de_limit, "H: gene-level log2FC")
  points(0.09, 0.408, pch = 16, cex = 0.72, col = "#F4C542")
  text(0.14, 0.408, "whole-transcriptome FDR < 0.05", adj = c(0, 0.5), cex = 0.62, col = "#3F484D")
  rect(0.075, 0.368, 0.115, 0.392, col = NA, border = "#111111", lwd = 1.0)
  text(0.14, 0.380, "representative gene", adj = c(0, 0.5), cex = 0.60, col = "#3F484D")
  rect(0.075, 0.328, 0.115, 0.352, col = "#DDE2E4", border = "#A9B2B6", lwd = 0.7)
  segments(0.078, 0.330, 0.112, 0.350, col = "#929DA2", lwd = 0.7)
  text(0.14, 0.340, "unmapped candidate OG", adj = c(0, 0.5), cex = 0.60, col = "#3F484D")

  text(0.04, 0.292, "SAMPLE BAND ORDER", adj = c(0, 0.5), cex = 0.72, font = 2, col = "#263036")
  sample_notes <- c(
    "P. massoniana R/S: inoculated then water, n = 6",
    "P. massoniana GxI: RI, RW, SI, SW, n = 12",
    "P. densiflora: PWN, nonpathogenic, water, n = 9",
    "Time series: 0, 2 and 4 weeks, n = 9"
  )
  y <- 0.258
  for (label in sample_notes) {
    text(0.06, y, label, adj = c(0, 0.5), cex = 0.57, col = "#4B565B")
    y <- y - 0.031
  }

  text(0.04, 0.124, "EVIDENCE BOUNDARIES", adj = c(0, 0.5), cex = 0.72, font = 2, col = "#263036")
  text(0.06, 0.092, "P. massoniana: exact candidate members", adj = c(0, 0.5), cex = 0.57, col = "#4B565B")
  text(0.06, 0.064, "P. densiflora: same-species homology", adj = c(0, 0.5), cex = 0.57, col = "#4B565B")
  text(0.06, 0.036, "P. strobus / P. thunbergii / hybrid: proxy mappings", adj = c(0, 0.5), cex = 0.57, col = "#4B565B")

  segments(0.61, 0.380, 0.72, 0.380, col = alpha_color("#C84C4C", 0.62), lwd = 4)
  text(0.74, 0.380, "hypothesis link", adj = c(0, 0.5), cex = 0.53, col = "#4B565B")
  segments(0.61, 0.340, 0.72, 0.340, col = alpha_color("#5078A5", 0.24), lwd = 4)
  text(0.74, 0.340, "same-OG link", adj = c(0, 0.5), cex = 0.53, col = "#4B565B")
}

draw_bottom_legend <- function() {
  par(mar = c(0.2, 0.8, 0.2, 0.8), family = "sans", xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  text(0.01, 0.88, "TRACKS", adj = c(0, 0.5), cex = 0.70, font = 2, col = "#263036")
  text(
    0.085, 0.88,
    "A sectors   G gene x sample z-score   H gene x contrast log2FC + FDR   B-E summary tracks   F priority",
    adj = c(0, 0.5), cex = 0.60, col = "#3F484D"
  )

  sample_palette <- vapply(seq(-sample_z_limit, sample_z_limit, length.out = 50), signed_color,
                           character(1), limit = sample_z_limit)
  de_palette <- vapply(seq(-de_limit, de_limit, length.out = 50), signed_color,
                       character(1), limit = de_limit)
  for (i in seq_along(sample_palette)) {
    left <- 0.01 + (i - 1) * 0.0036
    rect(left, 0.53, left + 0.0037, 0.66, col = sample_palette[[i]], border = NA)
    left2 <- 0.255 + (i - 1) * 0.0036
    rect(left2, 0.53, left2 + 0.0037, 0.66, col = de_palette[[i]], border = NA)
  }
  text(0.01, 0.43, "G within-gene expression z-score", adj = c(0, 0.5), cex = 0.55)
  text(0.255, 0.43, "H gene-level log2FC", adj = c(0, 0.5), cex = 0.55)
  points(0.475, 0.60, pch = 16, cex = 0.72, col = "#F4C542")
  text(0.49, 0.60, "whole-transcriptome FDR < 0.05", adj = c(0, 0.5), cex = 0.55)
  rect(0.475, 0.34, 0.490, 0.47, col = NA, border = "#111111", lwd = 0.9)
  text(0.50, 0.405, "representative gene", adj = c(0, 0.5), cex = 0.55)

  segments(0.65, 0.61, 0.705, 0.61, col = alpha_color("#5078A5", 0.42), lwd = 4)
  text(0.715, 0.61, "same-OG evidence continuity", adj = c(0, 0.5), cex = 0.55)
  segments(0.65, 0.39, 0.705, 0.39, col = alpha_color("#C84C4C", 0.62), lwd = 4)
  text(0.715, 0.39, "functional prior; not interaction", adj = c(0, 0.5), cex = 0.55)

  text(
    0.01, 0.10,
    "Expression colors are within-dataset summaries; P. densiflora uses same-species homology and P. strobus uses consensus projection.",
    adj = c(0, 0.5), cex = 0.53, col = "#4B565B"
  )
}

draw_figure <- function() {
  if (plot_mode == "main") {
    layout(matrix(c(1, 2), ncol = 1), heights = c(5.4, 0.85))
    draw_plot()
    draw_bottom_legend()
  } else {
    layout(matrix(c(1, 2), nrow = 1), widths = c(4.7, 1.45))
    draw_plot()
    draw_side_legend()
  }
}

dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)

png_width <- if (plot_mode == "main") 4200 else 4650
png_height <- if (plot_mode == "main") 3900 else 3750
pdf_width <- if (plot_mode == "main") 14.0 else 15.5
pdf_height <- if (plot_mode == "main") 13.0 else 12.5

png(png_path, width = png_width, height = png_height, res = 300, bg = "white")
dev.control(displaylist = "enable")
draw_figure()
recorded_figure <- recordPlot()
dev.off()

pdf(pdf_path, width = pdf_width, height = pdf_height, useDingbats = FALSE, bg = "white")
replayPlot(recorded_figure)
dev.off()

message(sprintf("Wrote %s and %s", pdf_path, png_path))
