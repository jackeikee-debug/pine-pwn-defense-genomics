#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(paste(
    "Usage: plot_manuscript_evidence_figures.R",
    "<candidate_summary.tsv> <apoplastic_audit.tsv>",
    "<fig4.pdf> <fig4.png> <fig6.pdf> <fig6.png>"
  ))
}

candidate_path <- args[[1]]
audit_path <- args[[2]]
fig4_pdf <- args[[3]]
fig4_png <- args[[4]]
fig6_pdf <- args[[5]]
fig6_png <- args[[6]]

candidates <- read.delim(candidate_path, stringsAsFactors = FALSE, check.names = FALSE)
audit <- read.delim(audit_path, stringsAsFactors = FALSE, check.names = FALSE)

if (nrow(candidates) != 12) {
  stop(sprintf("Expected 12 manuscript candidates, found %d", nrow(candidates)))
}
if (nrow(audit) != 11) {
  stop(sprintf("Expected 11 audited apoplastic candidates, found %d", nrow(audit)))
}

for (path in c(fig4_pdf, fig4_png, fig6_pdf, fig6_png)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

priority_label <- function(value) {
  ifelse(value == "mechanism_leading", "High", "Intermediate")
}

candidate_label <- function(row) {
  symbol <- row[["mapped_symbols"]]
  family <- row[["apoplastic_cell_wall_host_family"]]
  detail <- if (nzchar(symbol)) symbol else if (nzchar(family)) family else "unresolved homolog"
  sprintf("%s  %s", row[["orthogroup_id"]], detail)
}

draw_matrix <- function(values, cell_text, row_labels, column_labels, colors, main, subtitle) {
  nr <- nrow(values)
  nc <- ncol(values)
  plot.new()
  plot.window(xlim = c(0, nc + 5.6), ylim = c(0, nr + 2.2), xaxs = "i", yaxs = "i")
  title(main = main, line = 1.2, cex.main = 1.25, font.main = 2)
  text(0, nr + 1.25, subtitle, adj = 0, cex = 0.78, col = "#4A4A4A")

  for (i in seq_len(nr)) {
    y <- nr - i + 0.75
    text(0, y, row_labels[[i]], adj = 0, cex = 0.72)
    for (j in seq_len(nc)) {
      x <- 4.9 + j - 1
      value <- values[i, j]
      fill <- colors[[as.character(value)]]
      rect(x, y - 0.34, x + 0.82, y + 0.34, col = fill, border = "white", lwd = 1.5)
      label <- cell_text[i, j]
      if (nzchar(label)) {
        text(x + 0.41, y, label, cex = 0.60, col = "white", font = 2)
      }
    }
  }

  for (j in seq_len(nc)) {
    x <- 4.9 + j - 1 + 0.41
    text(x, nr + 0.65, column_labels[[j]], srt = 45, adj = 0, cex = 0.72, xpd = NA)
  }
  box(col = NA)
}

draw_fig4 <- function() {
  priority <- priority_label(candidates$overall_priority_class)
  values <- cbind(
    priority = ifelse(priority == "High", 3, 2),
    concordant = pmin(as.integer(candidates$nonproxy_concordant_species_count), 2) +
      ifelse(as.integer(candidates$nonproxy_concordant_species_count) > 0, 1, 0),
    discordant = pmin(as.integer(candidates$nonproxy_discordant_species_count), 2) +
      ifelse(as.integer(candidates$nonproxy_discordant_species_count) > 0, 1, 0),
    network = ifelse(candidates$experimental_or_database_supported == "true", 3, 1),
    structure = ifelse(as.integer(candidates$structurally_supported_effector_count) > 0, 3, 1),
    sequence = ifelse(candidates$apoplastic_cell_wall_sequence_support_class == "tier1_candidate_substrate", 3, 1)
  )
  values[, "concordant"] <- ifelse(
    as.integer(candidates$nonproxy_concordant_species_count) == 0,
    1,
    values[, "concordant"]
  )
  values[, "discordant"] <- ifelse(
    as.integer(candidates$nonproxy_discordant_species_count) == 0,
    1,
    values[, "discordant"]
  )
  cell_text <- matrix("", nrow = nrow(values), ncol = ncol(values))
  cell_text[, 1] <- ifelse(priority == "High", "H", "I")
  cell_text[, 2] <- ifelse(
    as.integer(candidates$nonproxy_concordant_species_count) > 0,
    candidates$nonproxy_concordant_species_count,
    ""
  )
  cell_text[, 3] <- ifelse(
    as.integer(candidates$nonproxy_discordant_species_count) > 0,
    candidates$nonproxy_discordant_species_count,
    ""
  )
  cell_text[, 4] <- ifelse(candidates$experimental_or_database_supported == "true", "yes", "")
  cell_text[, 5] <- ifelse(as.integer(candidates$structurally_supported_effector_count) > 0, "yes", "")
  cell_text[, 6] <- ifelse(candidates$apoplastic_cell_wall_sequence_support_class == "tier1_candidate_substrate", "yes", "")

  labels <- vapply(seq_len(nrow(candidates)), function(i) {
    candidate_label(candidates[i, , drop = FALSE])
  }, character(1))

  colors <- c("1" = "#ECECEC", "2" = "#D6A541", "3" = "#287D8E")
  draw_matrix(
    values,
    cell_text,
    labels,
    c("Priority", "Concordant species", "Discordant species", "Curated/experimental network", "Effector structure", "Sequence-resolved host candidate"),
    colors,
    "Integrated evidence for 12 prioritized host orthogroups",
    "Filled cells denote retained evidence; numbers indicate the count of directly mapped species."
  )
  legend(
    "bottomright",
    inset = c(0.01, 0.01),
    legend = c("Absent or unavailable", "Intermediate / caution", "High / supported"),
    fill = colors,
    border = NA,
    bty = "n",
    cex = 0.72,
    horiz = TRUE,
    xpd = NA
  )
}

draw_fig6 <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(3.2, 1.1))
  par(mar = c(2.5, 1.2, 4.7, 1.0), family = "sans")

  values <- cbind(
    module = ifelse(audit$manuscript_module_qualifying == "yes", 3, 1),
    concordant = ifelse(as.integer(audit$nonproxy_concordant_species_count) > 0, 3, 1),
    discordant = ifelse(as.integer(audit$nonproxy_discordant_species_count) > 0, 2, 1),
    network = ifelse(audit$network_supported == "true", 3, 1),
    sequence = ifelse(audit$sequence_support_class == "tier1_candidate_substrate", 3,
      ifelse(audit$sequence_support_class == "localization_or_sequence_caution", 2, 1))
  )
  cell_text <- matrix("", nrow = nrow(values), ncol = ncol(values))
  cell_text[, 1] <- ifelse(audit$manuscript_module_qualifying == "yes", "yes", "")
  cell_text[, 2] <- ifelse(
    as.integer(audit$nonproxy_concordant_species_count) > 0,
    audit$nonproxy_concordant_species_count,
    ""
  )
  cell_text[, 3] <- ifelse(
    as.integer(audit$nonproxy_discordant_species_count) > 0,
    audit$nonproxy_discordant_species_count,
    ""
  )
  cell_text[, 4] <- ifelse(audit$network_supported == "true", "yes", "")
  cell_text[, 5] <- ifelse(
    audit$sequence_support_class == "tier1_candidate_substrate",
    "yes",
    ifelse(audit$sequence_support_class == "localization_or_sequence_caution", "caution", "")
  )
  labels <- sprintf("%s  %s", audit$orthogroup_id, audit$host_family)
  colors <- c("1" = "#ECECEC", "2" = "#D6A541", "3" = "#287D8E")
  draw_matrix(
    values,
    cell_text,
    labels,
    c("Module criterion", "Concordant expression", "Discordant expression", "Network context", "Sequence resolution"),
    colors,
    "Independent apoplastic candidate audit",
    "Each orthogroup is counted once; pale cells indicate absent, unresolved, or unavailable evidence."
  )

  par(mar = c(5.0, 3.6, 4.7, 1.0), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title("Biochemical validation path", cex.main = 1.0, font.main = 2)
  path_labels <- c(
    "11 audited\northogroups",
    "4 module-supported\nfamilies",
    "2 sequence-resolved\nhost proteins",
    "Biochemical cleavage\ntest required"
  )
  path_colors <- c("#686868", "#287D8E", "#D6A541", "white")
  y <- c(0.82, 0.61, 0.40, 0.17)
  for (i in seq_along(y)) {
    rect(0.14, y[i] - 0.075, 0.86, y[i] + 0.075,
         col = path_colors[i], border = ifelse(i == 4, "#C7665B", NA),
         lwd = ifelse(i == 4, 1.8, 1), lty = ifelse(i == 4, 2, 1))
    text(0.50, y[i], path_labels[i], cex = 0.72,
         col = ifelse(i == 4, "#7A3D39", "white"), font = ifelse(i == 4, 2, 1))
    if (i < length(y)) arrows(0.50, y[i] - 0.085, 0.50, y[i + 1] + 0.085,
                               length = 0.06, lwd = 1.2, col = "#6F767A")
  }
  text(0.50, 0.02, "Proposed final step; not a zero-result\nscreening stage.", cex = 0.55, col = "#4A4A4A")
  layout(1)
}

render_pair <- function(pdf_path, png_path, draw_function, width, height) {
  pdf(pdf_path, width = width, height = height, family = "sans", useDingbats = FALSE)
  par(mar = c(2.5, 1.2, 4.7, 1.0), family = "sans")
  draw_function()
  dev.off()

  png(png_path, width = width, height = height, units = "in", res = 300, type = "cairo")
  par(mar = c(2.5, 1.2, 4.7, 1.0), family = "sans")
  draw_function()
  dev.off()
}

render_pair(fig4_pdf, fig4_png, draw_fig4, 12.0, 7.8)
render_pair(fig6_pdf, fig6_png, draw_fig6, 14.0, 7.8)

message("Wrote manuscript evidence figures: ", paste(c(fig4_pdf, fig4_png, fig6_pdf, fig6_png), collapse = ", "))
