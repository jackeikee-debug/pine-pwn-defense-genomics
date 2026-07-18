#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: plot_secretome_audit.R <stages.tsv> <combinations.tsv> <output.pdf> <output.png>")
}

stages <- read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
combinations <- read.delim(args[[2]], stringsAsFactors = FALSE, check.names = FALSE)
pdf_path <- args[[3]]
png_path <- args[[4]]

if (!identical(stages$protein_count, c(15860L, 2741L, 2164L, 1265L, 401L))) {
  stop("Secretome audit stage counts do not match the frozen analysis")
}
if (sum(combinations$protein_count) != 2164L) stop("UpSet combinations must total 2,164")

draw_figure <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.0, 1.35))
  palette <- c(ink = "#273238", teal = "#287D8E", gold = "#D7A83E", red = "#C45A52", pale = "#E8ECEC")

  par(mar = c(2.5, 1.0, 4.2, 1.2), family = "sans")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.02, 0.97, "A", adj = c(0, 1), font = 2, cex = 1.25)
  text(0.10, 0.97, "Sequence filter and independent support", adj = c(0, 1), font = 2, cex = 0.93)

  nested <- stages[stages$set_relation %in% c("starting_set", "nested_filter"), ]
  y <- c(0.78, 0.57, 0.36)
  widths <- sqrt(nested$protein_count / max(nested$protein_count)) * 0.82
  fills <- c(palette[["ink"]], palette[["teal"]], palette[["gold"]])
  for (i in seq_len(nrow(nested))) {
    rect(0.48 - widths[i] / 2, y[i] - 0.065, 0.48 + widths[i] / 2, y[i] + 0.065,
         col = fills[i], border = NA)
    text(0.48, y[i] + 0.012, format(nested$protein_count[i], big.mark = ","),
         col = "white", font = 2, cex = 0.88)
    text(0.48, y[i] - 0.032, nested$stage_label[i], col = "white", cex = 0.56)
    if (i < nrow(nested)) arrows(0.48, y[i] - 0.075, 0.48, y[i + 1] + 0.078,
                                  length = 0.06, lwd = 1.2, col = "#70787C")
  }
  support <- stages[stages$set_relation == "independent_support_within_2164", ]
  sx <- c(0.27, 0.70)
  scol <- c(palette[["teal"]], palette[["red"]])
  for (i in seq_len(nrow(support))) {
    rect(sx[i] - 0.18, 0.105, sx[i] + 0.18, 0.225, col = scol[i], border = NA)
    text(sx[i], 0.176, format(support$protein_count[i], big.mark = ","), col = "white", font = 2, cex = 0.85)
    text(sx[i], 0.135, support$stage_label[i], col = "white", cex = 0.50)
  }
  text(0.48, 0.275, "Independent evidence channels within the 2,164 candidates", cex = 0.58, col = "#4A5357")
  text(0.48, 0.025, "Support sets overlap and are not displayed as a false funnel.", cex = 0.55, col = "#565E62")

  par(mar = c(5.7, 4.5, 4.2, 1.0), family = "sans")
  combinations <- combinations[order(combinations$protein_count, decreasing = TRUE), ]
  x <- seq_len(nrow(combinations))
  ymax <- max(combinations$protein_count) * 1.12
  plot(NA, xlim = c(-0.45, length(x) + 0.55), ylim = c(-520, ymax), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i", yaxs = "i")
  text(0.50, ymax * 1.04, "B", adj = c(0, 1), font = 2, cex = 1.25, xpd = NA)
  text(1.0, ymax * 1.04, "Published-secretome set intersections", adj = c(0, 1), font = 2, cex = 0.93, xpd = NA)
  bars <- x
  bar_fill <- ifelse(combinations$published_support == "yes", palette[["red"]], "#B9C0C2")
  rect(bars - 0.34, 0, bars + 0.34, combinations$protein_count, col = bar_fill, border = NA)
  axis(2, las = 1, cex.axis = 0.68)
  mtext("Proteins", side = 2, line = 2.8, cex = 0.72)
  text(bars, combinations$protein_count + ymax * 0.025, format(combinations$protein_count, big.mark = ","), cex = 0.62, font = 2)

  set_y <- c(-120, -260, -400)
  set_names <- c("Shinya", "Cardoso", "Silva")
  membership <- rbind(combinations$shinya == "yes", combinations$cardoso == "yes", combinations$silva == "yes")
  for (j in seq_along(set_y)) text(0.48, set_y[j], set_names[j], adj = 1, cex = 0.67)
  for (i in seq_len(ncol(membership))) {
    active <- which(membership[, i])
    points(rep(bars[i], 3), set_y, pch = 21, bg = ifelse(membership[, i], palette[["ink"]], "white"), col = "#9AA1A4", cex = 0.95)
    if (length(active) > 1) segments(bars[i], min(set_y[active]), bars[i], max(set_y[active]), lwd = 1.6, col = palette[["ink"]])
  }
  abline(h = 0, col = "#BFC4C6")
  text(mean(range(bars)), -505, "All columns are subsets of the current 2,164 secreted-soluble candidates.", cex = 0.56, col = "#565E62")
  layout(1)
}

dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
pdf(pdf_path, width = 12.4, height = 7.1, family = "sans", useDingbats = FALSE); draw_figure(); dev.off()
png(png_path, width = 12.4, height = 7.1, units = "in", res = 300, type = "cairo"); draw_figure(); dev.off()
