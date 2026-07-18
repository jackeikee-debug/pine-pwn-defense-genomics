#!/usr/bin/env Rscript

parse_args <- function(args) {
  values <- list(dry_run = FALSE)
  i <- 1
  while (i <= length(args)) {
    token <- args[[i]]
    if (token == "--dry-run") {
      values$dry_run <- TRUE
      i <- i + 1
    } else if (startsWith(token, "--")) {
      if (i == length(args)) stop("Missing value for ", token)
      key <- gsub("-", "_", substring(token, 3), fixed = TRUE)
      values[[key]] <- args[[i + 1]]
      i <- i + 2
    } else {
      stop("Unexpected argument: ", token)
    }
  }
  required <- c("gene_values", "pathways", "coverage_output", "output_dir")
  missing <- required[!vapply(required, function(x) !is.null(values[[x]]), logical(1))]
  if (length(missing)) stop("Missing required arguments: ", paste(missing, collapse = ", "))
  values
}

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
             check.names = FALSE, quote = "", comment.char = "")
}

fetch_membership <- function(pathways) {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    stop("KEGGREST is required when --pathway-membership is not supplied")
  }
  rows <- lapply(pathways$pathway_id, function(pathway_id) {
    links <- KEGGREST::keggLink("ath", paste0("path:", pathway_id))
    data.frame(
      pathway_id = pathway_id,
      tair_locus = sub("^ath:", "", unname(links)),
      stringsAsFactors = FALSE
    )
  })
  unique(do.call(rbind, rows))
}

safe_name <- function(value) {
  gsub("[^A-Za-z0-9_.-]+", "_", value)
}

pathview_output_filename <- function(pathway_id, contrast_id) {
  paste0(pathway_id, ".", safe_name(contrast_id), ".png")
}

pathview_output_path <- function(output_dir, pathway_id, contrast_id) {
  normalizePath(
    file.path(output_dir, pathview_output_filename(pathway_id, contrast_id)),
    winslash = "/", mustWork = FALSE
  )
}

ensure_pathview_runtime <- function() {
  if (!requireNamespace("pathview", quietly = TRUE)) stop("pathview is not installed")
  # pathview 1.50.0 resolves its bundled `bods` table through the attached
  # package environment rather than its namespace when gene.idtype is used.
  suppressPackageStartupMessages(library("pathview", character.only = TRUE))
}

render_pathway <- function(pathway_id, contrast_id, values, output_dir) {
  ensure_pathview_runtime()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  suffix <- safe_name(contrast_id)
  pathway_number <- sub("^ath", "", pathway_id)
  output_path <- pathview_output_path(output_dir, pathway_id, contrast_id)
  gene_values <- values$log2_fold_change
  names(gene_values) <- values$tair_locus
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  setwd(output_dir)
  result <- suppressMessages(pathview::pathview(
    gene.data = gene_values,
    pathway.id = pathway_number,
    species = "ath",
    gene.idtype = "KEGG",
    out.suffix = suffix,
    kegg.native = TRUE,
    same.layer = TRUE,
    low = list(gene = "#2166AC", cpd = "#2166AC"),
    mid = list(gene = "#F7F7F7", cpd = "#F7F7F7"),
    high = list(gene = "#B2182B", cpd = "#B2182B"),
    limit = list(gene = 3, cpd = 1)
  ))
  pathway_result <- if (!is.null(result$plot.data.gene)) result else result[[1]]
  gene_nodes <- pathway_result$plot.data.gene
  rendered_node_count <- if (is.null(gene_nodes)) {
    0L
  } else {
    sum(!is.na(gene_nodes$all.mapped) & gene_nodes$all.mapped != "")
  }
  list(
    output_path = output_path,
    rendered_node_count = rendered_node_count
  )
}

build_coverage <- function(gene_values, pathways, membership, output_dir, dry_run) {
  contrasts <- sort(unique(gene_values$contrast_id))
  rows <- list()
  index <- 1
  for (pathway_index in seq_len(nrow(pathways))) {
    pathway <- pathways[pathway_index, , drop = FALSE]
    members <- unique(toupper(membership$tair_locus[membership$pathway_id == pathway$pathway_id]))
    minimum <- as.integer(pathway$minimum_mapped_genes)
    for (contrast_id in contrasts) {
      contrast_values <- gene_values[gene_values$contrast_id == contrast_id, , drop = FALSE]
      contrast_values <- contrast_values[toupper(contrast_values$tair_locus) %in% members, , drop = FALSE]
      mapped <- nrow(contrast_values)
      status <- if (mapped >= minimum) "eligible" else "insufficient_coverage"
      output_path <- ""
      rendered_node_count <- NA_integer_
      error <- ""
      if (!dry_run && status == "eligible") {
        result <- tryCatch(
          render_pathway(pathway$pathway_id, contrast_id, contrast_values, output_dir),
          error = function(condition) condition
        )
        if (inherits(result, "error")) {
          status <- "render_failed"
          error <- conditionMessage(result)
        } else {
          rendered_node_count <- result$rendered_node_count
          output_path <- result$output_path
          if (rendered_node_count > 0) {
            status <- "rendered"
          } else {
            status <- "render_failed"
            error <- "Pathview produced no mapped gene nodes"
          }
        }
      }
      rows[[index]] <- data.frame(
        pathway_id = pathway$pathway_id,
        pathway_name = pathway$pathway_name,
        mechanism_axis = pathway$mechanism_axis,
        contrast_id = contrast_id,
        pathway_gene_count = length(members),
        mapped_gene_count = mapped,
        rendered_node_count = rendered_node_count,
        minimum_mapped_genes = minimum,
        status = status,
        output_path = output_path,
        error = error,
        stringsAsFactors = FALSE
      )
      index <- index + 1
    }
  }
  do.call(rbind, rows)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  gene_values <- read_tsv(args$gene_values)
  gene_values$log2_fold_change <- as.numeric(gene_values$log2_fold_change)
  pathways <- read_tsv(args$pathways)
  membership <- if (!is.null(args$pathway_membership)) {
    read_tsv(args$pathway_membership)
  } else {
    fetched <- fetch_membership(pathways)
    if (!is.null(args$membership_output)) {
      dir.create(dirname(args$membership_output), recursive = TRUE, showWarnings = FALSE)
      write.table(fetched, args$membership_output, sep = "\t", row.names = FALSE, quote = FALSE)
    }
    fetched
  }
  coverage <- build_coverage(gene_values, pathways, membership, args$output_dir, args$dry_run)
  dir.create(dirname(args$coverage_output), recursive = TRUE, showWarnings = FALSE)
  write.table(coverage, args$coverage_output, sep = "\t", row.names = FALSE, quote = FALSE,
              na = "")
  message("Pathview audit: ", paste(names(table(coverage$status)), table(coverage$status),
                                    sep = "=", collapse = "; "))
}

if (sys.nframe() == 0) main()
