#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(tximport)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) %% 2 != 0) stop("Arguments must be --name value pairs", call. = FALSE)
options <- setNames(arguments[seq(2, length(arguments), 2)], sub("^--", "", arguments[seq(1, length(arguments), 2)]))
required <- c("sample-sheet", "quant-manifest", "output-two-week", "output-four-week")
if (!all(required %in% names(options))) stop(paste("Missing:", paste(setdiff(required, names(options)), collapse = ", ")), call. = FALSE)

normalize_path <- function(path) {
    path <- gsub("\\\\", "/", path)
    if (.Platform$OS.type == "windows") return(path)
    sub("^([A-Za-z]):/", "/mnt/\\L\\1/", path, perl = TRUE)
}

samples <- read.delim(options[["sample-sheet"]], stringsAsFactors = FALSE, check.names = FALSE)
manifest <- read.delim(options[["quant-manifest"]], stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample_id", "timepoint", "replicate") %in% names(samples))) stop("Sample sheet requires sample_id, timepoint, replicate", call. = FALSE)
if (!all(c("sample_id", "local_path") %in% names(manifest))) stop("Manifest requires sample_id and local_path", call. = FALSE)
expected_counts_text <- if ("expected-timepoint-counts" %in% names(options)) options[["expected-timepoint-counts"]] else "0w=3,2w=3,4w=3"
expected_parts <- strsplit(expected_counts_text, ",", fixed = TRUE)[[1]]
expected_pairs <- strsplit(expected_parts, "=", fixed = TRUE)
if (any(lengths(expected_pairs) != 2)) stop("Expected timepoint counts must use timepoint=count pairs", call. = FALSE)
expected_counts <- setNames(as.integer(vapply(expected_pairs, `[[`, character(1), 2)), vapply(expected_pairs, `[[`, character(1), 1))
if (any(is.na(expected_counts)) || any(expected_counts < 2)) stop("Each expected timepoint count must be an integer >= 2", call. = FALSE)
if (nrow(samples) != sum(expected_counts) || anyDuplicated(samples$sample_id)) stop("Sample count or sample IDs do not match the declared design", call. = FALSE)
timepoint_counts <- table(samples$timepoint)
if (!setequal(names(timepoint_counts), names(expected_counts)) || any(timepoint_counts[names(expected_counts)] != expected_counts)) stop("Observed timepoint counts do not match --expected-timepoint-counts", call. = FALSE)
if (!setequal(samples$sample_id, manifest$sample_id)) stop("Sample and manifest IDs differ", call. = FALSE)

manifest <- manifest[match(samples$sample_id, manifest$sample_id), ]
files <- setNames(vapply(manifest$local_path, normalize_path, character(1)), samples$sample_id)
if (!all(file.exists(files))) stop(paste("Missing quant files:", paste(names(files)[!file.exists(files)], collapse = ", ")), call. = FALSE)
transcripts <- read.delim(files[[1]], nrows = -1, stringsAsFactors = FALSE)$Name
tx2gene <- data.frame(TXNAME = transcripts, GENEID = transcripts, stringsAsFactors = FALSE)

samples$timepoint <- relevel(factor(samples$timepoint), ref = "0w")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM", dropInfReps = TRUE)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ timepoint)
dds <- DESeq(dds, quiet = TRUE)

write_contrast <- function(level, output) {
    result <- as.data.frame(results(dds, contrast = c("timepoint", level, "0w")))
    result$gene_id <- rownames(result)
    result <- result[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    write.table(result, output, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_contrast("2w", options[["output-two-week"]])
write_contrast("4w", options[["output-four-week"]])
message("Wrote 2w-vs-0w and 4w-vs-0w DESeq2 results for ", nrow(dds), " genes")
