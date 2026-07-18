#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(DESeq2)
    library(tximport)
})

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) %% 2 != 0) {
    stop("Arguments must be supplied as --name value pairs", call. = FALSE)
}
options <- setNames(arguments[seq(2, length(arguments), by = 2)], sub("^--", "", arguments[seq(1, length(arguments), by = 2)]))
required <- c("sample-sheet", "quant-manifest", "tx2gene", "output")
if (!all(required %in% names(options))) {
    stop(paste("Missing required arguments:", paste(setdiff(required, names(options)), collapse = ", ")), call. = FALSE)
}

normalize_path <- function(path) {
    path <- gsub("\\\\", "/", path)
    sub("^([A-Za-z]):/", "/mnt/\\L\\1/", path, perl = TRUE)
}

sample_sheet <- read.delim(options[["sample-sheet"]], check.names = FALSE, stringsAsFactors = FALSE)
manifest <- read.delim(options[["quant-manifest"]], check.names = FALSE, stringsAsFactors = FALSE)
tx2gene <- read.delim(options[["tx2gene"]], check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample_id", "phenotype", "inoculum") %in% names(sample_sheet))) {
    stop("Sample sheet must contain sample_id, phenotype, and inoculum", call. = FALSE)
}
if (!all(c("sample_id", "local_path") %in% names(manifest))) {
    stop("Quantification manifest must contain sample_id and local_path", call. = FALSE)
}
if (nrow(sample_sheet) != 12 || anyDuplicated(sample_sheet$sample_id)) {
    stop("Interaction model requires exactly 12 unique samples", call. = FALSE)
}
if (!setequal(sample_sheet$sample_id, manifest$sample_id)) {
    stop("Sample sheet and quantification manifest sample IDs differ", call. = FALSE)
}

manifest <- manifest[match(sample_sheet$sample_id, manifest$sample_id), ]
files <- setNames(vapply(manifest$local_path, normalize_path, character(1)), sample_sheet$sample_id)
if (!all(file.exists(files))) {
    stop(paste("Missing quant.sf file(s):", paste(names(files)[!file.exists(files)], collapse = ", ")), call. = FALSE)
}

sample_sheet$phenotype <- relevel(factor(sample_sheet$phenotype), ref = "resistant")
sample_sheet$inoculum <- relevel(factor(sample_sheet$inoculum), ref = "water_control")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
dds <- DESeqDataSetFromTximport(txi, colData = sample_sheet, design = ~ phenotype + inoculum + phenotype:inoculum)
dds <- DESeq(dds, quiet = TRUE)

interaction_candidates <- resultsNames(dds)[
    grepl("^phenotype.*susceptible.*inoculum.*b.xylophilus", resultsNames(dds))
]
if (length(interaction_candidates) != 1) {
    stop(paste("Could not identify a unique susceptible-by-PWN interaction coefficient:", paste(resultsNames(dds), collapse = ", ")), call. = FALSE)
}
coefficient <- interaction_candidates[[1]]
result <- as.data.frame(results(dds, name = coefficient))
result$feature_id <- rownames(result)
result$interaction_coefficient <- coefficient
result <- result[, c("feature_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "interaction_coefficient")]
names(result)[3:7] <- c("interaction_log2_fold_change", "interaction_lfc_se", "interaction_stat", "interaction_pvalue", "interaction_padj")

dir.create(dirname(options[["output"]]), recursive = TRUE, showWarnings = FALSE)
write.table(result, options[["output"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote ", nrow(result), " interaction results using coefficient ", coefficient)
