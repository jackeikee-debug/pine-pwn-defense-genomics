#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(fgsea))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) stop("Arguments must be --name value pairs", call. = FALSE)
opts <- setNames(args[seq(2, length(args), 2)], sub("^--", "", args[seq(1, length(args), 2)]))
required <- c("annotations", "mechanism-axes", "ranked-output", "ora-output")
if (!all(required %in% names(opts))) {
    stop(paste("Missing arguments:", paste(setdiff(required, names(opts)), collapse = ", ")), call. = FALSE)
}

annotations <- read.delim(opts[["annotations"]], check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
axes_config <- read.delim(opts[["mechanism-axes"]], check.names = FALSE, stringsAsFactors = FALSE, na.strings = character())

`%||%` <- function(x, y) if (is.null(x)) y else x

split_membership <- function(values) {
    result <- list()
    for (i in seq_along(values)) {
        if (is.na(values[[i]]) || values[[i]] == "") next
        for (set_id in strsplit(values[[i]], ";", fixed = TRUE)[[1]]) {
            result[[set_id]] <- unique(c(result[[set_id]], annotations$feature_id[[i]]))
        }
    }
    result
}

module_sets <- split_membership(annotations$module_ids)
axis_sets_observed <- split_membership(annotations$mechanism_axes)
configured_axes <- axes_config$mechanism_axis
axis_sets <- setNames(lapply(configured_axes, function(axis) axis_sets_observed[[axis]] %||% character()), configured_axes)

ranked_rows <- annotations[annotations$finite_wald == "yes" & is.finite(annotations$interaction_stat), ]
ranked_rows <- ranked_rows[order(ranked_rows$feature_id), ]
ranks <- setNames(ranked_rows$interaction_stat, ranked_rows$feature_id)
# fgsea requires a strict order. This deterministic perturbation is many orders
# below the reported Wald precision and only resolves exact ties by feature ID.
ranks <- ranks + (seq_along(ranks) - mean(seq_along(ranks))) * 1e-12
ranks <- sort(ranks, decreasing = TRUE)

run_ranked_level <- function(sets, level, max_size) {
    definitions <- data.frame(
        set_level = level,
        set_id = names(sets),
        set_size = vapply(sets, function(x) length(intersect(x, names(ranks))), integer(1)),
        stringsAsFactors = FALSE
    )
    definitions$nes <- NA_real_
    definitions$pvalue <- NA_real_
    definitions$padj <- NA_real_
    definitions$leading_edge_features <- ""
    definitions$competition_status <- "not_testable"
    testable <- definitions$set_size >= 10 & definitions$set_size <= max_size
    pathways <- sets[definitions$set_id[testable]]
    if (length(pathways) > 0) {
        result <- as.data.frame(fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 10, maxSize = max_size))
        for (i in seq_len(nrow(result))) {
            j <- match(result$pathway[[i]], definitions$set_id)
            definitions$nes[[j]] <- result$NES[[i]]
            definitions$pvalue[[j]] <- result$pval[[i]]
            definitions$padj[[j]] <- result$padj[[i]]
            definitions$leading_edge_features[[j]] <- paste(result$leadingEdge[[i]], collapse = ";")
            definitions$competition_status[[j]] <- if (!is.na(result$padj[[i]]) && result$padj[[i]] < 0.05) "significant" else "not_significant"
        }
    }
    definitions
}

ranked_module <- run_ranked_level(module_sets, "defense_module", 500)
ranked_axis <- run_ranked_level(axis_sets, "mechanism_axis", Inf)
ranked_output <- rbind(ranked_module, ranked_axis)

background <- annotations$feature_id[annotations$finite_padj == "yes"]
foreground <- annotations$feature_id[annotations$foreground_850 == "yes"]
if (!all(foreground %in% background)) stop("Every foreground feature must occur in the finite-P-value background", call. = FALSE)

run_ora_level <- function(sets, level, max_size) {
    rows <- lapply(names(sets), function(set_id) {
        members <- intersect(sets[[set_id]], background)
        a <- length(intersect(members, foreground))
        b <- length(setdiff(foreground, members))
        nonforeground <- setdiff(background, foreground)
        c_count <- length(intersect(members, nonforeground))
        d <- length(setdiff(nonforeground, members))
        testable <- length(members) >= 10 && length(members) <= max_size
        estimate <- pvalue <- ci_low <- ci_high <- NA_real_
        if (testable) {
            test <- fisher.test(matrix(c(a, b, c_count, d), nrow = 2), alternative = "two.sided")
            estimate <- unname(test$estimate)
            pvalue <- test$p.value
            ci_low <- test$conf.int[[1]]
            ci_high <- test$conf.int[[2]]
        }
        data.frame(
            set_level = level, set_id = set_id, set_size = length(members),
            foreground_in_set = a, foreground_not_in_set = b,
            background_nonforeground_in_set = c_count, background_nonforeground_not_in_set = d,
            odds_ratio = estimate, confidence_interval_low = ci_low, confidence_interval_high = ci_high,
            pvalue = pvalue, padj = NA_real_,
            competition_status = if (testable) "not_significant" else "not_testable",
            stringsAsFactors = FALSE
        )
    })
    output <- if (length(rows)) do.call(rbind, rows) else data.frame()
    if (nrow(output)) {
        testable <- is.finite(output$pvalue)
        output$padj[testable] <- p.adjust(output$pvalue[testable], method = "BH")
        output$competition_status[testable & output$padj < 0.05] <- "significant"
    }
    output
}

ora_module <- run_ora_level(module_sets, "defense_module", 500)
ora_axis <- run_ora_level(axis_sets, "mechanism_axis", Inf)
ora_output <- rbind(ora_module, ora_axis)

dir.create(dirname(opts[["ranked-output"]]), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opts[["ora-output"]]), recursive = TRUE, showWarnings = FALSE)
write.table(ranked_output, opts[["ranked-output"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(ora_output, opts[["ora-output"]], sep = "\t", quote = FALSE, row.names = FALSE, na = "")
message("Wrote ", nrow(ranked_output), " ranked and ", nrow(ora_output), " over-representation rows")
