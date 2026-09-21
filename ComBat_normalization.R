required_sva_version <- "3.54.0"
if (!requireNamespace("sva", quietly = TRUE)) stop("Install sva version ", required_sva_version, ".")
if (as.character(utils::packageVersion("sva")) != required_sva_version) stop("Required sva version: ", required_sva_version, "; installed: ", utils::packageVersion("sva"))
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
set.seed(42)

# ---- Settings ---------------------------------------------------------------
data_file <- "dat_untargeted_normalized_metab_df.csv"
batch_file <- "batch_collection_years_metab.csv"
group_file <- "biological_group_SGA_metab.csv"
output_dir <- "combat_output"
analysis <- "SGA" # Set to "APGAR" for the low 1-minute APGAR (<7) analysis.
apgar_score_column <- "APGAR_1min" # Amend if your metadata uses another column name.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Read and align data ----------------------------------------------------
expression_df <- utils::read.csv(data_file, row.names = 1, check.names = FALSE)
batch_info_df <- utils::read.csv(batch_file, row.names = 1, check.names = FALSE)
group_info_df <- utils::read.csv(group_file, row.names = 1, check.names = FALSE)
expression_matrix <- as.matrix(expression_df)
storage.mode(expression_matrix) <- "numeric"
if (any(!is.finite(expression_matrix))) stop("Expression matrix contains missing or non-finite values.")

sample_order <- colnames(expression_matrix)
if (length(setdiff(sample_order, rownames(batch_info_df)))) stop("Batch metadata is missing one or more samples.")
if (length(setdiff(sample_order, rownames(group_info_df)))) stop("Outcome metadata is missing one or more samples.")
batch_info_df <- batch_info_df[sample_order, , drop = FALSE]
group_info_df <- group_info_df[sample_order, , drop = FALSE]

collection_year <- suppressWarnings(as.numeric(as.character(batch_info_df$DateDel2)))
if (anyNA(collection_year)) stop("DateDel2 must be numeric and complete.")
batch_vector <- factor(ifelse(collection_year < 2015, "Before 2015", "2015 and after"))
if (nlevels(batch_vector) < 2) stop("ComBat requires at least two batches.")

# ---- Preserve the outcome relevant to this analysis ------------------------
if (analysis == "SGA") {
  if (!"SGA" %in% names(group_info_df)) stop("SGA column is absent.")
  outcome <- factor(group_info_df$SGA, levels = c(0, 1), labels = c("AGA", "SGA"))
  if (anyNA(outcome)) stop("SGA must be coded 0 (AGA) or 1 (SGA).")
  outcome_label <- "SGA status"
} else if (analysis == "APGAR") {
  if (!apgar_score_column %in% names(group_info_df)) stop("APGAR column is absent: ", apgar_score_column)
  apgar_score <- suppressWarnings(as.numeric(as.character(group_info_df[[apgar_score_column]])))
  if (anyNA(apgar_score)) stop("APGAR score must be numeric and complete.")
  outcome <- factor(ifelse(apgar_score < 7, "APGAR <7", "APGAR >=7"), levels = c("APGAR >=7", "APGAR <7"))
  outcome_label <- "Low 1-minute APGAR score (<7)"
} else stop("analysis must be 'SGA' or 'APGAR'.")
if (nlevels(droplevels(outcome)) < 2) stop("Both outcome groups must be present.")

# The design matrix protects outcome-related variation during correction.
model_matrix <- stats::model.matrix(~ outcome)
combat_corrected_matrix <- sva::ComBat(expression_matrix, batch_vector, model_matrix, par.prior = TRUE)

analysis_tag <- tolower(analysis)
utils::write.csv(combat_corrected_matrix, file.path(output_dir, paste0("combat_corrected_", analysis_tag, ".csv")))

# PCA plots provide a before/after quality-control record.
make_pca_plot <- function(matrix, stage) {
  pca <- stats::prcomp(t(matrix), scale. = FALSE)
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Batch = batch_vector, Outcome = outcome)
  variance <- summary(pca)$importance[2, 1:2] * 100
  ggplot2::ggplot(d, ggplot2::aes(PC1, PC2, color = Batch, shape = Outcome)) +
    ggplot2::geom_point(size = 3.5, alpha = .85) +
    ggplot2::scale_color_manual(values = c("Before 2015" = "#0072B2", "2015 and after" = "#D55E00")) +
    ggplot2::labs(title = stage, x = paste0("PC1: ", round(variance[1], 1), "% variance"), y = paste0("PC2: ", round(variance[2], 1), "% variance"), color = "Collection years", shape = outcome_label) +
    ggplot2::theme_classic(base_size = 14)
}
pca_before <- make_pca_plot(expression_matrix, "Before ComBat")
pca_after <- make_pca_plot(combat_corrected_matrix, "After ComBat")
ggplot2::ggsave(file.path(output_dir, paste0("pca_before_", analysis_tag, ".png")), pca_before, width = 7, height = 5)
ggplot2::ggsave(file.path(output_dir, paste0("pca_after_", analysis_tag, ".png")), pca_after, width = 7, height = 5)
writeLines(capture.output(sessionInfo()), file.path(output_dir, paste0("sessionInfo_", analysis_tag, ".txt")))
print(pca_before)
print(pca_after)
