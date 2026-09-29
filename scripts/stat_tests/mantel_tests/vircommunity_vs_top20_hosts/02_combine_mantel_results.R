# Usage: Rscript 02_combine_mantel_results.R prok_vs_caudo_mantel_commcomp/results or Rscript 02_combine_mantel_results.R euk_vs_mega_mantel_commcomp/results

args <- commandArgs(trailingOnly = TRUE)

results_dir <- args[1]

# Find all individual Mantel result files
files <- list.files(
  results_dir,
  pattern = "^mantel_[0-9]+\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No Mantel result files found in: ", results_dir)
}

# Read and combine
res <- do.call(
  rbind,
  lapply(files, read.csv)
)

# BH / FDR correction
res$p_adj_BH <- p.adjust(
  res$p_value,
  method = "BH"
)

# Add significance indicator
res$significant_FDR_0.05 <- res$p_adj_BH < 0.05

# Sort by adjusted p-value
res <- res[
  order(res$p_adj_BH, res$p_value),
]

# Write compiled result
outfile <- file.path(
  results_dir,
  "mantel_all_genera_BH.csv"
)

write.csv(
  res,
  outfile,
  row.names = FALSE
)

cat("Found", length(files), "Mantel result files\n")
cat("Compiled results written to:\n", outfile, "\n")
cat("Significant at FDR < 0.05:",
    sum(res$significant_FDR_0.05, na.rm = TRUE), "\n")

print(res)