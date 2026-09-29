# ============================================================
# 08_run_vcr_environment_mantel.R
#
# VCR ~ environmental parameters only
# For cellular-v-cellular only
#
# Usage:
#   Rscript 08_run_vcr_environment_mantel.R 999
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(vegan)
})

# ============================================================
# Settings
# ============================================================

base_dir <- "/path/to/base/dir"

args <- commandArgs(trailingOnly = TRUE)
permutations <- ifelse(length(args) >= 1, as.integer(args[1]), 999)

if (is.na(permutations)) {
  stop("permutations must be an integer.")
}

method <- "spearman"

run_label <- "cellular-v-cellular"
out_dir <- file.path(base_dir, paste0(run_label, "_mantel"))

input_file <- file.path(out_dir, "mantel_inputs_filtered.rds")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

cat("============================================================\n")
cat("Running VCR ~ environment Mantel tests\n")
cat("run_label: ", run_label, "\n", sep = "")
cat("permutations: ", permutations, "\n", sep = "")
cat("input_file: ", input_file, "\n", sep = "")
cat("out_dir: ", out_dir, "\n", sep = "")
cat("============================================================\n")

# ============================================================
# Load prepared inputs
# ============================================================

inputs <- readRDS(input_file)

env_by_node <- inputs$env_by_node
env_vars <- inputs$env_vars
vcr_nodes <- inputs$vcr_nodes

env_prok_filtered <- env_by_node$Prokaryotes
env_euk_filtered  <- env_by_node$Eukaryotes

cat("Loaded VCR nodes:\n")
print(names(vcr_nodes))

cat("Loaded env vars:\n")
print(env_vars)

# ============================================================
# Helpers
# ============================================================

make_vector_dist <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x) & is.finite(x)]
  
  if (length(x) < 3) return(NULL)
  if (length(unique(x)) < 2) return(NULL)
  
  d <- stats::dist(scale(x), method = "euclidean")
  
  if (any(is.na(d))) return(NULL)
  
  as.dist(as.matrix(d))
}

run_vcr_env_mantel <- function(
    vcr_nodes,
    env_table,
    env_vars,
    method = "spearman",
    permutations = 999,
    seed = 123
) {
  
  set.seed(seed)
  results <- list()
  
  for (vcr_name in names(vcr_nodes)) {
    
    vcr_df <- vcr_nodes[[vcr_name]]
    
    for (env_var in env_vars) {
      
      if (!env_var %in% colnames(env_table)) next
      
      tmp <- env_table %>%
        select(Run, env_value = all_of(env_var)) %>%
        filter(!is.na(env_value), is.finite(env_value)) %>%
        group_by(Run) %>%
        summarise(
          env_value = mean(env_value, na.rm = TRUE),
          .groups = "drop"
        )
      
      common_runs <- intersect(vcr_df$Run, tmp$Run)
      
      if (length(common_runs) < 3) next
      
      vcr_sub <- vcr_df %>%
        filter(Run %in% common_runs) %>%
        arrange(match(Run, common_runs))
      
      env_sub <- tmp %>%
        filter(Run %in% common_runs) %>%
        arrange(match(Run, common_runs))
      
      keep <- !is.na(vcr_sub$value) &
        is.finite(vcr_sub$value) &
        !is.na(env_sub$env_value) &
        is.finite(env_sub$env_value)
      
      vcr_sub <- vcr_sub[keep, , drop = FALSE]
      env_sub <- env_sub[keep, , drop = FALSE]
      
      if (nrow(vcr_sub) < 3) next
      if (length(unique(vcr_sub$value)) < 2) next
      if (length(unique(env_sub$env_value)) < 2) next
      
      vcr_dist <- make_vector_dist(vcr_sub$value)
      env_dist <- make_vector_dist(env_sub$env_value)
      
      if (is.null(vcr_dist) || is.null(env_dist)) next
      
      if (attr(vcr_dist, "Size") != attr(env_dist, "Size")) {
        stop(
          "Distance Size mismatch for ",
          vcr_name, " vs ", env_var,
          ": vcr_size=", attr(vcr_dist, "Size"),
          ", env_size=", attr(env_dist, "Size")
        )
      }
      
      message("Mantel: ", vcr_name, " vs ", env_var, " n=", attr(vcr_dist, "Size"))
      
      res <- vegan::mantel(
        vcr_dist,
        env_dist,
        method = method,
        permutations = permutations
      )
      
      results[[paste(vcr_name, env_var, sep = "__")]] <- tibble(
        dataset1 = vcr_name,
        dataset2 = env_var,
        comparison_type = "vcr_environment",
        mantel_r = unname(res$statistic),
        p_value = res$signif,
        n_samples = attr(vcr_dist, "Size"),
        permutations = permutations,
        method = method
      )
    }
  }
  
  bind_rows(results)
}

# ============================================================
# Run VCR ~ environment
# ============================================================

message("Running VCR_Caudo_Prok ~ environment using Prokaryotes env table")

vcr_env_results_caudo_prok <- run_vcr_env_mantel(
  vcr_nodes = list(
    VCR_Caudo_Prok = vcr_nodes$VCR_Caudo_Prok
  ),
  env_table = env_prok_filtered,
  env_vars = env_vars,
  method = method,
  permutations = permutations
)

message("Running VCR_Mega_Euk ~ environment using Eukaryotes env table")

vcr_env_results_mega_euk <- run_vcr_env_mantel(
  vcr_nodes = list(
    VCR_Mega_Euk = vcr_nodes$VCR_Mega_Euk
  ),
  env_table = env_euk_filtered,
  env_vars = env_vars,
  method = method,
  permutations = permutations
)

vcr_env_results <- bind_rows(
  vcr_env_results_caudo_prok,
  vcr_env_results_mega_euk
)

# ============================================================
# Save output
# ============================================================

outfile <- file.path(out_dir, "vcr_environment_mantel_08.tsv")

if (nrow(vcr_env_results) == 0) {
  message("No VCR-environment results produced. Writing empty placeholder.")
  
  vcr_env_results <- tibble(
    dataset1 = character(),
    dataset2 = character(),
    comparison_type = character(),
    mantel_r = numeric(),
    p_value = numeric(),
    n_samples = integer(),
    permutations = integer(),
    method = character()
  )
}

write_tsv(vcr_env_results, outfile)

message("Done.")
message("Wrote: ", outfile)