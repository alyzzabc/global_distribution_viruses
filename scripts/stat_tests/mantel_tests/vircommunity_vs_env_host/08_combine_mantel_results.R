# ============================================================
# 06_combine_mantel_results.R
#
# Combine ordinary + partial Mantel results into one Cytoscape edge table.
#
# Usage:
#   Rscript 06_combine_mantel_results.R virus-v-cellular
#   Rscript 06_combine_mantel_results.R cellular-v-cellular
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
})

# ============================================================
# Arguments and paths
# ============================================================

base_dir <- "/path/to/base/dir"

args <- commandArgs(trailingOnly = TRUE)

run_label <- ifelse(length(args) >= 1, args[1], "virus-v-cellular")

out_dir <- file.path(base_dir, paste0(run_label, "_mantel"))

if (!dir.exists(out_dir)) {
  stop("Output directory does not exist: ", out_dir)
}

message("Combining Mantel results from: ", out_dir)

# ============================================================
# Helper: read files safely
# ============================================================

read_result_files <- function(pattern, result_group) {
  
  files <- list.files(
    out_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    message("No files found for ", result_group, " using pattern: ", pattern)
    return(tibble())
  }
  
  message("Reading ", length(files), " files for ", result_group)
  
  purrr::map_dfr(files, function(f) {
    
    x <- readr::read_tsv(f, show_col_types = FALSE)
    
    if (nrow(x) == 0) {
      return(tibble())
    }
    
    x %>%
      mutate(
        source_file = basename(f),
        result_group = result_group
      )
  })
}

# ============================================================
# Read ordinary Mantel result files
# ============================================================

community_env_results <- read_result_files(
  pattern = "^community_env_.*\\.tsv$",
  result_group = "ordinary_community_environment"
)

community_community_results <- read_result_files(
  pattern = "^community_pair_.*\\.tsv$",
  result_group = "ordinary_community_community"
)

vcr_env_results <- read_result_files(
  pattern = "^vcr_environment_mantel\\.tsv$",
  result_group = "ordinary_vcr_environment"
)

vcr_community_results <- read_result_files(
  pattern = "^vcr_community_mantel\\.tsv$",
  result_group = "ordinary_vcr_community"
)

community_community_fraction_results <- read_result_files(
  pattern = "^fraction_mantel_.*\\.tsv$",
  result_group = "ordinary_community_community_fraction"
)

# ============================================================
# Read partial Mantel result files
# ============================================================

partial_results <- read_result_files(
  pattern = "^partial_mantel_.*\\.tsv$",
  result_group = "partial_mantel"
)

# ============================================================
# Combine all results
# ============================================================

mantel_results_all <- bind_rows(
  community_env_results,
  community_community_results,
  vcr_env_results,
  vcr_community_results,
  community_community_fraction_results,
  partial_results
)

if (nrow(mantel_results_all) == 0) {
  stop("No Mantel result rows found.")
}

# ============================================================
# Harmonize columns
# ============================================================

mantel_results_all <- mantel_results_all %>%
  mutate(
    dataset1 = as.character(dataset1),
    dataset2 = as.character(dataset2),
    comparison_type = as.character(comparison_type),
    mantel_r = as.numeric(mantel_r),
    p_value = as.numeric(p_value),
    n_samples = as.integer(n_samples),
    permutations = as.integer(permutations),
    method = as.character(method),
    
    partial_control = if ("partial_control" %in% colnames(.)) {
      as.character(partial_control)
    } else {
      NA_character_
    },
    
    host_node = if ("host_node" %in% colnames(.)) {
      as.character(host_node)
    } else {
      NA_character_
    },
    
    env_var = if ("env_var" %in% colnames(.)) {
      as.character(env_var)
    } else {
      NA_character_
    },
    
    mantel_type = if_else(
      str_detect(comparison_type, "partial") |
        str_detect(result_group, "partial"),
      "partial_mantel",
      "ordinary_mantel"
    )
  ) %>%
  filter(
    !is.na(dataset1),
    !is.na(dataset2),
    !is.na(mantel_r),
    !is.na(p_value)
  )

# ============================================================
# BH adjustment
# ============================================================
# I recommend one global BH correction across the combined network.
# If you want ordinary and partial adjusted separately, see below.

mantel_results_all_adj <- mantel_results_all %>%
  mutate(
    p_adj_BH_global = p.adjust(p_value, method = "BH"),
    significant_BH_0.05_global = p_adj_BH_global < 0.05,
    abs_mantel_r = abs(mantel_r),
    source = dataset1,
    target = dataset2,
    weight = abs_mantel_r,
    edge_sign = case_when(
      mantel_r > 0 ~ "positive",
      mantel_r < 0 ~ "negative",
      TRUE ~ "zero"
    )
  ) %>%
  group_by(mantel_type) %>%
  mutate(
    p_adj_BH_by_mantel_type = p.adjust(p_value, method = "BH"),
    significant_BH_0.05_by_mantel_type = p_adj_BH_by_mantel_type < 0.05
  ) %>%
  ungroup() %>%
  arrange(p_adj_BH_global, desc(abs_mantel_r))

# ============================================================
# Cytoscape edge table
# ============================================================

cytoscape_edges <- mantel_results_all_adj %>%
  transmute(
    source,
    target,
    interaction = mantel_type,
    comparison_type,
    mantel_type,
    partial_control,
    host_node,
    env_var,
    mantel_r,
    abs_mantel_r,
    weight,
    edge_sign,
    p_value,
    p_adj_BH_global,
    significant_BH_0.05_global,
    p_adj_BH_by_mantel_type,
    significant_BH_0.05_by_mantel_type,
    n_samples,
    permutations,
    method,
    result_group,
    source_file
  )

# ============================================================
# Save outputs
# ============================================================

all_results_file <- file.path(
  out_dir,
  "mantel_results_all_ordinary_plus_partial_adjusted_v4.tsv"
)

edges_file <- file.path(
  out_dir,
  "mantel_edges_cytoscape_ordinary_plus_partial_v4.tsv"
)

readr::write_tsv(
  mantel_results_all_adj,
  all_results_file
)

readr::write_tsv(
  cytoscape_edges,
  edges_file
)

message("Done.")
message("Wrote full results: ", all_results_file)
message("Wrote Cytoscape edges: ", edges_file)

message("Result counts:")
print(mantel_results_all_adj %>% count(mantel_type, comparison_type))