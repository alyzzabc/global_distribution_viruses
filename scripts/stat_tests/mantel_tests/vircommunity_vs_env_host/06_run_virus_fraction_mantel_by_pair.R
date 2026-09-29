# ============================================================
# 07_run_virus_fraction_mantel_by_pair.R
#
# Ordinary Mantel tests comparing viral-fraction viruses
# against cellular-fraction viruses.
#
# Tests:
#   1. All_viruses_known_tax viral fraction ~ cellular fraction
#   2. Caudoviricetes viral fraction ~ cellular fraction
#   3. Megaviricetes viral fraction ~ cellular fraction
#
# Usage:
#   Rscript 07_run_virus_fraction_mantel_by_pair.R 999
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(vegan)
})

# ============================================================
# 0. Command-line arguments and paths
# ============================================================

base_dir <- "/path/to/base/dir"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide pair_id, e.g. Rscript 06_run_virus_fraction_mantel_by_pair.R 1")
}

pair_id <- as.integer(args[1])
permutations <- ifelse(length(args) >= 2, as.integer(args[2]), 999)

if (is.na(pair_id)) {
  stop("pair_id must be an integer.")
}

if (is.na(permutations)) {
  stop("permutations must be an integer.")
}

method <- "spearman"

# Hard-coded comparison:
# viral-fraction viruses vs cellular-fraction viruses
viral_dir <- file.path(base_dir, "virus-v-cellular_mantel")
cellular_dir <- file.path(base_dir, "cellular-v-cellular_mantel")

viral_input_file <- file.path(viral_dir, "mantel_inputs_filtered.rds")
cellular_input_file <- file.path(cellular_dir, "mantel_inputs_filtered.rds")

# Put outputs in virus-v-cellular_mantel
output_dir <- viral_dir

if (!file.exists(viral_input_file)) {
  stop("Viral-fraction input file does not exist: ", viral_input_file)
}

if (!file.exists(cellular_input_file)) {
  stop("Cellular-fraction input file does not exist: ", cellular_input_file)
}

cat("============================================================\n")
cat("Running virus-fraction vs cellular-fraction Mantel test\n")
cat("pair_id: ", pair_id, "\n", sep = "")
cat("permutations: ", permutations, "\n", sep = "")
cat("viral_input_file: ", viral_input_file, "\n", sep = "")
cat("cellular_input_file: ", cellular_input_file, "\n", sep = "")
cat("output_dir: ", output_dir, "\n", sep = "")
cat("============================================================\n")

# ============================================================
# 1. Load prepared inputs
# ============================================================

viral_inputs <- readRDS(viral_input_file)
cellular_inputs <- readRDS(cellular_input_file)

viral_community_nodes <- viral_inputs$community_nodes
cellular_community_nodes <- cellular_inputs$community_nodes

viral_env_by_node <- viral_inputs$env_by_node
cellular_env_by_node <- cellular_inputs$env_by_node

match_cols <- viral_inputs$match_cols

if (is.null(match_cols)) {
  match_cols <- c("manual_longitude", "manual_depth", "manual_latitude", "manual_date")
}

cat("Loaded viral-fraction community nodes:\n")
print(names(viral_community_nodes))

cat("Loaded cellular-fraction community nodes:\n")
print(names(cellular_community_nodes))

# ============================================================
# 2. Define tests
# ============================================================

fraction_pair_table <- tibble::tibble(
  pair_id = 1:3,
  node_name = c(
    "All_viruses_known_tax",
    "Caudoviricetes",
    "Megaviricetes"
  )
)

pair_info <- fraction_pair_table %>%
  filter(pair_id == !!pair_id)

if (nrow(pair_info) != 1) {
  stop(
    "Invalid pair_id: ", pair_id,
    ". Valid pair_id values are: ",
    paste(fraction_pair_table$pair_id, collapse = ", ")
  )
}

node_name <- pair_info$node_name

if (!node_name %in% names(viral_community_nodes)) {
  stop("node_name not found in viral_community_nodes: ", node_name)
}

if (!node_name %in% names(cellular_community_nodes)) {
  stop("node_name not found in cellular_community_nodes: ", node_name)
}

if (!node_name %in% names(viral_env_by_node)) {
  stop("node_name not found in viral_env_by_node: ", node_name)
}

if (!node_name %in% names(cellular_env_by_node)) {
  stop("node_name not found in cellular_env_by_node: ", node_name)
}

cat("Selected test:\n")
cat("  ", node_name, "_viral_fraction ~ ", node_name, "_cellular_fraction\n", sep = "")

# ============================================================
# 3. Function: viral-fraction community ~ cellular-fraction community
# ============================================================

run_fraction_community_mantel <- function(
    viral_mat,
    viral_env,
    cellular_mat,
    cellular_env,
    node_name,
    match_cols = c("manual_longitude", "manual_depth", "manual_latitude", "manual_date"),
    method = "spearman",
    permutations = 999,
    seed = 123
) {
  
  set.seed(seed)
  
  viral_map <- viral_env %>%
    select(Run, all_of(match_cols)) %>%
    filter(if_all(all_of(match_cols), ~ !is.na(.x))) %>%
    group_by(Run) %>%
    summarise(
      across(all_of(match_cols), ~ dplyr::first(.x)),
      .groups = "drop"
    )
  
  cellular_map <- cellular_env %>%
    select(Run, all_of(match_cols)) %>%
    filter(if_all(all_of(match_cols), ~ !is.na(.x))) %>%
    group_by(Run) %>%
    summarise(
      across(all_of(match_cols), ~ dplyr::first(.x)),
      .groups = "drop"
    )
  
  matched <- viral_map %>%
    inner_join(
      cellular_map,
      by = match_cols,
      suffix = c("_viral_fraction", "_cellular_fraction"),
      relationship = "many-to-many"
    ) %>%
    filter(
      Run_viral_fraction %in% rownames(viral_mat),
      Run_cellular_fraction %in% rownames(cellular_mat)
    ) %>%
    distinct(Run_viral_fraction, Run_cellular_fraction, .keep_all = TRUE)
  
  if (nrow(matched) < 3) {
    message("  Skipping: fewer than 3 matched samples")
    return(NULL)
  }
  
  viral_sub <- viral_mat[matched$Run_viral_fraction, , drop = FALSE]
  cellular_sub <- cellular_mat[matched$Run_cellular_fraction, , drop = FALSE]
  
  # Remove empty features
  viral_sub <- viral_sub[, colSums(viral_sub > 0, na.rm = TRUE) > 0, drop = FALSE]
  cellular_sub <- cellular_sub[, colSums(cellular_sub > 0, na.rm = TRUE) > 0, drop = FALSE]
  
  # Remove matched pairs where either side is empty
  keep_pairs <- rowSums(viral_sub > 0, na.rm = TRUE) > 0 &
    rowSums(cellular_sub > 0, na.rm = TRUE) > 0
  
  matched <- matched[keep_pairs, , drop = FALSE]
  viral_sub <- viral_sub[keep_pairs, , drop = FALSE]
  cellular_sub <- cellular_sub[keep_pairs, , drop = FALSE]
  
  if (nrow(matched) < 3 || ncol(viral_sub) < 1 || ncol(cellular_sub) < 1) {
    message("  Skipping: fewer than 3 non-empty matched samples or no features")
    return(NULL)
  }
  
  rownames(viral_sub) <- seq_len(nrow(viral_sub))
  rownames(cellular_sub) <- seq_len(nrow(cellular_sub))
  
  message(
    "  Calculating distances: ",
    node_name, "_viral_fraction ~ ", node_name, "_cellular_fraction",
    "; n=", nrow(viral_sub),
    "; viral_features=", ncol(viral_sub),
    "; cellular_features=", ncol(cellular_sub)
  )
  
  viral_dist <- vegan::vegdist(viral_sub, method = "bray")
  cellular_dist <- vegan::vegdist(cellular_sub, method = "bray")
  
  if (any(is.na(viral_dist)) || any(is.na(cellular_dist))) {
    message("  Skipping: NA distances produced")
    return(NULL)
  }
  
  viral_dist <- as.dist(as.matrix(viral_dist))
  cellular_dist <- as.dist(as.matrix(cellular_dist))
  
  if (attr(viral_dist, "Size") != attr(cellular_dist, "Size")) {
    stop(
      "Distance Size mismatch: ",
      "viral_size=", attr(viral_dist, "Size"),
      ", cellular_size=", attr(cellular_dist, "Size")
    )
  }
  
  if (length(viral_dist) != length(cellular_dist)) {
    stop(
      "Distance length mismatch: ",
      "viral_dist=", length(viral_dist),
      ", cellular_dist=", length(cellular_dist)
    )
  }
  
  message(
    "  Mantel: ",
    node_name, "_viral_fraction ~ ", node_name, "_cellular_fraction",
    "; n=", attr(viral_dist, "Size"),
    "; length=", length(viral_dist)
  )
  
  res <- vegan::mantel(
    viral_dist,
    cellular_dist,
    method = method,
    permutations = permutations
  )
  
  tibble::tibble(
    dataset1 = paste0(node_name, "_viral_fraction"),
    dataset2 = paste0(node_name, "_cellular_fraction"),
    comparison_type = "viral_fraction_cellular_fraction_community",
    mantel_r = unname(res$statistic),
    p_value = res$signif,
    n_samples = attr(viral_dist, "Size"),
    permutations = permutations,
    method = method
  )
}

# ============================================================
# 4. Run selected Mantel test
# ============================================================

fraction_result <- run_fraction_community_mantel(
  viral_mat = viral_community_nodes[[node_name]],
  viral_env = viral_env_by_node[[node_name]],
  cellular_mat = cellular_community_nodes[[node_name]],
  cellular_env = cellular_env_by_node[[node_name]],
  node_name = node_name,
  match_cols = match_cols,
  method = method,
  permutations = permutations
)

# ============================================================
# 5. Save result
# ============================================================

safe_node <- gsub("[^A-Za-z0-9]+", "_", node_name)

outfile <- file.path(
  output_dir,
  paste0(
    "fraction_mantel_",
    pair_id,
    "_",
    safe_node,
    "_viral_fraction__cellular_fraction.tsv"
  )
)

if (is.null(fraction_result) || nrow(fraction_result) == 0) {
  
  message("No result produced. Writing empty placeholder: ", outfile)
  
  fraction_result <- tibble::tibble(
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

readr::write_tsv(fraction_result, outfile)

message("Done.")
message("Wrote: ", outfile)