# ============================================================
# 03_run_community_community_by_pair.R
#
# Usage:
#   Rscript 03_run_community_community_by_pair.R 1 virus-v-cellular 999
# ============================================================

library(dplyr)
library(tibble)
library(readr)
library(vegan)

base_dir <- "/path/to/base/dir"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide pair_id, e.g. Rscript 03_run_community_community_by_pair.R 1 virus-v-cellular 999")
}

pair_id <- as.integer(args[1])
run_label <- ifelse(length(args) >= 2, args[2], "virus-v-cellular")
permutations <- ifelse(length(args) >= 3, as.integer(args[3]), 999)

if (is.na(pair_id)) stop("pair_id must be an integer.")
if (is.na(permutations)) stop("permutations must be an integer.")

out_dir <- file.path(base_dir, paste0(run_label, "_mantel"))

input_file <- file.path(out_dir, "mantel_inputs_filtered.rds")
pair_table_file <- file.path(out_dir, "community_pair_table.tsv")

if (!file.exists(input_file)) {
  stop("Input RDS does not exist: ", input_file)
}

inputs <- readRDS(input_file)

community_nodes <- inputs$community_nodes
env_by_node <- inputs$env_by_node
match_cols <- inputs$match_cols
vir_category <- inputs$vir_category

# ============================================================
# Define source universe for each community node
# ============================================================

viral_node_source <- if (!is.null(vir_category) && vir_category == "cellular") {
  "cellular_fraction"
} else {
  "viral_fraction"
}

community_source <- c(
  Caudoviricetes = viral_node_source,
  Megaviricetes = viral_node_source,
  All_viruses_known_tax = viral_node_source,
  Prokaryotes = "cellular_fraction",
  Eukaryotes = "cellular_fraction"
)

community_source <- community_source[names(community_nodes)]

if (any(is.na(community_source))) {
  stop(
    "Some community nodes are missing source labels: ",
    paste(names(community_source)[is.na(community_source)], collapse = ", ")
  )
}

# ============================================================
# Create pair table if it does not exist
# ============================================================

if (!file.exists(pair_table_file)) {
  community_pairs <- combn(names(community_nodes), 2, simplify = FALSE)
  
  pair_table <- tibble::tibble(
    pair_id = seq_along(community_pairs),
    node1 = sapply(community_pairs, `[`, 1),
    node2 = sapply(community_pairs, `[`, 2)
  )
  
  readr::write_tsv(pair_table, pair_table_file)
} else {
  pair_table <- readr::read_tsv(pair_table_file, show_col_types = FALSE)
}

pair_info <- pair_table %>%
  filter(pair_id == !!pair_id)

if (nrow(pair_info) != 1) {
  stop(
    "Invalid pair_id: ", pair_id,
    ". Valid pair IDs: ",
    paste(pair_table$pair_id, collapse = ", ")
  )
}

name1 <- pair_info$node1
name2 <- pair_info$node2

source1 <- unname(community_source[[name1]])
source2 <- unname(community_source[[name2]])

message("Running community-community Mantel")
message("  pair_id: ", pair_id)
message("  pair: ", name1, " vs ", name2)
message("  source1: ", source1)
message("  source2: ", source2)
message("  run_label: ", run_label)
message("  permutations: ", permutations)
message("  out_dir: ", out_dir)

# ============================================================
# Function: make matched Run pairs
# ============================================================

make_community_pair_matches <- function(
    mat1,
    env1,
    mat2,
    env2,
    source1,
    source2,
    match_cols = c("manual_longitude", "manual_depth", "manual_latitude", "manual_date")
) {
  
  same_source <- identical(source1, source2)
  
  if (same_source) {
    
    common_runs <- intersect(rownames(mat1), rownames(mat2))
    
    matched <- tibble::tibble(
      Run_1 = common_runs,
      Run_2 = common_runs,
      matching_basis = "Run"
    )
    
  } else {
    
    map1 <- env1 %>%
      select(Run, all_of(match_cols)) %>%
      filter(if_all(all_of(match_cols), ~ !is.na(.x))) %>%
      group_by(Run) %>%
      summarise(
        across(all_of(match_cols), ~ dplyr::first(.x)),
        .groups = "drop"
      )
    
    map2 <- env2 %>%
      select(Run, all_of(match_cols)) %>%
      filter(if_all(all_of(match_cols), ~ !is.na(.x))) %>%
      group_by(Run) %>%
      summarise(
        across(all_of(match_cols), ~ dplyr::first(.x)),
        .groups = "drop"
      )
    
    matched <- map1 %>%
      inner_join(
        map2,
        by = match_cols,
        suffix = c("_1", "_2"),
        relationship = "many-to-many"
      ) %>%
      filter(
        Run_1 %in% rownames(mat1),
        Run_2 %in% rownames(mat2)
      ) %>%
      distinct(Run_1, Run_2, .keep_all = TRUE) %>%
      mutate(matching_basis = "metadata_match_cols")
  }
  
  matched
}

# ============================================================
# Function: one community-community Mantel pair
# ============================================================

run_single_community_pair_mantel <- function(
    mat1,
    env1,
    mat2,
    env2,
    name1,
    name2,
    source1,
    source2,
    match_cols = c("manual_longitude", "manual_depth", "manual_latitude", "manual_date"),
    method = "spearman",
    permutations = 99,
    seed = 123
) {
  
  set.seed(seed)
  
  matched <- make_community_pair_matches(
    mat1 = mat1,
    env1 = env1,
    mat2 = mat2,
    env2 = env2,
    source1 = source1,
    source2 = source2,
    match_cols = match_cols
  )
  
  if (nrow(matched) < 3) {
    message("  Skipping: fewer than 3 matched samples")
    return(NULL)
  }
  
  mat1_sub <- mat1[matched$Run_1, , drop = FALSE]
  mat2_sub <- mat2[matched$Run_2, , drop = FALSE]
  
  # Remove empty features
  mat1_sub <- mat1_sub[, colSums(mat1_sub > 0, na.rm = TRUE) > 0, drop = FALSE]
  mat2_sub <- mat2_sub[, colSums(mat2_sub > 0, na.rm = TRUE) > 0, drop = FALSE]
  
  # Remove matched pairs where either side is empty
  keep_pairs <- rowSums(mat1_sub > 0, na.rm = TRUE) > 0 &
    rowSums(mat2_sub > 0, na.rm = TRUE) > 0
  
  matched <- matched[keep_pairs, , drop = FALSE]
  mat1_sub <- mat1_sub[keep_pairs, , drop = FALSE]
  mat2_sub <- mat2_sub[keep_pairs, , drop = FALSE]
  
  if (nrow(matched) < 3 || ncol(mat1_sub) < 1 || ncol(mat2_sub) < 1) {
    message("  Skipping: fewer than 3 non-empty samples or no features")
    return(NULL)
  }
  
  n_samples <- nrow(matched)
  unique_runs_1 <- dplyr::n_distinct(matched$Run_1)
  unique_runs_2 <- dplyr::n_distinct(matched$Run_2)
  
  n_unique_match_keys <- if (all(match_cols %in% colnames(matched))) {
    matched %>%
      distinct(across(all_of(match_cols))) %>%
      nrow()
  } else {
    n_samples
  }
  
  duplicated_runs_1 <- n_samples - unique_runs_1
  duplicated_runs_2 <- n_samples - unique_runs_2
  
  matching_basis <- unique(matched$matching_basis)
  
  matching_type <- case_when(
    matching_basis == "Run" ~ "same_Run",
    matching_basis == "metadata_match_cols" &&
      duplicated_runs_1 == 0 &&
      duplicated_runs_2 == 0 ~ "one_to_one_metadata",
    matching_basis == "metadata_match_cols" &&
      (duplicated_runs_1 > 0 || duplicated_runs_2 > 0) ~ "many_to_many_metadata",
    TRUE ~ "unknown"
  )
  
  rownames(mat1_sub) <- seq_len(nrow(mat1_sub))
  rownames(mat2_sub) <- seq_len(nrow(mat2_sub))
  
  message(
    "  Matching: ", matching_basis,
    " / ", matching_type,
    " n=", n_samples,
    " unique_runs_1=", unique_runs_1,
    " unique_runs_2=", unique_runs_2
  )
  
  message(
    "  Calculating Bray-Curtis: ",
    name1, " features=", ncol(mat1_sub),
    "; ", name2, " features=", ncol(mat2_sub),
    "; n=", nrow(mat1_sub)
  )
  
  dist1 <- vegan::vegdist(mat1_sub, method = "bray")
  dist2 <- vegan::vegdist(mat2_sub, method = "bray")
  
  if (any(is.na(dist1)) || any(is.na(dist2))) {
    message("  Skipping: NA distances")
    return(NULL)
  }
  
  dist1 <- as.dist(as.matrix(dist1))
  dist2 <- as.dist(as.matrix(dist2))
  
  if (attr(dist1, "Size") != attr(dist2, "Size")) {
    stop(
      "Distance Size mismatch for ",
      name1, " vs ", name2,
      ": dist1=", attr(dist1, "Size"),
      ", dist2=", attr(dist2, "Size")
    )
  }
  
  if (length(dist1) != length(dist2)) {
    stop(
      "Distance length mismatch for ",
      name1, " vs ", name2,
      ": dist1=", length(dist1),
      ", dist2=", length(dist2)
    )
  }
  
  message(
    "  Mantel: ",
    name1, " vs ", name2,
    " n=", attr(dist1, "Size"),
    " length=", length(dist1)
  )
  
  res <- vegan::mantel(
    dist1,
    dist2,
    method = method,
    permutations = permutations
  )
  
  tibble::tibble(
    dataset1 = name1,
    dataset2 = name2,
    comparison_type = "community_community",
    source1 = source1,
    source2 = source2,
    matching_basis = matching_basis,
    matching_type = matching_type,
    mantel_r = unname(res$statistic),
    p_value = res$signif,
    n_samples = attr(dist1, "Size"),
    n_unique_match_keys = n_unique_match_keys,
    unique_runs_1 = unique_runs_1,
    unique_runs_2 = unique_runs_2,
    duplicated_runs_1 = duplicated_runs_1,
    duplicated_runs_2 = duplicated_runs_2,
    pair_inflation_factor = n_samples / n_unique_match_keys,
    permutations = permutations,
    method = method
  )
}

# ============================================================
# Run selected pair
# ============================================================

result <- run_single_community_pair_mantel(
  mat1 = community_nodes[[name1]],
  env1 = env_by_node[[name1]],
  mat2 = community_nodes[[name2]],
  env2 = env_by_node[[name2]],
  name1 = name1,
  name2 = name2,
  source1 = source1,
  source2 = source2,
  match_cols = match_cols,
  permutations = permutations
)

# ============================================================
# Save result
# ============================================================

safe_name1 <- gsub("[^A-Za-z0-9]+", "_", name1)
safe_name2 <- gsub("[^A-Za-z0-9]+", "_", name2)

outfile <- file.path(
  out_dir,
  paste0(
    "community_pair_",
    pair_id,
    "_",
    safe_name1,
    "__",
    safe_name2,
    ".tsv"
  )
)

if (is.null(result) || nrow(result) == 0) {
  result <- tibble::tibble(
    dataset1 = character(),
    dataset2 = character(),
    comparison_type = character(),
    source1 = character(),
    source2 = character(),
    matching_basis = character(),
    matching_type = character(),
    mantel_r = numeric(),
    p_value = numeric(),
    n_samples = integer(),
    n_unique_match_keys = integer(),
    unique_runs_1 = integer(),
    unique_runs_2 = integer(),
    duplicated_runs_1 = integer(),
    duplicated_runs_2 = integer(),
    pair_inflation_factor = numeric(),
    permutations = integer(),
    method = character()
  )
}

readr::write_tsv(result, outfile)

message("Done.")
message("Wrote: ", outfile)