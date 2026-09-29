# ============================================================
# 05_run_partial_mantel_by_pair.R
#
# Partial Mantel tests with three modes:
#
#   1. virus_env_partial_host
#      virus community ~ environmental variable | host community
#
#   2. host_env_partial_virus
#      host community ~ environmental variable | virus community
#
#   3. virus_host_partial_env
#      virus community ~ host community | environmental variable
#
# Usage:
#   Rscript 05_run_partial_mantel_by_pair.R 1 virus-v-cellular 999
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
  stop(
    "Please provide pair_id, e.g. ",
    "Rscript 05_run_partial_mantel_by_pair.R 1 virus-v-cellular 999"
  )
}

pair_id <- as.integer(args[1])
run_label <- ifelse(
  length(args) >= 2,
  args[2],
  "virus-v-cellular"
)
permutations <- ifelse(
  length(args) >= 3,
  as.integer(args[3]),
  999
)

if (is.na(pair_id)) {
  stop("pair_id must be an integer.")
}

if (is.na(permutations)) {
  stop("permutations must be an integer.")
}

method <- "spearman"

out_dir <- file.path(
  base_dir,
  paste0(run_label, "_mantel")
)

output_dir <- out_dir

input_file <- file.path(
  out_dir,
  "mantel_inputs_filtered.rds"
)

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

cat("============================================================\n")
cat("Running partial Mantel test\n")
cat("pair_id: ", pair_id, "\n", sep = "")
cat("run_label: ", run_label, "\n", sep = "")
cat("permutations: ", permutations, "\n", sep = "")
cat("out_dir: ", out_dir, "\n", sep = "")
cat("input_file: ", input_file, "\n", sep = "")
cat("============================================================\n")

# ============================================================
# 1. Load prepared inputs
# ============================================================

inputs <- readRDS(input_file)

community_nodes <- inputs$community_nodes
env_by_node <- inputs$env_by_node
match_cols <- inputs$match_cols

cat("Loaded community nodes:\n")
print(names(community_nodes))

# ============================================================
# 2. Define partial Mantel jobs
# ============================================================

partial_pair_table <- tibble::tibble(
  pair_id = 1:16,

  first_node = c(
    "Caudoviricetes", # 1
    "Caudoviricetes", # 2
    "Megaviricetes",  # 3
    "Megaviricetes",  # 4
    "Caudoviricetes", # 5
    "Prokaryotes",    # 6
    "Megaviricetes",  # 7
    "Eukaryotes",     # 8
    "Caudoviricetes", # 9
    "Megaviricetes",  # 10
    "Caudoviricetes", # 11
    "Prokaryotes",    # 12
    "Eukaryotes",     # 13
    "Prokaryotes",    # 14
    "Eukaryotes",     # 15
    "Megaviricetes"   # 16
  ),

  second_node = c(
    "Prokaryotes",    # 1
    "Prokaryotes",    # 2
    "Eukaryotes",     # 3
    "Eukaryotes",     # 4
    "Prokaryotes",    # 5
    "Caudoviricetes", # 6
    "Eukaryotes",     # 7
    "Megaviricetes",  # 8
    "Prokaryotes",    # 9
    "Eukaryotes",     # 10
    "Prokaryotes",    # 11
    "Caudoviricetes", # 12
    "Megaviricetes",  # 13
    "Caudoviricetes", # 14
    "Megaviricetes",  # 15
    "Eukaryotes"      # 16
  ),

  env_var = c(
    "env_temp",         # 1
    "env_abs_latitude", # 2
    "env_temp",         # 3
    "env_abs_latitude", # 4
    "env_PAR",          # 5
    "env_PAR",          # 6
    "env_PAR",          # 7
    "env_PAR",          # 8
    "env_temp",         # 9
    "env_temp",         # 10
    "env_O2",           # 11
    "env_temp",         # 12
    "env_temp",         # 13
    "env_abs_latitude", # 14
    "env_abs_latitude", # 15
    "env_O2"            # 16
  ),

  partial_mode = c(
    "virus_env_partial_host",  # 1
    "virus_env_partial_host",  # 2
    "virus_env_partial_host",  # 3
    "virus_env_partial_host",  # 4
    "virus_env_partial_host",  # 5
    "host_env_partial_virus",  # 6
    "virus_env_partial_host",  # 7
    "host_env_partial_virus",  # 8
    "virus_host_partial_env",  # 9
    "virus_host_partial_env",  # 10
    "virus_env_partial_host",  # 11
    "host_env_partial_virus",  # 12
    "host_env_partial_virus",  # 13
    "host_env_partial_virus",  # 14
    "host_env_partial_virus",  # 15
    "virus_env_partial_host"   # 16
  )
)

pair_info <- partial_pair_table %>%
  filter(pair_id == !!pair_id)

if (nrow(pair_info) != 1) {
  stop(
    "Invalid pair_id: ",
    pair_id,
    ". Valid pair_id values are: ",
    paste(partial_pair_table$pair_id, collapse = ", ")
  )
}

first_node <- pair_info$first_node
second_node <- pair_info$second_node
env_var <- pair_info$env_var
partial_mode <- pair_info$partial_mode

valid_modes <- c(
  "virus_env_partial_host",
  "host_env_partial_virus",
  "virus_host_partial_env"
)

if (!partial_mode %in% valid_modes) {
  stop("Unknown partial_mode: ", partial_mode)
}

if (!first_node %in% names(community_nodes)) {
  stop("first_node not found in community_nodes: ", first_node)
}

if (!second_node %in% names(community_nodes)) {
  stop("second_node not found in community_nodes: ", second_node)
}

if (!first_node %in% names(env_by_node)) {
  stop("first_node not found in env_by_node: ", first_node)
}

if (!second_node %in% names(env_by_node)) {
  stop("second_node not found in env_by_node: ", second_node)
}

cat("Selected test:\n")

if (partial_mode == "virus_host_partial_env") {
  cat(
    "  ",
    first_node,
    " ~ ",
    second_node,
    " | ",
    env_var,
    "\n",
    sep = ""
  )
} else {
  cat(
    "  ",
    first_node,
    " ~ ",
    env_var,
    " | ",
    second_node,
    "\n",
    sep = ""
  )
}

cat("  mode: ", partial_mode, "\n", sep = "")

# ============================================================
# 3. Function to run one partial Mantel test
# ============================================================

run_community_env_partial_mantel_control <- function(
    first_mat,
    first_env,
    second_mat,
    second_env,
    first_name,
    second_name,
    env_var,
    partial_mode,
    match_cols = c(
      "manual_longitude",
      "manual_depth",
      "manual_latitude",
      "manual_date"
    ),
    method = "spearman",
    permutations = 9999,
    seed = 123
) {

  set.seed(seed)

  if (!env_var %in% colnames(first_env)) {
    stop(
      "env_var not found in metadata for ",
      first_name,
      ": ",
      env_var
    )
  }

  # First community:
  # Run + matching metadata + environmental variable
  first_map <- first_env %>%
    select(
      Run,
      all_of(match_cols),
      env_value = all_of(env_var)
    ) %>%
    filter(
      if_all(all_of(match_cols), ~ !is.na(.x)),
      !is.na(env_value),
      is.finite(env_value)
    ) %>%
    group_by(Run) %>%
    summarise(
      across(
        all_of(match_cols),
        ~ dplyr::first(.x)
      ),
      env_value = mean(env_value, na.rm = TRUE),
      .groups = "drop"
    )

  # Second community:
  # Run + matching metadata
  second_map <- second_env %>%
    select(
      Run,
      all_of(match_cols)
    ) %>%
    filter(
      if_all(all_of(match_cols), ~ !is.na(.x))
    ) %>%
    group_by(Run) %>%
    summarise(
      across(
        all_of(match_cols),
        ~ dplyr::first(.x)
      ),
      .groups = "drop"
    )

  # Match samples from the two communities by metadata
  matched <- dplyr::inner_join(
    first_map,
    second_map,
    by = match_cols,
    suffix = c("_first", "_second"),
    relationship = "many-to-many"
  )

  matched <- matched %>%
    filter(
      Run_first %in% rownames(first_mat),
      Run_second %in% rownames(second_mat)
    ) %>%
    distinct(
      Run_first,
      Run_second,
      .keep_all = TRUE
    )

  if (nrow(matched) < 3) {
    message("  Skipping: fewer than 3 matched samples")
    return(NULL)
  }

  first_sub <- first_mat[
    matched$Run_first,
    ,
    drop = FALSE
  ]

  second_sub <- second_mat[
    matched$Run_second,
    ,
    drop = FALSE
  ]

  env_vec <- as.numeric(matched$env_value)

  # Remove features that are absent in all matched samples
  first_sub <- first_sub[
    ,
    colSums(first_sub > 0, na.rm = TRUE) > 0,
    drop = FALSE
  ]

  second_sub <- second_sub[
    ,
    colSums(second_sub > 0, na.rm = TRUE) > 0,
    drop = FALSE
  ]

  # Remove unusable matched pairs
  keep_pairs <-
    rowSums(first_sub > 0, na.rm = TRUE) > 0 &
    rowSums(second_sub > 0, na.rm = TRUE) > 0 &
    !is.na(env_vec) &
    is.finite(env_vec)

  matched <- matched[
    keep_pairs,
    ,
    drop = FALSE
  ]

  first_sub <- first_sub[
    keep_pairs,
    ,
    drop = FALSE
  ]

  second_sub <- second_sub[
    keep_pairs,
    ,
    drop = FALSE
  ]

  env_vec <- env_vec[keep_pairs]

  if (
    nrow(matched) < 3 ||
    ncol(first_sub) < 1 ||
    ncol(second_sub) < 1
  ) {
    message(
      "  Skipping: fewer than 3 non-empty matched samples ",
      "or no remaining features"
    )
    return(NULL)
  }

  if (length(unique(env_vec)) < 2) {
    message(
      "  Skipping: environmental variable has no variation"
    )
    return(NULL)
  }

  # Force identical anonymous sample order
  rownames(first_sub) <- seq_len(nrow(first_sub))
  rownames(second_sub) <- seq_len(nrow(second_sub))
  names(env_vec) <- seq_len(length(env_vec))

  message(
    "  Calculating distances: ",
    first_name,
    " and ",
    second_name,
    " with ",
    env_var,
    "; n=",
    nrow(first_sub),
    "; first_features=",
    ncol(first_sub),
    "; second_features=",
    ncol(second_sub)
  )

  first_dist <- vegan::vegdist(
    first_sub,
    method = "bray"
  )

  second_dist <- vegan::vegdist(
    second_sub,
    method = "bray"
  )

  env_dist <- stats::dist(
    scale(env_vec),
    method = "euclidean"
  )

  if (
    any(is.na(first_dist)) ||
    any(is.na(second_dist)) ||
    any(is.na(env_dist))
  ) {
    message("  Skipping: NA distances produced")
    return(NULL)
  }

  first_dist <- as.dist(as.matrix(first_dist))
  second_dist <- as.dist(as.matrix(second_dist))
  env_dist <- as.dist(as.matrix(env_dist))

  if (
    attr(first_dist, "Size") != attr(second_dist, "Size") ||
    attr(first_dist, "Size") != attr(env_dist, "Size")
  ) {
    stop(
      "Distance Size mismatch: ",
      "first_size=",
      attr(first_dist, "Size"),
      ", second_size=",
      attr(second_dist, "Size"),
      ", env_size=",
      attr(env_dist, "Size")
    )
  }

  if (
    length(first_dist) != length(second_dist) ||
    length(first_dist) != length(env_dist)
  ) {
    stop(
      "Distance length mismatch: ",
      "first_dist=",
      length(first_dist),
      ", second_dist=",
      length(second_dist),
      ", env_dist=",
      length(env_dist)
    )
  }

  message(
    "  Running partial Mantel test: ",
    partial_mode,
    "; n=",
    attr(first_dist, "Size"),
    "; distance_length=",
    length(first_dist)
  )

  if (partial_mode == "virus_env_partial_host") {

    # Virus community ~ environment | host community
    res <- vegan::mantel.partial(
      first_dist,
      env_dist,
      second_dist,
      method = method,
      permutations = permutations
    )

    out <- tibble::tibble(
      dataset1 = first_name,
      dataset2 = env_var,
      comparison_type = "virus_environment_partial_host",
      partial_control = second_name,
      host_node = second_name,
      env_var = env_var,
      mantel_r = unname(res$statistic),
      p_value = res$signif,
      n_samples = attr(first_dist, "Size"),
      permutations = permutations,
      method = method
    )

  } else if (partial_mode == "host_env_partial_virus") {

    # Host community ~ environment | virus community
    res <- vegan::mantel.partial(
      first_dist,
      env_dist,
      second_dist,
      method = method,
      permutations = permutations
    )

    out <- tibble::tibble(
      dataset1 = first_name,
      dataset2 = env_var,
      comparison_type = "host_environment_partial_virus",
      partial_control = second_name,
      host_node = first_name,
      env_var = env_var,
      mantel_r = unname(res$statistic),
      p_value = res$signif,
      n_samples = attr(first_dist, "Size"),
      permutations = permutations,
      method = method
    )

  } else if (partial_mode == "virus_host_partial_env") {

    # Virus community ~ host community | environment
    res <- vegan::mantel.partial(
      first_dist,
      second_dist,
      env_dist,
      method = method,
      permutations = permutations
    )

    out <- tibble::tibble(
      dataset1 = first_name,
      dataset2 = second_name,
      comparison_type = "virus_host_partial_environment",
      partial_control = env_var,
      host_node = second_name,
      env_var = env_var,
      mantel_r = unname(res$statistic),
      p_value = res$signif,
      n_samples = attr(first_dist, "Size"),
      permutations = permutations,
      method = method
    )

  } else {
    stop("Unknown partial_mode: ", partial_mode)
  }

  out
}

# ============================================================
# 4. Run selected partial Mantel test
# ============================================================

partial_result <- run_community_env_partial_mantel_control(
  first_mat = community_nodes[[first_node]],
  first_env = env_by_node[[first_node]],
  second_mat = community_nodes[[second_node]],
  second_env = env_by_node[[second_node]],
  first_name = first_node,
  second_name = second_node,
  env_var = env_var,
  partial_mode = partial_mode,
  match_cols = match_cols,
  method = method,
  permutations = permutations
)

# ============================================================
# 5. Save result
# ============================================================

safe_first <- gsub(
  "[^A-Za-z0-9]+",
  "_",
  first_node
)

safe_second <- gsub(
  "[^A-Za-z0-9]+",
  "_",
  second_node
)

safe_env <- gsub(
  "[^A-Za-z0-9]+",
  "_",
  env_var
)

safe_mode <- gsub(
  "[^A-Za-z0-9]+",
  "_",
  partial_mode
)

outfile <- file.path(
  output_dir,
  paste0(
    "partial_mantel_",
    pair_id,
    "_",
    safe_first,
    "__",
    safe_second,
    "__",
    safe_env,
    "__",
    safe_mode,
    ".tsv"
  )
)

if (
  is.null(partial_result) ||
  nrow(partial_result) == 0
) {

  message(
    "No result produced. Writing empty placeholder: ",
    outfile
  )

  partial_result <- tibble::tibble(
    dataset1 = character(),
    dataset2 = character(),
    comparison_type = character(),
    partial_control = character(),
    host_node = character(),
    env_var = character(),
    mantel_r = numeric(),
    p_value = numeric(),
    n_samples = integer(),
    permutations = integer(),
    method = character()
  )
}

readr::write_tsv(
  partial_result,
  outfile
)

message("Done.")
message("Wrote: ", outfile)