# ============================================================
# 04_run_vcr_mantel.R
#
# Run VCR-environment and VCR-community Mantel tests.
#
# Usage:
#   Rscript 04_run_vcr_mantel.R virus-v-cellular 999
#   Rscript 04_run_vcr_mantel.R cellular-v-cellular 999
#
# Arguments:
#   arg1 = run_label
#   arg2 = permutations
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

run_label <- ifelse(length(args) >= 1, args[1], "virus-v-cellular")
permutations <- ifelse(length(args) >= 2, as.integer(args[2]), 999)

if (is.na(permutations)) {
  stop("permutations must be an integer.")
}

method <- "spearman"

out_dir <- file.path(base_dir, paste0(run_label, "_mantel"))
output_dir <- out_dir

input_file <- file.path(out_dir, "mantel_inputs_filtered.rds")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

cat("============================================================\n")
cat("Running VCR Mantel tests\n")
cat("run_label: ", run_label, "\n", sep = "")
cat("permutations: ", permutations, "\n", sep = "")
cat("out_dir: ", out_dir, "\n", sep = "")
cat("input_file: ", input_file, "\n", sep = "")
cat("============================================================\n")

# ============================================================
# 1. Load prepared inputs from 01 script
# ============================================================

inputs <- readRDS(input_file)

community_nodes <- inputs$community_nodes
env_by_node <- inputs$env_by_node
env_vars <- inputs$env_vars
vcr_nodes <- inputs$vcr_nodes
match_cols <- inputs$match_cols

env_vir_filtered  <- env_by_node$Caudoviricetes
env_prok_filtered <- env_by_node$Prokaryotes
env_euk_filtered  <- env_by_node$Eukaryotes

cat("Loaded inputs:\n")
cat("Community nodes: ", paste(names(community_nodes), collapse = ", "), "\n", sep = "")
cat("VCR nodes: ", paste(names(vcr_nodes), collapse = ", "), "\n", sep = "")
cat("Env vars: ", paste(env_vars, collapse = ", "), "\n", sep = "")

# ============================================================
# 2. Helper functions
# ============================================================

clean_abundance_mat <- function(mat) {
  mat <- mat[, colSums(mat > 0, na.rm = TRUE) > 0, drop = FALSE]
  mat <- mat[rowSums(mat > 0, na.rm = TRUE) > 0, , drop = FALSE]
  mat
}

make_bray_dist <- function(mat) {
  mat <- clean_abundance_mat(mat)
  
  if (nrow(mat) < 3 || ncol(mat) < 1) {
    return(NULL)
  }
  
  d <- vegan::vegdist(mat, method = "bray")
  
  if (any(is.na(d))) {
    return(NULL)
  }
  
  as.dist(as.matrix(d))
}

make_vector_dist <- function(x) {
  x <- as.numeric(x)
  x <- x[!is.na(x) & is.finite(x)]
  
  if (length(x) < 3) {
    return(NULL)
  }
  
  if (length(unique(x)) < 2) {
    return(NULL)
  }
  
  d <- stats::dist(scale(x), method = "euclidean")
  
  if (any(is.na(d))) {
    return(NULL)
  }
  
  as.dist(as.matrix(d))
}

# ============================================================
# 3. Mantel: VCR-environment
# ============================================================

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
      
      if (!env_var %in% colnames(env_table)) {
        next
      }
      
      tmp <- env_table %>%
        select(Run, env_value = all_of(env_var)) %>%
        filter(!is.na(env_value), is.finite(env_value)) %>%
        group_by(Run) %>%
        summarise(
          env_value = mean(env_value, na.rm = TRUE),
          .groups = "drop"
        )
      
      common_runs <- intersect(vcr_df$Run, tmp$Run)
      
      if (length(common_runs) < 3) {
        next
      }
      
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
      
      if (nrow(vcr_sub) < 3) {
        next
      }
      
      if (length(unique(vcr_sub$value)) < 2) {
        next
      }
      
      if (length(unique(env_sub$env_value)) < 2) {
        next
      }
      
      vcr_dist <- make_vector_dist(vcr_sub$value)
      env_dist <- make_vector_dist(env_sub$env_value)
      
      if (is.null(vcr_dist) || is.null(env_dist)) {
        next
      }
      
      if (attr(vcr_dist, "Size") != attr(env_dist, "Size")) {
        stop(
          "Distance Size mismatch for ",
          vcr_name, " vs ", env_var,
          ": vcr_size=", attr(vcr_dist, "Size"),
          ", env_size=", attr(env_dist, "Size")
        )
      }
      
      if (length(vcr_dist) != length(env_dist)) {
        stop(
          "Distance length mismatch for ",
          vcr_name, " vs ", env_var,
          ": vcr_dist=", length(vcr_dist),
          ", env_dist=", length(env_dist)
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
# 4. Mantel: VCR-community
# ============================================================

run_vcr_community_mantel <- function(
    vcr_df,
    community_mat,
    vcr_env,
    community_env,
    vcr_name,
    community_name,
    vcr_source,
    community_source,
    match_cols = c("manual_longitude", "manual_depth", "manual_latitude", "manual_date"),
    method = "spearman",
    permutations = 999,
    seed = 123
) {
  
  set.seed(seed)
  
  same_source <- identical(vcr_source, community_source)
  
  if (same_source) {
    
    common_runs <- intersect(vcr_df$Run, rownames(community_mat))
    
    matched <- tibble::tibble(
      Run_1 = common_runs,
      Run_2 = common_runs,
      matching_basis = "Run"
    )
    
  } else {
    
    map1 <- vcr_env %>%
      select(Run, all_of(match_cols)) %>%
      filter(if_all(all_of(match_cols), ~ !is.na(.x))) %>%
      group_by(Run) %>%
      summarise(
        across(all_of(match_cols), ~ dplyr::first(.x)),
        .groups = "drop"
      )
    
    map2 <- community_env %>%
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
        Run_1 %in% vcr_df$Run,
        Run_2 %in% rownames(community_mat)
      ) %>%
      distinct(Run_1, Run_2, .keep_all = TRUE) %>%
      mutate(matching_basis = "metadata_match_cols")
  }
  
  if (nrow(matched) < 3) {
    message("Skipping ", vcr_name, " vs ", community_name, ": fewer than 3 matched samples")
    return(NULL)
  }
  
  vcr_sub <- vcr_df %>%
    rename(Run_1 = Run) %>%
    inner_join(matched, by = "Run_1") %>%
    arrange(seq_len(n()))
  
  mat_sub <- community_mat[vcr_sub$Run_2, , drop = FALSE]
  
  mat_sub <- mat_sub[, colSums(mat_sub > 0, na.rm = TRUE) > 0, drop = FALSE]
  
  keep_samples <- rowSums(mat_sub > 0, na.rm = TRUE) > 0 &
    !is.na(vcr_sub$value) &
    is.finite(vcr_sub$value)
  
  mat_sub <- mat_sub[keep_samples, , drop = FALSE]
  vcr_sub <- vcr_sub[keep_samples, , drop = FALSE]
  
  if (nrow(mat_sub) < 3 || ncol(mat_sub) < 1) {
    message("Skipping ", vcr_name, " vs ", community_name, ": not enough non-empty samples/features")
    return(NULL)
  }
  
  if (length(unique(vcr_sub$value)) < 2) {
    message("Skipping ", vcr_name, " vs ", community_name, ": VCR has no variation")
    return(NULL)
  }
  
  n_samples <- nrow(vcr_sub)
  unique_runs_1 <- dplyr::n_distinct(vcr_sub$Run_1)
  unique_runs_2 <- dplyr::n_distinct(vcr_sub$Run_2)
  
  n_unique_match_keys <- if (all(match_cols %in% colnames(vcr_sub))) {
    vcr_sub %>%
      distinct(across(all_of(match_cols))) %>%
      nrow()
  } else {
    n_samples
  }
  
  duplicated_runs_1 <- n_samples - unique_runs_1
  duplicated_runs_2 <- n_samples - unique_runs_2
  
  matching_basis <- unique(vcr_sub$matching_basis)
  
  matching_type <- case_when(
    matching_basis == "Run" ~ "same_Run",
    matching_basis == "metadata_match_cols" &
      duplicated_runs_1 == 0 &
      duplicated_runs_2 == 0 ~ "one_to_one_metadata",
    matching_basis == "metadata_match_cols" &
      (duplicated_runs_1 > 0 | duplicated_runs_2 > 0) ~ "many_to_many_metadata",
    TRUE ~ "unknown"
  )
  
  rownames(mat_sub) <- seq_len(nrow(mat_sub))
  
  community_dist <- make_bray_dist(mat_sub)
  vcr_dist <- make_vector_dist(vcr_sub$value)
  
  if (is.null(community_dist) || is.null(vcr_dist)) {
    message("Skipping ", vcr_name, " vs ", community_name, ": distance creation failed")
    return(NULL)
  }
  
  if (attr(community_dist, "Size") != attr(vcr_dist, "Size")) {
    stop(
      "Distance Size mismatch for ",
      vcr_name, " vs ", community_name,
      ": community_size=", attr(community_dist, "Size"),
      ", vcr_size=", attr(vcr_dist, "Size")
    )
  }
  
  message(
    "Mantel: ", vcr_name, " vs ", community_name,
    " n=", attr(community_dist, "Size"),
    " matching=", matching_basis,
    " / ", matching_type
  )
  
  res <- vegan::mantel(
    community_dist,
    vcr_dist,
    method = method,
    permutations = permutations
  )
  
  tibble(
    dataset1 = vcr_name,
    dataset2 = community_name,
    comparison_type = "vcr_community",
    vcr_source = vcr_source,
    community_source = community_source,
    matching_basis = matching_basis,
    matching_type = matching_type,
    mantel_r = unname(res$statistic),
    p_value = res$signif,
    n_samples = attr(community_dist, "Size"),
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
# 5. Run VCR-environment tests
# ============================================================

message("Running Mantel tests: VCR ~ environment")

message("Overlap checks:")
message("  VCR_Caudo_Prok with Prok env: ", sum(vcr_nodes$VCR_Caudo_Prok$Run %in% env_prok_filtered$Run))
message("  VCR_Caudo_Prok with Viral env: ", sum(vcr_nodes$VCR_Caudo_Prok$Run %in% env_vir_filtered$Run))
message("  VCR_Mega_Euk with Euk env: ", sum(vcr_nodes$VCR_Mega_Euk$Run %in% env_euk_filtered$Run))
message("  VCR_Mega_Euk with Viral env: ", sum(vcr_nodes$VCR_Mega_Euk$Run %in% env_vir_filtered$Run))

vcr_env_results_caudo_prok <- run_vcr_env_mantel(
  vcr_nodes = list(
    VCR_Caudo_Prok = vcr_nodes$VCR_Caudo_Prok
  ),
  env_table = env_prok_filtered,
  env_vars = env_vars,
  method = method,
  permutations = permutations
)

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
# 6. Run VCR-community tests
# ============================================================

message("Running Mantel tests: VCR ~ community")

viral_node_source <- if (!is.null(inputs$vir_category) && inputs$vir_category == "cellular") {
  "cellular_fraction"
} else {
  "viral_fraction"
}

vcr_community_results <- bind_rows(
  run_vcr_community_mantel(
    vcr_df = vcr_nodes$VCR_Caudo_Prok,
    community_mat = community_nodes$Caudoviricetes,
    vcr_env = env_prok_filtered,
    community_env = env_by_node$Caudoviricetes,
    vcr_name = "VCR_Caudo_Prok",
    community_name = "Caudoviricetes",
    vcr_source = "cellular_fraction",
    community_source = viral_node_source,
    match_cols = match_cols,
    method = method,
    permutations = permutations
  ),
  
  run_vcr_community_mantel(
    vcr_df = vcr_nodes$VCR_Caudo_Prok,
    community_mat = community_nodes$Prokaryotes,
    vcr_env = env_prok_filtered,
    community_env = env_by_node$Prokaryotes,
    vcr_name = "VCR_Caudo_Prok",
    community_name = "Prokaryotes",
    vcr_source = "cellular_fraction",
    community_source = "cellular_fraction",
    match_cols = match_cols,
    method = method,
    permutations = permutations
  ),
  
  run_vcr_community_mantel(
    vcr_df = vcr_nodes$VCR_Mega_Euk,
    community_mat = community_nodes$Megaviricetes,
    vcr_env = env_euk_filtered,
    community_env = env_by_node$Megaviricetes,
    vcr_name = "VCR_Mega_Euk",
    community_name = "Megaviricetes",
    vcr_source = "cellular_fraction",
    community_source = viral_node_source,
    match_cols = match_cols,
    method = method,
    permutations = permutations
  ),
  
  run_vcr_community_mantel(
    vcr_df = vcr_nodes$VCR_Mega_Euk,
    community_mat = community_nodes$Eukaryotes,
    vcr_env = env_euk_filtered,
    community_env = env_by_node$Eukaryotes,
    vcr_name = "VCR_Mega_Euk",
    community_name = "Eukaryotes",
    vcr_source = "cellular_fraction",
    community_source = "cellular_fraction",
    match_cols = match_cols,
    method = method,
    permutations = permutations
  )
)

# ============================================================
# 7. Save outputs
# ============================================================

write_tsv(
  vcr_env_results,
  file.path(output_dir, "vcr_environment_mantel.tsv")
)

write_tsv(
  vcr_community_results,
  file.path(output_dir, "vcr_community_mantel.tsv")
)

vcr_mantel_results <- bind_rows(
  vcr_env_results,
  vcr_community_results
)

write_tsv(
  vcr_mantel_results,
  file.path(output_dir, "vcr_mantel_all.tsv")
)

message("Done.")
message("Wrote: ", file.path(output_dir, "vcr_environment_mantel.tsv"))
message("Wrote: ", file.path(output_dir, "vcr_community_mantel.tsv"))
message("Wrote: ", file.path(output_dir, "vcr_mantel_all.tsv"))