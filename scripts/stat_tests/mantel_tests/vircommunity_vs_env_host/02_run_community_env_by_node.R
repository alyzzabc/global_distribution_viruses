library(dplyr)
library(tibble)
library(readr)
library(vegan)

# ============================================================
# 02_run_community_env_by_node.R
#
# Usage:
#   Rscript 02_run_community_env_by_node.R Caudoviricetes virus-v-cellular 999
#   Rscript 02_run_community_env_by_node.R Prokaryotes cellular-v-cellular 999
# ============================================================

base_dir <- "/path/to/base/dir"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide a node name, e.g. Caudoviricetes")
}

node_to_run <- args[1]
run_label <- ifelse(length(args) >= 2, args[2], "virus-v-cellular")
permutations <- ifelse(length(args) >= 3, as.integer(args[3]), 999)

if (is.na(permutations)) {
  stop("permutations must be an integer.")
}

out_dir <- file.path(base_dir, paste0(run_label, "_mantel"))
input_file <- file.path(out_dir, "mantel_inputs_filtered.rds")

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

message("Running community-env Mantel")
message("  node_to_run: ", node_to_run)
message("  run_label: ", run_label)
message("  permutations: ", permutations)
message("  out_dir: ", out_dir)
message("  input_file: ", input_file)

inputs <- readRDS(input_file)

community_nodes <- inputs$community_nodes
env_by_node <- inputs$env_by_node
env_vars <- inputs$env_vars
match_cols <- inputs$match_cols

if (!node_to_run %in% names(community_nodes)) {
  stop(
    "node_to_run not found in community_nodes: ",
    node_to_run,
    "\nAvailable nodes: ",
    paste(names(community_nodes), collapse = ", ")
  )
}

if (!node_to_run %in% names(env_by_node)) {
  stop(
    "node_to_run not found in env_by_node: ",
    node_to_run,
    "\nAvailable nodes: ",
    paste(names(env_by_node), collapse = ", ")
  )
}

# Keep only the requested node
community_nodes <- community_nodes[node_to_run]
env_by_node <- env_by_node[node_to_run]

message("CONFIRMED only running node: ", paste(names(community_nodes), collapse = ", "))

# ============================================================
# Mantel: Community–environment
# ============================================================

run_community_env_mantel <- function(
    community_nodes,
    env_by_node,
    env_vars,
    output_file = "community_env_mantel_checkpoint.tsv",
    method = "spearman",
    permutations = 999,
    seed = 123
) {
  
  set.seed(seed)
  
  if (file.exists(output_file)) {
    completed <- readr::read_tsv(output_file, show_col_types = FALSE)
    completed_keys <- paste(completed$dataset1, completed$dataset2, sep = "__")
    message("Loaded existing checkpoint with ", nrow(completed), " completed tests")
  } else {
    completed_keys <- character()
  }
  
  for (node_name in names(community_nodes)) {
    
    message("Processing community node: ", node_name)
    
    mat <- community_nodes[[node_name]]
    env_df <- env_by_node[[node_name]]
    
    for (env_var in env_vars) {
      
      key <- paste(node_name, env_var, sep = "__")
      
      if (key %in% completed_keys) {
        message("  Skipping already completed: ", key)
        next
      }
      
      if (!env_var %in% colnames(env_df)) next
      
      # One environmental value per Run
      tmp <- env_df %>%
        select(Run, env_value = all_of(env_var)) %>%
        filter(!is.na(env_value), is.finite(env_value)) %>%
        group_by(Run) %>%
        summarise(
          env_value = mean(env_value, na.rm = TRUE),
          .groups = "drop"
        )
      
      common_runs <- intersect(rownames(mat), tmp$Run)
      
      if (length(common_runs) < 3) next
      
      # Put environmental table in the exact order used for the matrix
      tmp_sub <- tmp %>%
        filter(Run %in% common_runs) %>%
        arrange(match(Run, common_runs))
      
      mat_sub <- mat[tmp_sub$Run, , drop = FALSE]
      
      # Remove empty features
      mat_sub <- mat_sub[, colSums(mat_sub > 0, na.rm = TRUE) > 0, drop = FALSE]
      
      # Remove empty samples and matching environmental values
      keep_samples <- rowSums(mat_sub > 0, na.rm = TRUE) > 0 &
        !is.na(tmp_sub$env_value) &
        is.finite(tmp_sub$env_value)
      
      mat_sub <- mat_sub[keep_samples, , drop = FALSE]
      tmp_sub <- tmp_sub[keep_samples, , drop = FALSE]
      
      if (nrow(mat_sub) < 3 || ncol(mat_sub) < 1) next
      if (length(unique(tmp_sub$env_value)) < 2) next
      
      # Anonymous shared sample order
      rownames(mat_sub) <- seq_len(nrow(mat_sub))
      
      env_vec <- as.numeric(tmp_sub$env_value)
      names(env_vec) <- seq_len(length(env_vec))
      
      comm_dist <- vegan::vegdist(mat_sub, method = "bray")
      env_dist <- stats::dist(scale(env_vec), method = "euclidean")
      
      if (any(is.na(comm_dist)) || any(is.na(env_dist))) next
      
      # Force clean dist objects
      comm_dist <- as.dist(as.matrix(comm_dist))
      env_dist <- as.dist(as.matrix(env_dist))
      
      if (attr(comm_dist, "Size") != attr(env_dist, "Size")) {
        stop(
          "Distance Size mismatch for ",
          node_name, " vs ", env_var,
          ": comm_size=", attr(comm_dist, "Size"),
          ", env_size=", attr(env_dist, "Size")
        )
      }
      
      if (length(comm_dist) != length(env_dist)) {
        stop(
          "Distance length mismatch for ",
          node_name, " vs ", env_var,
          ": comm_dist=", length(comm_dist),
          ", env_dist=", length(env_dist)
        )
      }
      
      # Extra diagnostic: if this works, mantel should work
      test_cor <- suppressWarnings(
        cor(
          as.vector(comm_dist),
          as.vector(env_dist),
          method = method,
          use = "complete.obs"
        )
      )
      
      if (is.na(test_cor)) {
        message("  Skipping ", node_name, " vs ", env_var, ": cor() returned NA")
        next
      }
      
      message(
        "  Mantel: ", node_name, " vs ", env_var,
        " n=", attr(comm_dist, "Size"),
        " length=", length(comm_dist)
      )
      
      res <- vegan::mantel(
        comm_dist,
        env_dist,
        method = method,
        permutations = permutations
      )
      
      one_result <- tibble::tibble(
        dataset1 = node_name,
        dataset2 = env_var,
        comparison_type = "community_environment",
        mantel_r = unname(res$statistic),
        p_value = res$signif,
        n_samples = attr(comm_dist, "Size"),
        permutations = permutations,
        method = method
      )
      
      if (file.exists(output_file)) {
        readr::write_tsv(
          one_result,
          output_file,
          append = TRUE,
          col_names = FALSE
        )
      } else {
        readr::write_tsv(
          one_result,
          output_file
        )
      }
      
      completed_keys <- c(completed_keys, key)
    }
  }
  
  readr::read_tsv(output_file, show_col_types = FALSE)
}


community_env_results <- run_community_env_mantel(
  community_nodes = community_nodes,
  env_by_node = env_by_node,
  env_vars = env_vars,
  output_file = file.path(out_dir, paste0("community_env_", node_to_run, ".tsv")),
  permutations = permutations
)