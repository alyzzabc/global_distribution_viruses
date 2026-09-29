#!/usr/bin/env Rscript

library(vegan)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript 01_run_virus_v_genus.R ",
    "<top20_table.rds> <virus_comm_table.rds> <results_dir>"
  )
}

top20_file <- args[1]
virus_file <- args[2]
results_dir <- args[3]

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Array index
# ------------------------------------------------------------

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(task_id)) {
  stop("SLURM_ARRAY_TASK_ID is not set.")
}

# ------------------------------------------------------------
# Read tables
#
# Expected:
#   top20_table: Run + 20 genus columns
#   virus_table: Run + viral taxon/community columns
# ------------------------------------------------------------

top20_table <- readRDS(top20_file)
virus_table <- readRDS(virus_file)

# ------------------------------------------------------------
# Make sure Runs match and are in same order
# ------------------------------------------------------------

shared_runs <- intersect(
  top20_table$Run,
  virus_table$Run
)

top20_table <- top20_table[
  top20_table$Run %in% shared_runs,
]

virus_table <- virus_table[
  virus_table$Run %in% shared_runs,
]

top20_table <- top20_table[
  match(shared_runs, top20_table$Run),
]

virus_table <- virus_table[
  match(shared_runs, virus_table$Run),
]

stopifnot(
  identical(
    as.character(top20_table$Run),
    as.character(virus_table$Run)
  )
)

# ------------------------------------------------------------
# Pick genus according to array index 1-20
# ------------------------------------------------------------

genus_names <- setdiff(
  colnames(top20_table),
  "Run"
)

if (task_id < 1 || task_id > length(genus_names)) {
  stop(
    "Array task ",
    task_id,
    " is outside available genus columns (1-",
    length(genus_names),
    ")."
  )
}

genus <- genus_names[task_id]

cat("Task:", task_id, "\n")
cat("Genus:", genus, "\n")
cat("Shared Runs:", length(shared_runs), "\n")

# ------------------------------------------------------------
# Genus abundance vector
# ------------------------------------------------------------

genus_abund <- as.numeric(
  top20_table[[genus]]
)

# Need variation for a meaningful distance matrix
if (length(unique(genus_abund)) < 2) {
  stop(
    "Genus ", genus,
    " has no variation across shared Runs."
  )
}

genus_dist <- dist(
  genus_abund,
  method = "euclidean"
)

# ------------------------------------------------------------
# Viral community matrix
# ------------------------------------------------------------

virus_mat <- virus_table[
  ,
  setdiff(colnames(virus_table), "Run"),
  drop = FALSE
]

virus_mat <- as.matrix(virus_mat)

storage.mode(virus_mat) <- "numeric"

# Remove viral features that are zero across every shared Run
keep_virus <- colSums(
  virus_mat,
  na.rm = TRUE
) > 0

virus_mat <- virus_mat[, keep_virus, drop = FALSE]

cat(
  "Viral features retained:",
  ncol(virus_mat),
  "\n"
)

if (ncol(virus_mat) == 0) {
  stop("No non-zero viral features remain.")
}

# Bray-Curtis community distance
virus_dist <- vegdist(
  virus_mat,
  method = "bray"
)

# ------------------------------------------------------------
# Mantel test
# ------------------------------------------------------------

set.seed(12345 + task_id)

mantel_result <- vegan::mantel(
  genus_dist,
  virus_dist,
  method = "spearman",
  permutations = 999
)

# ------------------------------------------------------------
# Save compact result
# ------------------------------------------------------------

result_df <- data.frame(
  task_id = task_id,
  genus = genus,
  n_runs = length(shared_runs),
  n_virus_features = ncol(virus_mat),
  mantel_r = unname(mantel_result$statistic),
  p_value = mantel_result$signif
)

outfile <- file.path(
  results_dir,
  paste0("mantel_", task_id, ".csv")
)

write.csv(
  result_df,
  outfile,
  row.names = FALSE
)

cat("\nResult:\n")
print(result_df)

cat("\nSaved:", outfile, "\n")