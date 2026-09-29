#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(vegan)
  library(readr)
})

base_dir <- "/path/to/base/dir"
vp_dir <- file.path(base_dir, "permanova")
result_dir <- file.path(vp_dir, "results")
permanova_dir <- file.path(result_dir, "permanova_tests")
dir.create(permanova_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide test_id 1-6.")
}

test_id <- as.integer(args[1])
permutations <- ifelse(length(args) >= 2, as.integer(args[2]), 999)

test_table <- tibble::tibble(
  test_id = 1:6,
  group = c(
    "Caudoviricetes_Prokaryotes",
    "Caudoviricetes_Prokaryotes",
    "Caudoviricetes_Prokaryotes",
    "Megaviricetes_Eukaryotes",
    "Megaviricetes_Eukaryotes",
    "Megaviricetes_Eukaryotes"
  ),
  vp_file = c(
    "vp_caudo_temp_abs_lat_prok.rds",
    "vp_caudo_temp_abs_lat_prok.rds",
    "vp_caudo_temp_abs_lat_prok.rds",
    "vp_mega_temp_abs_lat_euk.rds",
    "vp_mega_temp_abs_lat_euk.rds",
    "vp_mega_temp_abs_lat_euk.rds"
  ),
  test_name = c(
    "temp_given_abs_lat_host",
    "abs_lat_given_temp_host",
    "host_given_temp_abs_lat",
    "temp_given_abs_lat_host",
    "abs_lat_given_temp_host",
    "host_given_temp_abs_lat"
  ),
  fraction = c(
    "temp | abs_lat + host",
    "abs_lat | temp + host",
    "host | temp + abs_lat",
    "temp | abs_lat + host",
    "abs_lat | temp + host",
    "host | temp + abs_lat"
  )
)

if (!test_id %in% test_table$test_id) {
  stop("Invalid test_id: ", test_id)
}

this_test <- test_table %>%
  filter(test_id == !!test_id)

vp <- readRDS(file.path(result_dir, this_test$vp_file))

message("Running PERMANOVA test_id: ", test_id)
message("Group: ", this_test$group)
message("Test: ", this_test$fraction)
message("Permutations: ", permutations)

# ------------------------------------------------------------
# Build viral Bray-Curtis distance matrix
# ------------------------------------------------------------

virus_dist <- vegan::vegdist(vp$virus_sub, method = "bray")

# ------------------------------------------------------------
# Build predictor table
# ------------------------------------------------------------

adonis_df <- bind_cols(
  vp$temp_scaled,
  vp$abs_lat_scaled,
  vp$host_axes
)

host_terms <- colnames(vp$host_axes)
host_formula <- paste(host_terms, collapse = " + ")

# ------------------------------------------------------------
# Partial PERMANOVA models
# ------------------------------------------------------------

if (this_test$test_name == "temp_given_abs_lat_host") {
  
  formula_use <- as.formula(
    paste0(
      "virus_dist ~ env_abs_latitude + ",
      host_formula,
      " + env_temp"
    )
  )
  
  target_term <- "env_temp"
  
} else if (this_test$test_name == "abs_lat_given_temp_host") {
  
  formula_use <- as.formula(
    paste0(
      "virus_dist ~ env_temp + ",
      host_formula,
      " + env_abs_latitude"
    )
  )
  
  target_term <- "env_abs_latitude"
  
} else if (this_test$test_name == "host_given_temp_abs_lat") {
  
  formula_use <- as.formula(
    paste0(
      "virus_dist ~ env_temp + env_abs_latitude + ",
      host_formula
    )
  )
  
  target_term <- "host_axes_block"
  
} else {
  stop("Unknown test_name: ", this_test$test_name)
}

message("Formula: ", deparse(formula_use))

set.seed(123)

adonis_res <- vegan::adonis2(
  formula_use,
  data = adonis_df,
  permutations = permutations,
  by = "terms"
)

adonis_df_out <- as.data.frame(adonis_res) %>%
  rownames_to_column("term") %>%
  as_tibble() %>%
  mutate(
    test_id = test_id,
    group = this_test$group,
    test_name = this_test$test_name,
    fraction = this_test$fraction,
    permutations = permutations,
    n_samples = vp$n_samples,
    n_virus_features = vp$n_virus_features,
    n_host_features = vp$n_host_features,
    n_host_axes = vp$n_host_axes
  ) %>%
  relocate(
    test_id,
    group,
    test_name,
    fraction,
    permutations,
    n_samples,
    n_virus_features,
    n_host_features,
    n_host_axes
  )

safe_group <- gsub("[^A-Za-z0-9_]+", "_", this_test$group)
safe_test <- gsub("[^A-Za-z0-9_]+", "_", this_test$test_name)

outfile <- file.path(
  permanova_dir,
  paste0("permanova_", test_id, "_", safe_group, "_", safe_test, ".tsv")
)

readr::write_tsv(adonis_df_out, outfile)

txtfile <- sub("\\.tsv$", ".txt", outfile)
writeLines(capture.output(adonis_res), txtfile)

message("Done.")
message("Wrote: ", outfile)