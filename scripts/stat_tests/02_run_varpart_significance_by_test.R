#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(vegan)
  library(readr)
})

base_dir <- "/path/to/base/dir"
vp_dir <- file.path(base_dir, "variation_partitioning")
result_dir <- file.path(vp_dir, "results")
sig_dir <- file.path(result_dir, "significance_tests")
dir.create(sig_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide test_id 1-6.")
}

test_id <- as.integer(args[1])
permutations <- ifelse(length(args) >= 2, as.integer(args[2]), 99)

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

message("Running test_id: ", test_id)
message("Group: ", this_test$group)
message("Test: ", this_test$fraction)
message("Permutations: ", permutations)

if (this_test$test_name == "temp_given_abs_lat_host") {
  
  mod <- vegan::rda(
    vp$virus_hell,
    vp$temp_scaled,
    cbind(vp$abs_lat_scaled, vp$host_axes)
  )
  
} else if (this_test$test_name == "abs_lat_given_temp_host") {
  
  mod <- vegan::rda(
    vp$virus_hell,
    vp$abs_lat_scaled,
    cbind(vp$temp_scaled, vp$host_axes)
  )
  
} else if (this_test$test_name == "host_given_temp_abs_lat") {
  
  mod <- vegan::rda(
    vp$virus_hell,
    vp$host_axes,
    cbind(vp$temp_scaled, vp$abs_lat_scaled)
  )
  
} else {
  stop("Unknown test_name: ", this_test$test_name)
}

set.seed(123)

anova_res <- vegan::anova.cca(
  mod,
  permutations = permutations
)

anova_df <- as.data.frame(anova_res) %>%
  rownames_to_column("term") %>%
  as_tibble()

out <- anova_df %>%
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
  sig_dir,
  paste0("partial_rda_significance_", test_id, "_", safe_group, "_", safe_test, ".tsv")
)

readr::write_tsv(out, outfile)

txtfile <- sub("\\.tsv$", ".txt", outfile)
writeLines(capture.output(anova_res), txtfile)

message("Done.")
message("Wrote: ", outfile)