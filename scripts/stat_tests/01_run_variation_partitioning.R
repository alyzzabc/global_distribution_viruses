#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(vegan)
  library(readr)
})

base_dir <- "/path/to/base/dir"
vp_dir <- file.path(base_dir, "variation_partitioning")

input_file <- file.path(vp_dir, "variation_partitioning_input_tables.rds")
out_dir <- file.path(vp_dir, "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

inputs <- readRDS(input_file)

caudo_mat_filt <- inputs$caudo_mat_filt
prok_mat_filt <- inputs$prok_mat_filt
env_prok_filtered <- inputs$env_prok_filtered

mega_mat_filt <- inputs$mega_mat_filt
euk_mat_filt <- inputs$euk_mat_filt
env_euk_filtered <- inputs$env_euk_filtered

run_quick_varpart_3sets <- function(
    virus_mat,
    host_mat,
    virus_env,
    host_env,
    env_var1 = "env_temp",
    env_var2 = "env_abs_latitude",
    match_cols = c("manual_longitude", "manual_depth", "manual_latitude", "manual_date"),
    n_host_axes = 5
) {
  
  env_vars_use <- c(env_var1, env_var2)
  
  virus_map <- virus_env %>%
    select(Run, all_of(match_cols), all_of(env_vars_use)) %>%
    filter(if_all(all_of(match_cols), ~ !is.na(.x))) %>%
    filter(if_all(all_of(env_vars_use), ~ !is.na(.x))) %>%
    group_by(Run) %>%
    summarise(
      across(all_of(c(match_cols, env_vars_use)), ~ dplyr::first(.x)),
      .groups = "drop"
    )
  
  host_map <- host_env %>%
    select(Run, all_of(match_cols)) %>%
    filter(if_all(all_of(match_cols), ~ !is.na(.x))) %>%
    group_by(Run) %>%
    summarise(
      across(all_of(match_cols), ~ dplyr::first(.x)),
      .groups = "drop"
    )
  
  matched <- virus_map %>%
    inner_join(
      host_map,
      by = match_cols,
      suffix = c("_virus", "_host"),
      relationship = "many-to-many"
    ) %>%
    filter(
      Run_virus %in% rownames(virus_mat),
      Run_host %in% rownames(host_mat)
    ) %>%
    distinct(Run_virus, Run_host, .keep_all = TRUE)
  
  virus_sub <- virus_mat[matched$Run_virus, , drop = FALSE]
  host_sub <- host_mat[matched$Run_host, , drop = FALSE]
  
  env_sub <- matched %>%
    select(all_of(env_vars_use)) %>%
    mutate(across(everything(), as.numeric))
  
  keep <- rowSums(virus_sub, na.rm = TRUE) > 0 &
    rowSums(host_sub, na.rm = TRUE) > 0 &
    stats::complete.cases(env_sub)
  
  virus_sub <- virus_sub[keep, , drop = FALSE]
  host_sub <- host_sub[keep, , drop = FALSE]
  env_sub <- env_sub[keep, , drop = FALSE]
  matched <- matched[keep, , drop = FALSE]
  
  virus_sub <- virus_sub[, colSums(virus_sub > 0, na.rm = TRUE) > 0, drop = FALSE]
  host_sub <- host_sub[, colSums(host_sub > 0, na.rm = TRUE) > 0, drop = FALSE]
  
  message("Matched samples after filtering: ", nrow(host_sub))
  message("Virus features: ", ncol(virus_sub))
  message("Host features: ", ncol(host_sub))
  
  k_use <- min(n_host_axes, nrow(host_sub) - 1)
  
  if (k_use < 1) {
    stop("Too few matched samples after filtering for host PCoA.")
  }
  
  host_dist <- vegan::vegdist(host_sub, method = "bray")
  host_pcoa <- cmdscale(host_dist, k = k_use, eig = TRUE)
  
  host_axes <- as.data.frame(host_pcoa$points)
  colnames(host_axes) <- paste0("host_axis", seq_len(ncol(host_axes)))
  
  temp_scaled <- as.data.frame(scale(env_sub[[env_var1]]))
  colnames(temp_scaled) <- env_var1
  
  abs_lat_scaled <- as.data.frame(scale(env_sub[[env_var2]]))
  colnames(abs_lat_scaled) <- env_var2
  
  virus_hell <- vegan::decostand(virus_sub, method = "hellinger")
  
  vp <- vegan::varpart(
    virus_hell,
    temp_scaled,
    abs_lat_scaled,
    host_axes
  )
  
  host_axis_variance <- host_pcoa$eig / sum(host_pcoa$eig[host_pcoa$eig > 0])
  
  list(
    varpart = vp,
    matched = matched,
    virus_sub = virus_sub,
    host_sub = host_sub,
    virus_hell = virus_hell,
    temp_scaled = temp_scaled,
    abs_lat_scaled = abs_lat_scaled,
    host_axes = host_axes,
    host_pcoa_eig = host_pcoa$eig,
    host_axis_variance = host_axis_variance,
    n_samples = nrow(host_sub),
    n_virus_features = ncol(virus_sub),
    n_host_features = ncol(host_sub),
    n_host_axes = k_use
  )
}

message("Running Caudoviricetes ~ temp + abs latitude + Prokaryotes")

vp_caudo_3 <- run_quick_varpart_3sets(
  virus_mat = caudo_mat_filt,
  host_mat = prok_mat_filt,
  virus_env = env_prok_filtered,
  host_env = env_prok_filtered,
  env_var1 = "env_temp",
  env_var2 = "env_abs_latitude",
  n_host_axes = 5
)

saveRDS(
  vp_caudo_3,
  file.path(out_dir, "vp_caudo_temp_abs_lat_prok.rds")
)

writeLines(
  capture.output(vp_caudo_3$varpart),
  file.path(out_dir, "vp_caudo_temp_abs_lat_prok_summary.txt")
)

pdf(file.path(out_dir, "vp_caudo_temp_abs_lat_prok_plot.pdf"))
plot(vp_caudo_3$varpart)
dev.off()

readr::write_tsv(
  tibble(
    group = "Caudoviricetes_Prokaryotes",
    n_samples = vp_caudo_3$n_samples,
    n_virus_features = vp_caudo_3$n_virus_features,
    n_host_features = vp_caudo_3$n_host_features,
    n_host_axes = vp_caudo_3$n_host_axes,
    host_axis_variance_first5 = sum(vp_caudo_3$host_axis_variance[seq_len(vp_caudo_3$n_host_axes)], na.rm = TRUE)
  ),
  file.path(out_dir, "vp_caudo_diagnostics.tsv")
)

message("Running Megaviricetes ~ temp + abs latitude + Eukaryotes")

vp_mega_3 <- run_quick_varpart_3sets(
  virus_mat = mega_mat_filt,
  host_mat = euk_mat_filt,
  virus_env = env_euk_filtered,
  host_env = env_euk_filtered,
  env_var1 = "env_temp",
  env_var2 = "env_abs_latitude",
  n_host_axes = 5
)

saveRDS(
  vp_mega_3,
  file.path(out_dir, "vp_mega_temp_abs_lat_euk.rds")
)

writeLines(
  capture.output(vp_mega_3$varpart),
  file.path(out_dir, "vp_mega_temp_abs_lat_euk_summary.txt")
)

pdf(file.path(out_dir, "vp_mega_temp_abs_lat_euk_plot.pdf"))
plot(vp_mega_3$varpart)
dev.off()

readr::write_tsv(
  tibble(
    group = "Megaviricetes_Eukaryotes",
    n_samples = vp_mega_3$n_samples,
    n_virus_features = vp_mega_3$n_virus_features,
    n_host_features = vp_mega_3$n_host_features,
    n_host_axes = vp_mega_3$n_host_axes,
    host_axis_variance_first5 = sum(vp_mega_3$host_axis_variance[seq_len(vp_mega_3$n_host_axes)], na.rm = TRUE)
  ),
  file.path(out_dir, "vp_mega_diagnostics.tsv")
)

message("Done variation partitioning.")
message("Results written to: ", out_dir)(base)