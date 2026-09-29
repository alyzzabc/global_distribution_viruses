# Inputs from mantel_tests/01_prepare_mantel_inputs.R

# ============================================================
# Spearman correlations with audit table + full heatmap
# Same source = match by Run
# Different source = match by metadata
# ============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(readr)

# ============================================================
# 1. Total abundance per Run
# ============================================================

make_total_abund_df <- function(mat, node_name) {
  tibble::tibble(
    Run = rownames(mat),
    !!node_name := rowSums(mat, na.rm = TRUE)
  )
}

# ============================================================
# 2. Define node values
# ============================================================

node_values <- list(
  Caudoviricetes = make_total_abund_df(
    caudo_mat_filt,
    "Caudoviricetes"
  ),
  
  Megaviricetes = make_total_abund_df(
    mega_mat_filt,
    "Megaviricetes"
  ),
  
  All_viruses_known_tax = make_total_abund_df(
    all_viruses_known_mat_filt,
    "All_viruses_known_tax"
  ),
  
  Prokaryotes = make_total_abund_df(
    prok_mat_filt,
    "Prokaryotes"
  ),
  
  Eukaryotes = make_total_abund_df(
    euk_mat_filt,
    "Eukaryotes"
  ),
  
  VCR_Caudo_Prok = vcr_nodes$VCR_Caudo_Prok %>%
    rename(VCR_Caudo_Prok = value),
  
  VCR_Mega_Euk = vcr_nodes$VCR_Mega_Euk %>%
    rename(VCR_Mega_Euk = value)
)

# ============================================================
# 3. Environmental metadata table for each node
# ============================================================

env_by_node_spearman <- list(
  Caudoviricetes = env_vir_filtered,
  Megaviricetes = env_vir_filtered,
  All_viruses_known_tax = env_vir_filtered,
  Prokaryotes = env_prok_filtered,
  Eukaryotes = env_euk_filtered,
  VCR_Caudo_Prok = env_prok_filtered,
  VCR_Mega_Euk = env_euk_filtered
)

# ============================================================
# 4. Define source universe for each node
# ============================================================

viral_node_source <- if (exists("vir_category") && vir_category == "cellular") {
  "cellular_fraction"
} else {
  "viral_fraction"
}

node_source <- c(
  Caudoviricetes = viral_node_source,
  Megaviricetes = viral_node_source,
  All_viruses_known_tax = viral_node_source,
  Prokaryotes = "cellular_fraction",
  Eukaryotes = "cellular_fraction",
  VCR_Caudo_Prok = "cellular_fraction",
  VCR_Mega_Euk = "cellular_fraction"
)

stopifnot(all(names(node_values) %in% names(node_source)))

# ============================================================
# 5. Matching columns for cross-source comparisons
# ============================================================

match_cols <- c(
  "manual_longitude",
  "manual_depth",
  "manual_latitude",
  "manual_date"
)

# ============================================================
# 6. Helper functions
# ============================================================

get_value_col <- function(df) {
  value_cols <- setdiff(colnames(df), "Run")
  
  if (length(value_cols) != 1) {
    stop(
      "Expected exactly one value column besides Run, found: ",
      paste(value_cols, collapse = ", ")
    )
  }
  
  value_cols[1]
}

# ============================================================
# 7. Node-node Spearman with audit table
# ============================================================

run_pairwise_matched_spearman <- function(
    node_values,
    env_by_node,
    node_source,
    match_cols = c("manual_longitude", "manual_depth", "manual_latitude", "manual_date"),
    method = "spearman",
    min_samples = 3
) {
  
  node_pairs <- combn(names(node_values), 2, simplify = FALSE)
  
  results <- list()
  pairwise_inputs <- list()
  
  for (pair in node_pairs) {
    
    name1 <- pair[1]
    name2 <- pair[2]
    
    df1 <- node_values[[name1]]
    df2 <- node_values[[name2]]
    
    env1 <- env_by_node[[name1]]
    env2 <- env_by_node[[name2]]
    
    source1 <- unname(node_source[[name1]])
    source2 <- unname(node_source[[name2]])
    
    same_source <- identical(source1, source2)
    
    if (is.null(df1) || is.null(df2)) next
    if (!same_source && (is.null(env1) || is.null(env2))) next
    
    value_col1 <- get_value_col(df1)
    value_col2 <- get_value_col(df2)
    
    matched <- if (same_source) {
      
      df1 %>%
        rename(value1 = all_of(value_col1)) %>%
        inner_join(
          df2 %>%
            rename(value2 = all_of(value_col2)),
          by = "Run"
        ) %>%
        mutate(
          Run_1 = Run,
          Run_2 = Run,
          matching_basis = "Run"
        ) %>%
        select(
          Run_1,
          Run_2,
          matching_basis,
          value1,
          value2
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
      
      map1 %>%
        inner_join(
          map2,
          by = match_cols,
          suffix = c("_1", "_2"),
          relationship = "many-to-many"
        ) %>%
        inner_join(
          df1 %>%
            rename(Run_1 = Run, value1 = all_of(value_col1)),
          by = "Run_1"
        ) %>%
        inner_join(
          df2 %>%
            rename(Run_2 = Run, value2 = all_of(value_col2)),
          by = "Run_2"
        ) %>%
        distinct(Run_1, Run_2, .keep_all = TRUE) %>%
        mutate(matching_basis = "metadata_match_cols")
    }
    
    matched <- matched %>%
      mutate(
        dataset1 = name1,
        dataset2 = name2,
        source1 = source1,
        source2 = source2
      ) %>%
      relocate(dataset1, dataset2, source1, source2, matching_basis) %>%
      filter(
        !is.na(value1),
        !is.na(value2),
        is.finite(value1),
        is.finite(value2)
      )
    
    pair_key <- paste(name1, name2, sep = "__")
    pairwise_inputs[[pair_key]] <- matched
    
    n_samples <- nrow(matched)
    
    if (n_samples < min_samples) next
    if (length(unique(matched$value1)) < 2) next
    if (length(unique(matched$value2)) < 2) next
    
    test <- suppressWarnings(
      cor.test(
        matched$value1,
        matched$value2,
        method = method,
        exact = FALSE
      )
    )
    
    n_unique_match_keys <- if (all(match_cols %in% colnames(matched))) {
      matched %>%
        distinct(across(all_of(match_cols))) %>%
        nrow()
    } else {
      n_samples
    }
    
    unique_runs_1 <- dplyr::n_distinct(matched$Run_1)
    unique_runs_2 <- dplyr::n_distinct(matched$Run_2)
    
    results[[pair_key]] <- tibble(
      dataset1 = name1,
      dataset2 = name2,
      comparison_type = "matched_node_node",
      matching_basis = unique(matched$matching_basis),
      source1 = source1,
      source2 = source2,
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      n_samples = n_samples,
      n_unique_match_keys = n_unique_match_keys,
      unique_runs_1 = unique_runs_1,
      unique_runs_2 = unique_runs_2,
      duplicated_runs_1 = n_samples - unique_runs_1,
      duplicated_runs_2 = n_samples - unique_runs_2,
      method = method
    )
  }
  
  list(
    results = bind_rows(results),
    inputs = bind_rows(pairwise_inputs)
  )
}

spearman_node_node_out <- run_pairwise_matched_spearman(
  node_values = node_values,
  env_by_node = env_by_node_spearman,
  node_source = node_source,
  match_cols = match_cols,
  method = "spearman",
  min_samples = 3
)

spearman_node_node <- spearman_node_node_out$results
spearman_node_node_inputs <- spearman_node_node_out$inputs

# ============================================================
# 8. Node-environment Spearman with audit table
# ============================================================

run_node_env_spearman <- function(
    node_values,
    env_by_node,
    env_vars,
    node_source,
    method = "spearman",
    min_samples = 3
) {
  
  results <- list()
  pairwise_inputs <- list()
  
  for (node_name in names(node_values)) {
    
    df <- node_values[[node_name]]
    env_df <- env_by_node[[node_name]]
    
    if (is.null(df) || is.null(env_df)) next
    
    value_col <- get_value_col(df)
    
    for (env_var in env_vars) {
      
      if (!env_var %in% colnames(env_df)) next
      
      tmp <- df %>%
        rename(value1 = all_of(value_col)) %>%
        inner_join(
          env_df %>%
            select(Run, value2 = all_of(env_var)) %>%
            distinct(),
          by = "Run"
        ) %>%
        mutate(
          dataset1 = node_name,
          dataset2 = env_var,
          source1 = unname(node_source[[node_name]]),
          source2 = "environment",
          matching_basis = "Run",
          Run_1 = Run,
          Run_2 = Run
        ) %>%
        select(
          dataset1,
          dataset2,
          source1,
          source2,
          matching_basis,
          Run_1,
          Run_2,
          value1,
          value2
        ) %>%
        filter(
          !is.na(value1),
          !is.na(value2),
          is.finite(value1),
          is.finite(value2)
        )
      
      pair_key <- paste(node_name, env_var, sep = "__")
      pairwise_inputs[[pair_key]] <- tmp
      
      n_samples <- nrow(tmp)
      
      if (n_samples < min_samples) next
      if (length(unique(tmp$value1)) < 2) next
      if (length(unique(tmp$value2)) < 2) next
      
      test <- suppressWarnings(
        cor.test(
          tmp$value1,
          tmp$value2,
          method = method,
          exact = FALSE
        )
      )
      
      results[[pair_key]] <- tibble(
        dataset1 = node_name,
        dataset2 = env_var,
        comparison_type = "node_environment",
        matching_basis = "Run",
        source1 = unname(node_source[[node_name]]),
        source2 = "environment",
        spearman_rho = unname(test$estimate),
        p_value = test$p.value,
        n_samples = n_samples,
        n_unique_match_keys = n_samples,
        unique_runs_1 = dplyr::n_distinct(tmp$Run_1),
        unique_runs_2 = dplyr::n_distinct(tmp$Run_2),
        duplicated_runs_1 = n_samples - dplyr::n_distinct(tmp$Run_1),
        duplicated_runs_2 = 0,
        method = method
      )
    }
  }
  
  list(
    results = bind_rows(results),
    inputs = bind_rows(pairwise_inputs)
  )
}

spearman_node_env_out <- run_node_env_spearman(
  node_values = node_values,
  env_by_node = env_by_node_spearman,
  env_vars = env_vars,
  node_source = node_source,
  method = "spearman",
  min_samples = 3
)

spearman_node_env <- spearman_node_env_out$results
spearman_node_env_inputs <- spearman_node_env_out$inputs

# ============================================================
# 9. Environment-environment Spearman with audit table
# ============================================================

env_all_for_spearman <- bind_rows(
  env_vir_filtered %>% mutate(env_source = "viral_fraction"),
  env_prok_filtered %>% mutate(env_source = "cellular_prok"),
  env_euk_filtered %>% mutate(env_source = "cellular_euk")
) %>%
  select(any_of(c("Run", "env_source", match_cols, env_vars))) %>%
  distinct()

run_env_env_spearman <- function(
    env_df,
    env_vars,
    method = "spearman",
    min_samples = 3
) {
  
  env_vars_use <- env_vars[env_vars %in% colnames(env_df)]
  env_pairs <- combn(env_vars_use, 2, simplify = FALSE)
  
  results <- list()
  pairwise_inputs <- list()
  
  for (pair in env_pairs) {
    
    v1 <- pair[1]
    v2 <- pair[2]
    
    tmp <- env_df %>%
      select(
        Run,
        any_of(c("env_source", match_cols)),
        value1 = all_of(v1),
        value2 = all_of(v2)
      ) %>%
      mutate(
        dataset1 = v1,
        dataset2 = v2,
        source1 = "environment",
        source2 = "environment",
        matching_basis = "environment_table",
        Run_1 = Run,
        Run_2 = Run
      ) %>%
      select(
        dataset1,
        dataset2,
        source1,
        source2,
        matching_basis,
        Run_1,
        Run_2,
        any_of(c("env_source", match_cols)),
        value1,
        value2
      ) %>%
      filter(
        !is.na(value1),
        !is.na(value2),
        is.finite(value1),
        is.finite(value2)
      )
    
    pair_key <- paste(v1, v2, sep = "__")
    pairwise_inputs[[pair_key]] <- tmp
    
    n_samples <- nrow(tmp)
    
    if (n_samples < min_samples) next
    if (length(unique(tmp$value1)) < 2) next
    if (length(unique(tmp$value2)) < 2) next
    
    test <- suppressWarnings(
      cor.test(
        tmp$value1,
        tmp$value2,
        method = method,
        exact = FALSE
      )
    )
    
    results[[pair_key]] <- tibble(
      dataset1 = v1,
      dataset2 = v2,
      comparison_type = "environment_environment",
      matching_basis = "environment_table",
      source1 = "environment",
      source2 = "environment",
      spearman_rho = unname(test$estimate),
      p_value = test$p.value,
      n_samples = n_samples,
      n_unique_match_keys = n_samples,
      unique_runs_1 = dplyr::n_distinct(tmp$Run_1),
      unique_runs_2 = dplyr::n_distinct(tmp$Run_2),
      duplicated_runs_1 = n_samples - dplyr::n_distinct(tmp$Run_1),
      duplicated_runs_2 = 0,
      method = method
    )
  }
  
  list(
    results = bind_rows(results),
    inputs = bind_rows(pairwise_inputs)
  )
}

spearman_env_env_out <- run_env_env_spearman(
  env_df = env_all_for_spearman,
  env_vars = env_vars,
  method = "spearman",
  min_samples = 3
)

spearman_env_env <- spearman_env_env_out$results
spearman_env_env_inputs <- spearman_env_env_out$inputs

# ============================================================
# 10. Combine all audit inputs
# ============================================================

spearman_pairwise_input_df <- bind_rows(
  spearman_node_node_inputs,
  spearman_node_env_inputs,
  spearman_env_env_inputs
)

# ============================================================
# 11. Combine all Spearman results and adjust p-values
# ============================================================

spearman_results_df_adj <- bind_rows(
  spearman_node_node,
  spearman_node_env,
  spearman_env_env
) %>%
  filter(!is.na(p_value), !is.na(spearman_rho)) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significant_BH_0.05 = p_adj_BH < 0.05,
    abs_spearman_rho = abs(spearman_rho),
    source = dataset1,
    target = dataset2,
    weight = abs_spearman_rho,
    matching_type = case_when(
      matching_basis == "Run" ~ "same_Run",
      matching_basis == "metadata_match_cols" &
        duplicated_runs_1 == 0 &
        duplicated_runs_2 == 0 ~ "one_to_one_metadata",
      matching_basis == "metadata_match_cols" &
        (duplicated_runs_1 > 0 | duplicated_runs_2 > 0) ~ "many_to_many_metadata",
      matching_basis == "environment_table" ~ "environment_environment",
      TRUE ~ "unknown"
    ),
    pair_inflation_factor = n_samples / n_unique_match_keys
  ) %>%
  arrange(p_adj_BH, desc(abs_spearman_rho))



# ============================================================
# Heatmaps from spearman_results_df_adj
# Produces:
#   1. p_full
#   2. p_half
#   3. p_clustered_full
#   4. p_clustered_half
# ============================================================

library(dplyr)
library(tibble)
library(ggplot2)

# ------------------------------------------------------------
# 1. Node order
# ------------------------------------------------------------

node_order <- c(
  "Caudoviricetes",
  "Megaviricetes",
  "All_viruses_known_tax",
  "Prokaryotes",
  "Eukaryotes",
  "VCR_Caudo_Prok",
  "VCR_Mega_Euk",
  "env_latitude",
  "env_abs_latitude",
  "env_longitude",
  "env_depth",
  "env_temp",
  "env_chl",
  "env_O2",
  "env_nitrate",
  "env_salinity",
  "env_mld",
  "env_sea_surface_height",
  "env_phosphate",
  "env_silicate",
  "env_npp",
  "env_chl_0m_castant",
  "env_PAR_0m",
  "env_PAR"
)

all_nodes <- unique(c(
  spearman_results_df_adj$dataset1,
  spearman_results_df_adj$dataset2
))

node_order <- node_order[node_order %in% all_nodes]

# ------------------------------------------------------------
# 2. Build full symmetric rho matrix
# ------------------------------------------------------------

corr_mat <- matrix(
  NA_real_,
  nrow = length(node_order),
  ncol = length(node_order),
  dimnames = list(node_order, node_order)
)

diag(corr_mat) <- 1

for (i in seq_len(nrow(spearman_results_df_adj))) {
  a <- spearman_results_df_adj$dataset1[i]
  b <- spearman_results_df_adj$dataset2[i]
  r <- spearman_results_df_adj$spearman_rho[i]
  
  if (a %in% node_order && b %in% node_order) {
    corr_mat[a, b] <- r
    corr_mat[b, a] <- r
  }
}

# ------------------------------------------------------------
# 3. Build full symmetric adjusted p-value matrix
# ------------------------------------------------------------

p_mat <- matrix(
  NA_real_,
  nrow = length(node_order),
  ncol = length(node_order),
  dimnames = list(node_order, node_order)
)

for (i in seq_len(nrow(spearman_results_df_adj))) {
  a <- spearman_results_df_adj$dataset1[i]
  b <- spearman_results_df_adj$dataset2[i]
  p <- spearman_results_df_adj$p_adj_BH[i]
  
  if (a %in% node_order && b %in% node_order) {
    p_mat[a, b] <- p
    p_mat[b, a] <- p
  }
}

# ------------------------------------------------------------
# 4. Convert matrices to plotting table
# ------------------------------------------------------------

make_heatmap_df <- function(corr_mat, p_mat, plot_order) {
  
  corr_df <- as.data.frame(as.table(corr_mat)) %>%
    rename(
      dataset1 = Var1,
      dataset2 = Var2,
      spearman_rho = Freq
    )
  
  p_df <- as.data.frame(as.table(p_mat)) %>%
    rename(
      dataset1 = Var1,
      dataset2 = Var2,
      p_adj_BH = Freq
    )
  
  corr_df %>%
    left_join(p_df, by = c("dataset1", "dataset2")) %>%
    mutate(
      dataset1 = as.character(dataset1),
      dataset2 = as.character(dataset2),
      row_id = match(dataset1, plot_order),
      col_id = match(dataset2, plot_order),
      sig_label = case_when(
        dataset1 == dataset2 ~ "",
        is.na(p_adj_BH) ~ "",
        p_adj_BH < 0.001 ~ "***",
        p_adj_BH < 0.01 ~ "**",
        p_adj_BH < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
}

spearman_heatmap_full_df <- make_heatmap_df(
  corr_mat = corr_mat,
  p_mat = p_mat,
  plot_order = node_order
)

spearman_heatmap_half_df <- spearman_heatmap_full_df %>%
  filter(row_id >= col_id)

# ============================================================
# 5. FULL heatmap, fixed order
# ============================================================

p_full <- ggplot(
  spearman_heatmap_full_df,
  aes(
    x = factor(dataset2, levels = node_order),
    y = factor(dataset1, levels = rev(node_order)),
    fill = spearman_rho
  )
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho",
    na.value = "grey90"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Pairwise Spearman correlations",
    subtitle = "Full matrix; stars show BH-adjusted significance"
  )

p_full

# ============================================================
# 6. HALF heatmap, fixed order
# ============================================================

p_half <- ggplot(
  spearman_heatmap_half_df,
  aes(
    x = factor(dataset2, levels = node_order),
    y = factor(dataset1, levels = rev(node_order)),
    fill = spearman_rho
  )
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho",
    na.value = "grey90"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Pairwise Spearman correlations",
    subtitle = "Lower triangle; stars show BH-adjusted significance"
  )

p_half

# ============================================================
# 7. Clustered order
# ============================================================

corr_for_clustering <- corr_mat
corr_for_clustering[is.na(corr_for_clustering)] <- 0
diag(corr_for_clustering) <- 1

cluster_dist <- as.dist(1 - corr_for_clustering)
hc <- hclust(cluster_dist, method = "average")

clustered_order <- hc$labels[hc$order]

corr_mat_clustered <- corr_mat[clustered_order, clustered_order, drop = FALSE]
p_mat_clustered <- p_mat[clustered_order, clustered_order, drop = FALSE]

spearman_heatmap_clustered_full_df <- make_heatmap_df(
  corr_mat = corr_mat_clustered,
  p_mat = p_mat_clustered,
  plot_order = clustered_order
)

spearman_heatmap_clustered_half_df <- spearman_heatmap_clustered_full_df %>%
  filter(row_id >= col_id)

# ============================================================
# 8. CLUSTERED FULL heatmap
# ============================================================

p_clustered_full <- ggplot(
  spearman_heatmap_clustered_full_df,
  aes(
    x = factor(dataset2, levels = clustered_order),
    y = factor(dataset1, levels = rev(clustered_order)),
    fill = spearman_rho
  )
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho",
    na.value = "grey90"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Clustered pairwise Spearman correlations",
    subtitle = "Full matrix; variables clustered by Spearman correlation"
  )

p_clustered_full

# ============================================================
# 9. CLUSTERED HALF heatmap
# ============================================================

p_clustered_half <- ggplot(
  spearman_heatmap_clustered_half_df,
  aes(
    x = factor(dataset2, levels = clustered_order),
    y = factor(dataset1, levels = rev(clustered_order)),
    fill = spearman_rho
  )
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho",
    na.value = "grey90"
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Clustered pairwise Spearman correlations",
    subtitle = "Lower triangle; variables clustered by Spearman correlation"
  )

p_clustered_half


# ============================================================
# Useful checks
# ============================================================

# Check exact rows used for a pair:
 spearman_pairwise_input_df %>%
   filter(dataset1 == "Megaviricetes", dataset2 == "Eukaryotes") %>%
   select(dataset1, dataset2, matching_basis, Run_1, Run_2,
          any_of(match_cols), value1, value2) %>%
   head(20)

# Manually verify one correlation:
 tmp <- spearman_pairwise_input_df %>%
   filter(dataset1 == "Megaviricetes", dataset2 == "Eukaryotes")

 cor(tmp$value1, tmp$value2, method = "spearman", use = "complete.obs")