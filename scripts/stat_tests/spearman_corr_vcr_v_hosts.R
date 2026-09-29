library(dplyr)
library(stringr)
library(tidyr)

# All epipelagic samples
epi_runs_proks <- prok_cov %>% # prok_cov is the coverage file of prok otus with taxonomy and metadata
  filter(layer == "epi" & manual_env == "seawater") %>%
  distinct(Run) %>%
  pull(Run)


genus_by_run_epi_proks <- prok_cov %>%
  filter(layer == "epi") %>%
  
  # exclude rows containing multiple taxonomy paths
  filter(str_count(taxonomy, fixed("Root;")) == 1) %>%
  
  mutate(
    genus = str_extract(taxonomy, "(?<=^|;\\s)g__[^;]+"),
    genus = str_remove(genus, "^g__")
  ) %>%
  filter(!is.na(genus)) %>%
  
  group_by(Run, genus) %>%
  summarise(
    genus_norm_cov = sum(mean_norm_coverage, na.rm = TRUE),
    .groups = "drop"
  )

prok_prev <- genus_by_run_epi_proks %>%
  group_by(genus) %>%
  summarise(
    n_positive_runs = n_distinct(Run[genus_norm_cov > 0]),
    .groups = "drop"
  )

top20_genera_epi_proks <- genus_by_run_epi_proks %>%
  complete(
    Run = epi_runs_proks,
    genus,
    fill = list(genus_norm_cov = 0)
  ) %>%
  group_by(genus) %>%
  summarise(
    mean_norm_cov = mean(genus_norm_cov),
    .groups = "drop"
  ) %>%
  left_join(prok_prev, by = "genus") %>%
  filter(n_positive_runs >= 20) %>%
  arrange(desc(mean_norm_cov)) %>%
  slice_head(n = 20)

top20_genera_epi_proks

# Top 20 genus names
top20_prok_names <- top20_genera_epi_proks$genus

# Create table inputs for Mantel test
prok_comm_table <- genus_by_run_epi_proks %>%
  filter(
    genus %in% top20_prok_names
  ) %>%
  complete(
    Run = epi_runs_proks,
    genus = top20_prok_names,
    fill = list(genus_norm_cov = 0)
  ) %>%
  pivot_wider(
    names_from = genus,
    values_from = genus_norm_cov,
    values_fill = 0
  ) %>%
  arrange(Run)

caudo_comm_long <- votu_cov_sw_epi %>% # votu_cov_sw_epi is the coverage file for votus with taxonomy and metadata, filtered for epipelagic seawater
  filter(
    Run %in% epi_runs_proks,
    class == "Caudoviricetes"
  ) %>%
  group_by(Run, virus_genome) %>%
  summarise(
    virus_norm_cov = sum(mean_norm_coverage, na.rm = TRUE),
    .groups = "drop"
  )

caudo_runs <- caudo_comm_long %>%
  distinct(Run) %>%
  pull(Run)

shared_runs_caudo_comm <- intersect(
  epi_runs_proks,
  caudo_runs
)

prok_comm_table <- prok_comm_table %>%
  filter(Run %in% shared_runs_caudo_comm) %>%
  arrange(Run)

caudo_taxa <- caudo_comm_long %>%
  distinct(virus_genome) %>%
  pull(virus_genome)

caudo_comm_table <- caudo_comm_long %>%
  filter(Run %in% shared_runs_caudo_comm) %>%
  complete(
    Run = shared_runs_caudo_comm,
    virus_genome = caudo_taxa,
    fill = list(virus_norm_cov = 0)
  ) %>%
  pivot_wider(
    names_from = virus_genome,
    values_from = virus_norm_cov,
    values_fill = 0
  ) %>%
  arrange(Run)

identical(
  prok_comm_table$Run,
  caudo_comm_table$Run
)

dim(prok_comm_table)
dim(caudo_comm_table)

###################################

library(dplyr)
library(stringr)
library(tidyr)

# All epipelagic samples, BEFORE any taxonomic filtering
epi_runs_euks <- euk_cov %>% # euk_cov is the coverage file of euk otus with taxonomy and metadata
  filter(layer == "epi" & manual_env == "seawater") %>%
  distinct(Run) %>%
  pull(Run)


# Genus abundance within each epipelagic Run
genus_by_run_epi_euks <- euk_cov %>%
  filter(layer == "epi") %>%
  
  # Remove multiple taxonomy paths
  filter(!str_detect(taxon, fixed(","))) %>%
  
  # Require a pipe: entries without one are unknown eukaryotes
  filter(str_detect(taxon, fixed("|"))) %>%
  
  # Extract genus
  mutate(
    genus = taxon %>%
      str_remove("^.*\\|") %>%     # remove everything before/including pipe
      str_trim() %>%               # remove leading/trailing whitespace
      word(1) %>%                  # first word after pipe
      str_remove("_.*$")           # remove underscore and everything after
  ) %>%
  
  # Keep only sensible genus names
  filter(
    !is.na(genus),
    genus != "",
    str_detect(genus, "^[A-Za-z]")
  ) %>%
  
  # Collapse species/sequences belonging to the same genus
  # within each sample
  group_by(Run, genus) %>%
  summarise(
    genus_norm_cov = sum(mean_norm_coverage, na.rm = TRUE),
    .groups = "drop"
  )

euk_prev <- genus_by_run_epi_euks %>%
  group_by(genus) %>%
  summarise(
    n_positive_runs = n_distinct(Run[genus_norm_cov > 0]),
    .groups = "drop"
  )

top20_genera_epi_euks <- genus_by_run_epi_euks %>%
  complete(
    Run = epi_runs_euks,
    genus,
    fill = list(genus_norm_cov = 0)
  ) %>%
  group_by(genus) %>%
  summarise(
    mean_norm_cov = mean(genus_norm_cov),
    .groups = "drop"
  ) %>%
  left_join(euk_prev, by = "genus") %>%
  filter(n_positive_runs >= 20) %>%
  arrange(desc(mean_norm_cov)) %>%
  slice_head(n = 20)


top20_genera_epi_euks

# Top 20 genus names
top20_euk_names <- top20_genera_epi_euks$genus

# Create table inputs for Mantel test
euk_comm_table <- genus_by_run_epi_euks %>%
  filter(
    genus %in% top20_euk_names
  ) %>%
  complete(
    Run = epi_runs_euks,
    genus = top20_euk_names,
    fill = list(genus_norm_cov = 0)
  ) %>%
  pivot_wider(
    names_from = genus,
    values_from = genus_norm_cov,
    values_fill = 0
  ) %>%
  arrange(Run)

mega_comm_long <- votu_cov_sw_epi %>%
  filter(
    Run %in% epi_runs_euks,
    class == "Megaviricetes"
  ) %>%
  group_by(Run, virus_genome) %>%
  summarise(
    virus_norm_cov = sum(mean_norm_coverage, na.rm = TRUE),
    .groups = "drop"
  )

mega_runs <- mega_comm_long %>%
  distinct(Run) %>%
  pull(Run)

shared_runs_mega_comm <- intersect(
  epi_runs_euks,
  mega_runs
)

euk_comm_table <- euk_comm_table %>%
  filter(Run %in% shared_runs_mega_comm) %>%
  arrange(Run)

mega_taxa <- mega_comm_long %>%
  distinct(virus_genome) %>%
  pull(virus_genome)

mega_comm_table <- mega_comm_long %>%
  filter(Run %in% shared_runs_mega_comm) %>%
  complete(
    Run = shared_runs_mega_comm,
    virus_genome = mega_taxa,
    fill = list(virus_norm_cov = 0)
  ) %>%
  pivot_wider(
    names_from = virus_genome,
    values_from = virus_norm_cov,
    values_fill = 0
  ) %>%
  arrange(Run)

identical(
  euk_comm_table$Run,
  mega_comm_table$Run
)

dim(euk_comm_table)
dim(mega_comm_table)


##### Make spearman correlations
make_heatmap_df <- function(corr_mat, p_mat, plot_order) {

  # ----------------------------------------------------------
  # BH correction using unique pairwise tests only
  # ----------------------------------------------------------

  p_adj_mat <- matrix(
    NA_real_,
    nrow = nrow(p_mat),
    ncol = ncol(p_mat),
    dimnames = dimnames(p_mat)
  )

  # Upper triangle = unique comparisons
  upper_idx <- upper.tri(p_mat)

  # BH-adjust those p-values
  p_adj_mat[upper_idx] <- p.adjust(
    p_mat[upper_idx],
    method = "BH"
  )

  # Mirror upper triangle onto lower triangle
  p_adj_mat[lower.tri(p_adj_mat)] <- t(p_adj_mat)[lower.tri(p_adj_mat)]

  # Diagonal isn't tested
  diag(p_adj_mat) <- NA_real_

  # ----------------------------------------------------------
  # Convert matrices to long format
  # ----------------------------------------------------------

  corr_df <- as.data.frame(as.table(corr_mat)) %>%
    rename(
      dataset1 = Var1,
      dataset2 = Var2,
      spearman_rho = Freq
    )

  p_df <- as.data.frame(as.table(p_adj_mat)) %>%
    rename(
      dataset1 = Var1,
      dataset2 = Var2,
      p_adj_BH = Freq
    )

  # ----------------------------------------------------------
  # Combine + significance stars
  # ----------------------------------------------------------

  out <- corr_df %>%
    left_join(
      p_df,
      by = c("dataset1", "dataset2")
    ) %>%
    mutate(
      dataset1 = as.character(dataset1),
      dataset2 = as.character(dataset2),

      sig_label = case_when(
        dataset1 == dataset2 ~ "",
        p_adj_BH < 0.001 ~ "***",
        p_adj_BH < 0.01  ~ "**",
        p_adj_BH < 0.05  ~ "*",
        TRUE ~ ""
      ),

      row_id = match(dataset1, plot_order),
      col_id = match(dataset2, plot_order)
    )

  out
}

prok_cor_data <- prok_mantel_table %>%
  inner_join(caudo_ratio_table, by = "Run") %>%
  select(-Run)

names(prok_cor_data)[
  names(prok_cor_data) == "virus_prok_ratio"
] <- "Caudo_VCR"

prok_rcorr <- Hmisc::rcorr(
  as.matrix(prok_cor_data),
  type = "spearman"
)

corr_mat <- prok_rcorr$r
p_mat    <- prok_rcorr$P

corr_for_clustering <- corr_mat
corr_for_clustering[is.na(corr_for_clustering)] <- 0
diag(corr_for_clustering) <- 1

cluster_dist <- as.dist(1 - corr_for_clustering)
hc <- hclust(cluster_dist, method = "average")

clustered_order <- hc$labels[hc$order]

corr_mat_clustered <- corr_mat[
  clustered_order,
  clustered_order,
  drop = FALSE
]

p_mat_clustered <- p_mat[
  clustered_order,
  clustered_order,
  drop = FALSE
]

prok_heatmap_df <- make_heatmap_df(
  corr_mat = corr_mat_clustered,
  p_mat = p_mat_clustered,
  plot_order = clustered_order
)

p_prok <- ggplot(
  prok_heatmap_df,
  aes(
    x = factor(dataset2, levels = clustered_order),
    y = factor(dataset1, levels = rev(clustered_order)),
    fill = spearman_rho
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sig_label),
    color = "black",
    size = 3
  ) +
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
    title = "Caudoviricetes V:C vs prokaryotic genera",
    subtitle = "Clustered pairwise Spearman correlations"
  )

p_prok

euk_cor_data <- euk_mantel_table %>%
  inner_join(mega_ratio_table, by = "Run") %>%
  select(-Run)

names(euk_cor_data)[
  names(euk_cor_data) == "virus_euk_ratio"
] <- "Mega_VCR"

euk_rcorr <- Hmisc::rcorr(
  as.matrix(euk_cor_data),
  type = "spearman"
)

corr_mat <- euk_rcorr$r
p_mat    <- euk_rcorr$P

corr_for_clustering <- corr_mat
corr_for_clustering[is.na(corr_for_clustering)] <- 0
diag(corr_for_clustering) <- 1

cluster_dist <- as.dist(1 - corr_for_clustering)
hc <- hclust(cluster_dist, method = "average")

clustered_order <- hc$labels[hc$order]

corr_mat_clustered <- corr_mat[
  clustered_order,
  clustered_order,
  drop = FALSE
]

p_mat_clustered <- p_mat[
  clustered_order,
  clustered_order,
  drop = FALSE
]

euk_heatmap_df <- make_heatmap_df(
  corr_mat = corr_mat_clustered,
  p_mat = p_mat_clustered,
  plot_order = clustered_order
)

p_euk <- ggplot(
  euk_heatmap_df,
  aes(
    x = factor(dataset2, levels = clustered_order),
    y = factor(dataset1, levels = rev(clustered_order)),
    fill = spearman_rho
  )
) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = sig_label),
    color = "black",
    size = 3
  ) +
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
    title = "Megaviricetes V:C vs eukaryotic genera",
    subtitle = "Clustered pairwise Spearman correlations"
  )

p_euk
