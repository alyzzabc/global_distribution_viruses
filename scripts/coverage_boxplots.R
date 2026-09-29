library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)
library(scales)

### Desired plotting order
desired_order <- c("virus", "girus", "prok", "cellular", "large", "verylarge")

### Groups to plot
all_groups <- c(
  "Caudoviricetes",
  "Megaviricetes",
  "Prokaryotes",
  "Eukaryotes"
)

### Colors
group_cols <- c(
  "Caudoviricetes" = "mediumorchid",
  "Megaviricetes"  = "orange",
  "Prokaryotes"    = "steelblue4",
  "Eukaryotes"     = "chartreuse3"
)

### Get all Runs
all_runs <- sort(unique(metadata$Run))

### Use only desired fractions
all_fractions <- desired_order

### Complete grid
run_categories <- metadata %>%
  distinct(Run, category) %>%
  filter(category %in% desired_order)

complete_grid <- run_categories %>%
  tidyr::crossing(group = all_groups)

### Viral groups: sum mean_norm_coverage per Run, category, class
df_virus_sum <- votu_cov_tax_sw %>% # coverage file for votus with taxonomy, filtered only for seawater
  filter(
    class %in% c("Caudoviricetes", "Megaviricetes"),
    category %in% desired_order
  ) %>%
  group_by(Run, category, group = class) %>%
  summarize(
    total_norm_coverage = sum(mean_norm_coverage, na.rm = TRUE),
    .groups = "drop"
  )

### Prokaryotes: sum mean_norm_coverage per Run, category
df_prok_sum <- prok_cov_sw %>% # coverage file for prok otus, filtered only for seawater
  filter(category %in% desired_order) %>%
  group_by(Run, category) %>%
  summarize(
    total_norm_coverage = sum(mean_norm_coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(group = "Prokaryotes")

### Eukaryotes: sum mean_norm_coverage per Run, category
df_euk_sum <- euk_cov_sw %>%  coverage file for euk otus, filtered only for seawater
  filter(category %in% desired_order) %>%
  group_by(Run, category) %>%
  summarize(
    total_norm_coverage = sum(mean_norm_coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(group = "Eukaryotes")

### Combine all groups
df_sum_all <- bind_rows(
  df_virus_sum,
  df_prok_sum,
  df_euk_sum
)

### Join into complete grid and fill missing combinations with 0
df_complete <- complete_grid %>%
  left_join(
    df_sum_all,
    by = c("Run", "category", "group")
  ) %>%
  mutate(
    total_norm_coverage = replace_na(total_norm_coverage, 0),
    category = factor(category, levels = desired_order),
    group = factor(group, levels = all_groups)
  )

### Keep only nonzero values for log10 plotting
df_nonzero_filt <- df_complete %>%
  filter(total_norm_coverage > 0)


#### Separate box plots ####

group_cols <- c(
  "Caudoviricetes" = "mediumorchid",
  "Megaviricetes"  = "orange",
  "Prokaryotes"    = "cyan2",
  "Eukaryotes"     = "chartreuse2"
)

outline_point_cols <- c(
  # boxplot outlines
  "Caudoviricetes_box" = alpha("black", 0.7),
  "Megaviricetes_box"  = alpha("black", 0.7),
  "Prokaryotes_box"    = alpha("cyan2", 0.65),
  "Eukaryotes_box"     = alpha("chartreuse2", 0.65),
  
  # point colors
  "Caudoviricetes_pt"  = "mediumorchid",
  "Megaviricetes_pt"   = "orange",
  "Prokaryotes_pt"     = "cyan2",
  "Eukaryotes_pt"      = "chartreuse2"
)

df_nonzero_filt <- df_nonzero_filt %>%
  mutate(
    category = factor(category, levels = desired_order),
    
    group = factor(
      group,
      levels = c(
        "Caudoviricetes",
        "Megaviricetes",
        "Prokaryotes",
        "Eukaryotes"
      )
    ),
    
    broad_type = case_when(
      group %in% c("Caudoviricetes", "Megaviricetes") ~ "Viruses",
      group %in% c("Prokaryotes", "Eukaryotes") ~ "Hosts",
      TRUE ~ NA_character_
    ),
    
    broad_type = factor(
      broad_type,
      levels = c("Viruses", "Hosts")
    ),
    
    box_col = factor(
      paste0(as.character(group), "_box"),
      levels = c(
        "Caudoviricetes_box",
        "Megaviricetes_box",
        "Prokaryotes_box",
        "Eukaryotes_box"
      )
    ),
    
    point_col = factor(
      paste0(as.character(group), "_pt"),
      levels = c(
        "Caudoviricetes_pt",
        "Megaviricetes_pt",
        "Prokaryotes_pt",
        "Eukaryotes_pt"
      )
    )
  )

ggplot(df_nonzero_filt, aes(x = category, y = total_norm_coverage)) +
  geom_boxplot(
    aes(
      group = interaction(category, group),
      fill = group,
      color = box_col,
      alpha = broad_type
    ),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  geom_jitter(
    aes(
      group = interaction(category, group),
      color = point_col,
      alpha = broad_type
    ),
    position = position_jitterdodge(
      jitter.width = 0.12,
      jitter.height = 0,
      dodge.width = 0.8
    ),
    size = 0.55
  ) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_fill_manual(
    values = group_cols,
    drop = FALSE
  ) +
  scale_color_manual(
    values = outline_point_cols,
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c(
      "Viruses" = 0.75,
      "Hosts"   = 0.35
    ),
    guide = "none"
  ) +
  labs(
    title = NULL,
    x = "Size fraction",
    y = "Total normalized coverage",
    fill = "Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 9)
  )


## Add bar plot for zero-coverage samples
viral_run_fraction <- votu_cov_tax_sw %>%
  filter(category %in% desired_order) %>%
  distinct(Run, category)

prok_run_fraction <- prok_cov_sw %>%
  filter(category %in% desired_order) %>%
  distinct(Run, category)

euk_run_fraction <- euk_cov_sw %>%
  filter(category %in% desired_order) %>%
  distinct(Run, category)

complete_grid_virus <- viral_run_fraction %>%
  tidyr::crossing(
    group = c("Caudoviricetes", "Megaviricetes")
  )

complete_grid_prok <- prok_run_fraction %>%
  mutate(group = "Prokaryotes")

complete_grid_euk <- euk_run_fraction %>%
  mutate(group = "Eukaryotes")

complete_grid <- bind_rows(
  complete_grid_virus,
  complete_grid_prok,
  complete_grid_euk
)

zero_summary <- df_complete %>%
  group_by(category, group) %>%
  summarize(
    n_zero = sum(total_norm_coverage == 0),
    n_total = n(),
    prop_zero = n_zero / n_total,
    .groups = "drop"
  ) %>%
  mutate(
    category = factor(category, levels = desired_order),
    group = factor(
      group,
      levels = c(
        "Caudoviricetes",
        "Megaviricetes",
        "Prokaryotes",
        "Eukaryotes"
      )
    )
  )

df_complete <- complete_grid %>%
  left_join(df_sum_all, by = c("Run", "category", "group")) %>%
  mutate(
    total_norm_coverage = replace_na(total_norm_coverage, 0),
    category = factor(category, levels = desired_order),
    group = factor(group, levels = all_groups)
  )

ggplot(zero_summary, aes(x = category, y = prop_zero, fill = group)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.72,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = group_cols,
    drop = FALSE
  ) +
  labs(
    x = "Size fraction",
    y = "Number of zeroes",
    fill = "Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1, size = 9)
  )