# Needs votu_lat_dist with columns: virus_genome, final_tag (biogeographic classification)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpattern)
library(viridisLite)
library(grid)

# ============================================================
# 3) Original mono-/bihemispherical counts
# ============================================================

geo_plot_df <- votu_lat_dist %>%
  filter(
    !is.na(source),
    !is.na(final_geo),
    !is.na(excl_status),
    !is.na(class),
    final_geo != "discordant",
    class %in% c("Caudoviricetes", "Megaviricetes")
  ) %>%
  mutate(
    source = as.character(source),
    final_geo = as.character(final_geo),
    excl_status = str_to_lower(as.character(excl_status)),
    class = as.character(class)
  ) %>%
  filter(
    excl_status %in% c("exclusive", "preferential")
  ) %>%
  count(
    class,
    source,
    final_geo,
    excl_status,
    name = "n"
  ) %>%
  group_by(class, source, final_geo) %>%
  complete(
    excl_status = c("exclusive", "preferential"),
    fill = list(n = 0)
  ) %>%
  ungroup()


# ============================================================
# 4) Cosmopolitan counts
#
# Cosmopolitan IDs are NOT removed from the geographic data.
# They are counted separately here.
# ============================================================
tax_lookup <- taxonomy_split %>%
  transmute(
    virus_genome = seq_name,
    class = class
  ) %>%
  distinct(virus_genome, .keep_all = TRUE)

cosmo_plot_df <- tibble(
  virus_genome = unique(cosmopolitan_detection_ids)
) %>%
  left_join(tax_lookup, by = "virus_genome") %>%
  filter(
    !is.na(class),
    class %in% c("Caudoviricetes", "Megaviricetes")
  ) %>%
  distinct(virus_genome, class) %>%
  count(class, name = "n") %>%
  complete(
    class = c("Caudoviricetes", "Megaviricetes"),
    fill = list(n = 0)
  ) %>%
  mutate(
    source = "cosmopolitan",
    final_geo = "cosmopolitan",
    excl_status = "cosmopolitan"
  )


# ============================================================
# 5) Combine
# ============================================================

plot_df <- bind_rows(
  geo_plot_df,
  cosmo_plot_df
)

plot_df <- plot_df %>%
  mutate(
    source = if_else(
      final_geo == "tropical",
      "monohemispherical",
      as.character(source)
    )
  )

# ============================================================
# 6) Order the ACTUAL labels already present in final_geo
#
# No geographic labels are renamed.
# ============================================================

plot_df <- plot_df %>%
  mutate(
    source_rank = case_when(
      source == "monohemispherical" ~ 1,
      source == "bihemispherical"   ~ 2,
      source == "cosmopolitan"      ~ 3,
      TRUE                          ~ 99
    ),

    geo_rank = case_when(
      # Monohemispherical: requested top-to-bottom order
      source == "monohemispherical" &
        final_geo == "polar_north" ~ 1,

      source == "monohemispherical" &
        final_geo == "subpolar_north" ~ 2,

      source == "monohemispherical" &
        final_geo == "high_latitude_north" ~ 3,

      source == "monohemispherical" &
        final_geo == "subtropical_north" ~ 4,  

      source == "monohemispherical" &
        final_geo == "tropical" ~ 5,

      source == "monohemispherical" &
        final_geo == "subtropical_south" ~ 6,

      source == "monohemispherical" &
        final_geo == "high_latitude_south" ~ 7,

      source == "monohemispherical" &
        final_geo == "subpolar_south" ~ 8,

      source == "monohemispherical" &
        final_geo == "polar_south" ~ 9,

      # Bihemispherical:
      # works with either bipolar/bisubpolar/etc.
      # or bipolar/subpolar/high_latitude/subtropical
      source == "bihemispherical" &
        str_detect(final_geo, "^(bipolar|polar)$") ~ 1,

      source == "bihemispherical" &
        str_detect(final_geo, "^(bisubpolar|subpolar)$") ~ 2,

      source == "bihemispherical" &
        str_detect(
          final_geo,
          "^(bi_?high_latitude|high_latitude)$"
        ) ~ 3,

      source == "bihemispherical" &
        str_detect(
          final_geo,
          "^(bisubtropical|subtropical)$"
        ) ~ 4,

      source == "cosmopolitan" &
        final_geo == "cosmopolitan" ~ 1,

      # Preserve any unexpected categories rather than turning them into NA
      TRUE ~ 99
    )
  )


# Derive factor levels from the exact labels in the table
ordered_levels <- plot_df %>%
  distinct(source, final_geo, source_rank, geo_rank) %>%
  arrange(source_rank, geo_rank, final_geo) %>%
  pull(final_geo) %>%
  unique()

# Reverse because the first discrete y level is plotted at the bottom
plot_df <- plot_df %>%
  mutate(
    final_geo = factor(
      final_geo,
      levels = rev(ordered_levels)
    ),

    source = factor(
      source,
      levels = c(
        "monohemispherical",
        "bihemispherical",
        "cosmopolitan"
      )
    ),

    class = factor(
      class,
      levels = c(
        "Caudoviricetes",
        "Megaviricetes"
      )
    ),

    excl_status = factor(
      excl_status,
      levels = c(
        "exclusive",
        "preferential",
        "cosmopolitan"
      )
    )
  )


# Optional diagnostic: should be empty
plot_df %>%
  filter(
    is.na(source) |
      is.na(class) |
      is.na(final_geo)
  )


# ============================================================
# 7) Colors derived from the actual labels
# ============================================================

mono_geos <- plot_df %>%
  filter(source == "monohemispherical") %>%
  distinct(final_geo, geo_rank) %>%
  arrange(geo_rank, final_geo) %>%
  pull(final_geo) %>%
  as.character()

bi_geos <- plot_df %>%
  filter(source == "bihemispherical") %>%
  distinct(final_geo, geo_rank) %>%
  arrange(geo_rank, final_geo) %>%
  pull(final_geo) %>%
  as.character()

mono_cols <- viridisLite::plasma(
  length(mono_geos),
  begin = 0.10,
  end = 0.90
)

bi_cols <- viridisLite::plasma(
  length(bi_geos),
  begin = 0.20,
  end = 1.00
)

pal <- c(
  setNames(mono_cols, mono_geos),
  setNames(bi_cols, bi_geos),
  cosmopolitan = "purple4"
)


# ============================================================
# 8) Plot
#
# Columns = viral classes
# Rows = geographic groups
#
# facet_grid(scales = "free") means:
# - one shared x-range for all Caudoviricetes panels
# - one shared x-range for all Megaviricetes panels
# - a separate y-category scale for each source row
# ============================================================

p <- ggplot(
  plot_df,
  aes(
    x = n,
    y = final_geo,
    fill = final_geo,
    pattern = excl_status
  )
) +
  ggpattern::geom_col_pattern(
    position = position_dodge2(
      width = 0.8,
      preserve = "single"
    ),
    width = 0.7,
    color = "black",
    linewidth = 0.2,
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.25,
    pattern_spacing = 0.03,
    pattern_size = 0.008
  ) +
  facet_grid(
    rows = vars(source),
    cols = vars(class),
    scales = "free",
    space = "free_y",
    labeller = labeller(
      source = c(
        monohemispherical = "Monohemispherical",
        bihemispherical   = "Bihemispherical",
        cosmopolitan      = "Cosmopolitan"
      )
    )
  ) +
  scale_fill_manual(
    values = pal,
    guide = "none"
  ) +
  scale_pattern_manual(
    values = c(
      exclusive    = "none",
      preferential = "circle",
      cosmopolitan = "none"
    ),
    breaks = c(
      "exclusive",
      "preferential"
    ),
    labels = c(
      exclusive    = "Exclusive",
      preferential = "Preferential"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_y_discrete(
    drop = TRUE
  ) +
  theme_bw() +
  labs(
    x = "Number of taxa",
    y = NULL,
    pattern = NULL
  ) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.spacing.x = unit(1, "lines"),
    panel.spacing.y = unit(1, "lines")
  )

p
