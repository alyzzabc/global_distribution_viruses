# Needs votu_normcov_metadata_epi with columns: sample_id, virus_genome, class, norm_coverage, latitude, and longitude from epipelagic seawater samples
# Needs votu_lat_dist with columns: virus_genome, final_tag (biogeographic classification)

plot_otu_group_map <- function(derep_df,
                               tag_df,
                               final_tag_value,
                               coverage_col = "norm_coverage",
                               otu_col_derep = "virus_genome",
                               otu_col_tag = "virus_genome",
                               tag_col = "final_tag",
                               run_col = "sample_id",
                               lat_col = "latitude",
                               lon_col = "longitude",
                               class_col = "class",
                               class_value = NULL,
                               aggregate_fun = sum,
                               world_fill = "grey95",
                               world_line = "grey70",
                               ring_color = "grey30",
                               point_color = "steelblue",
                               ring_size = 1,
                               show_zeroes = FALSE,
                               size_limits = NULL,
                               size_breaks = NULL,
                               size_range = c(1, 8)) {
  
  library(dplyr)
  library(ggplot2)
  
  selected_ids <- tag_df %>%
    dplyr::filter(.data[[tag_col]] == final_tag_value) %>%
    dplyr::pull(all_of(otu_col_tag)) %>%
    unique()
  
  if (length(selected_ids) == 0) {
    stop("No taxon_id found for final_tag = ", final_tag_value)
  }
  
  derep_df_filt <- derep_df
  if (!is.null(class_value)) {
    derep_df_filt <- derep_df_filt %>%
      dplyr::filter(.data[[class_col]] %in% class_value)
  }
  
  all_runs <- derep_df_filt %>%
    dplyr::select(all_of(c(run_col, lat_col, lon_col))) %>%
    dplyr::rename(
      sample_id = all_of(run_col),
      latitude = all_of(lat_col),
      longitude = all_of(lon_col)
    ) %>%
    dplyr::distinct() %>%
    dplyr::filter(!is.na(latitude), !is.na(longitude))
  
  matched_df <- derep_df_filt %>%
    dplyr::filter(.data[[otu_col_derep]] %in% selected_ids) %>%
    dplyr::filter(!is.na(.data[[lat_col]]), !is.na(.data[[lon_col]]))
  
  detected_runs <- matched_df %>%
    dplyr::group_by(
      sample_id = .data[[run_col]],
      latitude = .data[[lat_col]],
      longitude = .data[[lon_col]]
    ) %>%
    dplyr::summarise(
      plot_coverage = aggregate_fun(.data[[coverage_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(plot_coverage > 0) %>%
    dplyr::arrange(dplyr::desc(plot_coverage))
  
  absent_runs <- all_runs %>%
    dplyr::anti_join(
      detected_runs %>% dplyr::select(sample_id),
      by = "sample_id"
    )
  
  world <- ggplot2::map_data("world")
  
  p <- ggplot() +
    geom_polygon(
      data = world,
      aes(x = long, y = lat, group = group),
      fill = world_fill,
      color = world_line,
      linewidth = 0.2
    ) +
    geom_hline(
      yintercept = c(-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90),
      linetype = "solid",
      color = "grey90",
      linewidth = 0.3
    ) +
    geom_vline(
      xintercept = c(-120, -60, 0, 60, 120),
      linetype = "solid",
      color = "grey90",
      linewidth = 0.3
    ) +
    geom_point(
      data = detected_runs,
      aes(
        x = longitude,
        y = latitude,
        size = plot_coverage
      ),
      color = point_color,
      alpha = 0.5,
      stroke = 0
    ) +
    coord_quickmap(
      xlim = c(-180.5, 180.5),
      ylim = c(-90.5, 90.5),
      expand = TRUE
    ) +
    scale_x_continuous(breaks = c(-120, 120)) +
    scale_y_continuous(breaks = c(-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90)) +
    scale_size_continuous(
      name = coverage_col,
      limits = size_limits,
      breaks = size_breaks,
      range = size_range
    ) +
    scale_y_continuous(breaks = c(-90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90)) +
    labs(
      title = if (is.null(class_value)) {
        paste0("Distribution of ", final_tag_value, " (all classes)")
      } else {
        paste0("Distribution of ", final_tag_value, " | class: ", paste(class_value, collapse = ", "))
      },
      x = "Longitude",
      y = "Latitude"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      panel.border = element_blank()
    )
  
  if (show_zeroes) {
    p <- p +
      geom_point(
        data = absent_runs,
        aes(
          x = longitude,
          y = latitude
        ),
        shape = 1,
        color = ring_color,
        stroke = 0.3,
        size = ring_size
      )
  }
  
  p
}

get_shared_size_limits <- function(derep_df,
                                   tag_df,
                                   final_tag_values,
                                   coverage_col = "norm_coverage",
                                   otu_col_derep = "virus_genome",
                                   otu_col_tag = "virus_genome",
                                   tag_col = "final_tag",
                                   run_col = "sample_id",
                                   lat_col = "latitude",
                                   lon_col = "longitude",
                                   class_col = "class",
                                   class_value = NULL,
                                   aggregate_fun = sum) {
  
  library(dplyr)
  
  selected_ids <- tag_df %>%
    dplyr::filter(.data[[tag_col]] %in% final_tag_values) %>%
    dplyr::pull(all_of(otu_col_tag)) %>%
    unique()
  
  if (length(selected_ids) == 0) {
    stop("No OTU IDs found for requested final_tag_values")
  }
  
  derep_df_filt <- derep_df
  if (!is.null(class_value)) {
    derep_df_filt <- derep_df_filt %>%
      dplyr::filter(.data[[class_col]] %in% class_value)
  }
  
  matched_df <- derep_df_filt %>%
    dplyr::filter(.data[[otu_col_derep]] %in% selected_ids) %>%
    dplyr::filter(!is.na(.data[[lat_col]]), !is.na(.data[[lon_col]]))
  
  detected_runs <- matched_df %>%
    dplyr::group_by(
      sample_id = .data[[run_col]],
      latitude = .data[[lat_col]],
      longitude = .data[[lon_col]]
    ) %>%
    dplyr::summarise(
      plot_coverage = aggregate_fun(.data[[coverage_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(plot_coverage > 0)
  
  max_cov <- max(detected_runs$plot_coverage, na.rm = TRUE)
  
  c(0, max_cov)
}

caudo_limits <- get_shared_size_limits(
  derep_df = votu_normcov_metadata_epi,
  tag_df = votu_lat_dist,
  final_tag_values = c("preferential_tropical", "preferential_high_latitude"),
  class_value = "Caudoviricetes"
)

caudo_breaks <- pretty(caudo_limits, n = 4)
#caudo_breaks <- caudo_breaks[caudo_breaks > 0]

p_caudo_trop <- plot_otu_group_map(
  derep_df = votu_normcov_metadata_epi,
  tag_df = votu_lat_dist,
  final_tag_value = "preferential_tropical",
  class_value = "Caudoviricetes",
  size_limits = caudo_limits,
  size_breaks = caudo_breaks,
  point_color = "purple",
  show_zeroes = TRUE
)

p_caudo_high <- plot_otu_group_map(
  derep_df = votu_normcov_metadata_epi,
  tag_df = votu_lat_dist,
  final_tag_value = "preferential_high_latitude",
  class_value = "Caudoviricetes",
  size_limits = caudo_limits,
  size_breaks = caudo_breaks,
  point_color = "purple",
  show_zeroes = TRUE
)

#p_caudo_trop
#p_caudo_high

library(patchwork)

p_caudo_trop / p_caudo_high

mega_limits <- get_shared_size_limits(
  derep_df = votu_normcov_metadata_epi,
  tag_df = votu_lat_dist,
  final_tag_values = c("preferential_tropical", "preferential_high_latitude"),
  class_value = "Megaviricetes"
)

mega_breaks <- pretty(mega_limits, n = 4)
#mega_breaks <- mega_breaks[mega_breaks > 0]

p_mega_trop <- plot_otu_group_map(
  derep_df = votu_normcov_metadata_epi,
  tag_df = votu_lat_dist,
  final_tag_value = "preferential_tropical",
  class_value = "Megaviricetes",
  size_limits = mega_limits,
  size_breaks = mega_breaks,
  point_color = "orange",
  show_zeroes = TRUE
)

p_mega_high <- plot_otu_group_map(
  derep_df = votu_normcov_metadata_epi,
  tag_df = votu_lat_dist,
  final_tag_value = "preferential_high_latitude",
  class_value = "Megaviricetes",
  size_limits = mega_limits,
  size_breaks = mega_breaks,
  point_color = "orange",
  show_zeroes = TRUE
)

#p_mega_trop
#p_mega_high

p_mega_trop / p_mega_high
