# Needs metadata and gp (for proks) and ge (for euks) with columns: sample_id, otu_id, norm_coverage, size_fraction, temp, env, latitude, and longitude from all seawater samples)
# Needs copernicus (modeled physicochemical parameters)

# ----------------------------------------------------------------------
# UMAP for prokaryotic and eukaryotic OTUs
# Color = Latitude Zone, Symbol = Depth
# ----------------------------------------------------------------------

## For prokaryotes
run_map <- metadata[, .(sample_id, depth, size_fraction, latitude, temp, env)] # <- start here

# Attach sample_id to prok table
gp_meta <- run_map[
  gp,
  on = .(sample_id = prok_sample_id),
  nomatch = 0L
]

gp_meta[, depth := suppressWarnings(as.numeric(depth))]
gp_meta[, latitude := suppressWarnings(as.numeric(latitude))]
gp_meta[, temp := suppressWarnings(as.numeric(temp))]

library(data.table)
library(vegan)
library(uwot)
library(ggplot2)
library(patchwork)

# Merge with Copernicus
gp_dt <- gp_meta

gp_dt <- as.data.table(gp_dt)
copernicus <- as.data.table(copernicus)

# check current missingness
sum(is.na(gp_dt$temp))

if (!"temp" %in% colnames(gp_dt)) {
  gp_dt[, temp := NA_real_]
}

cop_temp <- unique(
  copernicus[, .(sample_id, copernicus_Temperature_C_Daily)]
)

gp_dt <- merge(
  gp_dt,
  cop_temp,
  by = "sample_id",
  all.x = TRUE
)

gp_dt[
  is.na(temp),
  temp := copernicus_Temperature_C_Daily
]

gp_dt[, copernicus_Temperature_C_Daily := NULL]

sum(is.na(gp_dt$temp))

zone_of <- function(lat) {
  data.table::fcase(
    is.na(lat), NA_character_,
    lat <= -60, "south_polar",
    lat <  -45, "south_subpolar",
    lat <  -15, "south_subtropical",
    lat <   15, "tropical",
    lat <   45, "north_subtropical",
    lat <   60, "north_subpolar",
    lat >=  60, "north_polar",
    default = NA_character_
  )
}

# ----------------------------------------------------------------------
# Clean final table for rerunning UMAPs
# Keep only seawater samples and remove rows still missing temperature
# ----------------------------------------------------------------------

gp_dt <- gp_dt[
    !is.na(sample_id) &
    !is.na(sequence) &
    !is.na(norm_coverage) &
    !is.na(latitude) &
    !is.na(depth) &
    !is.na(temp) &
    env == "seawater"
]
gp_dt[, depth := suppressWarnings(as.numeric(depth))]
gp_dt[, latitude := suppressWarnings(as.numeric(latitude))]
gp_dt[, temp := suppressWarnings(as.numeric(temp))]

gp_dt[, lat_zone := zone_of(latitude)]
gp_dt[, lat_zone := factor(
  lat_zone,
  levels = c(
    "south_polar",
    "south_subpolar",
    "south_subtropical",
    "tropical",
    "north_subtropical",
    "north_subpolar",
    "north_polar"
  )
)]

# define epi / apho from depth
gp_dt[, layer := ifelse(depth > 200, "apho", "epi")]
gp_dt[, layer := factor(layer, levels = c("epi", "apho"))]

# quick checks
cat("\n=== FINAL FILTERED TABLE ===\n")
cat("Rows:", nrow(gp_dt), "\n")
cat("Unique sample_id:", uniqueN(gp_dt$sample_id), "\n")
cat("Unique OTUs:", uniqueN(gp_dt$sequence), "\n")
cat("Missing temperatures:", sum(is.na(gp_dt$temp)), "\n")
cat("env values:\n")
print(table(gp_dt$env, useNA = "ifany"))

## Now for eukaryotes
run_map <- metadata[, .(sample_id, depth, size_fraction, latitude, temp, env)] # <- start here

# Attach sample_id to euk table
ge_meta <- run_map[
  ge,
  on = .(sample_id = euk_sample_id),
  nomatch = 0L
]

ge_meta[, depth := suppressWarnings(as.numeric(depth))]
ge_meta[, latitude := suppressWarnings(as.numeric(latitude))]
ge_meta[, temp := suppressWarnings(as.numeric(temp))]

# Merge with Copernicus
ge_dt <- ge_meta

ge_dt <- as.data.table(ge_dt)
copernicus <- as.data.table(copernicus)

# check current missingness
sum(is.na(ge_dt$temp))

if (!"temp" %in% colnames(ge_dt)) {
  ge_dt[, temp := NA_real_]
}

cop_temp <- unique(
  copernicus[, .(sample_id, copernicus_Temperature_C_Daily)]
)

ge_dt <- merge(
  ge_dt,
  cop_temp,
  by = "sample_id",
  all.x = TRUE
)

ge_dt[
  is.na(temp),
  temp := copernicus_Temperature_C_Daily
]

ge_dt[, copernicus_Temperature_C_Daily := NULL]

sum(is.na(ge_dt$temp))

# ----------------------------------------------------------------------
# Clean final table for rerunning UMAPs
# Keep only seawater samples and remove rows still missing temperature
# ----------------------------------------------------------------------

ge_dt <- ge_dt[
    !is.na(sample_id) &
    !is.na(sequence) &
    !is.na(norm_coverage) &
    !is.na(latitude) &
    !is.na(depth) &
    !is.na(temp) &
    env == "seawater"
]
ge_dt[, depth := suppressWarnings(as.numeric(depth))]
ge_dt[, latitude := suppressWarnings(as.numeric(latitude))]
ge_dt[, temp := suppressWarnings(as.numeric(temp))]

ge_dt[, lat_zone := zone_of(latitude)]
ge_dt[, lat_zone := factor(
  lat_zone,
  levels = c(
    "south_polar",
    "south_subpolar",
    "south_subtropical",
    "tropical",
    "north_subtropical",
    "north_subpolar",
    "north_polar"
  )
)]

# define epi / apho from depth
ge_dt[, layer := ifelse(depth > 200, "apho", "epi")]
ge_dt[, layer := factor(layer, levels = c("epi", "apho"))]

# quick checks
cat("\n=== FINAL FILTERED TABLE ===\n")
cat("Rows:", nrow(ge_dt), "\n")
cat("Unique sample_id:", uniqueN(ge_dt$sample_id), "\n")
cat("Unique OTUs:", uniqueN(ge_dt$taxon), "\n")
cat("Missing temperatures:", sum(is.na(ge_dt$temp)), "\n")
cat("env values:\n")
print(table(ge_dt$env, useNA = "ifany"))


# Function to make the UMAP 
make_hellinger_umap_repro <- function(dat,
                                      target_class = NULL,
                                      top_features = 10000,
                                      coverage_quantile = 0.10,
                                      n_neighbors = 100,
                                      min_dist = 0.75,
                                      seed = 123,
                                      feature_col = "virus_genome",
                                      abundance_col = "norm_coverage") {

  d <- data.table::copy(data.table::as.data.table(dat))

  # Optionally filter by class
  if (!is.null(target_class)) {
    if (!"class" %in% names(d)) {
      stop("target_class was supplied, but dat has no 'class' column.")
    }

    message("-----")
    message("Class: ", target_class)

    d <- d[class == target_class]
  } else {
    message("-----")
    message("Using all observations; no class filter applied.")
  }

  if (nrow(d) == 0L) {
    stop("No rows remain after filtering.")
  }

  # Check required columns
  required_cols <- c(
    "sample_id",
    feature_col,
    abundance_col,
    "latitude",
    "depth",
    "temp",
    "layer",
    "size_fraction",
    "lat_zone"
  )

  missing_cols <- setdiff(required_cols, names(d))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Standardize feature and abundance column names internally
  d[, feature_id := get(feature_col)]
  d[, abundance_value := get(abundance_col)]

  # Deterministic sample metadata
  data.table::setorder(d, sample_id, feature_id)

  sample_meta <- d[, .(
    latitude = first(latitude),
    depth = first(depth),
    temp = first(temp),
    layer = first(layer),
    size_fraction = first(size_fraction),
    lat_zone = first(lat_zone)
  ), by = sample_id]

  # Aggregate duplicate sample-feature rows
  d <- d[, .(
    abundance = sum(abundance_value, na.rm = TRUE)
  ), by = .(sample_id, feature_id)]

  data.table::setorder(d, sample_id, feature_id)

  message("Rows after aggregation: ", nrow(d))
  message("Unique samples before filtering: ", uniqueN(d$sample_id))
  message("Unique features before filtering: ", uniqueN(d$feature_id))

  # Keep top shared features
  feature_counts <- d[
    abundance > 0,
    .(sharing_count = uniqueN(sample_id)),
    by = feature_id
  ]

  data.table::setorder(feature_counts, -sharing_count, feature_id)

  n_keep <- min(top_features, nrow(feature_counts))
  keep_features <- feature_counts$feature_id[seq_len(n_keep)]

  d <- d[feature_id %in% keep_features]

  message("Features retained after top-shared filter: ", length(keep_features))

  # Filter low-total-abundance samples
  run_totals <- d[, .(
    total_abundance = sum(abundance, na.rm = TRUE)
  ), by = sample_id]

  data.table::setorder(run_totals, sample_id)

  threshold <- stats::quantile(
    run_totals$total_abundance,
    coverage_quantile,
    na.rm = TRUE,
    names = FALSE
  )

  keep_runs <- run_totals[
    total_abundance > threshold,
    sample_id
  ]

  d <- d[sample_id %in% keep_runs]
  sample_meta <- sample_meta[sample_id %in% keep_runs]

  message("Samples retained after abundance filter: ", length(keep_runs))

  if (nrow(d) == 0L) {
    stop("No observations remain after the abundance filter.")
  }

  # Wide matrix
  abund <- data.table::dcast(
    d,
    sample_id ~ feature_id,
    value.var = "abundance",
    fun.aggregate = sum,
    fill = 0
  )

  data.table::setorder(abund, sample_id)

  run_ids <- abund$sample_id
  abund[, sample_id := NULL]

  mat <- as.matrix(abund)
  storage.mode(mat) <- "numeric"

  # Remove features found in only one sample
  keep_cols <- colSums(mat > 0) > 1
  mat <- mat[, keep_cols, drop = FALSE]

  if (ncol(mat) == 0L) {
    stop("No features occur in more than one retained sample.")
  }

  # Remove empty samples
  keep_rows <- rowSums(mat) > 0
  mat <- mat[keep_rows, , drop = FALSE]
  run_ids <- run_ids[keep_rows]

  sample_meta <- sample_meta[match(run_ids, sample_meta$sample_id)]

  rownames(mat) <- run_ids

  message(
    "Final matrix: ",
    nrow(mat),
    " samples x ",
    ncol(mat),
    " features"
  )

  if (nrow(mat) < 3L) {
    stop("At least three non-empty samples are needed for UMAP.")
  }

  # n_neighbors must be smaller than the number of observations
  effective_neighbors <- min(n_neighbors, nrow(mat) - 1L)

  if (effective_neighbors != n_neighbors) {
    message(
      "Reducing n_neighbors from ",
      n_neighbors,
      " to ",
      effective_neighbors,
      " because only ",
      nrow(mat),
      " samples remain."
    )
  }

  # Hellinger transformation
  mat_hell <- vegan::decostand(
    mat,
    method = "hellinger"
  )

  # Reproducible UMAP
  set.seed(seed)

  umap_res <- uwot::umap(
    mat_hell,
    n_neighbors = effective_neighbors,
    min_dist = min_dist,
    metric = "euclidean",
    scale = FALSE,
    verbose = TRUE,
    ret_model = FALSE
  )

  out <- data.table::data.table(
    sample_id = run_ids,
    UMAP1 = umap_res[, 1],
    UMAP2 = umap_res[, 2]
  )

  out <- merge(
    out,
    sample_meta,
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  if (!is.null(target_class)) {
    out[, class := target_class]
  }

  data.table::setorder(out, sample_id)

  list(
    umap = out,
    matrix = mat,
    matrix_hellinger = mat_hell,
    params = list(
      target_class = target_class,
      feature_col = feature_col,
      abundance_col = abundance_col,
      top_features = top_features,
      coverage_quantile = coverage_quantile,
      n_neighbors_requested = n_neighbors,
      n_neighbors_used = effective_neighbors,
      min_dist = min_dist,
      seed = seed
    )
  )
}

# Plot
plot_umap_lat_outline <- function(res_obj, title_text, alpha_apho = 0.7) {
  
  umap_df <- as.data.frame(res_obj$umap)
  
  ggplot(umap_df, aes(UMAP1, UMAP2)) +
    # epipelagic first: opaque circles
    geom_point(
      data = umap_df[umap_df$layer == "epi", ],
      aes(color = lat_zone),
      shape = 16,
      size = 1,
      alpha = 1
    ) +
    # aphotic second: semi-transparent triangles with black outline
    geom_point(
      data = umap_df[umap_df$layer == "apho", ],
      aes(fill = lat_zone),
      shape = 24,
      size = 1.5,
      alpha = alpha_apho,
      color = "black",
      stroke = 0.35
    ) +
    scale_color_manual(
      values = c(
        south_polar = "#5b72a5",
        south_subpolar = "#14B8A6",
        south_subtropical = "#A5D6A7",
        tropical = "#D73027",
        north_subtropical = "#388E3C",
        north_subpolar = "#019299",
        north_polar = "#1E3A8A"
      ),
      drop = FALSE
    ) +
    scale_fill_manual(
      values = c(
        south_polar = "#5b72a5",
        south_subpolar = "#14B8A6",
        south_subtropical = "#A5D6A7",
        tropical = "#D73027",
        north_subtropical = "#388E3C",
        north_subpolar = "#019299",
        north_polar = "#1E3A8A"
      ),
      drop = FALSE
    ) +
    labs(
      title = title_text,
      x = "UMAP1",
      y = "UMAP2",
      color = "Latitude zone",
      fill = "Latitude zone"
    ) +
    guides(fill = "none") +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

plot_umap_temp_outline <- function(res_obj, title_text, alpha_apho = 0.7) {
  
  umap_df <- as.data.frame(res_obj$umap)
  
  ggplot(umap_df, aes(UMAP1, UMAP2)) +
    # epipelagic first: opaque circles
    geom_point(
      data = umap_df[umap_df$layer == "epi", ],
      aes(color = temp),
      shape = 16,
      size = 1,
      alpha = 1
    ) +
    # aphotic second: semi-transparent triangles with black outline
    geom_point(
      data = umap_df[umap_df$layer == "apho", ],
      aes(fill = temp),
      shape = 24,
      size = 1.5,
      alpha = alpha_apho,
      color = "black",
      stroke = 0.35
    ) +
    scale_color_gradientn(
      colors = c("darkblue", "deepskyblue3", "yellow", "orange", "red"),
      na.value = "grey80"
    ) +
    scale_fill_gradientn(
      colors = c("darkblue", "deepskyblue3", "yellow", "orange", "red"),
      na.value = "grey80"
    ) +
    labs(
      title = title_text,
      x = "UMAP1",
      y = "UMAP2",
      color = "Temperature",
      fill = "Temperature"
    ) +
    guides(fill = "none") +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

# Usage

prok_umap <- make_hellinger_umap_repro(
  dat = gp_dt,
  feature_col = "sequence",
  abundance_col = "norm_coverage",
  top_features = 10000,
  coverage_quantile = 0.10,
  n_neighbors = 100,
  min_dist = 0.75,
  seed = 123
)

prok_umap_lat <- plot_umap_lat_outline(prok_flip, "Prokaryotes", alpha_apho = 0.7)
prok_umap_temp <- plot_umap_temp_outline(prok_flip, "Prokaryotes", alpha_apho = 0.7)
patchwork::wrap_plots(prok_umap_lat, prok_umap_temp, ncol = 2)



euk_umap <- make_hellinger_umap_repro(
  dat = ge_dt,
  feature_col = "taxon",
  abundance_col = "norm_coverage",
  top_features = 10000,
  coverage_quantile = 0.10,
  n_neighbors = 100,
  min_dist = 0.75,
  seed = 123
)

euk_umap_lat <- plot_umap_lat_outline(euk_flip, "Eukaryotes", alpha_apho = 0.7)
euk_umap_temp <- plot_umap_temp_outline(euk_flip, "Eukaryotes", alpha_apho = 0.7)
patchwork::wrap_plots(euk_umap_lat, euk_umap_temp, ncol = 2)


prok_flip <- prok_umap
prok_flip$umap[, UMAP1 := -UMAP1]
#prok_flip$umap[, UMAP2 := -UMAP2]

euk_flip <- euk_umap
euk_flip$umap[, UMAP1 := -UMAP1]

# Save files
#saveRDS(prok_umap, "~/Documents/global_metagenomes/figures/prok_hellinger_umap_coords_30jul2026.rds")
#saveRDS(euk_umap, "~/Documents/global_metagenomes/figures/euk_hellinger_umap_coords_30jul2026.rds")



umaps_all <- wrap_plots(
  p_caudo_lat, p_caudo_temp,
  prok_umap_lat, prok_umap_temp,
  p_mega_lat, p_mega_temp,
  euk_umap_lat, euk_umap_temp,
  ncol = 4
)

