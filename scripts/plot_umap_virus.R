# Needs metadata and gv with columns: sample_id, virus_genome, class, norm_coverage, size_fraction, temp, env, latitude, and longitude from all seawater samples)
# Needs copernicus (modeled physicochemical parameters)

# ----------------------------------------------------------------------
# UMAP for viral OTUs
# Color = Latitude Zone, Symbol = Depth
# ----------------------------------------------------------------------

library(data.table)
library(fst)
library(vegan)
library(uwot)


## --- Import data ---
run_map <- metadata[, .(sample_id, depth, size_fraction, latitude, temp, env)] # <- start here

# Attach sample_id to viral table
gv_meta <- run_map[
  gv,
  on = .(sample_id = virus_sample_id),
  nomatch = 0L
]

gv_meta[, depth := suppressWarnings(as.numeric(depth))]
gv_meta[, latitude := suppressWarnings(as.numeric(latitude))]
gv_meta[, temp := suppressWarnings(as.numeric(temp))]

setDT(taxonomy_split)
tax_class <- taxonomy_split[, .(seq_name, class)]

gv_meta_tax <- tax_class[
  gv_meta,
  on = .(seq_name = virus_genome)
]

library(data.table)
library(vegan)
library(uwot)
library(ggplot2)
library(patchwork)

gv_dt <- as.data.table(gv_meta_tax)
setnames(gv_dt,"seq_name","virus_genome")

# Merge with Copernicus
library(data.table)

gv_dt <- as.data.table(gv_dt)
copernicus <- as.data.table(copernicus)

# check current missingness
sum(is.na(gv_dt$temp))

if (!"temp" %in% colnames(gv_dt)) {
  gv_dt[, temp := NA_real_]
}

cop_temp <- unique(
  copernicus[, .(sample_id, copernicus_Temperature_C_Daily)]
)

gv_dt <- merge(
  gv_dt,
  cop_temp,
  by = "sample_id",
  all.x = TRUE
)

gv_dt[
  is.na(temp),
  temp := copernicus_Temperature_C_Daily
]

gv_dt[, copernicus_Temperature_C_Daily := NULL]

sum(is.na(gv_dt$temp))

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

gv_dt <- gv_dt[
  class %in% c("Caudoviricetes", "Megaviricetes") &
    !is.na(sample_id) &
    !is.na(virus_genome) &
    !is.na(norm_coverage) &
    !is.na(latitude) &
    !is.na(depth) &
    !is.na(temp) &
    env == "seawater"
]

gv_dt[, depth := suppressWarnings(as.numeric(depth))]
gv_dt[, latitude := suppressWarnings(as.numeric(latitude))]
gv_dt[, temp := suppressWarnings(as.numeric(temp))]

gv_dt[, lat_zone := zone_of(latitude)]
gv_dt[, lat_zone := factor(
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
gv_dt[, layer := ifelse(depth > 200, "apho", "epi")]
gv_dt[, layer := factor(layer, levels = c("epi", "apho"))]

# quick checks
cat("\n=== FINAL FILTERED TABLE ===\n")
cat("Rows:", nrow(gv_dt), "\n")
cat("Unique sample_id:", uniqueN(gv_dt$sample_id), "\n")
cat("Unique vOTUs:", uniqueN(gv_dt$virus_genome), "\n")
cat("Missing temperatures:", sum(is.na(gv_dt$temp)), "\n")
cat("env values:\n")
print(table(gv_dt$env, useNA = "ifany"))

make_hellinger_umap_repro <- function(dat,
                                      target_class,
                                      top_votus = 10000,
                                      coverage_quantile = 0.10,
                                      n_neighbors = 100,
                                      min_dist = 0.75,
                                      seed = 123) {
  
  message("-----")
  message("Class: ", target_class)
  
  d <- copy(dat[class == target_class])
  
  # deterministic sample metadata
  sample_meta <- d[order(sample_id, virus_genome), .(
    latitude = first(latitude),
    depth = first(depth),
    temp = first(temp),
    layer = first(layer),
    size_fraction = first(size_fraction),
    lat_zone = first(lat_zone)
  ), by = sample_id]
  
  # aggregate duplicate sample-vOTU rows
  d <- d[, .(
    abundance = sum(norm_coverage, na.rm = TRUE)
  ), by = .(sample_id, virus_genome)]
  
  # deterministic ordering
  setorder(d, sample_id, virus_genome)
  
  message("Rows after aggregation: ", nrow(d))
  message("Unique samples before filtering: ", uniqueN(d$sample_id))
  message("Unique vOTUs before filtering: ", uniqueN(d$virus_genome))
  
  # keep top shared vOTUs
  votu_counts <- d[abundance > 0,
                   .(sharing_count = uniqueN(sample_id)),
                   by = virus_genome]
  setorder(votu_counts, -sharing_count, virus_genome)
  
  keep_votus <- votu_counts$virus_genome[1:min(top_votus, nrow(votu_counts))]
  d <- d[virus_genome %in% keep_votus]
  
  message("vOTUs retained after top-shared filter: ", length(keep_votus))
  
  # filter low-total-coverage samples
  run_totals <- d[, .(total_coverage = sum(abundance, na.rm = TRUE)), by = sample_id]
  setorder(run_totals, sample_id)
  
  threshold <- quantile(run_totals$total_coverage, coverage_quantile, na.rm = TRUE)
  keep_runs <- run_totals[total_coverage > threshold, sample_id]
  
  d <- d[sample_id %in% keep_runs]
  sample_meta <- sample_meta[sample_id %in% keep_runs]
  
  message("Samples retained after total-coverage filter: ", length(keep_runs))
  
  # wide matrix
  abund <- dcast(
    d,
    sample_id ~ virus_genome,
    value.var = "abundance",
    fun.aggregate = sum,
    fill = 0
  )
  
  # deterministic ordering again
  setorder(abund, sample_id)
  
  run_ids <- abund$sample_id
  abund[, sample_id := NULL]
  mat <- as.matrix(abund)
  
  # remove ultra-rare columns and empty rows
  keep_cols <- colSums(mat > 0) > 1
  mat <- mat[, keep_cols, drop = FALSE]
  
  keep_rows <- rowSums(mat) > 0
  mat <- mat[keep_rows, , drop = FALSE]
  run_ids <- run_ids[keep_rows]
  
  sample_meta <- sample_meta[match(run_ids, sample_meta$sample_id)]
  
  # final deterministic rownames
  rownames(mat) <- run_ids
  
  message("Final matrix: ", nrow(mat), " samples x ", ncol(mat), " vOTUs")
  
  # Hellinger transform
  mat_hell <- vegan::decostand(mat, method = "hellinger")
  
  # reproducible UMAP
  set.seed(seed)
  umap_res <- uwot::umap(
    mat_hell,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = "euclidean",
    scale = FALSE,
    verbose = TRUE,
    ret_model = FALSE
  )
  
  out <- data.table(
    sample_id = run_ids,
    UMAP1 = umap_res[, 1],
    UMAP2 = umap_res[, 2]
  )
  
  out <- merge(out, sample_meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  out[, class := target_class]
  
  # keep deterministic order in output too
  setorder(out, sample_id)
  
  list(
    umap = out,
    matrix = mat,
    matrix_hellinger = mat_hell,
    params = list(
      top_votus = top_votus,
      coverage_quantile = coverage_quantile,
      n_neighbors = n_neighbors,
      min_dist = min_dist,
      seed = seed
    )
  )
}

res_caudo_hell_norm <- make_hellinger_umap_repro(
  gv_dt,
  target_class = "Caudoviricetes",
  top_votus = 10000,
  coverage_quantile = 0.10,
  n_neighbors = 100,
  min_dist = 0.85,
  seed = 123
)

res_mega_hell_norm <- make_hellinger_umap_repro(
  gv_dt,
  target_class = "Megaviricetes",
  top_votus = 10000,
  coverage_quantile = 0.10,
  n_neighbors = 100,
  min_dist = 0.85,
  seed = 123
)

mega_far <- res_mega_hell_norm$umap[UMAP1 <= -40]

mega_no_far <- res_mega_hell_norm$umap[!sample_id %in% mega_far$sample_id]
mega_no_far_res <- list(umap = mega_no_far)


## Plotting
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

p_caudo_lat <- plot_umap_lat_outline(res_caudo_hell_norm, "Caudoviricetes", alpha_apho = 0.7)
p_mega_lat  <- plot_umap_lat_outline(mega_flip_res, "Megaviricetes", alpha_apho = 0.7)

patchwork::wrap_plots(p_caudo_lat, p_mega_lat, ncol = 2)

p_caudo_temp <- plot_umap_temp_outline(res_caudo_hell_norm, "Caudoviricetes", alpha_apho = 0.7)
p_mega_temp  <- plot_umap_temp_outline(mega_flip_res, "Megaviricetes", alpha_apho = 0.7)

patchwork::wrap_plots(p_caudo_temp, p_mega_temp, ncol = 2)

## Flipping both axes if necessary
mega_flip_res <- mega_no_far_res

mega_flip_res$umap[, UMAP1 := -UMAP1]
mega_flip_res$umap[, UMAP2 := -UMAP2]


