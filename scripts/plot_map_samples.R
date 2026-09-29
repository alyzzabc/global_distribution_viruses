library(data.table)
library(ggplot2)
library(sf)
library(rnaturalearth)

metadata <- fread(
  "/path/to/metadata",
  select = c(
    "Run", "manual_date", "manual_latitude", "manual_longitude",
    "manual_depth", "manual_env", "category", "manual_temp",
    "manual_chl", "manual_O2", "manual_nitrate", "study",
    "manual_size_text_strings"
  )
)

setDT(metadata)

# For this map, only seawater samples with coordinates are relevant
metadata <- metadata[
  manual_env == "seawater" &
    !is.na(manual_latitude) &
    !is.na(manual_longitude)
]

# Make size strings UTF-8 safe
metadata[, manual_size_text_strings := iconv(
  manual_size_text_strings,
  from = "",
  to = "UTF-8",
  sub = ""
)]

# Normalize common placeholders to NA
metadata[, manual_size_text_strings := fifelse(
  is.na(manual_size_text_strings) |
    trimws(manual_size_text_strings) %chin% c("", "missing", "#N/A", "NA", "N/A"),
  NA_character_,
  trimws(manual_size_text_strings)
)]

# Missing size values stay unique, so they are not assumed to be duplicates
metadata[, size_group := fcoalesce(
  manual_size_text_strings,
  paste0("__no_size__", .I)
)]

# Known dates collapse across time; missing dates stay unique
metadata[, time_group := fifelse(
  is.na(manual_date),
  paste0("__no_date__", .I),
  "__has_date__"
)]

# Deduplicate seawater metadata for the map
dedup_keys <- c(
  "manual_latitude",
  "manual_longitude",
  "manual_depth",
  "size_group",
  "time_group"
)

metadata[, rep_Run := Run[order(Run)][1], by = dedup_keys]

# replicate counts per representative Run
rep_counts <- metadata[, .(n_reps = .N), by = rep_Run]

# one row per representative Run
dedup_runs <- unique(metadata[Run == rep_Run])

# attach replicate counts
dedup_runs[rep_counts, n_reps := i.n_reps, on = .(rep_Run)]

# Convert depth and remove unknown/non-numeric depth before plotting
dedup_runs[, depth_m := suppressWarnings(as.numeric(manual_depth))]
dedup_runs <- dedup_runs[!is.na(depth_m)]

# Assign layer
dedup_runs[, layer := fifelse(depth_m < 200, "epi", "apho")]

env_check <- metadata[
  , .(
    n_rows = .N,
    n_env = uniqueN(manual_env),
    manual_env_values = paste(sort(unique(manual_env)), collapse = "; ")
  ),
  by = rep_Run
][n_env > 1]

env_check

# Aggregate per location/layer for plotting
pts <- dedup_runs[
  , .(n_runs_layer = sum(n_reps)),
  by = .(manual_longitude, manual_latitude, layer)
]


## Plot ##
max_cex <- 7
maxN <- max(pts$n_runs_layer, na.rm = TRUE)

pts[, size_plot := max_cex * log1p(n_runs_layer) / log1p(maxN)]

# draw big first so small are plotted on top
setorder(pts, -size_plot)

breaks_n <- c(1, 50, 100, 300)
breaks_n <- breaks_n[breaks_n <= maxN]

breaks_size <- max_cex * log1p(breaks_n) / log1p(maxN)

# keep only breaks that are not larger than the observed maximum
breaks_n <- breaks_n[breaks_n <= maxN]

breaks_size <- max_cex * log1p(breaks_n) / log1p(maxN)

# 5) Minimal graticule + labeled axes
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

lon_breaks <- c(-120, -60, 0, 60, 120)
lat_breaks <- c(-60, -45, -30, -15, 0, 15, 30, 45, 60)

label_lon <- function(x) {
  ifelse(x == 0, "0°",
         paste0(abs(x), "°", ifelse(x < 0, "W", "E")))
}

label_lat <- function(x) {
  ifelse(x == 0, "0°",
         paste0(abs(x), "°", ifelse(x < 0, "S", "N")))
}

grat <- st_graticule(
  lon = lon_breaks,
  lat = lat_breaks,
  crs = st_crs(world)
)

ggplot() +
  geom_sf(data = world, fill = "grey90", color = "black", linewidth = 0.2) +
  geom_sf(data = grat, color = "grey80", linewidth = 0.25) +
  geom_point(
    data = pts,
    aes(
      x = manual_longitude,
      y = manual_latitude,
      color = layer,
      size = size_plot
    ),
    shape = 16
  ) +
  scale_color_manual(
    values = c(epi = "blue", apho = "red"),
    name = NULL
  ) +
  scale_size_identity(
    name = "Samples per location/layer",
    breaks = breaks_size,
    labels = breaks_n,
    guide = "legend"
  ) +
  scale_x_continuous(breaks = lon_breaks, labels = label_lon) +
  scale_y_continuous(breaks = lat_breaks, labels = label_lat) +
  coord_sf(xlim = c(-180, 180), ylim = c(-80, 85), expand = FALSE, clip = "off") +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid  = element_blank(),
    axis.title  = element_text(color = "black"),
    axis.text   = element_text(color = "black"),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(10, 20, 10, 20)
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    size  = guide_legend(override.aes = list(color = "grey20"))
  )

## Make the pies

# category presence per deduplicated sample
# one row per original run with its category and deduplicated group
run_cats <- unique(metadata[, .(Run, rep_Run, category)])

# representative lat/layer info from the deduplicated table
rep_info <- dedup_runs[, .(rep_Run, manual_latitude, layer)]

# attach lat/layer to each original run
dt_pie <- rep_info[run_cats, on = "rep_Run"]

# keep usable rows
dt_pie <- dt_pie[
  !is.na(manual_latitude) &
    !is.na(layer) &
    !is.na(category)
]

lat_levels <- c(
  "south_polar",
  "south_subpolar",
  "south_subtropical",
  "equatorial",
  "north_subtropical",
  "north_subpolar",
  "north_polar"
)

dt_pie[, lat_bin := fcase(
  abs(manual_latitude) < 15,                         "equatorial",
  
  manual_latitude >= 15 & manual_latitude < 45,      "north_subtropical",
  manual_latitude <= -15 & manual_latitude > -45,    "south_subtropical",
  
  manual_latitude >= 45 & manual_latitude < 60,      "north_subpolar",
  manual_latitude <= -45 & manual_latitude > -60,    "south_subpolar",
  
  manual_latitude >= 60,                             "north_polar",
  manual_latitude <= -60,                            "south_polar",
  
  default = NA_character_
)]

dt_pie[, lat_bin := factor(lat_bin, levels = lat_levels)]

dt_pie[, cat4 := fcase(
  trimws(tolower(category)) == "cellular", "cellular",
  trimws(tolower(category)) == "virus",    "virus",
  trimws(tolower(category)) == "prok",     "prok",
  trimws(tolower(category)) == "large",    "large",
  default = "other"
)]

dt_pie[, cat4 := factor(
  cat4,
  levels = c("cellular", "prok", "virus", "large", "other")
)]

dt_pie[, layer := factor(layer, levels = c("epi", "apho"))]

bar_dt <- dt_pie[
  , .(n_runs = .N),
  by = .(lat_bin, layer, cat4)
]

pie_dt <- copy(bar_dt)
pie_dt[, prop := n_runs / sum(n_runs), by = .(lat_bin, layer)]

library(RColorBrewer)

cols <- c(
  cellular = brewer.pal(8, "Dark2")[1],
  prok     = brewer.pal(8, "Dark2")[4],
  virus    = brewer.pal(8, "Dark2")[6],
  large    = brewer.pal(8, "Dark2")[3],
  other    = "grey70"
)

ggplot(pie_dt, aes(x = "", y = prop, fill = cat4)) +
  geom_col(width = 1, color = "white", linewidth = 0.2) +
  coord_polar(theta = "y") +
  facet_grid(lat_bin ~ layer, drop = FALSE) +
  scale_fill_manual(values = cols, drop = FALSE) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_void(base_size = 11) +
  theme(
    aspect.ratio = 1,
    strip.background = element_rect(fill = "grey90", color = "grey50"),
    strip.text = element_text(color = "black"),
    panel.spacing = unit(0.6, "lines"),
    legend.position = "right"
  )