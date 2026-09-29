library(data.table)
library(readr)

# ============================================================
# 1) Read prokaryote coverage table and normalize coverage per Gb
# ============================================================

global_proks <- read_tsv(
  "path/to/coverage/file",
  col_names = TRUE
)

setDT(global_proks)

# Expected columns:
# Run, sequence, coverage, taxonomy

# Read sequencing-depth information
bases <- fread(
  "/path/to/seq/lengths/file"
)

bases[, total_length_bases := as.numeric(total_length_bases)]

# Keep only runs with at least 1 Mbp sequenced
bases_filtered <- bases[total_length_bases >= 1e6]

total_bases_vec <- setNames(
  bases_filtered$total_length_bases,
  bases_filtered$sample_name
)

# Keep only prokaryote coverage rows from runs with enough bases
global_proks <- global_proks[
  Run %chin% bases_filtered$sample_name
]

# Normalize coverage per Gb
global_proks[, norm_coverage := coverage / (total_bases_vec[Run] / 1e9)]

# Rename Run in prok table to avoid confusion after joins
gp <- copy(global_proks)
setnames(gp, "Run", "prok_Run")

# gp now contains:
# prok_Run, sequence, coverage, taxonomy, norm_coverage


# ============================================================
# 2) Read metadata
# ============================================================

metadata <- fread(
  "path/to/metadata",
  select = c(
    "Run",
    "manual_date",
    "manual_latitude",
    "manual_longitude",
    "manual_depth",
    "manual_env",
    "category",
    "manual_temp",
    "manual_chl",
    "manual_O2",
    "manual_nitrate",
    "study",
    "manual_size_text_strings"
  )
)

setDT(metadata)

# Add stable row ID.
# This is used to make missing deduplication fields unique per row.
metadata[, row_id_for_dedup := .I]


# ============================================================
# 3) Clean size-fraction strings
# ============================================================

metadata[, manual_size_text_strings := iconv(
  manual_size_text_strings,
  from = "",
  to = "UTF-8",
  sub = ""
)]

metadata[, manual_size_text_strings := fifelse(
  is.na(manual_size_text_strings) |
    trimws(manual_size_text_strings) %chin% c("", "missing", "#N/A", "NA", "N/A"),
  NA_character_,
  trimws(manual_size_text_strings)
)]

# Missing size values are made unique per row.
# This prevents unrelated samples with unknown size fraction from collapsing together.
metadata[, size_group := fcoalesce(
  manual_size_text_strings,
  paste0("__no_size__", row_id_for_dedup)
)]


# ============================================================
# 4) Create NA-safe deduplication keys
# ============================================================

# For technical-replicate deduplication, missing values should NOT group together.
# So every missing deduplication field gets a row-specific placeholder.

make_na_safe_key <- function(x, prefix, row_id) {
  x_chr <- trimws(as.character(x))
  
  fifelse(
    is.na(x) |
      is.na(x_chr) |
      x_chr %chin% c("", "missing", "#N/A", "NA", "N/A"),
    paste0("__no_", prefix, "__", row_id),
    x_chr
  )
}

metadata[, key_latitude := make_na_safe_key(
  manual_latitude,
  "latitude",
  row_id_for_dedup
)]

metadata[, key_longitude := make_na_safe_key(
  manual_longitude,
  "longitude",
  row_id_for_dedup
)]

metadata[, key_date := make_na_safe_key(
  manual_date,
  "date",
  row_id_for_dedup
)]

metadata[, key_depth := make_na_safe_key(
  manual_depth,
  "depth",
  row_id_for_dedup
)]

metadata[, key_env := make_na_safe_key(
  manual_env,
  "env",
  row_id_for_dedup
)]

metadata[, key_category := make_na_safe_key(
  category,
  "category",
  row_id_for_dedup
)]


# ============================================================
# 5) Deduplicate technical replicates only
# ============================================================

# IMPORTANT:
# For prokaryote coverage, we do NOT collapse time-series samples.
# Therefore, actual date is part of the deduplication key.
#
# Technical replicates collapse only if they share:
# - latitude
# - longitude
# - date
# - depth
# - environment
# - category
# - size fraction
#
# Missing values in these fields were made unique above,
# so missing metadata cannot accidentally cause deduplication.

dedup_keys <- c(
  "key_latitude",
  "key_longitude",
  "key_date",
  "key_depth",
  "key_env",
  "key_category",
  "size_group"
)

metadata[, rep_Run := Run[order(Run)][1], by = dedup_keys]

# Number of original runs collapsed into each representative Run
rep_counts <- metadata[, .(n_reps = .N), by = rep_Run]

# One metadata row per representative Run
dedup_metadata <- unique(metadata[Run == rep_Run])

# Attach replicate counts
dedup_metadata[rep_counts, n_reps := i.n_reps, on = "rep_Run"]


# ============================================================
# 6) Sanity checks for deduplication
# ============================================================

# Largest technical-replicate groups
metadata[, .N, by = rep_Run][order(-N)][1:20]

# These should be empty because env and category are part of dedup_keys
env_check <- metadata[
  , .(
    n_env = uniqueN(manual_env),
    envs = paste(sort(unique(manual_env)), collapse = "; ")
  ),
  by = rep_Run
][n_env > 1]

cat_check <- metadata[
  , .(
    n_category = uniqueN(category),
    categories = paste(sort(unique(category)), collapse = "; ")
  ),
  by = rep_Run
][n_category > 1]

env_check
cat_check

# Prok runs missing metadata
gp[!prok_Run %chin% metadata$Run, unique(prok_Run)]


# ============================================================
# 7) Attach rep_Run to prokaryote coverage rows
# ============================================================

# Map every original Run to its representative Run
run_map <- metadata[, .(Run, rep_Run)]

# Attach rep_Run to prokaryote table
gp_meta <- run_map[
  gp,
  on = .(Run = prok_Run),
  nomatch = 0L
]

# gp_meta now contains:
# Run, rep_Run, sequence, coverage, taxonomy, norm_coverage


# ============================================================
# 8) Average prokaryote coverage only within technical replicate groups
# ============================================================

# Important sparse-table issue:
# The prokaryote coverage table may be sparse.
# It may not contain every possible Run x sequence combination.
#
# We avoid making the full global Run x sequence matrix.
# Instead, for each technical-replicate group, we only expand:
#   all Runs in that rep_Run group
#   x
#   sequences observed in at least one Run of that rep_Run group
#
# If a sequence is present in one technical replicate but absent/missing
# in another replicate from the same group, the missing value is treated as zero.
#
# This gives a technical-replicate mean without collapsing time-series samples.

# All original runs per representative group
rep_runs <- unique(run_map[, .(rep_Run, Run)])

# Prokaryote sequences observed at least once within each representative group
rep_proks <- unique(gp_meta[, .(rep_Run, sequence)])

# Build group-local grid:
# for each rep_Run, all original runs x sequences observed in that group
rep_grid <- rep_runs[
  rep_proks,
  on = "rep_Run",
  allow.cartesian = TRUE
]

# Attach observed coverage to the group-local grid
rep_cov_prok <- gp_meta[
  rep_grid,
  on = .(rep_Run, Run, sequence)
]

# Missing coverage within a technical-replicate group is treated as zero
rep_cov_prok[is.na(coverage), coverage := 0]
rep_cov_prok[is.na(norm_coverage), norm_coverage := 0]

# Preserve taxonomy for each sequence.
# This assumes each sequence has one taxonomy.
# If a sequence has multiple taxonomy strings, they are collapsed with "; ".
seq_taxonomy <- gp_meta[
  !is.na(taxonomy),
  .(taxonomy = paste(sort(unique(taxonomy)), collapse = "; ")),
  by = sequence
]

# Average coverage across technical replicates only
avg_prok <- rep_cov_prok[
  , .(
    mean_norm_coverage = mean(norm_coverage),
    mean_coverage = mean(coverage),
    n_runs_collapsed = uniqueN(Run),
    n_runs_detected = uniqueN(Run[norm_coverage > 0])
  ),
  by = .(rep_Run, sequence)
]

# Attach taxonomy
avg_prok <- seq_taxonomy[
  avg_prok,
  on = "sequence"
]


# ============================================================
# 9) Attach representative metadata to averaged prokaryote coverage
# ============================================================

avg_prok_meta <- dedup_metadata[
  avg_prok,
  on = "rep_Run",
  nomatch = 0L
]

# Convert coordinates/depth to numeric
avg_prok_meta[, manual_longitude := suppressWarnings(as.numeric(manual_longitude))]
avg_prok_meta[, manual_latitude := suppressWarnings(as.numeric(manual_latitude))]
avg_prok_meta[, depth_m := suppressWarnings(as.numeric(manual_depth))]

# Assign layer.
# Unknown/non-numeric depth stays NA.
avg_prok_meta[, layer := fifelse(
  is.na(depth_m),
  NA_character_,
  fifelse(depth_m < 200, "epi", "apho")
)]


# ============================================================
# 10) Quick inspection
# ============================================================

avg_prok_meta[
  1:5,
  .(
    rep_Run,
    manual_env,
    category,
    manual_latitude,
    manual_longitude,
    manual_date,
    manual_depth,
    manual_size_text_strings,
    sequence,
    taxonomy,
    mean_norm_coverage,
    mean_coverage,
    n_runs_collapsed,
    n_runs_detected,
    layer
  )
]

# Check largest technical-replicate groups represented in prokaryote table
avg_prok_meta[
  , .(n_sequences = uniqueN(sequence)),
  by = .(rep_Run, n_runs_collapsed)
][order(-n_runs_collapsed)][1:20]


## Actually prepare GAM input
library(dplyr)
library(ggplot2)

# 1) Restrict to epipelagic seawater only
avg_prok_meta_epi_sw <- avg_prok_meta %>%
  filter(layer == "epi", manual_env == "seawater")

# Quick sanity check
avg_prok_meta_epi_sw %>% count(layer, manual_env, sort = TRUE)

# 2) OTU-level summary in epi + seawater
avg_prok_meta_epi_sw_renamed <- avg_prok_meta_epi_sw %>%
  filter(layer == "epi") %>%
  #select(Run, manual_latitude, sequence, mean_norm_coverage) %>%
  rename(
    sample_name = Run,
    latitude = manual_latitude,
    taxon_id = otu_id,
    norm_coverage = mean_norm_coverage
  )

avg_prok_meta_epi_sw_summary <- avg_prok_meta_epi_sw_renamed %>%
  group_by(taxon_id) %>%
  summarise(
    mean_norm_coverage_epi = if (all(is.na(norm_coverage))) NA_real_ else mean(norm_coverage, na.rm = TRUE),
    median_norm_coverage_epi = if (all(is.na(norm_coverage))) NA_real_ else median(norm_coverage, na.rm = TRUE),
    max_norm_coverage_epi = if (all(is.na(norm_coverage))) NA_real_ else max(norm_coverage, na.rm = TRUE),
    n_runs_epi = n(),
    n_nonmissing_cov_epi = sum(!is.na(norm_coverage)),
    n_detected_epi = sum(norm_coverage > 0, na.rm = TRUE),
    n_studies_epi = n_distinct(study[!is.na(study)]),
    lat_range_epi = {
      lat_non_na <- latitude[!is.na(latitude)]
      if (length(lat_non_na) == 0) NA_real_ else diff(range(lat_non_na))
    },
    .groups = "drop"
  )

ggplot(avg_prok_meta_epi_sw_summary, aes(x = n_detected_epi)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  labs(x = "Number of detected epipelagic seawater runs per OTU",
       y = "Number of OTUs")

ggplot(avg_prok_meta_epi_sw_summary, aes(x = log10(mean_norm_coverage_epi + 1e-6))) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(x = "log10(mean_norm_coverage_epi + 1e-6)",
       y = "Number of OTUs")

avg_prok_meta_epi_sw_summary %>%
  mutate(pass = n_detected_epi >= 5 & mean_norm_coverage_epi > 0) %>%
  ggplot(aes(x = n_detected_epi,
             y = log10(mean_norm_coverage_epi + 1e-6),
             color = pass)) +
  geom_point(alpha = 0.2, size = 0.7) +
  theme_bw() +
  labs(x = "n_detected_epi",
       y = "log10(mean_norm_coverage_epi + 1e-6)")

# Filter OTUs for those with at least 5 detections
otus_for_bigeog_v1_sw <- avg_prok_meta_epi_sw_summary %>%
  mutate(pass_v1 = n_detected_epi >= 5 & mean_norm_coverage_epi > 0) %>%
  filter(pass_v1)

n_retained_sw <- nrow(otus_for_bigeog_v1_sw)
n_retained_sw

# Now make the table for GAM
otu_ids_bigeog_v1_sw <- otus_for_bigeog_v1_sw$taxon_id

derep_OTU_bigeog_v1_epi_sw <- avg_prok_meta %>%
  filter(layer == "epi",
         manual_env == "seawater",
         otu_id %in% otu_ids_bigeog_v1_sw)

# sanity check
c(
  rows = nrow(derep_OTU_bigeog_v1_epi_sw),
  unique_otus = n_distinct(derep_OTU_bigeog_v1_epi_sw$otu_id)
)

# GAM prok + cellular fractions
derep_OTU_bigeog_v1_epi_sw_prok_cell <- derep_OTU_bigeog_v1_epi_sw %>%
  filter(category %in% c("prok", "cellular"))

forGAM_OTU_cellular_prokfracs <- derep_OTU_bigeog_v1_epi_sw_prok_cell %>%
  transmute(
    sample_name = rep_Run,
    latitude = manual_latitude,
    taxon_id = otu_id,
    norm_coverage = mean_norm_coverage
  )


# GAM all fractions
forGAM_OTU_allfracs <- derep_OTU_bigeog_v1_epi_sw %>%
  transmute(
    sample_name = rep_Run,
    latitude = manual_latitude,
    taxon_id = otu_id,
    norm_coverage = mean_norm_coverage
  )

