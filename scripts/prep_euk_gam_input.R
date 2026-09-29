library(data.table)
library(readr)

# ============================================================
# 1) Read eukaryote coverage table and normalize coverage per Gb
# ============================================================

global_euks <- read_tsv(
  "/path/to/coverage/file",
  col_names = TRUE
)

setDT(global_euks)

# Expected useful columns include:
# sample, taxon, coverage, cpm, taxon_num_reads, taxon_num_alignments, etc.

# Rename sample column to Run
setnames(global_euks, "sample", "Run")

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

# Keep only eukaryote coverage rows from runs with enough bases
global_euks <- global_euks[
  Run %chin% bases_filtered$sample_name
]

# Normalize coverage per Gb
global_euks[, norm_coverage := coverage / (total_bases_vec[Run] / 1e9)]

# Rename Run in euk table to avoid confusion after joins
ge <- copy(global_euks)
setnames(ge, "Run", "euk_Run")

# ge now contains:
# euk_Run, taxon, coverage, cpm, taxon_num_reads, ..., norm_coverage


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
# Every missing deduplication field gets a row-specific placeholder.

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
# For eukaryote coverage, as with viral/prokaryote coverage,
# we collapse only technical replicates.
# We do NOT collapse time-series samples.
#
# Therefore, actual manual_date is part of the deduplication key.
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

# Euk runs missing metadata
ge[!euk_Run %chin% metadata$Run, unique(euk_Run)]


# ============================================================
# 7) Attach rep_Run to eukaryote coverage rows
# ============================================================

# Map every original Run to its representative Run
run_map <- metadata[, .(Run, rep_Run)]

# Attach rep_Run to eukaryote table
ge_meta <- run_map[
  ge,
  on = .(Run = euk_Run),
  nomatch = 0L
]

# ge_meta now contains:
# Run, rep_Run, taxon, coverage, cpm, taxon_num_reads, ..., norm_coverage


# ============================================================
# 8) Average eukaryote coverage only within technical replicate groups
# ============================================================

# Important sparse-table issue:
# The eukaryote coverage table may be sparse.
# It may not contain every possible Run x taxon combination.
#
# We avoid making the full global Run x taxon matrix.
# Instead, for each technical-replicate group, we only expand:
#   all Runs in that rep_Run group
#   x
#   taxa observed in at least one Run of that rep_Run group
#
# If a taxon is present in one technical replicate but absent/missing
# in another replicate from the same group, the missing value is treated as zero.
#
# This gives a technical-replicate mean without collapsing time-series samples.

# All original runs per representative group
rep_runs <- unique(run_map[, .(rep_Run, Run)])

# Eukaryote taxa observed at least once within each representative group
rep_euks <- unique(ge_meta[, .(rep_Run, taxon)])

# Build group-local grid:
# for each rep_Run, all original runs x taxa observed in that group
rep_grid <- rep_runs[
  rep_euks,
  on = "rep_Run",
  allow.cartesian = TRUE
]

# Attach observed coverage to the group-local grid
rep_cov_euk <- ge_meta[
  rep_grid,
  on = .(rep_Run, Run, taxon)
]

# Missing coverage within a technical-replicate group is treated as zero
rep_cov_euk[is.na(coverage), coverage := 0]
rep_cov_euk[is.na(norm_coverage), norm_coverage := 0]

# Average coverage across technical replicates only
avg_euk <- rep_cov_euk[
  , .(
    mean_norm_coverage = mean(norm_coverage),
    mean_coverage = mean(coverage),
    mean_cpm = mean(fifelse(is.na(cpm), 0, cpm)),
    mean_taxon_num_reads = mean(fifelse(is.na(taxon_num_reads), 0, taxon_num_reads)),
    mean_taxon_num_alignments = mean(fifelse(is.na(taxon_num_alignments), 0, taxon_num_alignments)),
    n_runs_collapsed = uniqueN(Run),
    n_runs_detected = uniqueN(Run[norm_coverage > 0])
  ),
  by = .(rep_Run, taxon)
]


# ============================================================
# 9) Attach representative metadata to averaged eukaryote coverage
# ============================================================

avg_euk_meta <- dedup_metadata[
  avg_euk,
  on = "rep_Run",
  nomatch = 0L
]

# Convert coordinates/depth to numeric
avg_euk_meta[, manual_longitude := suppressWarnings(as.numeric(manual_longitude))]
avg_euk_meta[, manual_latitude := suppressWarnings(as.numeric(manual_latitude))]
avg_euk_meta[, depth_m := suppressWarnings(as.numeric(manual_depth))]

# Assign layer.
# Unknown/non-numeric depth stays NA.
avg_euk_meta[, layer := fifelse(
  is.na(depth_m),
  NA_character_,
  fifelse(depth_m < 200, "epi", "apho")
)]


# ============================================================
# 10) Quick inspection
# ============================================================

avg_euk_meta[
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
    taxon,
    mean_norm_coverage,
    mean_coverage,
    mean_cpm,
    n_runs_collapsed,
    n_runs_detected,
    layer
  )
]

# Check largest technical-replicate groups represented in eukaryote table
avg_euk_meta[
  , .(n_taxa = uniqueN(taxon)),
  by = .(rep_Run, n_runs_collapsed)
][order(-n_runs_collapsed)][1:20]


