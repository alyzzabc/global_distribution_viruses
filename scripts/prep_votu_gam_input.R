library(data.table)
library(dplyr)
library(readr)

# ============================================================
# 1) Read viral coverage table and normalize coverage per Gb
# ============================================================

global_viruses <- read_csv(
  "path/to/coverage/file"
)

setDT(global_viruses)

setnames(
  global_viruses,
  c("Sample_name", "Contig", "Value"),
  c("Run", "virus_genome", "coverage")
)

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

# Keep only viral coverage rows from runs with enough bases
global_viruses <- global_viruses[
  Run %chin% bases_filtered$sample_name
]

# Normalize coverage per Gb
global_viruses[, norm_coverage := coverage / (total_bases_vec[Run] / 1e9)]

# Rename Run in viral table to avoid confusion after joins
gv <- copy(global_viruses)
setnames(gv, "Run", "virus_Run")

# gv now contains:
# virus_Run, virus_genome, coverage, norm_coverage


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

# Make size strings UTF-8 safe
metadata[, manual_size_text_strings := iconv(
  manual_size_text_strings,
  from = "",
  to = "UTF-8",
  sub = ""
)]

# Normalize common missing-like placeholders to real NA
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

# Helper function:
# For known values, use the actual value.
# For missing / placeholder values, create a unique value per row.
# This prevents samples from being deduplicated just because both have NA.
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
# For viral coverage, we do NOT collapse time-series samples.
# Therefore, the actual date is part of the deduplication key.
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
# Any missing value in those fields was made unique above,
# so missing values cannot accidentally cause deduplication.

dedup_keys <- c(
  "key_latitude",
  "key_longitude",
  "key_date",
  "key_depth",
  "key_env",
  "key_category",
  "size_group"
)

# Representative Run for each technical-replicate group
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

# Check largest technical-replicate groups
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

# Viral runs missing metadata
gv[!virus_Run %chin% metadata$Run, unique(virus_Run)]


# ============================================================
# 7) Attach rep_Run to viral coverage rows
# ============================================================

# Map every original Run to its representative Run
run_map <- metadata[, .(Run, rep_Run)]

# Attach rep_Run to viral table
gv_meta <- run_map[
  gv,
  on = .(Run = virus_Run),
  nomatch = 0L
]

# gv_meta now contains:
# Run, rep_Run, virus_genome, coverage, norm_coverage


# ============================================================
# 8) Average coverage only within technical replicate groups
# ============================================================

# Important sparse-table issue:
# The viral coverage table is sparse.
# It does not contain every possible Run x virus_genome combination.
#
# We avoid creating the full global Run x virus grid.
# Instead, for each technical-replicate group, we only expand:
#   all Runs in that rep_Run group
#   x
#   viruses observed in at least one Run of that rep_Run group
#
# If a virus is present in one technical replicate but absent/missing
# in another replicate from the same group, the missing value is treated as zero.
#
# This gives a technical-replicate mean without exploding to the full
# global Run x virus matrix.

# All original runs per representative group
rep_runs <- unique(run_map[, .(rep_Run, Run)])

# Viruses observed at least once within each representative group
rep_viruses <- unique(gv_meta[, .(rep_Run, virus_genome)])

# Build group-local grid:
# for each rep_Run, all original runs x viruses observed in that group
rep_grid <- rep_runs[
  rep_viruses,
  on = "rep_Run",
  allow.cartesian = TRUE
]

# Attach observed coverage to the group-local grid
rep_cov <- gv_meta[
  rep_grid,
  on = .(rep_Run, Run, virus_genome)
]

# Missing within a technical-replicate group is treated as zero
rep_cov[is.na(coverage), coverage := 0]
rep_cov[is.na(norm_coverage), norm_coverage := 0]

# Average coverage across technical replicates only
avg_cov <- rep_cov[
  , .(
    mean_norm_coverage = mean(norm_coverage),
    mean_coverage = mean(coverage),
    n_runs_collapsed = uniqueN(Run),
    n_runs_detected = uniqueN(Run[norm_coverage > 0])
  ),
  by = .(rep_Run, virus_genome)
]


# ============================================================
# 9) Attach representative metadata to averaged viral coverage
# ============================================================

avg_cov_meta <- dedup_metadata[
  avg_cov,
  on = "rep_Run",
  nomatch = 0L
]

# Convert coordinates/depth to numeric
avg_cov_meta[, manual_longitude := suppressWarnings(as.numeric(manual_longitude))]
avg_cov_meta[, manual_latitude := suppressWarnings(as.numeric(manual_latitude))]
avg_cov_meta[, depth_m := suppressWarnings(as.numeric(manual_depth))]

# Assign layer.
# Unknown/non-numeric depth stays NA.
avg_cov_meta[, layer := fifelse(
  is.na(depth_m),
  NA_character_,
  fifelse(depth_m < 200, "epi", "apho")
)]


# ============================================================
# 10) Quick inspection
# ============================================================

avg_cov_meta[
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
    virus_genome,
    mean_norm_coverage,
    mean_coverage,
    n_runs_collapsed,
    n_runs_detected,
    layer
  )
]

# Check largest technical-replicate groups represented in viral table
avg_cov_meta[
  , .N,
  by = .(rep_Run, n_runs_collapsed)
][order(-n_runs_collapsed)][1:20]

############################
# Sanity checks
# How many original Runs collapse into each rep_Run?
metadata[, .N, by = rep_Run][order(-N)][1:10]

# Any virus rows missing metadata?
global_viruses[!Run %chin% metadata$Run, unique(Run)]
############################

# Convert columns to numeric
avg_cov_meta[, manual_longitude := as.numeric(manual_longitude)]
avg_cov_meta[, manual_latitude := as.numeric(manual_latitude)]
avg_cov_meta[, manual_depth := as.numeric(manual_depth)]

# Classify layer
avg_cov_meta[, layer := fifelse(manual_depth < 200, "epi",
                                fifelse(manual_depth >= 200, "apho", NA_character_))]



###################################################################################
################################ Prepare GAM inputs ###############################
###################################################################################

avg_cov_meta_sw <- avg_cov_meta %>% filter (manual_env=="seawater")

forGAM_vOTU_epi_allfracs <- avg_cov_meta_sw %>%
  filter(layer == "epi") %>%
  select(Run, manual_latitude, virus_genome, mean_norm_coverage) %>%
  rename(
    sample_name = Run,
    latitude = manual_latitude,
    taxon_id = virus_genome,
    norm_coverage = mean_norm_coverage
  )

forGAM_vOTU_epi_cellular <- avg_cov_meta_sw %>% 
  filter(layer=="epi" & category %in% c("prok", "cellular"))  %>%
           select(Run, manual_latitude, virus_genome, mean_norm_coverage) %>%
           rename(
             sample_name = Run,
             latitude = manual_latitude,
             taxon_id = virus_genome,
             norm_coverage = mean_norm_coverage
           )

forGAM_vOTU_epi_viruses <- avg_cov_meta_sw %>% 
  filter(layer=="epi" & category=="virus")  %>%
  select(Run, manual_latitude, virus_genome, mean_norm_coverage) %>%
  rename(
    sample_name = Run,
    latitude = manual_latitude,
    taxon_id = virus_genome,
    norm_coverage = mean_norm_coverage
  )

### Link avg_cov_meta with taxonomy ##
taxonomy <- read_tsv("path/to/taxonomy/file")
taxonomy_split <- taxonomy %>%
  separate(taxonomy, into = c("domain", "realm", "kingdom", "phylum", "class", "order", "family"), sep = ";", fill = "right", extra = "merge")

avg_cov_meta_tax <- avg_cov_meta %>%
  left_join(
    taxonomy_split %>% select(seq_name, class),
    by = c("virus_genome" = "seq_name")
  )
sum(!unique(avg_cov_meta$virus_genome) %in% taxonomy_split$seq_name)

avg_cov_meta_sw_tax <- avg_cov_meta_sw %>%
  left_join(
    taxonomy_split %>% select(seq_name, class),
    by = c("virus_genome" = "seq_name")
  )
sum(!unique(avg_cov_meta_sw$virus_genome) %in% taxonomy_split$seq_name)