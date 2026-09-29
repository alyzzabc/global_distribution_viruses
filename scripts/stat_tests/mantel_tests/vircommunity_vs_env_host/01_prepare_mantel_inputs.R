# ============================================================
# 01_prepare_mantel_inputs.R
#
# Usage:
#   Rscript 01_prepare_mantel_inputs.R virus 50000 virus-v-cellular
#   Rscript 01_prepare_mantel_inputs.R cellular 50000 cellular-v-cellular
#
# Arguments:
#   arg1 = vir_category: "virus" or "cellular"
#   arg2 = top_n_to_use: number of top abundant features to retain
#   arg3 = run_label: output folder prefix/name
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(vegan)
  library(readr)
  library(fst)
})

# ============================================================
# 0. Command-line arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

vir_category <- ifelse(length(args) >= 1, args[1], "virus")
top_n_to_use <- ifelse(length(args) >= 2, as.integer(args[2]), 50000)
run_label <- ifelse(length(args) >= 3, args[3], paste0("", vir_category, "-v-cellular"))

if (!vir_category %in% c("virus", "cellular")) {
  stop("vir_category must be 'virus' or 'cellular'. Got: ", vir_category)
}

if (is.na(top_n_to_use)) {
  stop("top_n_to_use must be an integer.")
}

cat("============================================================\n")
cat("Preparing Mantel inputs\n")
cat("vir_category: ", vir_category, "\n", sep = "")
cat("top_n_to_use: ", top_n_to_use, "\n", sep = "")
cat("run_label: ", run_label, "\n", sep = "")
cat("working directory: ", getwd(), "\n", sep = "")
cat("============================================================\n")

# ============================================================
# 1. Paths
# ============================================================

base_dir <- "/path/to/base/dir"

derep_prok_file <- file.path(base_dir, "derep_avg_cov_OTU_with_tax_6may2026.fst")
derep_euk_file  <- file.path(base_dir, "derep_avg_cov_euks_with_tax_6may2026.fst")
derep_vir_file  <- file.path(base_dir, "derep_avg_normcov_vOTU_with_metadata_and_tax_6may26.fst")
copernicus_file <- file.path(base_dir, "copernicus_castant_1july2026.txt")

vcr_caudo_file <- file.path(base_dir, "vcr_res_cellular_caudo_vp_28may2026.rds")
vcr_mega_file  <- file.path(base_dir, "vcr_res_cellular_mega_ve_28may2026.rds")

out_dir <- file.path(base_dir, paste0(run_label, "_mantel"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Output directory: ", out_dir, "\n", sep = "")

# ============================================================
# 2. Matching columns
# ============================================================

match_cols <- c(
  "manual_longitude",
  "manual_depth",
  "manual_latitude",
  "manual_date"
)

# ============================================================
# 3. Read input files
# ============================================================

cat("Reading input files...\n")

derep_prok <- read_fst(derep_prok_file)
derep_euk  <- read_fst(derep_euk_file)
derep_vir  <- read_fst(derep_vir_file)

copernicus <- read_tsv(
  copernicus_file,
  show_col_types = FALSE
)

res_cellular_caudo_vp <- readRDS(vcr_caudo_file)
res_cellular_mega_ve  <- readRDS(vcr_mega_file)

cat("Finished reading input files.\n")

# ============================================================
# 4. Helper: standardize matching columns
# ============================================================

standardize_match_cols <- function(df) {
  df %>%
    mutate(
      manual_longitude = suppressWarnings(as.numeric(manual_longitude)),
      manual_latitude  = suppressWarnings(as.numeric(manual_latitude)),
      manual_depth     = suppressWarnings(as.numeric(manual_depth)),
      manual_date      = as.character(manual_date)
    )
}

# ============================================================
# 5. Helper: make hybrid environmental table
# ============================================================
# Manual values are preferred where available.
# Copernicus values are used as fallback for selected variables.
# Velocity variables are intentionally excluded.

make_hybrid_env_table <- function(metadata_source, copernicus) {
  
  copernicus_clean <- copernicus %>%
    distinct(Run, .keep_all = TRUE) %>%
    mutate(
      across(
        starts_with("copernicus_"),
        ~ suppressWarnings(as.numeric(na_if(as.character(.x), "nan")))
      ),
      Castant_PAR_0m = suppressWarnings(as.numeric(Castant_PAR_0m)),
      Castant_PAR = suppressWarnings(as.numeric(Castant_PAR)),
      Castant_Chl_0m = suppressWarnings(as.numeric(Castant_Chl_0m))
    )
  
  metadata_source %>%
    standardize_match_cols() %>%
    distinct(Run, .keep_all = TRUE) %>%
    left_join(
      copernicus_clean,
      by = "Run",
      relationship = "many-to-one"
    ) %>%
    transmute(
      Run,
      manual_longitude,
      manual_latitude,
      manual_depth,
      manual_date,
      manual_env,
      
      env_latitude = as.numeric(manual_latitude),
      env_abs_latitude = abs(as.numeric(manual_latitude)),
      env_longitude = as.numeric(manual_longitude),
      env_depth = as.numeric(manual_depth),
      
      env_temp = coalesce(
        suppressWarnings(as.numeric(manual_temp)),
        suppressWarnings(as.numeric(copernicus_Temperature_C_Daily))
      ),
      
      env_chl = coalesce(
        suppressWarnings(as.numeric(manual_chl)),
        suppressWarnings(as.numeric(copernicus_Chlorophyll_mg_m3_Daily))
      ),
      
      env_O2 = coalesce(
        suppressWarnings(as.numeric(manual_O2)),
        suppressWarnings(as.numeric(copernicus_Oxygen_mmol_m3_Daily))
      ),
      
      env_nitrate = coalesce(
        suppressWarnings(as.numeric(manual_nitrate)),
        suppressWarnings(as.numeric(copernicus_Nitrate_mmol_m3_Daily))
      ),
      
      env_salinity = suppressWarnings(as.numeric(copernicus_Salinity_psu_Daily)),
      env_mld = suppressWarnings(as.numeric(copernicus_MLD_m_Daily)),
      env_sea_surface_height = suppressWarnings(as.numeric(copernicus_Sea_Surface_Height_m_Daily)),
      env_phosphate = suppressWarnings(as.numeric(copernicus_Phosphate_mmol_m3_Daily)),
      env_silicate = suppressWarnings(as.numeric(copernicus_Silicate_mmol_m3_Daily)),
      env_npp = suppressWarnings(as.numeric(copernicus_NPP_mgC_m3_day_Daily)),
      
      env_PAR_0m = suppressWarnings(as.numeric(Castant_PAR_0m)),
      env_PAR = suppressWarnings(as.numeric(Castant_PAR)),
      env_chl_0m_castant = suppressWarnings(as.numeric(Castant_Chl_0m))
    )
}

# ============================================================
# 6. Helper: make abundance matrix
# ============================================================
# Rows = Run
# Columns = features
# Values = mean_norm_coverage

make_abundance_matrix <- function(
    df,
    feature_col,
    category_filter = NULL,
    class_filter = NULL
) {
  
  out <- df
  
  if (!is.null(category_filter)) {
    out <- out %>%
      filter(category == category_filter)
  }
  
  if (!is.null(class_filter)) {
    out <- out %>%
      filter(class %in% class_filter)
  }
  
  out <- out %>%
    select(
      Run,
      feature = all_of(feature_col),
      mean_norm_coverage
    ) %>%
    filter(!is.na(feature), feature != "")
  
  if (nrow(out) == 0) {
    warning("No rows left after filtering for feature_col = ", feature_col)
    return(data.frame(row.names = character()))
  }
  
  out %>%
    group_by(Run, feature) %>%
    summarise(
      mean_norm_coverage = mean(mean_norm_coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = feature,
      values_from = mean_norm_coverage,
      values_fill = 0
    ) %>%
    as.data.frame() %>%
    column_to_rownames("Run")
}

# ============================================================
# 7. Build viral community matrices
# ============================================================

cat("Building viral abundance matrices from category: ", vir_category, "\n", sep = "")

caudo_mat <- make_abundance_matrix(
  df = derep_vir,
  feature_col = "virus_genome",
  category_filter = vir_category,
  class_filter = "Caudoviricetes"
)

mega_mat <- make_abundance_matrix(
  df = derep_vir,
  feature_col = "virus_genome",
  category_filter = vir_category,
  class_filter = "Megaviricetes"
)

all_viruses_known_mat <- derep_vir %>%
  filter(
    category == vir_category,
    !is.na(class),
    class != ""
  ) %>%
  make_abundance_matrix(
    feature_col = "virus_genome"
  )

# ============================================================
# 8. Build host community matrices
# ============================================================

cat("Building host abundance matrices...\n")

prok_mat <- make_abundance_matrix(
  df = derep_prok,
  feature_col = "sequence",
  category_filter = vir_category
)

euk_mat <- make_abundance_matrix(
  df = derep_euk,
  feature_col = "taxon",
  category_filter = vir_category
)

cat("Raw matrix dimensions:\n")
print(dim(caudo_mat))
print(dim(mega_mat))
print(dim(all_viruses_known_mat))
print(dim(prok_mat))
print(dim(euk_mat))

# ============================================================
# 9. Build environmental tables
# ============================================================

cat("Building hybrid environmental tables...\n")

env_vir  <- make_hybrid_env_table(derep_vir, copernicus)
env_prok <- make_hybrid_env_table(derep_prok, copernicus)
env_euk  <- make_hybrid_env_table(derep_euk, copernicus)

# ============================================================
# 10. Filter environmental scope
# ============================================================

filter_env_scope <- function(
    env_df,
    env_type = "seawater",
    depth_zone = "all",
    epi_max_depth = 200
) {
  
  depth_zone <- match.arg(depth_zone, c("all", "epipelagic", "aphotic"))
  
  out <- env_df %>%
    mutate(
      manual_env_clean = tolower(trimws(manual_env)),
      depth_num = suppressWarnings(as.numeric(manual_depth))
    )
  
  if (!is.null(env_type)) {
    out <- out %>%
      filter(manual_env_clean == tolower(env_type))
  }
  
  if (depth_zone == "epipelagic") {
    out <- out %>%
      filter(!is.na(depth_num), depth_num <= epi_max_depth)
  }
  
  if (depth_zone == "aphotic") {
    out <- out %>%
      filter(!is.na(depth_num), depth_num > epi_max_depth)
  }
  
  out %>%
    select(-manual_env_clean, -depth_num)
}

env_type_to_use <- "seawater"
depth_zone_to_use <- "epipelagic"
epi_max_depth <- 200

env_vir_filtered <- filter_env_scope(
  env_vir,
  env_type = env_type_to_use,
  depth_zone = depth_zone_to_use,
  epi_max_depth = epi_max_depth
)

env_prok_filtered <- filter_env_scope(
  env_prok,
  env_type = env_type_to_use,
  depth_zone = depth_zone_to_use,
  epi_max_depth = epi_max_depth
)

env_euk_filtered <- filter_env_scope(
  env_euk,
  env_type = env_type_to_use,
  depth_zone = depth_zone_to_use,
  epi_max_depth = epi_max_depth
)

stopifnot(exists("env_vir_filtered"))
stopifnot(exists("env_prok_filtered"))
stopifnot(exists("env_euk_filtered"))

cat("Environmental table row counts:\n")
cat("env_vir: ", nrow(env_vir), " -> ", nrow(env_vir_filtered), "\n", sep = "")
cat("env_prok: ", nrow(env_prok), " -> ", nrow(env_prok_filtered), "\n", sep = "")
cat("env_euk: ", nrow(env_euk), " -> ", nrow(env_euk_filtered), "\n", sep = "")

# ============================================================
# 11. Build VCR nodes
# ============================================================

cat("Building VCR nodes...\n")

vcr_mega_euk <- res_cellular_mega_ve$plot_data %>%
  transmute(
    Run,
    VCR_Mega_Euk = as.numeric(virus_euk_ratio)
  ) %>%
  filter(!is.na(VCR_Mega_Euk), is.finite(VCR_Mega_Euk)) %>%
  distinct()

vcr_caudo_prok <- res_cellular_caudo_vp$plot_data %>%
  transmute(
    Run,
    VCR_Caudo_Prok = as.numeric(virus_prok_ratio)
  ) %>%
  filter(!is.na(VCR_Caudo_Prok), is.finite(VCR_Caudo_Prok)) %>%
  distinct()

cat("VCR overlap checks:\n")
cat("VCR_Caudo_Prok with prok env: ", sum(vcr_caudo_prok$Run %in% env_prok_filtered$Run), "\n", sep = "")
cat("VCR_Caudo_Prok with viral env: ", sum(vcr_caudo_prok$Run %in% env_vir_filtered$Run), "\n", sep = "")
cat("VCR_Mega_Euk with euk env: ", sum(vcr_mega_euk$Run %in% env_euk_filtered$Run), "\n", sep = "")
cat("VCR_Mega_Euk with viral env: ", sum(vcr_mega_euk$Run %in% env_vir_filtered$Run), "\n", sep = "")

# ============================================================
# 12. Helper: filter community matrices
# ============================================================

filter_community_mat <- function(
    mat,
    min_prevalence = 3,
    top_n = 50000
) {
  
  mat <- mat[, colSums(mat > 0, na.rm = TRUE) > 0, drop = FALSE]
  mat <- mat[rowSums(mat > 0, na.rm = TRUE) > 0, , drop = FALSE]
  
  prev <- colSums(mat > 0, na.rm = TRUE)
  mat <- mat[, prev >= min_prevalence, drop = FALSE]
  
  total_abund <- colSums(mat, na.rm = TRUE)
  
  if (ncol(mat) > top_n) {
    keep <- names(sort(total_abund, decreasing = TRUE))[1:top_n]
    mat <- mat[, keep, drop = FALSE]
  }
  
  mat
}

# ============================================================
# 13. Filter matrices
# ============================================================

cat("Filtering community matrices with top_n = ", top_n_to_use, "\n", sep = "")

caudo_mat_filt <- filter_community_mat(
  caudo_mat,
  min_prevalence = 3,
  top_n = top_n_to_use
)

mega_mat_filt <- filter_community_mat(
  mega_mat,
  min_prevalence = 3,
  top_n = top_n_to_use
)

all_viruses_known_mat_filt <- filter_community_mat(
  all_viruses_known_mat,
  min_prevalence = 3,
  top_n = top_n_to_use
)

prok_mat_filt <- filter_community_mat(
  prok_mat,
  min_prevalence = 3,
  top_n = top_n_to_use
)

euk_mat_filt <- filter_community_mat(
  euk_mat,
  min_prevalence = 3,
  top_n = top_n_to_use
)

cat("Filtered matrix dimensions:\n")
cat("Caudoviricetes raw/filtered:\n")
print(dim(caudo_mat))
print(dim(caudo_mat_filt))

cat("Megaviricetes raw/filtered:\n")
print(dim(mega_mat))
print(dim(mega_mat_filt))

cat("All viruses known raw/filtered:\n")
print(dim(all_viruses_known_mat))
print(dim(all_viruses_known_mat_filt))

cat("Prokaryotes raw/filtered:\n")
print(dim(prok_mat))
print(dim(prok_mat_filt))

cat("Eukaryotes raw/filtered:\n")
print(dim(euk_mat))
print(dim(euk_mat_filt))

# ============================================================
# 14. Define nodes and variables
# ============================================================

community_nodes <- list(
  Caudoviricetes = caudo_mat_filt,
  Megaviricetes = mega_mat_filt,
  All_viruses_known_tax = all_viruses_known_mat_filt,
  Prokaryotes = prok_mat_filt,
  Eukaryotes = euk_mat_filt
)

env_by_node <- list(
  Caudoviricetes = env_vir_filtered,
  Megaviricetes = env_vir_filtered,
  All_viruses_known_tax = env_vir_filtered,
  Prokaryotes = env_prok_filtered,
  Eukaryotes = env_euk_filtered
)

env_vars <- c(
  "env_latitude",
  "env_abs_latitude",
  "env_longitude",
  "env_depth",
  "env_temp",
  "env_chl",
  "env_O2",
  "env_nitrate",
  "env_salinity",
  "env_mld",
  "env_sea_surface_height",
  "env_phosphate",
  "env_silicate",
  "env_npp",
  "env_chl_0m_castant",
  "env_PAR_0m",
  "env_PAR"
)

vcr_nodes <- list(
  VCR_Caudo_Prok = vcr_caudo_prok %>%
    select(Run, value = VCR_Caudo_Prok),
  
  VCR_Mega_Euk = vcr_mega_euk %>%
    select(Run, value = VCR_Mega_Euk)
)

# ============================================================
# 15. Sanity checks
# ============================================================

cat("Sanity checks:\n")

cat("Community node dimensions:\n")
print(sapply(community_nodes, dim))

cat("Environmental table rows by node:\n")
print(sapply(env_by_node, nrow))

cat("VCR node rows:\n")
print(sapply(vcr_nodes, nrow))

cat("Community nodes missing env tables:\n")
print(setdiff(names(community_nodes), names(env_by_node)))

cat("Env vars missing in viral env table:\n")
print(setdiff(env_vars, colnames(env_vir_filtered)))

cat("Env vars missing in prok env table:\n")
print(setdiff(env_vars, colnames(env_prok_filtered)))

cat("Env vars missing in euk env table:\n")
print(setdiff(env_vars, colnames(env_euk_filtered)))

if (length(setdiff(names(community_nodes), names(env_by_node))) > 0) {
  stop("Some community nodes are missing env_by_node entries.")
}

if (length(setdiff(env_vars, colnames(env_vir_filtered))) > 0) {
  stop("Some env_vars are missing in env_vir_filtered.")
}

if (length(setdiff(env_vars, colnames(env_prok_filtered))) > 0) {
  stop("Some env_vars are missing in env_prok_filtered.")
}

if (length(setdiff(env_vars, colnames(env_euk_filtered))) > 0) {
  stop("Some env_vars are missing in env_euk_filtered.")
}

# ============================================================
# 16. Save RDS
# ============================================================

output_file <- file.path(out_dir, "mantel_inputs_filtered.rds")

saveRDS(
  list(
    community_nodes = community_nodes,
    env_by_node = env_by_node,
    env_vars = env_vars,
    vcr_nodes = vcr_nodes,
    match_cols = match_cols,
    vir_category = vir_category,
    top_n_to_use = top_n_to_use,
    run_label = run_label,
    env_type_to_use = env_type_to_use,
    depth_zone_to_use = depth_zone_to_use,
    epi_max_depth = epi_max_depth
  ),
  output_file
)

cat("============================================================\n")
cat("Saved Mantel input RDS\n")
cat("File: ", output_file, "\n", sep = "")
cat("vir_category: ", vir_category, "\n", sep = "")
cat("top_n_to_use: ", top_n_to_use, "\n", sep = "")
cat("run_label: ", run_label, "\n", sep = "")
cat("============================================================\n")