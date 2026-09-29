library(data.table)
library(vegan)

# Start from your virus coverage file with taxonomy
# merged_table has: Run, virus_genome, norm_coverage, layer, etc.

# Step 1: Define presence/absence threshold
presence_threshold <- 0.01

# Step 2: Make a presence/absence table per layer
make_presence_matrix <- function(dt_layer) {
  dt_layer[, presence := as.integer(mean_norm_coverage > presence_threshold)]
  mat <- dcast(dt_layer, rep_Run ~ virus_genome, value.var = "presence", fill = 0)
  # Set sample names as rownames and remove the 'Run' column
  mat_out <- as.data.frame(mat)
  rownames(mat_out) <- mat_out$Run
  mat_out$rep_Run <- NULL
  return(mat_out)
}

# Step 3: Split data by layer
merged_table <- avg_cov_meta_tax_sw
setDT(merged_table)
epi_dt  <- merged_table[layer == "epi"]
deep_dt <- merged_table[layer == "deep"]

# Step 4: Make community matrices
epi_mat  <- make_presence_matrix(epi_dt)
deep_mat <- make_presence_matrix(deep_dt)

# Optional: combined all samples
combined_dt  <- merged_table
combined_mat <- make_presence_matrix(combined_dt)

# Step 5: Compute accumulation curves with random permutations
set.seed(42)
acc_epi  <- specaccum(epi_mat, method = "random", permutations = 100)
acc_deep <- specaccum(deep_mat, method = "random", permutations = 100)
acc_comb <- specaccum(combined_mat, method = "random", permutations = 100)

# Step 6: Plot the curves
plot(acc_epi,  ci.type = "poly", ci.col = rgb(0,0,1,0.2), col = "blue", lwd = 2,
     xlab = "Number of samples", ylab = "Number of vOTUs detected",
     main = "Viral richness accumulation curves by depth layer")
plot(acc_deep, ci.type = "poly", ci.col = rgb(1,0,0,0.2), col = "red", lwd = 2, add = TRUE)
plot(acc_comb, ci.type = "poly", ci.col = rgb(0,0,0,0.1), col = "black", lwd = 2, add = TRUE)

legend("bottomright", legend = c("Epipelagic", "Deep", "All samples"),
       col = c("blue", "red", "black"), lwd = 2, bty = "n")
grid()

#### ggplot version
library(ggplot2)
library(dplyr)
library(tidyr)

# Convert specaccum objects into data frames
acc_to_df <- function(acc_obj, layer_name) {
  data.frame(
    samples = acc_obj$sites,
    richness = acc_obj$richness,
    sd = acc_obj$sd,
    layer = layer_name
  )
}

# Combine your three curves
df_epi  <- acc_to_df(acc_epi,  "Epipelagic")
df_deep <- acc_to_df(acc_deep, "Deep")
df_comb <- acc_to_df(acc_comb, "Combined")

acc_df <- bind_rows(df_epi, df_deep, df_comb)

# Compute upper and lower confidence limits
acc_df <- acc_df %>%
  mutate(
    ymin = richness - sd,
    ymax = richness + sd
  )

# Plot
ggplot(acc_df, aes(x = samples, y = richness, color = layer, fill = layer)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.25, color = NA) +
  scale_color_manual(values = c("Epipelagic" = "blue",
                                "Deep" = "red",
                                "Combined" = "#333333")) +
  scale_fill_manual(values = c("Epipelagic" = "blue",
                               "Deep" = "red",
                               "Combined" = "#333333")) +
  labs(
    x = "Number of samples sequenced",
    y = "Number of vOTUs detected",
    title = "Viral richness accumulation curves by depth layer"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )