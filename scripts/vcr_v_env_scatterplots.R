# Load VCR data from compute_vcr.R
res_cellular_caudo_vp
res_cellular_mega_ve

# Load environmental data from prepare_for_mantel.R
# env_prok_filtered
# env_euk_filtered


library(dplyr)
library(ggplot2)
library(rlang)

plot_vcr_env <- function(
    vcr_type,
    env_variable,
    log_y = TRUE,
    smooth_method = "lm",
    correlation_method = "spearman"
) {
  
  vcr_type <- match.arg(
    vcr_type,
    c("caudo_prok", "mega_euk")
  )
  
  # Select the corresponding VCR and environmental datasets
  if (vcr_type == "caudo_prok") {
    
    vcr_data <- res_cellular_caudo_vp$data %>%
      transmute(
        Run,
        vcr = virus_prok_ratio
      )
    
    env_data <- env_prok_filtered
    vcr_label <- "Caudovirus-to-prokaryote ratio"
    
  } else {
    
    vcr_data <- res_cellular_mega_ve$data %>%
      transmute(
        Run,
        vcr = virus_euk_ratio
      )
    
    env_data <- env_euk_filtered
    vcr_label <- "Megavirus-to-eukaryote ratio"
  }
  
  # Join by Run and remove incomplete, non-finite and zero VCR values
  data <- env_data %>%
    select(
      Run,
      all_of(env_variable),
      env_temp
    ) %>%
    inner_join(
      vcr_data,
      by = "Run"
    ) %>%
    filter(
      !is.na(.data[[env_variable]]),
      !is.na(vcr),
      is.finite(.data[[env_variable]]),
      is.finite(vcr),
      vcr > 0
    )
  
  n_pairs <- nrow(data)
  
  if (n_pairs < 3) {
    stop(
      "Fewer than three complete positive-VCR pairs remain for ",
      env_variable,
      "."
    )
  }
  
  # Correlation uses exactly the same observations shown in the plot
  correlation <- cor.test(
    x = data[[env_variable]],
    y = data$vcr,
    method = correlation_method,
    exact = FALSE
  )
  
  estimate <- unname(correlation$estimate)
  
  estimate_symbol <- switch(
    correlation_method,
    spearman = "\u03c1",
    pearson = "r",
    kendall = "\u03c4",
    "estimate"
  )
  
  annotation <- paste0(
    estimate_symbol, " = ", round(estimate, 3),
    "\np = ", format.pval(
      correlation$p.value,
      digits = 3
    ),
    "\nn = ", n_pairs
  )
  
  p <- ggplot(
    data,
    aes(
      x = .data[[env_variable]],
      y = vcr
    )
  ) +
    geom_point(
      aes(color = env_temp),
      alpha = 0.6
    ) +
    scale_color_gradientn(
      colours = c(
        "darkblue",
        "deepskyblue3",
        "yellow",
        "orange",
        "red"
      ),
      name = "Temperature"
    ) +
    geom_smooth(
      method = smooth_method,
      se = TRUE,
      color = "black"
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = annotation,
      hjust = 1.1,
      vjust = 1.1
    ) +
    labs(
      x = env_variable,
      y = vcr_label
    ) +
    theme_bw()
  
  if (log_y) {
    p <- p +
      scale_y_log10() +
      labs(
        y = paste0(vcr_label, " (log10 scale)")
      )
  }
  
  correlation_summary <- tibble(
    vcr_type = vcr_type,
    environmental_variable = env_variable,
    method = correlation$method,
    estimate = estimate,
    p_value = correlation$p.value,
    n_pairs = n_pairs
  )
  
  list(
    plot = p,
    correlation = correlation,
    correlation_summary = correlation_summary,
    n_pairs = n_pairs,
    data = data
  )
}

library(purrr)
library(patchwork)

env_variables <- c(
  "env_temp",
  "env_phosphate",
  "env_PAR",
  "env_chl_0m_castant"
)

results <- c(
  map(
    env_variables,
    ~ plot_vcr_env(
      "caudo_prok",
      .x,
      log_y = TRUE
    )
  ),
  map(
    env_variables,
    ~ plot_vcr_env(
      "mega_euk",
      .x,
      log_y = TRUE
    )
  )
)

plots <- map(results, "plot")

combined_plot <- wrap_plots(
  plots,
  ncol = 4
) +
  plot_annotation(
    title = "Environmental associations with virus-to-cell ratios",
    tag_levels = "A"
  )

combined_plot