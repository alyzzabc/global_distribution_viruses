## Viruses
#run_map_vir <- metadata[, .(sample_id, latitude, depth, size_frac, env)]
#run_map_vir <- run_map_vir[env == "seawater"]
#run_map_vir[, depth_num := as.numeric(depth)]

#run_map_vir[, layer := fifelse(
#  is.na(depth_num),
#  NA_character_,
#  fifelse(depth_num < 200, "epi", "deep")
#)]

# Attach sample_id to viral table
gv_meta <- run_map_vir[
  gv,
  on = .(sample_id = virus_sample_id),
  nomatch = 0L
]

# Load taxonomy
taxonomy <- read_tsv("/path/to/virus/taxonomy/file")
taxonomy_split <- taxonomy %>%
  separate(taxonomy, into = c("domain", "realm", "kingdom", "phylum", "class", "order", "family"), sep = ";", fill = "right", extra = "merge")

gv_meta_tax <- gv_meta %>%
  left_join(
    taxonomy_split %>% select(seq_name, class),
    by = c("virus_genome" = "seq_name")
  )


sum(!unique(gv_meta_tax$virus_genome) %in% taxonomy_split$seq_name)


## Prokaryotes
#run_map_prok <- metadata[, .(sample_id, size_frac, env, depth)]
#run_map_prok <- run_map_prok[env == "seawater"]
#run_map_prok[, depth_num := as.numeric(depth)]

#run_map_prok[, layer := fifelse(
#  is.na(depth_num),
#  NA_character_,
#  fifelse(depth_num < 200, "epi", "deep")
#)]

# Attach sample_id to prokaryote table
gp_meta <- run_map_prok[
  gp,
  on = .(sample_id = prok_sample_id),
  nomatch = 0L
]


## Eukaryotes
#run_map_euk <- metadata[, .(Run, sample_id, size_frac, env, depth)]
#run_map_euk <- run_map_euk[env == "seawater"]
#run_map_euk[, depth_num := as.numeric(depth)]

#run_map_euk[, layer := fifelse(
#    is.na(depth_num),
#    NA_character_,
#    fifelse(depth_num < 200, "epi", "deep")
#)]


# Attach sample_id to eukaryote table
ge_meta <- run_map_euk[
  ge,
  on = .(sample_id = euk_sample_id),
  nomatch = 0L
]

## Plot
library(data.table)
library(ggplot2)

plot_prok_virus_ratio <- function(
    gp_meta,
    gv_meta_tax,
    target_size_frac = NULL,
    virus_class = NULL,
    target_layer = NULL,
    max_depth = 200,
    ratio_direction = c("prok_per_virus", "virus_per_prok"),
    use_pseudocount = FALSE,
    pseudocount = 1e-6,
    log_y = TRUE,
    loess_span = 0.75,
    cap_quantile = NULL,
    smooth_lat_min = NULL,
    smooth_lat_max = NULL,
    point_alpha = 0.5,
    point_size = 2
) {
  ratio_direction <- match.arg(ratio_direction)
  
  gp <- copy(as.data.table(gp_meta))
  gv <- copy(as.data.table(gv_meta_tax))
  
  # Optionally restrict analysis to selected layer(s)
  if (!is.null(target_layer)) {
    gp <- gp[layer %in% target_layer]
    gv <- gv[layer %in% target_layer]
  }
  
  # Make sure depth is numeric
  if (!"depth_num" %in% names(gv)) {
    gv[, depth_num := as.numeric(depth)]
  }
  
  # If no size_frac are specified, use all size_frac shared by both tables
  if (is.null(target_size_frac)) {
    target_size_frac <- sort(intersect(unique(gp$size_frac), unique(gv$size_frac)))
    size_frac_label <- "all size_frac"
  } else {
    size_frac_label <- paste(target_size_frac, collapse = " + ")
  }
  
  # 1. Define sampled epipelagic Runs using selected size_frac
  # Do NOT filter by virus_class here.
  # This keeps samples with no target virus class as true zeroes.
  sample_grid <- gv[
    depth_num <= max_depth &
      size_frac %in% target_size_frac,
    .(
      latitude = first(na.omit(latitude)),
      depth = first(na.omit(depth)),
      depth_num = first(na.omit(depth_num))
    ),
    by = Run
  ]
  
  # 2. Sum prokaryote norm_coverage across selected size_frac per Run
  prok_sum <- gp[
    size_frac %in% target_size_frac,
    .(
      prok_norm_cov = sum(norm_coverage, na.rm = TRUE)
    ),
    by = Run
  ]
  
  # 3. Filter virus table by depth/size_frac/class
  gv_selected <- gv[
    depth_num <= max_depth &
      size_frac %in% target_size_frac
  ]
  
  if (!is.null(virus_class)) {
    gv_selected <- gv_selected[
      class %in% virus_class
    ]
  }
  
  # 4. Sum viral norm_coverage across selected size_frac per Run
  virus_sum <- gv_selected[
    ,
    .(
      virus_norm_cov = sum(norm_coverage, na.rm = TRUE)
    ),
    by = Run
  ]
  
  # 5. Join onto sampled Run scaffold
  ratio_dt <- merge(
    sample_grid,
    prok_sum,
    by = "sample_id",
    all.x = TRUE
  )
  
  ratio_dt <- merge(
    ratio_dt,
    virus_sum,
    by = "sample_id",
    all.x = TRUE
  )
  
  # 6. Fill missing detections as true zeroes
  ratio_dt[is.na(prok_norm_cov), prok_norm_cov := 0]
  ratio_dt[is.na(virus_norm_cov), virus_norm_cov := 0]
  
  # 7. Calculate ratio
  if (use_pseudocount) {
    
    if (ratio_direction == "prok_per_virus") {
      ratio_dt[
        ,
        ratio := (prok_norm_cov + pseudocount) /
          (virus_norm_cov + pseudocount)
      ]
      ratio_label <- "Prokaryote:virus ratio"
      ratio_column <- "prok_virus_ratio"
      
    } else if (ratio_direction == "virus_per_prok") {
      ratio_dt[
        ,
        ratio := (virus_norm_cov + pseudocount) /
          (prok_norm_cov + pseudocount)
      ]
      ratio_label <- "Virus:prokaryote ratio"
      ratio_column <- "virus_prok_ratio"
    }
    
  } else {
    
    if (ratio_direction == "prok_per_virus") {
      ratio_dt[
        virus_norm_cov > 0,
        ratio := prok_norm_cov / virus_norm_cov
      ]
      ratio_dt[
        virus_norm_cov == 0,
        ratio := NA_real_
      ]
      ratio_label <- "Prokaryote:virus ratio"
      ratio_column <- "prok_virus_ratio"
      
    } else if (ratio_direction == "virus_per_prok") {
      ratio_dt[
        prok_norm_cov > 0,
        ratio := virus_norm_cov / prok_norm_cov
      ]
      ratio_dt[
        prok_norm_cov == 0,
        ratio := NA_real_
      ]
      ratio_label <- "Virus:prokaryote ratio"
      ratio_column <- "virus_prok_ratio"
    }
  }
  
  # Keep direction-specific ratio column too
  ratio_dt[, (ratio_column) := ratio]
  
  # 8. Optional capping
  if (!is.null(cap_quantile)) {
    cap_value <- quantile(
      ratio_dt$ratio,
      probs = cap_quantile,
      na.rm = TRUE
    )
    
    ratio_dt[
      !is.na(ratio) & ratio > cap_value,
      ratio_capped := cap_value
    ]
    
    ratio_dt[
      !is.na(ratio) & ratio <= cap_value,
      ratio_capped := ratio
    ]
    
    y_col <- "ratio_capped"
    cap_label <- paste0(", capped at q", cap_quantile)
    
  } else {
    y_col <- "ratio"
    cap_value <- NA_real_
    cap_label <- ""
  }
  
  # 9. Plotting data
  plot_dt <- ratio_dt[
    !is.na(get(y_col)) &
      get(y_col) > 0
  ]
  
  # 10. Smoothing data, optionally restricted by latitude
  smooth_dt <- copy(plot_dt)
  
  if (!is.null(smooth_lat_min)) {
    smooth_dt <- smooth_dt[
      latitude >= smooth_lat_min
    ]
  }
  
  if (!is.null(smooth_lat_max)) {
    smooth_dt <- smooth_dt[
      latitude <= smooth_lat_max
    ]
  }
  
  virus_label <- ifelse(
    is.null(virus_class),
    "all viruses",
    paste(virus_class, collapse = " + ")
  )
  
  pseudocount_label <- ifelse(
    use_pseudocount,
    paste0(", pseudocount = ", pseudocount),
    ""
  )
  
  smooth_label <- ""
  if (!is.null(smooth_lat_min) | !is.null(smooth_lat_max)) {
    smooth_label <- paste0(
      ", LOESS fit ",
      ifelse(is.null(smooth_lat_min), "-Inf", smooth_lat_min),
      " to ",
      ifelse(is.null(smooth_lat_max), "Inf", smooth_lat_max),
      "°"
    )
  }
  
  p <- ggplot(
    plot_dt,
    aes(x = latitude, y = get(y_col))
  ) +
    geom_point(alpha = point_alpha, size = point_size) +
    geom_smooth(
      data = smooth_dt,
      aes(x = latitude, y = get(y_col)),
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      span = loess_span,
      color = "red",
      fill = "red",
      alpha = 0.2,
      linewidth = 1
    ) +
    theme_bw() +
    labs(
      x = "Latitude",
      y = ratio_label,
      title = paste0(
        ratio_label,
        " across latitude; ",
        size_frac_label,
        "; ",
        virus_label,
        "; ≤",
        max_depth,
        " m",
        pseudocount_label,
        cap_label,
        smooth_label
      )
    )
  
  if (log_y) {
    p <- p +
      scale_y_log10() +
      labs(y = paste0(ratio_label, ", log10 scale"))
  }
  
  return(
    list(
      data = ratio_dt,
      plot_data = plot_dt,
      smooth_data = smooth_dt,
      plot = p,
      size_frac_used = target_size_frac,
      virus_class_used = virus_class,
      ratio_direction = ratio_direction,
      ratio_column = ratio_column,
      y_column_plotted = y_col,
      use_pseudocount = use_pseudocount,
      pseudocount = pseudocount,
      cap_quantile = cap_quantile,
      cap_value = cap_value,
      smooth_lat_min = smooth_lat_min,
      smooth_lat_max = smooth_lat_max
    )
  )
}


plot_euk_virus_ratio <- function(
    ge_meta,
    gv_meta_tax,
    target_size_frac = NULL,
    virus_class = NULL,
    target_layer = NULL,
    max_depth = 200,
    ratio_direction = c("euk_per_virus", "virus_per_euk"),
    use_pseudocount = FALSE,
    pseudocount = 1e-6,
    log_y = TRUE,
    loess_span = 0.75,
    cap_quantile = NULL,
    smooth_lat_min = NULL,
    smooth_lat_max = NULL,
    point_alpha = 0.5,
    point_size = 2
) {
  ratio_direction <- match.arg(ratio_direction)
  
  ge <- copy(as.data.table(ge_meta))
  gv <- copy(as.data.table(gv_meta_tax))
  
  # Optionally restrict analysis to selected layer(s)
  if (!is.null(target_layer)) {
    gp <- gp[layer %in% target_layer]
    gv <- gv[layer %in% target_layer]
  }
  
  # Make sure depth is numeric
  if (!"depth_num" %in% names(gv)) {
    gv[, depth_num := as.numeric(depth)]
  }
  
  # If no size_frac are specified, use all size_frac shared by both tables
  if (is.null(target_size_frac)) {
    target_size_frac <- sort(intersect(unique(ge$size_frac), unique(gv$size_frac)))
    size_frac_label <- "all size_frac"
  } else {
    size_frac_label <- paste(target_size_frac, collapse = " + ")
  }
  
  # 1. Define sampled epipelagic Runs using selected size_frac.
  # Do NOT filter by virus_class here.
  # This keeps samples with no target virus class as true zeroes.
  sample_grid <- gv[
    depth_num <= max_depth &
      size_frac %in% target_size_frac,
    .(
      latitude = first(na.omit(latitude)),
      depth = first(na.omit(depth)),
      depth_num = first(na.omit(depth_num))
    ),
    by = Run
  ]
  
  # 2. Sum eukaryote norm_coverage across selected size_frac per Run
  euk_sum <- ge[
    size_frac %in% target_size_frac,
    .(
      euk_norm_cov = sum(norm_coverage, na.rm = TRUE)
    ),
    by = Run
  ]
  
  # 3. Filter virus table by depth/size_frac/class
  gv_selected <- gv[
    depth_num <= max_depth &
      size_frac %in% target_size_frac
  ]
  
  if (!is.null(virus_class)) {
    gv_selected <- gv_selected[
      class %in% virus_class
    ]
  }
  
  # 4. Sum viral norm_coverage across selected size_frac per Run
  virus_sum <- gv_selected[
    ,
    .(
      virus_norm_cov = sum(norm_coverage, na.rm = TRUE)
    ),
    by = Run
  ]
  
  # 5. Join onto sampled Run scaffold
  ratio_dt <- merge(
    sample_grid,
    euk_sum,
    by = "sample_id",
    all.x = TRUE
  )
  
  ratio_dt <- merge(
    ratio_dt,
    virus_sum,
    by = "sample_id",
    all.x = TRUE
  )
  
  # 6. Fill missing detections as true zeroes
  ratio_dt[is.na(euk_norm_cov), euk_norm_cov := 0]
  ratio_dt[is.na(virus_norm_cov), virus_norm_cov := 0]
  
  # 7. Calculate ratio
  if (use_pseudocount) {
    
    if (ratio_direction == "euk_per_virus") {
      ratio_dt[
        ,
        ratio := (euk_norm_cov + pseudocount) /
          (virus_norm_cov + pseudocount)
      ]
      ratio_label <- "Eukaryote:virus ratio"
      ratio_column <- "euk_virus_ratio"
      
    } else if (ratio_direction == "virus_per_euk") {
      ratio_dt[
        ,
        ratio := (virus_norm_cov + pseudocount) /
          (euk_norm_cov + pseudocount)
      ]
      ratio_label <- "Virus:eukaryote ratio"
      ratio_column <- "virus_euk_ratio"
    }
    
  } else {
    
    if (ratio_direction == "euk_per_virus") {
      ratio_dt[
        virus_norm_cov > 0,
        ratio := euk_norm_cov / virus_norm_cov
      ]
      ratio_dt[
        virus_norm_cov == 0,
        ratio := NA_real_
      ]
      ratio_label <- "Eukaryote:virus ratio"
      ratio_column <- "euk_virus_ratio"
      
    } else if (ratio_direction == "virus_per_euk") {
      ratio_dt[
        euk_norm_cov > 0,
        ratio := virus_norm_cov / euk_norm_cov
      ]
      ratio_dt[
        euk_norm_cov == 0,
        ratio := NA_real_
      ]
      ratio_label <- "Virus:eukaryote ratio"
      ratio_column <- "virus_euk_ratio"
    }
  }
  
  # Keep direction-specific ratio column too
  ratio_dt[, (ratio_column) := ratio]
  
  # 8. Optional capping
  if (!is.null(cap_quantile)) {
    cap_value <- quantile(
      ratio_dt$ratio,
      probs = cap_quantile,
      na.rm = TRUE
    )
    
    ratio_dt[
      !is.na(ratio) & ratio > cap_value,
      ratio_capped := cap_value
    ]
    
    ratio_dt[
      !is.na(ratio) & ratio <= cap_value,
      ratio_capped := ratio
    ]
    
    y_col <- "ratio_capped"
    cap_label <- paste0(", capped at q", cap_quantile)
    
  } else {
    y_col <- "ratio"
    cap_value <- NA_real_
    cap_label <- ""
  }
  
  # 9. Plotting data
  plot_dt <- ratio_dt[
    !is.na(get(y_col)) &
      get(y_col) > 0
  ]
  
  # 10. Smoothing data, optionally restricted by latitude
  smooth_dt <- copy(plot_dt)
  
  if (!is.null(smooth_lat_min)) {
    smooth_dt <- smooth_dt[
      latitude >= smooth_lat_min
    ]
  }
  
  if (!is.null(smooth_lat_max)) {
    smooth_dt <- smooth_dt[
      latitude <= smooth_lat_max
    ]
  }
  
  virus_label <- ifelse(
    is.null(virus_class),
    "all viruses",
    paste(virus_class, collapse = " + ")
  )
  
  pseudocount_label <- ifelse(
    use_pseudocount,
    paste0(", pseudocount = ", pseudocount),
    ""
  )
  
  smooth_label <- ""
  if (!is.null(smooth_lat_min) | !is.null(smooth_lat_max)) {
    smooth_label <- paste0(
      ", LOESS fit ",
      ifelse(is.null(smooth_lat_min), "-Inf", smooth_lat_min),
      " to ",
      ifelse(is.null(smooth_lat_max), "Inf", smooth_lat_max),
      "°"
    )
  }
  
  p <- ggplot(
    plot_dt,
    aes(x = latitude, y = get(y_col))
  ) +
    geom_point(alpha = point_alpha, size = point_size) +
    geom_smooth(
      data = smooth_dt,
      aes(x = latitude, y = get(y_col)),
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      span = loess_span,
      color = "red",
      fill = "red",
      alpha = 0.2,
      linewidth = 1
    ) +
    theme_bw() +
    labs(
      x = "Latitude",
      y = ratio_label,
      title = paste0(
        ratio_label,
        " across latitude; ",
        size_frac_label,
        "; ",
        virus_label,
        "; ≤",
        max_depth,
        " m",
        pseudocount_label,
        cap_label,
        smooth_label
      )
    )
  
  if (log_y) {
    p <- p +
      scale_y_log10() +
      labs(y = paste0(ratio_label, ", log10 scale"))
  }
  
  return(
    list(
      data = ratio_dt,
      plot_data = plot_dt,
      smooth_data = smooth_dt,
      plot = p,
      size_frac_used = target_size_frac,
      virus_class_used = virus_class,
      ratio_direction = ratio_direction,
      ratio_column = ratio_column,
      y_column_plotted = y_col,
      use_pseudocount = use_pseudocount,
      pseudocount = pseudocount,
      cap_quantile = cap_quantile,
      cap_value = cap_value,
      smooth_lat_min = smooth_lat_min,
      smooth_lat_max = smooth_lat_max
    )
  )
}

## Virus per prokaryote
# Smoother only until -75
res_cellular_caudo_vp <- plot_prok_virus_ratio(
  gp_meta = gp_meta,
  gv_meta_tax = gv_meta_tax,
  target_size_frac = "cellular",
  virus_class = "Caudoviricetes",
  ratio_direction = "virus_per_prok",
  max_depth = 200,
  use_pseudocount = FALSE,
  log_y = TRUE,
  loess_span = 0.75,
  #cap_quantile = 0.95,
  smooth_lat_min = -75
)

## Virus per eukaryote
# Smoother only until -75
res_cellular_mega_ve <- plot_euk_virus_ratio(
  ge_meta = ge_meta,
  gv_meta_tax = gv_meta_tax,
  target_size_frac = "cellular",
  virus_class = "Megaviricetes",
  ratio_direction = "virus_per_euk",
  max_depth = 200,
  use_pseudocount = FALSE,
  log_y = TRUE,
  loess_span = 0.75,
  #cap_quantile = 0.95,
  smooth_lat_min = -75
)

res_cellular_caudo_vp$plot
res_cellular_mega_ve$plot



plot_ratio_latitude_colored <- function(
    res,
    env_data,
    point_alpha = 0.7,
    point_size = 2,
    smooth_color = "black",
    smooth_fill = "grey40",
    smooth_alpha = 0.25
) {
  
  if (!is.list(res)) {
    stop("res must be a result list from plot_prok_virus_ratio() or plot_euk_virus_ratio().")
  }
  
  required_fields <- c("plot", "plot_data", "y_column_plotted")
  missing_fields <- setdiff(required_fields, names(res))
  
  if (length(missing_fields) > 0) {
    stop(
      "res is missing required fields: ",
      paste(missing_fields, collapse = ", ")
    )
  }
  
  # Copy plotting data
  plot_data <- data.table::copy(
    data.table::as.data.table(res$plot_data)
  )
  
  y_col <- res$y_column_plotted
  
  if (!"Run" %in% names(plot_data)) {
    stop("plot_data must contain Run so temperature can be joined.")
  }
  
  if (!"latitude" %in% names(plot_data)) {
    stop("plot_data must contain latitude.")
  }
  
  if (!y_col %in% names(plot_data)) {
    stop("plot_data does not contain the plotted y-column: ", y_col)
  }
  
  if (!all(c("Run", "env_temp") %in% names(env_data))) {
    stop("env_data must contain Run and env_temp.")
  }
  
  # Add temperature by matching Run
  temp_data <- env_data %>%
    dplyr::select(
      Run,
      env_temp
    ) %>%
    dplyr::distinct(Run, .keep_all = TRUE)
  
  plot_data <- plot_data %>%
    dplyr::left_join(
      temp_data,
      by = "Run"
    )
  
  # Remove rows without valid temperature
  plot_data <- plot_data %>%
    dplyr::filter(
      !is.na(env_temp),
      is.finite(env_temp)
    )
  
  # Start from original plot
  p <- res$plot
  
  # Keep your existing smoother styling
  smooth_layers <- vapply(
    p$layers,
    function(x) inherits(x$geom, "GeomSmooth"),
    logical(1)
  )
  
  for (i in which(smooth_layers)) {
    p$layers[[i]]$aes_params$colour <- smooth_color
    p$layers[[i]]$aes_params$fill <- smooth_fill
    p$layers[[i]]$aes_params$alpha <- smooth_alpha
  }
  
  # Remove original point layer
  point_layers <- vapply(
    p$layers,
    function(x) inherits(x$geom, "GeomPoint"),
    logical(1)
  )
  
  p$layers <- p$layers[!point_layers]
  
  # Add points colored by temperature
  p <- p +
    ggplot2::geom_point(
      data = plot_data,
      mapping = ggplot2::aes(
        x = latitude,
        y = .data[[y_col]],
        color = env_temp
      ),
      inherit.aes = FALSE,
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::scale_color_gradientn(
      colours = c(
        "darkblue",
        "deepskyblue3",
        "yellow",
        "orange",
        "red"
      ),
      name = "Temperature"
    )
  
  return(p)
}

caudo_colored <- plot_ratio_latitude_colored(
  res = res_cellular_caudo_vp,
  env_data = env_prok_filtered # <--- from prepare_for_mantel.R
)

caudo_colored

mega_colored <- plot_ratio_latitude_colored(
  res = res_cellular_mega_ve,
  env_data = env_euk_filtered # <--- from prepare_for_mantel.R
)

mega_colored