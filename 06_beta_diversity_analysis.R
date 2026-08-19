library(tidyverse)
library(haven)
library(vegan)
library(tools)

# ------------------------------------------------------------
# Beta-diversity analysis
# WARNING: Computationally intensive.
# Each beta-diversity analysis may take approximately 1 hour to complete,
# depending on available hardware. PERMANOVA and beta-dispersion tests
# are performed using 999 permutations for each distance metric.
# ------------------------------------------------------------

run_beta_analysis <- function(
  beta_file,
  metadata_file,
  metadata_col = "depression_dsm_flag",
  permutations = 999,
  output_dir = "data/processed/publication_images"
) {

  # ------------------------------------------------------------
  # Load data
  # ------------------------------------------------------------

  df_beta <- read.table(beta_file, header = TRUE)

  df_metadata <- read_xpt(metadata_file) %>%
    select(SEQN, all_of(metadata_col)) %>%
    mutate(SEQN = as.character(SEQN))

  # ------------------------------------------------------------
  # Prepare beta-diversity matrix
  # ------------------------------------------------------------

  rownames(df_beta) <- as.character(df_beta$SEQN)
  df_beta$SEQN <- NULL
  colnames(df_beta) <- sub("^X", "", colnames(df_beta))

  df_metadata_clean <- df_metadata %>%
    filter(SEQN %in% rownames(df_beta)) %>%
    filter(!is.na(.data[[metadata_col]])) %>%
    mutate(bucket = factor(.data[[metadata_col]], levels = c(0, 1)))

  matched_ids <- intersect(
    df_metadata_clean$SEQN,
    intersect(rownames(df_beta), colnames(df_beta))
  )

  df_beta_clean <- df_beta[
    matched_ids,
    matched_ids,
    drop = FALSE
  ]

  df_metadata_clean <- df_metadata_clean[
    match(matched_ids, df_metadata_clean$SEQN),
    ,
    drop = FALSE
  ]

  beta_dist_clean <- as.dist(df_beta_clean)

  # ------------------------------------------------------------
  # Beta-diversity metric label
  # ------------------------------------------------------------

  beta_name <- file_path_sans_ext(basename(beta_file))

  beta_label <- case_when(
    str_detect(str_to_lower(beta_name), "bray") ~ "Bray–Curtis dissimilarity",
    str_detect(str_to_lower(beta_name), "unwunifrac|unweighted") ~ "Unweighted UniFrac distance",
    str_detect(str_to_lower(beta_name), "wunifrac|weighted") ~ "Weighted UniFrac distance",
    TRUE ~ beta_name
  )

  # ------------------------------------------------------------
  # PERMANOVA
  # ------------------------------------------------------------

  set.seed(123)

  permanova_result <- adonis2(
    beta_dist_clean ~ bucket,
    data = df_metadata_clean,
    permutations = permutations
  )

  permanova_f <- permanova_result$F[1]
  permanova_r2 <- permanova_result$R2[1]
  permanova_p <- permanova_result$`Pr(>F)`[1]

  # ------------------------------------------------------------
  # Beta dispersion
  # ------------------------------------------------------------

  disp <- betadisper(
    beta_dist_clean,
    df_metadata_clean$bucket
  )

  set.seed(123)

  disp_test <- permutest(
    disp,
    permutations = permutations
  )

  disp_f <- disp_test$tab[1, "F"]
  disp_p <- disp_test$tab[1, "Pr(>F)"]

  # ------------------------------------------------------------
  # Beta-dispersion figure
  # ------------------------------------------------------------

  disp_df <- data.frame(
    SEQN = names(disp$distances),
    Distance = disp$distances
  ) %>%
    left_join(df_metadata_clean, by = "SEQN")

  p_label <- ifelse(
    disp_p < 0.001,
    "p < 0.001",
    paste0("p = ", formatC(disp_p, format = "f", digits = 3))
  )

  p_disp <- ggplot(
    disp_df,
    aes(x = bucket, y = Distance, fill = bucket)
  ) +
    geom_boxplot(alpha = 0.7, linewidth = 0.7) +
    scale_fill_manual(
      values = c(
        "0" = "#FA8072",
        "1" = "#20B2AA"
      )
    ) +
    scale_x_discrete(
      labels = c(
        "0" = "Control",
        "1" = "Depression"
      )
    ) +
    labs(
      x = "Depression status",
      y = "Distance to centroid"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(disp_df$Distance, na.rm = TRUE) * 1.08,
      label = p_label
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(2, 2, 2, 2),
      legend.position = "none"
    )

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  ggsave(
    filename = file.path(
      output_dir,
      paste0(beta_name, "_betadispersion.png")
    ),
    plot = p_disp,
    width = 120,
    height = 90,
    units = "mm",
    dpi = 600,
    bg = "white"
  )

  # ------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------

  tibble(
    beta_metric = beta_label,
    n = nrow(df_metadata_clean),
    permanova_F = permanova_f,
    permanova_R2 = permanova_r2,
    permanova_p = permanova_p,
    betadisper_F = disp_f,
    betadisper_p = disp_p
  )
}

# ------------------------------------------------------------
# Run all beta-diversity analyses
# ------------------------------------------------------------

beta_files <- c(
  "data/beta/dada2rsv-braycurtis-beta.txt",
  "data/beta/dada2rsv-unwunifrac-beta.txt",
  "data/beta/dada2rsv-wunifrac-beta.txt"
)

beta_results <- map_dfr(
  beta_files,
  ~ run_beta_analysis(
    beta_file = .x,
    metadata_file = "data/metadata/processed/df_depression_with_score.xpt"
  )
)

write_csv(
  beta_results,
  "data/processed/beta_diversity_results.csv"
)

print(beta_results)
