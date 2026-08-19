library(tidyverse)
library(haven)

# ------------------------------------------------------------
# Alpha-diversity analysis
# ------------------------------------------------------------

run_alpha_analysis <- function(
  alpha_file,
  metadata_file,
  alpha_col,
  metadata_col = "depression_dsm_flag",
  output_dir = "data/processed/publication_images"
) {

  # Load processed alpha-diversity and depression data
  df_alpha <- read_xpt(alpha_file) %>%
    select(SEQN, all_of(alpha_col))

  df_metadata <- read_xpt(metadata_file) %>%
    select(SEQN, all_of(metadata_col))

  # Merge datasets
  stratified_df <- inner_join(df_alpha, df_metadata, by = "SEQN") %>%
    filter(!is.na(.data[[alpha_col]]), !is.na(.data[[metadata_col]])) %>%
    mutate(Strata = factor(.data[[metadata_col]], levels = c(0, 1)))

  # Publication-friendly metric name
  alpha_label <- case_when(
    alpha_col == "RSV_ObservedOTUs_6000" ~ "Observed OTUs",
    alpha_col == "RSV_FaPhyloDiv_6000" ~ "Faith's phylogenetic diversity",
    alpha_col == "RSV_ShanWienDiv_6000" ~ "Shannon diversity",
    alpha_col == "RSV_InverseSimpson_6000" ~ "Inverse Simpson diversity",
    TRUE ~ alpha_col
  )

  # ------------------------------------------------------------
  # One-way ANOVA
  # ------------------------------------------------------------

  anova_result <- aov(
    reformulate("Strata", response = alpha_col),
    data = stratified_df
  )

  p_value <- summary(anova_result)[[1]][1, "Pr(>F)"]

  # ------------------------------------------------------------
  # Cohen's d
  # Depression - Control
  # ------------------------------------------------------------

  control <- stratified_df %>%
    filter(Strata == "0") %>%
    pull(all_of(alpha_col))

  depression <- stratified_df %>%
    filter(Strata == "1") %>%
    pull(all_of(alpha_col))

  n_control <- length(control)
  n_depression <- length(depression)

  pooled_sd <- sqrt(
    ((n_control - 1) * var(control) +
       (n_depression - 1) * var(depression)) /
      (n_control + n_depression - 2)
  )

  cohens_d <- (mean(depression) - mean(control)) / pooled_sd

  # ------------------------------------------------------------
  # Group summary statistics
  # ------------------------------------------------------------

  group_summary <- stratified_df %>%
    group_by(Strata) %>%
    summarise(
      n = n(),
      mean = mean(.data[[alpha_col]], na.rm = TRUE),
      sd = sd(.data[[alpha_col]], na.rm = TRUE),
      median = median(.data[[alpha_col]], na.rm = TRUE),
      IQR = IQR(.data[[alpha_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(group = recode(as.character(Strata), "0" = "Control", "1" = "Depression"))

  # ------------------------------------------------------------
  # Publication figure
  # ------------------------------------------------------------

  p_label <- ifelse(
    p_value < 0.001,
    "p < 0.001",
    paste0("p = ", formatC(p_value, format = "f", digits = 3))
  )

  p_alpha <- ggplot(stratified_df, aes(x = Strata, y = .data[[alpha_col]], fill = Strata)) +
    geom_boxplot(alpha = 0.7, linewidth = 0.7) +
    scale_fill_manual(values = c("0" = "#FA8072", "1" = "#20B2AA")) +
    scale_x_discrete(labels = c("0" = "Control", "1" = "Depression")) +
    labs(x = "Depression status", y = alpha_label) +
    annotate(
      "text",
      x = 1.5,
      y = max(stratified_df[[alpha_col]], na.rm = TRUE) * 1.08,
      label = p_label
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(2, 2, 2, 2),
      legend.position = "none"
    )

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  ggsave(
    filename = file.path(output_dir, paste0(alpha_col, "_alpha.png")),
    plot = p_alpha,
    width = 120,
    height = 90,
    units = "mm",
    dpi = 600,
    bg = "white"
  )

  # ------------------------------------------------------------
  # Return manuscript results
  # ------------------------------------------------------------

  tibble(
    alpha_metric = alpha_label,
    anova_p = p_value,
    cohens_d = cohens_d,
    control_n = group_summary$n[group_summary$group == "Control"],
    control_mean = group_summary$mean[group_summary$group == "Control"],
    control_sd = group_summary$sd[group_summary$group == "Control"],
    depression_n = group_summary$n[group_summary$group == "Depression"],
    depression_mean = group_summary$mean[group_summary$group == "Depression"],
    depression_sd = group_summary$sd[group_summary$group == "Depression"]
  )
}

# ------------------------------------------------------------
# Run alpha-diversity analyses
# ------------------------------------------------------------

alpha_vars <- c(
  "RSV_ObservedOTUs_6000",
  "RSV_FaPhyloDiv_6000",
  "RSV_ShanWienDiv_6000",
  "RSV_InverseSimpson_6000"
)

alpha_results <- map_dfr(
  alpha_vars,
  ~ run_alpha_analysis(
    alpha_file = "data/alpha/processed/average_alpha.xpt",
    metadata_file = "data/metadata/processed/df_depression_with_score.xpt",
    alpha_col = .x
  )
)

# Save results
write_csv(
  alpha_results,
  "data/processed/alpha_diversity_results.csv"
)

# Display results
print(alpha_results)
