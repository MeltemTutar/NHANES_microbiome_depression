library(tidyverse)
library(broom)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

df <- read_csv(
  "data/processed/nhanes_analysis_dataset_filtered.csv",
  show_col_types = FALSE
)

var_dict <- read_csv(
  "data/processed/nhanes_variable_dictionary_corrected.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Define variables
# ------------------------------------------------------------

alpha_vars <- c(
  "RSV_ObservedOTUs_6000",
  "RSV_FaPhyloDiv_6000",
  "RSV_ShanWienDiv_6000",
  "RSV_InverseSimpson_6000"
)

outcome <- "depression_dsm_flag"

covariates <- c(
  "RIAGENDR",
  "RIDAGEYR",
  "RIDRETH1",
  "DMDEDUC2",
  "INDFMPIR",
  "SMQ020",
  "BMXBMI",
  "PAQ605",
  "PAQ635",
  "PAQ650",
  "OHQ845"
)

categorical_vars <- c(
  "RIAGENDR",
  "RIDRETH1",
  "DMDEDUC2",
  "SMQ020",
  "PAQ605",
  "PAQ635",
  "PAQ650",
  "OHQ845"
)

continuous_vars <- c(
  "RIDAGEYR",
  "INDFMPIR",
  "BMXBMI"
)

# ------------------------------------------------------------
# Convert NHANES missing-value codes to NA
# ------------------------------------------------------------

invalid_codes <- c(
  7, 9, 77, 99, 777, 999, 9999
)

df <- df %>%
  mutate(
    across(
      all_of(categorical_vars),
      ~ replace(.x, .x %in% invalid_codes, NA)
    )
  )

# ------------------------------------------------------------
# Impute missing covariates
# ------------------------------------------------------------

# Continuous covariates: mean imputation
df_imputed <- df %>%
  mutate(
    across(
      all_of(continuous_vars),
      ~ replace(.x, is.na(.x), mean(.x, na.rm = TRUE))
    )
  )

# Categorical covariates: mode imputation
mode_impute <- function(x) {
  observed <- x[!is.na(x)]
  names(sort(table(observed), decreasing = TRUE))[1]
}

df_imputed <- df_imputed %>%
  mutate(
    across(
      all_of(categorical_vars),
      ~ replace(
        as.character(.x),
        is.na(.x),
        mode_impute(.x)
      )
    )
  ) %>%
  mutate(
    across(
      all_of(categorical_vars),
      as.factor
    )
  )

# ------------------------------------------------------------
# Fit one logistic regression model per alpha-diversity metric
# ------------------------------------------------------------

results_list <- list()

for (alpha in alpha_vars) {

  model_formula <- reformulate(
    c(alpha, covariates),
    response = outcome
  )

  model <- glm(
    model_formula,
    data = df_imputed,
    family = binomial
  )

  results_list[[alpha]] <- tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    mutate(alpha_metric = alpha)
}

# ------------------------------------------------------------
# Combine model results
# ------------------------------------------------------------

final_results <- bind_rows(results_list) %>%
  select(
    alpha_metric,
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  ) %>%
  rename(
    `Alpha Metric` = alpha_metric,
    Variable = term,
    `Odds Ratio` = estimate,
    `CI Lower` = conf.low,
    `CI Upper` = conf.high,
    `p-value` = p.value
  )

write_csv(
  final_results,
  "data/processed/results_logistic_regression_alpha_metrics.csv"
)
