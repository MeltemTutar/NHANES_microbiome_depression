library(tidyverse)

# ------------------------------------------------------------
# Load analysis dataset and variable dictionary
# ------------------------------------------------------------

df <- read_csv(
  "data/processed/nhanes_analysis_dataset_filtered.csv",
  show_col_types = FALSE
)

dict <- read_csv(
  "data/processed/nhanes_variable_dictionary_corrected.csv",
  show_col_types = FALSE
) %>%
  filter(`Variable Name` != "SEQN") %>%
  rename(
    Variable = `Variable Name`,
    Description = Description
  )

continuous_vars <- dict %>%
  filter(Type == "continuous") %>%
  pull(Variable)

categorical_vars <- dict %>%
  filter(Type == "categorical") %>%
  pull(Variable)

# ------------------------------------------------------------
# Continuous variables
# ------------------------------------------------------------

cont_summary <- df %>%
  select(all_of(continuous_vars)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  group_by(Variable) %>%
  summarise(
    N = sum(!is.na(Value)),
    Missing = sum(is.na(Value)),
    `Missing (%)` = round(mean(is.na(Value)) * 100, 1),
    Mean = round(mean(Value, na.rm = TRUE), 2),
    SD = round(sd(Value, na.rm = TRUE), 2),
    Min = round(min(Value, na.rm = TRUE), 2),
    Max = round(max(Value, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  left_join(dict, by = "Variable") %>%
  select(Variable, Description, everything())

# ------------------------------------------------------------
# Categorical variables
# ------------------------------------------------------------

cat_summary <- map_dfr(
  categorical_vars,
  function(variable) {

    x <- df[[variable]]

    missing_n <- sum(is.na(x))
    missing_pct <- round(mean(is.na(x)) * 100, 1)
    valid_n <- sum(!is.na(x))

    tibble(Category = as.character(x)) %>%
      filter(!is.na(Category)) %>%
      count(Category, name = "Count") %>%
      mutate(
        Variable = variable,
        Percent = round(100 * Count / valid_n, 1),
        Missing = missing_n,
        `Missing (%)` = missing_pct
      )
  }
) %>%
  left_join(dict, by = "Variable") %>%
  select(
    Variable,
    Description,
    Category,
    Count,
    Percent,
    Missing,
    `Missing (%)`
  ) %>%
  arrange(Variable, desc(Count))

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------

write_csv(
  cont_summary,
  "data/processed/summary_continuous.csv"
)

write_csv(
  cat_summary,
  "data/processed/summary_categorical.csv"
)
