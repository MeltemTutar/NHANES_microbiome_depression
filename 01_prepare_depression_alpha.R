library(tidyverse)
library(haven)

# ------------------------------------------------------------
# Depression data: NHANES 2009–2012
# ------------------------------------------------------------

df_depression <- bind_rows(
  read_xpt("data/metadata/DPQ_F.xpt"),
  read_xpt("data/metadata/DPQ_G.xpt")
)

# Require complete responses to DPQ010–DPQ090
dpq_cols_1_to_9 <- paste0("DPQ", sprintf("%03d", seq(10, 90, by = 10)))

df_depression <- df_depression[
  complete.cases(df_depression[, dpq_cols_1_to_9]),
]

# PHQ-9 summed score
df_depression$depression_summed_score <- rowSums(
  df_depression[, dpq_cols_1_to_9]
)

df_depression$depression_summed_flag <-
  df_depression$depression_summed_score >= 10

# ------------------------------------------------------------
# DSM-IV depression algorithm
# ------------------------------------------------------------

# At least one core symptom:
# DPQ010 (depressed mood) or DPQ020 (anhedonia) scored >= 2
core_symptom <- df_depression$DPQ010 >= 2 |
                df_depression$DPQ020 >= 2

# DPQ010–DPQ080 count when scored >= 2
dpq_cols_1_to_8 <- paste0(
  "DPQ",
  sprintf("%03d", seq(10, 80, by = 10))
)

symptom_flags <- df_depression[, dpq_cols_1_to_8] >= 2

# DPQ090 (suicidal ideation) counts when scored >= 1
dpq090_flag <- df_depression$DPQ090 >= 1

qualifying_symptom_count <-
  rowSums(symptom_flags) + dpq090_flag

# DSM-IV depression flag
df_depression$depression_dsm_flag <-
  core_symptom & qualifying_symptom_count >= 5

dir.create("data/metadata/processed", recursive = TRUE, showWarnings = FALSE)

write_xpt(
  df_depression,
  "data/metadata/processed/df_depression_with_score.xpt"
)

# ------------------------------------------------------------
# Alpha-diversity data
# ------------------------------------------------------------

df_alpha <- read.table(
  "data/alpha/dada2rsv-alpha.txt",
  header = TRUE
)

create_averaged_alpha <- function(df_alpha) {

  prefixes <- unique(
    sub("_\\d+$", "", colnames(df_alpha)[-1])
  )

  averaged_df <- data.frame(SEQN = df_alpha$SEQN)

  for (prefix in prefixes) {

    cols <- grep(
      paste0("^", prefix, "_\\d+$"),
      colnames(df_alpha),
      value = TRUE
    )

    numeric_data <- data.frame(
      lapply(
        df_alpha[cols],
        function(x) as.numeric(as.character(x))
      )
    )

    averaged_df[[prefix]] <- rowMeans(
      numeric_data,
      na.rm = TRUE
    )
  }

  averaged_df
}

df_alpha_average <- create_averaged_alpha(df_alpha)

dir.create("data/alpha/processed", recursive = TRUE, showWarnings = FALSE)

write_xpt(
  df_alpha_average,
  "data/alpha/processed/average_alpha.xpt"
)
