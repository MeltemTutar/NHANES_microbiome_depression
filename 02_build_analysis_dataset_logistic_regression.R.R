library(tidyverse)
library(haven)

# ------------------------------------------------------------
# Load processed depression and alpha-diversity data
# ------------------------------------------------------------

df_depression <- read_xpt(
  "data/metadata/processed/df_depression_with_score.xpt"
) %>%
  select(SEQN, depression_dsm_flag)

df_alpha <- read_xpt(
  "data/alpha/processed/average_alpha.xpt"
) %>%
  select(
    SEQN,
    RSV_ObservedOTUs_6000,
    RSV_FaPhyloDiv_6000,
    RSV_ShanWienDiv_6000,
    RSV_InverseSimpson_6000
  )

# ------------------------------------------------------------
# Helper for loading NHANES cycles
# ------------------------------------------------------------

read_and_tag <- function(path, cycle) {
  read_xpt(path) %>%
    mutate(cycle = cycle)
}

# ------------------------------------------------------------
# Load covariates from NHANES 2009–2012
# ------------------------------------------------------------

demo <- bind_rows(
  read_and_tag("data/metadata/DEMO_F.xpt", "2009-2010"),
  read_and_tag("data/metadata/DEMO_G.xpt", "2011-2012")
)

smq <- bind_rows(
  read_and_tag("data/metadata/SMQ_F.xpt", "2009-2010"),
  read_and_tag("data/metadata/SMQ_G.xpt", "2011-2012")
)

ohq <- bind_rows(
  read_and_tag("data/metadata/OHQ_F.xpt", "2009-2010"),
  read_and_tag("data/metadata/OHQ_G.xpt", "2011-2012")
)

paq <- bind_rows(
  read_and_tag("data/metadata/PAQ_F.xpt", "2009-2010"),
  read_and_tag("data/metadata/PAQ_G.xpt", "2011-2012")
)

bmx <- bind_rows(
  read_and_tag("data/metadata/BMX_F.xpt", "2009-2010"),
  read_and_tag("data/metadata/BMX_G.xpt", "2011-2012")
)

# ------------------------------------------------------------
# Merge datasets by participant identifier
# ------------------------------------------------------------

nhanes_full <- df_depression %>%
  left_join(df_alpha, by = "SEQN") %>%
  left_join(demo, by = "SEQN") %>%
  left_join(paq, by = "SEQN") %>%
  left_join(smq, by = "SEQN") %>%
  left_join(ohq, by = "SEQN") %>%
  left_join(bmx, by = "SEQN")

# ------------------------------------------------------------
# Restrict to participants with oral microbiome data
# ------------------------------------------------------------

analysis_df <- nhanes_full %>%
  filter(!is.na(RSV_ObservedOTUs_6000)) %>%
  select(
    SEQN,
    depression_dsm_flag,

    RSV_ObservedOTUs_6000,
    RSV_FaPhyloDiv_6000,
    RSV_ShanWienDiv_6000,
    RSV_InverseSimpson_6000,

    RIAGENDR,
    RIDAGEYR,
    RIDRETH1,
    DMDEDUC2,
    INDFMPIR,
    SMQ020,
    BMXBMI,
    PAQ605,
    PAQ635,
    PAQ650,
    OHQ845
  )

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write_csv(
  analysis_df,
  "data/processed/nhanes_analysis_dataset_filtered.csv"
)
