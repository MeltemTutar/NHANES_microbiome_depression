# Oral Microbiome Diversity and Depression in NHANES

This repository contains the R code used to analyze associations between oral microbiome diversity and depression among adults participating in the 2009–2012 National Health and Nutrition Examination Survey (NHANES).

## Study Overview

The analysis evaluates associations between depression and oral microbiome diversity using publicly available NHANES 2009–2012 data.

Alpha diversity was assessed using:
- Observed OTUs
- Faith's phylogenetic diversity
- Shannon diversity
- Inverse Simpson diversity

Beta diversity was assessed using:
- Bray–Curtis dissimilarity
- Unweighted UniFrac distance
- Weighted UniFrac distance

Depression status was derived from responses to the NHANES Depression Screener Questionnaire using DSM-IV criteria.

## Repository Structure

The analysis scripts are organized in the order in which they should be run:

1. `01_prepare_depression_alpha.R`
   - Constructs the DSM-IV depression variable.
   - Prepares averaged alpha-diversity measures.

2. `02_build_analysis_dataset.R`
   - Combines microbiome, depression, demographic, behavioral, and oral health data.
   - Creates the primary analysis dataset.

3. `03_descriptive_statistics.R`
   - Generates descriptive statistics for categorical and continuous variables.

4. `04_logistic_regression.R`
   - Performs multivariable logistic regression analyses for each alpha-diversity metric.

5. `05_alpha_diversity_analysis.R`
   - Compares alpha diversity by depression status using one-way ANOVA.
   - Calculates Cohen's d effect sizes.
   - Generates alpha-diversity boxplots.

6. `06_beta_diversity_analysis.R`
   - Performs PERMANOVA for each beta-diversity distance metric.
   - Evaluates within-group dispersion using `betadisper`.
   - Generates beta-dispersion boxplots.

## Data

NHANES demographic, questionnaire, examination, and depression data are publicly available from the National Center for Health Statistics.

Raw data files are not included in this repository. The scripts assume that the required NHANES and oral microbiome files have been downloaded and placed in the corresponding directories under `data/`.

## Statistical Analysis

Unadjusted differences in alpha diversity between participants with and without depression were evaluated using one-way analysis of variance (ANOVA), with Cohen's d used to quantify effect sizes.

Associations between alpha diversity and depression were further evaluated using multivariable logistic regression adjusting for demographic, socioeconomic, behavioral, and oral health covariates.

Differences in overall microbial community composition were evaluated using PERMANOVA with 999 permutations. Differences in within-group microbial dispersion were evaluated using `betadisper` with permutation testing (999 permutations).

Missing continuous covariates were imputed using the mean, and missing categorical covariates were imputed using the mode.

## Runtime

Beta-diversity analyses are computationally intensive. On the computer used for the original analysis, each beta-diversity metric required approximately 1 hour to complete. Runtime will vary depending on available hardware.

## Software

Analyses were conducted in R using packages including:

- `tidyverse`
- `haven`
- `vegan`
- `broom`
- `ggplot2`

## Reproducibility

Random seeds are set before permutation-based analyses to improve reproducibility. Results may vary if the number of permutations, random seed, software versions, or input data are changed.
