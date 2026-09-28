############################################################
# Project:
# Temporal and Geographic Patterns in Antimicrobial
# Resistance Reported to WHO GLASS, 2020–2023
#
# Script:
# 04_statistical_analysis.R
#
# Purpose:
#   Conduct the main inferential analysis of temporal AMR
#   patterns using country-aware binomial regression with
#   cluster-robust inference.
#
# Important interpretation:
#   The data contain counts of antimicrobial susceptibility
#   test (AST) results.
#
#   Models therefore describe temporal associations in
#   reported AST-level resistance among contributing
#   surveillance systems.
#
#   They do NOT estimate population-representative global
#   AMR prevalence.
#
# Input:
#   data/processed/AMR_GLASS_clean.csv
#
# Outputs:
#   Statistical model results and coverage diagnostics
#   saved in outputs/
############################################################


# ----------------------------------------------------------
# 1. Load packages
# ----------------------------------------------------------

library(tidyverse)
library(here)
library(sandwich)
library(lmtest)


# ----------------------------------------------------------
# 2. Load cleaned dataset
# ----------------------------------------------------------

amr_data <- read_csv(
  here(
    "data",
    "processed",
    "AMR_GLASS_clean.csv"
  ),
  show_col_types = FALSE
)


# ----------------------------------------------------------
# 3. Create output directory
# ----------------------------------------------------------

dir.create(
  here("outputs"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ----------------------------------------------------------
# 4. Confirm required variables
# ----------------------------------------------------------

required_variables <- c(
  "PathogenName",
  "AntibioticName",
  "Iso3",
  "CountryTerritoryArea",
  "WHORegionName",
  "InterpretableAST",
  "Resistant",
  "Year"
)


missing_variables <- setdiff(
  required_variables,
  names(amr_data)
)


if (length(missing_variables) > 0) {
  
  stop(
    "Required variables are missing: ",
    paste(
      missing_variables,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------
# 5. Prepare analysis dataset
# ----------------------------------------------------------
#
# NonResistant is required for grouped binomial regression.
#
# Year is centred at 2020:
#
#   2020 -> 0
#   2021 -> 1
#   2022 -> 2
#   2023 -> 3
#
# This does NOT change the temporal slope.
#
# It makes the intercept interpretable as the modelled
# baseline at 2020 and improves numerical stability.
# ----------------------------------------------------------

analysis_data <- amr_data %>%
  mutate(
    
    NonResistant =
      InterpretableAST - Resistant,
    
    YearCentered =
      Year - 2020,
    
    PathogenAntibiotic =
      interaction(
        PathogenName,
        AntibioticName,
        drop = TRUE,
        sep = " | "
      )
  )


############################################################
# PART A — DATA VALIDATION
############################################################


# ----------------------------------------------------------
# 6. Basic validation
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "STATISTICAL ANALYSIS\n",
  "========================================\n"
)


cat(
  "\nRows:",
  nrow(analysis_data),
  "\n"
)


cat(
  "Countries/territories:",
  n_distinct(
    analysis_data$Iso3
  ),
  "\n"
)


cat(
  "Years:",
  paste(
    sort(
      unique(
        analysis_data$Year
      )
    ),
    collapse = ", "
  ),
  "\n"
)


cat(
  "Pathogens:",
  n_distinct(
    analysis_data$PathogenName
  ),
  "\n"
)


# ----------------------------------------------------------
# 7. Validate study years
# ----------------------------------------------------------

expected_years <- c(
  2020,
  2021,
  2022,
  2023
)


observed_years <- sort(
  unique(
    analysis_data$Year
  )
)


if (
  length(observed_years) != length(expected_years) ||
  any(observed_years != expected_years)
) {
  
  stop(
    "Unexpected study years detected. Found: ",
    paste(
      observed_years,
      collapse = ", "
    )
  )
}

# ----------------------------------------------------------
# 8. Check impossible values
# ----------------------------------------------------------

invalid_count_rows <- analysis_data %>%
  filter(
    is.na(Resistant) |
      is.na(InterpretableAST) |
      Resistant < 0 |
      InterpretableAST <= 0 |
      Resistant > InterpretableAST |
      NonResistant < 0
  )


cat(
  "\nInvalid count rows:",
  nrow(invalid_count_rows),
  "\n"
)


if (nrow(invalid_count_rows) > 0) {
  
  stop(
    "Invalid resistant/non-resistant counts detected."
  )
}


############################################################
# PART B — PATHOGEN–ANTIBIOTIC COVERAGE
############################################################


# ----------------------------------------------------------
# 9. Coverage by combination and year
# ----------------------------------------------------------

combination_year_coverage <- analysis_data %>%
  group_by(
    PathogenName,
    AntibioticName,
    Year
  ) %>%
  summarise(
    
    Countries =
      n_distinct(
        Iso3
      ),
    
    InterpretableAST =
      sum(
        InterpretableAST,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ----------------------------------------------------------
# 10. Overall coverage for each combination
# ----------------------------------------------------------

combination_overall_coverage <- analysis_data %>%
  group_by(
    PathogenName,
    AntibioticName
  ) %>%
  summarise(
    
    NumberOfYears =
      n_distinct(
        Year
      ),
    
    OverallCountries =
      n_distinct(
        Iso3
      ),
    
    TotalInterpretableAST =
      sum(
        InterpretableAST,
        na.rm = TRUE
      ),
    
    Has2020 =
      any(
        Year == 2020
      ),
    
    Has2023 =
      any(
        Year == 2023
      ),
    
    .groups = "drop"
  )


# ----------------------------------------------------------
# 11. Year-specific country coverage
# ----------------------------------------------------------

combination_year_summary <- combination_year_coverage %>%
  group_by(
    PathogenName,
    AntibioticName
  ) %>%
  summarise(
    
    MinCountriesPerYear =
      min(
        Countries
      ),
    
    MedianCountriesPerYear =
      median(
        Countries
      ),
    
    MaxCountriesPerYear =
      max(
        Countries
      ),
    
    MinASTPerYear =
      min(
        InterpretableAST
      ),
    
    .groups = "drop"
  )


# ----------------------------------------------------------
# 12. Define coverage criteria
# ----------------------------------------------------------
#
# Primary analysis:
#
#   - observations in all four years
#   - 2020 and 2023 represented
#   - at least 20 contributing countries in EVERY year
#
# Sensitivity definition:
#
#   - same temporal requirements
#   - at least 10 contributing countries in every year
#
# The sensitivity threshold will be examined later.
#
# These thresholds concern longitudinal surveillance
# coverage, not AST sample size alone.
# ----------------------------------------------------------

PRIMARY_MIN_COUNTRIES_PER_YEAR <- 20

SENSITIVITY_MIN_COUNTRIES_PER_YEAR <- 10


combination_coverage <- combination_overall_coverage %>%
  left_join(
    combination_year_summary,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  mutate(
    
    PrimaryEligible =
      NumberOfYears == 4 &
      Has2020 &
      Has2023 &
      MinCountriesPerYear >=
      PRIMARY_MIN_COUNTRIES_PER_YEAR,
    
    SensitivityEligible10 =
      NumberOfYears == 4 &
      Has2020 &
      Has2023 &
      MinCountriesPerYear >=
      SENSITIVITY_MIN_COUNTRIES_PER_YEAR,
    
    PrimaryExclusionReason =
      case_when(
        
        NumberOfYears != 4 ~
          "Not represented in all four study years",
        
        !Has2020 | !Has2023 ~
          "Missing 2020 or 2023 endpoint",
        
        MinCountriesPerYear <
          PRIMARY_MIN_COUNTRIES_PER_YEAR ~
          "Fewer than 20 countries in at least one year",
        
        TRUE ~
          "Eligible"
      )
  ) %>%
  arrange(
    desc(
      PrimaryEligible
    ),
    MinCountriesPerYear
  )


write_csv(
  combination_coverage,
  here(
    "outputs",
    "pathogen_antibiotic_coverage_for_inference.csv"
  )
)


# ----------------------------------------------------------
# 13. Primary combinations
# ----------------------------------------------------------

primary_combinations <- combination_coverage %>%
  filter(
    PrimaryEligible
  )


sensitivity_combinations_10 <- combination_coverage %>%
  filter(
    SensitivityEligible10
  )


cat(
  "\n========================================\n",
  "PATHOGEN–ANTIBIOTIC ELIGIBILITY\n",
  "========================================\n"
)


cat(
  "\nTotal combinations:",
  nrow(combination_coverage),
  "\n"
)


cat(
  "Primary combinations (>=20 countries/year):",
  nrow(primary_combinations),
  "\n"
)


cat(
  "Sensitivity combinations (>=10 countries/year):",
  nrow(sensitivity_combinations_10),
  "\n"
)


cat(
  "\nPrimary exclusions:\n"
)


print(
  combination_coverage %>%
    filter(
      !PrimaryEligible
    ) %>%
    select(
      PathogenName,
      AntibioticName,
      OverallCountries,
      MinCountriesPerYear,
      TotalInterpretableAST,
      PrimaryExclusionReason
    )
)


############################################################
# PART C — CREATE PRIMARY INFERENTIAL DATASET
############################################################


# ----------------------------------------------------------
# 14. Restrict inferential analysis to combinations with
#     adequate longitudinal country coverage
# ----------------------------------------------------------

primary_analysis_data <- analysis_data %>%
  semi_join(
    primary_combinations %>%
      select(
        PathogenName,
        AntibioticName
      ),
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  mutate(
    PathogenAntibiotic =
      droplevels(
        PathogenAntibiotic
      )
  )


cat(
  "\nPrimary inferential dataset rows:",
  nrow(primary_analysis_data),
  "\n"
)


cat(
  "Countries represented:",
  n_distinct(
    primary_analysis_data$Iso3
  ),
  "\n"
)


cat(
  "Pathogen-antibiotic combinations:",
  n_distinct(
    primary_analysis_data$PathogenAntibiotic
  ),
  "\n"
)


############################################################
# PART D — CLUSTER-ROBUST INFERENCE FUNCTION
############################################################


# ----------------------------------------------------------
# 15. Function to extract YearCentered effect
# ----------------------------------------------------------
#
# This helper:
#
#   1. calculates country-clustered covariance
#   2. applies an HC1 finite-sample adjustment
#   3. uses country-cluster degrees of freedom
#   4. extracts the temporal coefficient
#   5. converts log-odds to odds ratios
#   6. calculates a Pearson dispersion diagnostic
#
# The OR per year assumes a linear trend in log odds.
#
# OR_2023_vs_2020 represents the implied three-year
# contrast under that linear trend assumption.
# ----------------------------------------------------------

extract_cluster_robust_year <- function(
    model,
    model_data,
    model_label
) {
  
  number_of_clusters <-
    n_distinct(
      model_data$Iso3
    )
  
  
  if (number_of_clusters < 2) {
    
    stop(
      "Fewer than two country clusters in model: ",
      model_label
    )
  }
  
  
  robust_vcov <- vcovCL(
    model,
    cluster = model_data$Iso3,
    type = "HC1",
    cadjust = TRUE
  )
  
  
  robust_test <- coeftest(
    model,
    vcov. = robust_vcov,
    df = number_of_clusters - 1
  )
  
  
  if (!"YearCentered" %in%
      rownames(
        robust_test
      )) {
    
    stop(
      "YearCentered coefficient missing in model: ",
      model_label
    )
  }
  
  
  estimate <-
    unname(
      robust_test[
        "YearCentered",
        1
      ]
    )
  
  
  robust_se <-
    unname(
      robust_test[
        "YearCentered",
        2
      ]
    )
  
  
  test_statistic <-
    unname(
      robust_test[
        "YearCentered",
        3
      ]
    )
  
  
  p_value <-
    unname(
      robust_test[
        "YearCentered",
        4
      ]
    )
  
  
  degrees_freedom <-
    number_of_clusters - 1
  
  
  critical_value <-
    qt(
      0.975,
      df = degrees_freedom
    )
  
  
  ci_log_low <-
    estimate -
    critical_value *
    robust_se
  
  
  ci_log_high <-
    estimate +
    critical_value *
    robust_se
  
  
  pearson_dispersion <-
    sum(
      residuals(
        model,
        type = "pearson"
      )^2
    ) /
    df.residual(
      model
    )
  
  
  tibble(
    
    Model =
      model_label,
    
    Countries =
      number_of_clusters,
    
    EstimateLogOddsPerYear =
      estimate,
    
    RobustSE =
      robust_se,
    
    DegreesFreedom =
      degrees_freedom,
    
    TestStatistic =
      test_statistic,
    
    P_value =
      p_value,
    
    OddsRatioPerYear =
      exp(
        estimate
      ),
    
    CI_low =
      exp(
        ci_log_low
      ),
    
    CI_high =
      exp(
        ci_log_high
      ),
    
    OR_2023_vs_2020 =
      exp(
        estimate * 3
      ),
    
    OR_2023_vs_2020_CI_low =
      exp(
        3 * ci_log_low
      ),
    
    OR_2023_vs_2020_CI_high =
      exp(
        3 * ci_log_high
      ),
    
    PearsonDispersion =
      pearson_dispersion,
    
    ModelConverged =
      model$converged
  )
}


############################################################
# PART E — ADJUSTED OVERALL TEMPORAL MODEL
############################################################


# ----------------------------------------------------------
# 16. Overall adjusted model
# ----------------------------------------------------------
#
# This is NOT interpreted as "global AMR prevalence".
#
# The model estimates a common temporal association among
# the adequately represented pathogen-antibiotic
# combinations.
#
# It adjusts for:
#
#   - pathogen-antibiotic combination
#   - time-invariant country differences
#
# Country-clustered SE account for correlation of repeated
# observations from the same country.
#
# The result remains AST-volume weighted.
# ----------------------------------------------------------

overall_adjusted_model <- glm(
  
  cbind(
    Resistant,
    NonResistant
  ) ~
    
    YearCentered +
    
    PathogenAntibiotic +
    
    factor(
      Iso3
    ),
  
  data =
    primary_analysis_data,
  
  family =
    binomial
)


if (!overall_adjusted_model$converged) {
  
  stop(
    "Overall adjusted model did not converge."
  )
}


overall_adjusted_result <-
  extract_cluster_robust_year(
    model =
      overall_adjusted_model,
    model_data =
      primary_analysis_data,
    model_label =
      "Overall adjusted AST-level temporal model"
  )


cat(
  "\n========================================\n",
  "OVERALL ADJUSTED TEMPORAL MODEL\n",
  "========================================\n"
)


print(
  overall_adjusted_result
)


write_csv(
  overall_adjusted_result,
  here(
    "outputs",
    "overall_adjusted_temporal_model_robust.csv"
  )
)


############################################################
# PART F — PATHOGEN-SPECIFIC ADJUSTED MODELS
############################################################


# ----------------------------------------------------------
# 17. Pathogen-specific models
# ----------------------------------------------------------
#
# Within each pathogen, the model adjusts for:
#
#   - antibiotic identity
#   - time-invariant country differences
#
# This is an improvement over pooling antibiotics without
# accounting for which antibiotic generated each AST result.
# ----------------------------------------------------------

pathogen_list <- sort(
  unique(
    primary_analysis_data$PathogenName
  )
)


pathogen_results_list <- list()


for (pathogen in pathogen_list) {
  
  
  pathogen_data <- primary_analysis_data %>%
    filter(
      PathogenName == pathogen
    )
  
  
  pathogen_model <- glm(
    
    cbind(
      Resistant,
      NonResistant
    ) ~
      
      YearCentered +
      
      AntibioticName +
      
      factor(
        Iso3
      ),
    
    data =
      pathogen_data,
    
    family =
      binomial
  )
  
  
  if (!pathogen_model$converged) {
    
    stop(
      "Pathogen model did not converge: ",
      pathogen
    )
  }
  
  
  pathogen_result <-
    extract_cluster_robust_year(
      model =
        pathogen_model,
      model_data =
        pathogen_data,
      model_label =
        pathogen
    ) %>%
    mutate(
      
      PathogenName =
        pathogen,
      
      Antibiotics =
        n_distinct(
          pathogen_data$AntibioticName
        ),
      
      Observations =
        nrow(
          pathogen_data
        ),
      
      .before =
        Model
    )
  
  
  pathogen_results_list[[pathogen]] <-
    pathogen_result
}


pathogen_results <- bind_rows(
  pathogen_results_list
) %>%
  mutate(
    
    FDR_P_value =
      p.adjust(
        P_value,
        method = "BH"
      )
  ) %>%
  arrange(
    FDR_P_value
  )


cat(
  "\n========================================\n",
  "PATHOGEN-SPECIFIC ADJUSTED MODELS\n",
  "========================================\n"
)


print(
  pathogen_results
)


write_csv(
  pathogen_results,
  here(
    "outputs",
    "pathogen_adjusted_temporal_models_robust.csv"
  )
)


############################################################
# PART G — PATHOGEN–ANTIBIOTIC MODELS
############################################################


# ----------------------------------------------------------
# 18. Fit primary pathogen-antibiotic models
# ----------------------------------------------------------
#
# Each model contains only one pathogen-antibiotic
# combination.
#
# The model therefore adjusts for country baseline using
# country fixed effects:
#
#   resistance ~ year + country
#
# Country-clustered robust inference is then used.
#
# These models provide the clearest inferential estimates
# of specific resistance phenotypes in Script 04.
# ----------------------------------------------------------

combination_results_list <- list()


for (
  i in seq_len(
    nrow(
      primary_combinations
    )
  )
) {
  
  
  pathogen <-
    primary_combinations$PathogenName[i]
  
  
  antibiotic <-
    primary_combinations$AntibioticName[i]
  
  
  combination_data <- primary_analysis_data %>%
    filter(
      PathogenName == pathogen,
      AntibioticName == antibiotic
    )
  
  
  combination_model <- glm(
    
    cbind(
      Resistant,
      NonResistant
    ) ~
      
      YearCentered +
      
      factor(
        Iso3
      ),
    
    data =
      combination_data,
    
    family =
      binomial
  )
  
  
  if (!combination_model$converged) {
    
    stop(
      "Combination model did not converge: ",
      pathogen,
      " / ",
      antibiotic
    )
  }
  
  
  result <-
    extract_cluster_robust_year(
      model =
        combination_model,
      model_data =
        combination_data,
      model_label =
        paste(
          pathogen,
          antibiotic,
          sep = " | "
        )
    ) %>%
    mutate(
      
      PathogenName =
        pathogen,
      
      AntibioticName =
        antibiotic,
      
      Observations =
        nrow(
          combination_data
        ),
      
      .before =
        Model
    )
  
  
  combination_results_list[[i]] <-
    result
}


combination_trends <- bind_rows(
  combination_results_list
)


############################################################
# PART H — MULTIPLE TESTING
############################################################


# ----------------------------------------------------------
# 19. Benjamini-Hochberg FDR correction
# ----------------------------------------------------------
#
# A separate statistical test is performed for each
# eligible pathogen-antibiotic combination.
#
# Benjamini-Hochberg adjustment controls the expected false
# discovery rate across this family of tests.
#
# Statistical significance alone will NOT determine
# scientific interpretation.
# ----------------------------------------------------------

combination_trends <- combination_trends %>%
  mutate(
    
    FDR_P_value =
      p.adjust(
        P_value,
        method = "BH"
      ),
    
    FDR_below_0_05 =
      FDR_P_value < 0.05
  ) %>%
  arrange(
    FDR_P_value
  )


write_csv(
  combination_trends,
  here(
    "outputs",
    "pathogen_antibiotic_primary_temporal_models_robust.csv"
  )
)


############################################################
# PART I — DESCRIPTIVE 2020 VS 2023 VALUES
############################################################


# ----------------------------------------------------------
# 20. Annual descriptive resistance within primary
#     pathogen-antibiotic combinations
# ----------------------------------------------------------

primary_keys <- primary_combinations %>%
  select(
    PathogenName,
    AntibioticName
  )


combination_yearly <- analysis_data %>%
  semi_join(
    primary_keys,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  group_by(
    PathogenName,
    AntibioticName,
    Year
  ) %>%
  summarise(
    
    Countries =
      n_distinct(
        Iso3
      ),
    
    Resistant =
      sum(
        Resistant,
        na.rm = TRUE
      ),
    
    InterpretableAST =
      sum(
        InterpretableAST,
        na.rm = TRUE
      ),
    
    ResistancePercentage =
      Resistant /
      InterpretableAST *
      100,
    
    .groups = "drop"
  )


# ----------------------------------------------------------
# 21. Create 2020–2023 descriptive comparison
# ----------------------------------------------------------

combination_change <- combination_yearly %>%
  filter(
    Year %in% c(
      2020,
      2023
    )
  ) %>%
  select(
    PathogenName,
    AntibioticName,
    Year,
    Countries,
    InterpretableAST,
    ResistancePercentage
  ) %>%
  pivot_wider(
    
    names_from =
      Year,
    
    values_from = c(
      Countries,
      InterpretableAST,
      ResistancePercentage
    ),
    
    names_glue =
      "{.value}_{Year}"
  ) %>%
  mutate(
    
    AbsoluteChangePercentagePoints =
      ResistancePercentage_2023 -
      ResistancePercentage_2020
  )


############################################################
# PART J — FINAL PRIMARY RESULTS TABLE
############################################################


# ----------------------------------------------------------
# 22. Add coverage and descriptive values to model results
# ----------------------------------------------------------

final_trend_results <- combination_trends %>%
  left_join(
    combination_change,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  left_join(
    primary_combinations %>%
      select(
        PathogenName,
        AntibioticName,
        OverallCountries,
        MinCountriesPerYear,
        MedianCountriesPerYear,
        MaxCountriesPerYear,
        TotalInterpretableAST
      ),
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  arrange(
    FDR_P_value
  )


write_csv(
  final_trend_results,
  here(
    "outputs",
    "final_temporal_trend_results_primary.csv"
  )
)


# ----------------------------------------------------------
# 23. Convenience table:
#     combinations with BH-adjusted p < 0.05
# ----------------------------------------------------------
#
# This table is for convenient inspection only.
#
# Inclusion here does NOT automatically imply that a trend
# is robust, causal, clinically important, or suitable for
# highlighting in the final interpretation.
# ----------------------------------------------------------

fdr_below_005 <- final_trend_results %>%
  filter(
    FDR_P_value < 0.05
  )


write_csv(
  fdr_below_005,
  here(
    "outputs",
    "fdr_below_0_05_primary.csv"
  )
)


############################################################
# PART K — DESCRIPTIVE EFFECT DIRECTION
############################################################


# ----------------------------------------------------------
# 24. Direction of descriptive 2020–2023 changes
# ----------------------------------------------------------

direction_summary <- final_trend_results %>%
  mutate(
    
    DescriptiveDirection =
      case_when(
        
        AbsoluteChangePercentagePoints > 0 ~
          "Increase",
        
        AbsoluteChangePercentagePoints < 0 ~
          "Decrease",
        
        AbsoluteChangePercentagePoints == 0 ~
          "No change",
        
        TRUE ~
          NA_character_
      )
  ) %>%
  count(
    DescriptiveDirection,
    name = "Combinations"
  )


write_csv(
  direction_summary,
  here(
    "outputs",
    "primary_trend_direction_summary.csv"
  )
)

############################################################
# PART L — PRIMARY RESULTS INSPECTION AND DIAGNOSTICS
############################################################


# ----------------------------------------------------------
# 25. Detailed results for combinations passing BH FDR
# ----------------------------------------------------------
#
# Statistical significance is not sufficient on its own.
#
# For combinations passing FDR < 0.05, inspect:
#
#   - effect size
#   - confidence interval
#   - absolute descriptive change
#   - country coverage
#   - model dispersion
#
# This prevents interpretation based only on p-values.
# ----------------------------------------------------------

fdr_detailed_results <- final_trend_results %>%
  filter(
    FDR_P_value < 0.05
  ) %>%
  select(
    PathogenName,
    AntibioticName,
    Countries,
    MinCountriesPerYear,
    OddsRatioPerYear,
    CI_low,
    CI_high,
    OR_2023_vs_2020,
    OR_2023_vs_2020_CI_low,
    OR_2023_vs_2020_CI_high,
    P_value,
    FDR_P_value,
    ResistancePercentage_2020,
    ResistancePercentage_2023,
    AbsoluteChangePercentagePoints,
    PearsonDispersion
  )


cat(
  "\n========================================\n",
  "FDR-SIGNIFICANT COMBINATIONS\n",
  "========================================\n"
)


print(
  fdr_detailed_results,
  n = Inf,
  width = Inf
)


write_csv(
  fdr_detailed_results,
  here(
    "outputs",
    "fdr_significant_detailed_results.csv"
  )
)


# ----------------------------------------------------------
# 26. Distribution of Pearson dispersion
# ----------------------------------------------------------
#
# A simple grouped-binomial model expects a dispersion
# statistic approximately around 1.
#
# Values substantially above 1 indicate extra-binomial
# variation: the observations vary more than the simple
# binomial variance assumption predicts.
#
# Cluster-robust inference helps protect the standard errors
# against variance/correlation misspecification, but this
# diagnostic is still important because it motivates the
# mixed-effects and sensitivity analyses in Script 05.
# ----------------------------------------------------------

dispersion_summary <- final_trend_results %>%
  summarise(
    
    Minimum =
      min(
        PearsonDispersion,
        na.rm = TRUE
      ),
    
    FirstQuartile =
      quantile(
        PearsonDispersion,
        0.25,
        na.rm = TRUE
      ),
    
    Median =
      median(
        PearsonDispersion,
        na.rm = TRUE
      ),
    
    Mean =
      mean(
        PearsonDispersion,
        na.rm = TRUE
      ),
    
    ThirdQuartile =
      quantile(
        PearsonDispersion,
        0.75,
        na.rm = TRUE
      ),
    
    Maximum =
      max(
        PearsonDispersion,
        na.rm = TRUE
      )
  )


cat(
  "\n========================================\n",
  "PEARSON DISPERSION SUMMARY\n",
  "========================================\n"
)


print(
  dispersion_summary
)


write_csv(
  dispersion_summary,
  here(
    "outputs",
    "pathogen_antibiotic_dispersion_summary.csv"
  )
)


# ----------------------------------------------------------
# 27. Ten strongest pathogen-antibiotic signals
# ----------------------------------------------------------
#
# These are ordered by BH-adjusted p-value for diagnostic
# inspection.
#
# Being in this table does NOT imply that a result is
# statistically significant, clinically important, or
# suitable for highlighting in the final report.
# ----------------------------------------------------------

top_10_signals <- final_trend_results %>%
  arrange(
    FDR_P_value
  ) %>%
  select(
    PathogenName,
    AntibioticName,
    OddsRatioPerYear,
    CI_low,
    CI_high,
    P_value,
    FDR_P_value,
    PearsonDispersion
  ) %>%
  slice_head(
    n = 10
  )


cat(
  "\n========================================\n",
  "TEN STRONGEST STATISTICAL SIGNALS\n",
  "========================================\n"
)


print(
  top_10_signals,
  n = 10,
  width = Inf
)


write_csv(
  top_10_signals,
  here(
    "outputs",
    "top_10_pathogen_antibiotic_signals.csv"
  )
)

############################################################
# PART L — FINAL CONSOLE SUMMARY
############################################################


cat(
  "\n========================================\n",
  "MODEL SUMMARY\n",
  "========================================\n"
)


cat(
  "\nPrimary pathogen-antibiotic combinations:",
  nrow(primary_combinations),
  "\n"
)


cat(
  "Sensitivity combinations at >=10 countries/year:",
  nrow(sensitivity_combinations_10),
  "\n"
)


cat(
  "Combinations with BH-adjusted p < 0.05:",
  nrow(fdr_below_005),
  "\n"
)


cat(
  "\nOverall adjusted temporal model:\n"
)


print(
  overall_adjusted_result %>%
    select(
      Countries,
      OddsRatioPerYear,
      CI_low,
      CI_high,
      P_value,
      PearsonDispersion
    )
)


cat(
  "\nPathogen-specific adjusted temporal models:\n"
)


print(
  pathogen_results %>%
    select(
      PathogenName,
      Countries,
      Antibiotics,
      OddsRatioPerYear,
      CI_low,
      CI_high,
      P_value,
      FDR_P_value,
      PearsonDispersion
    )
)


cat(
  "\nDirection of descriptive 2020–2023 changes:\n"
)


print(
  direction_summary
)


############################################################
# PART M — FINAL MESSAGE
############################################################


cat(
  "\n========================================\n",
  "STATISTICAL ANALYSIS COMPLETE\n",
  "========================================\n"
)


cat(
  "\nResults saved to:\n",
  here("outputs"),
  "\n"
)


cat(
  "\nPrimary inferential files:\n",
  "- pathogen_antibiotic_coverage_for_inference.csv\n",
  "- overall_adjusted_temporal_model_robust.csv\n",
  "- pathogen_adjusted_temporal_models_robust.csv\n",
  "- pathogen_antibiotic_primary_temporal_models_robust.csv\n",
  "- final_temporal_trend_results_primary.csv\n",
  "- fdr_below_0_05_primary.csv\n"
)


cat(
  "\nIMPORTANT:\n",
  paste(
    "These models estimate temporal associations in reported",
    "AST-level resistance among contributing WHO GLASS",
    "surveillance systems. They should not be interpreted as",
    "population-representative global AMR trends or causal effects."
  ),
  "\n"
)