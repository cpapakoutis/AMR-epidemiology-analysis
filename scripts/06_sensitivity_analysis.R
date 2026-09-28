############################################################
# Project:
# Temporal and Geographic Patterns in Antimicrobial
# Resistance Reported to WHO GLASS, 2020–2023
#
# Script:
# 06_sensitivity_analysis.R
#
# Purpose:
#   Test whether the main temporal findings are robust to
#   alternative reasonable analytical assumptions.
#
# Sensitivity analyses:
#
#   1. Beta-binomial mixed models
#   2. Stable-country analysis
#   3. Categorical-year analysis
#   4. Equal-country-weighted descriptive analysis
#
# Primary reference analysis:
#   Script 04
#
# Additional robustness model:
#   Script 05
#
# Important:
#   These analyses remain observational analyses of
#   surveillance data.
#
#   They do NOT establish causal effects or
#   population-representative global AMR prevalence.
############################################################


# ----------------------------------------------------------
# 1. Load packages
# ----------------------------------------------------------

library(tidyverse)
library(here)
library(glmmTMB)
library(sandwich)
library(lmtest)


# ----------------------------------------------------------
# 2. Create output directory
# ----------------------------------------------------------

dir.create(
  here("outputs"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ----------------------------------------------------------
# 3. Load cleaned data
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


if (
  length(
    missing_variables
  ) > 0
) {
  
  stop(
    "Required variables are missing: ",
    paste(
      missing_variables,
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------
# 5. Prepare analysis data
# ----------------------------------------------------------

analysis_data <- amr_data %>%
  mutate(
    
    NonResistant =
      InterpretableAST -
      Resistant,
    
    YearCentered =
      Year - 2020,
    
    YearFactor =
      factor(
        Year,
        levels = c(
          2020,
          2021,
          2022,
          2023
        )
      ),
    
    Iso3 =
      factor(
        Iso3
      ),
    
    PathogenName =
      factor(
        PathogenName
      ),
    
    AntibioticName =
      factor(
        AntibioticName
      )
  )


############################################################
# PART A — VALIDATION
############################################################


# ----------------------------------------------------------
# 6. Validate input data
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "SENSITIVITY ANALYSES\n",
  "========================================\n"
)


cat(
  "\nRows:",
  nrow(
    analysis_data
  ),
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


invalid_rows <- analysis_data %>%
  filter(
    
    is.na(
      Resistant
    ) |
      
      is.na(
        InterpretableAST
      ) |
      
      Resistant < 0 |
      
      InterpretableAST <= 0 |
      
      Resistant >
      InterpretableAST |
      
      NonResistant < 0
  )


cat(
  "Invalid count rows:",
  nrow(
    invalid_rows
  ),
  "\n"
)


if (
  nrow(
    invalid_rows
  ) > 0
) {
  
  stop(
    "Invalid resistant/non-resistant counts detected."
  )
}


############################################################
# PART B — REPRODUCE PRIMARY 37-COMBINATION SET
############################################################


# ----------------------------------------------------------
# 7. Calculate year-specific coverage
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
# 8. Summarise combination coverage
# ----------------------------------------------------------

combination_coverage <- combination_year_coverage %>%
  group_by(
    PathogenName,
    AntibioticName
  ) %>%
  summarise(
    
    NumberOfYears =
      n_distinct(
        Year
      ),
    
    MinCountriesPerYear =
      min(
        Countries
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
# 9. Apply same primary eligibility rule as Scripts 04/05
# ----------------------------------------------------------

PRIMARY_MIN_COUNTRIES_PER_YEAR <- 20


primary_combinations <- combination_coverage %>%
  filter(
    
    NumberOfYears == 4,
    
    Has2020,
    
    Has2023,
    
    MinCountriesPerYear >=
      PRIMARY_MIN_COUNTRIES_PER_YEAR
  )


cat(
  "\nPrimary pathogen-antibiotic combinations:",
  nrow(
    primary_combinations
  ),
  "\n"
)


if (
  nrow(
    primary_combinations
  ) != 37
) {
  
  stop(
    "Expected 37 primary combinations but found ",
    nrow(
      primary_combinations
    ),
    ". Check upstream data or eligibility definitions."
  )
}


# ----------------------------------------------------------
# 10. Restrict to primary analysis combinations
# ----------------------------------------------------------

primary_analysis_data <- analysis_data %>%
  semi_join(
    primary_combinations,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  droplevels()


############################################################
# PART C — HELPER FUNCTIONS
############################################################


# ----------------------------------------------------------
# 11. Cluster-robust continuous-year extractor
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
  
  
  robust_vcov <- vcovCL(
    model,
    cluster =
      model_data$Iso3,
    type =
      "HC1",
    cadjust =
      TRUE
  )
  
  
  robust_test <- coeftest(
    model,
    vcov. =
      robust_vcov,
    df =
      number_of_clusters - 1
  )
  
  
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
  
  
  statistic <-
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
  
  
  critical_value <-
    qt(
      0.975,
      df =
        number_of_clusters - 1
    )
  
  
  ci_log_low <-
    estimate -
    critical_value *
    robust_se
  
  
  ci_log_high <-
    estimate +
    critical_value *
    robust_se
  
  
  tibble(
    
    Model =
      model_label,
    
    Countries =
      number_of_clusters,
    
    EstimateLogOddsPerYear =
      estimate,
    
    RobustSE =
      robust_se,
    
    TestStatistic =
      statistic,
    
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
        ci_log_low * 3
      ),
    
    OR_2023_vs_2020_CI_high =
      exp(
        ci_log_high * 3
      )
  )
}


# ----------------------------------------------------------
# 12. Cluster-robust categorical-year extractor
# ----------------------------------------------------------
#
# YearFactor has 2020 as its reference category.
#
# Therefore YearFactor2023 directly compares:
#
#   2023 versus 2020
#
# without assuming a linear trend through 2021 and 2022.
# ----------------------------------------------------------

extract_categorical_2023 <- function(
    model,
    model_data,
    model_label
) {
  
  number_of_clusters <-
    n_distinct(
      model_data$Iso3
    )
  
  
  robust_vcov <- vcovCL(
    model,
    cluster =
      model_data$Iso3,
    type =
      "HC1",
    cadjust =
      TRUE
  )
  
  
  robust_test <- coeftest(
    model,
    vcov. =
      robust_vcov,
    df =
      number_of_clusters - 1
  )
  
  
  target_term <-
    "YearFactor2023"
  
  
  if (
    !target_term %in%
    rownames(
      robust_test
    )
  ) {
    
    stop(
      "2023 categorical-year term missing for: ",
      model_label
    )
  }
  
  
  estimate <-
    unname(
      robust_test[
        target_term,
        1
      ]
    )
  
  
  robust_se <-
    unname(
      robust_test[
        target_term,
        2
      ]
    )
  
  
  statistic <-
    unname(
      robust_test[
        target_term,
        3
      ]
    )
  
  
  p_value <-
    unname(
      robust_test[
        target_term,
        4
      ]
    )
  
  
  critical_value <-
    qt(
      0.975,
      df =
        number_of_clusters - 1
    )
  
  
  ci_log_low <-
    estimate -
    critical_value *
    robust_se
  
  
  ci_log_high <-
    estimate +
    critical_value *
    robust_se
  
  
  tibble(
    
    Model =
      model_label,
    
    Countries =
      number_of_clusters,
    
    EstimateLogOdds2023vs2020 =
      estimate,
    
    RobustSE =
      robust_se,
    
    TestStatistic =
      statistic,
    
    P_value =
      p_value,
    
    OR_2023_vs_2020 =
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
      )
  )
}


# ----------------------------------------------------------
# 13. Extract beta-binomial temporal effect
# ----------------------------------------------------------

extract_betabinomial_year <- function(
    model,
    model_data,
    model_label
) {
  
  coefficient_table <-
    summary(
      model
    )$coefficients$cond
  
  
  if (
    !"YearCentered" %in%
    rownames(
      coefficient_table
    )
  ) {
    
    stop(
      "YearCentered coefficient missing for: ",
      model_label
    )
  }
  
  
  estimate <-
    coefficient_table[
      "YearCentered",
      "Estimate"
    ]
  
  
  standard_error <-
    coefficient_table[
      "YearCentered",
      "Std. Error"
    ]
  
  
  z_value <-
    coefficient_table[
      "YearCentered",
      "z value"
    ]
  
  
  p_value <-
    coefficient_table[
      "YearCentered",
      "Pr(>|z|)"
    ]
  
  
  ci_log_low <-
    estimate -
    1.96 *
    standard_error
  
  
  ci_log_high <-
    estimate +
    1.96 *
    standard_error
  
  
  convergence_code <-
    model$fit$convergence
  
  
  positive_definite_hessian <-
    isTRUE(
      model$sdr$pdHess
    )
  
  
  convergence_message <-
    model$fit$message
  
  
  random_effect_sd <-
    as.numeric(
      attr(
        VarCorr(
          model
        )$cond$Iso3,
        "stddev"
      )[1]
    )
  
  
  beta_binomial_phi <-
    sigma(
      model
    )
  
  
  tibble(
    
    Model =
      model_label,
    
    Countries =
      n_distinct(
        model_data$Iso3
      ),
    
    Observations =
      nrow(
        model_data
      ),
    
    EstimateLogOddsPerYear =
      estimate,
    
    StandardError =
      standard_error,
    
    Z =
      z_value,
    
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
        ci_log_low * 3
      ),
    
    OR_2023_vs_2020_CI_high =
      exp(
        ci_log_high * 3
      ),
    
    CountryRandomInterceptSD =
      random_effect_sd,
    
    BetaBinomialPhi =
      beta_binomial_phi,
    
    ConvergenceCode =
      convergence_code,
    
    PositiveDefiniteHessian =
      positive_definite_hessian,
    
    ConvergenceMessage =
      convergence_message
  )
}


############################################################
# PART D — SENSITIVITY 1:
# BETA-BINOMIAL MIXED MODELS
############################################################


# ----------------------------------------------------------
# 14. Why beta-binomial?
# ----------------------------------------------------------
#
# Script 05 showed substantial residual overdispersion in
# many ordinary binomial GLMMs.
#
# A beta-binomial distribution allows greater variability
# between grouped observations than a simple binomial model.
#
# This analysis therefore tests whether temporal signals
# remain when extra-binomial variability is modelled
# explicitly.
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "SENSITIVITY 1: BETA-BINOMIAL MODELS\n",
  "========================================\n"
)


betabinomial_results_list <-
  list()


betabinomial_failures <-
  list()


for (
  i in seq_len(
    nrow(
      primary_combinations
    )
  )
) {
  
  
  pathogen <-
    as.character(
      primary_combinations$PathogenName[i]
    )
  
  
  antibiotic <-
    as.character(
      primary_combinations$AntibioticName[i]
    )
  
  
  combination_data <- primary_analysis_data %>%
    filter(
      
      PathogenName ==
        pathogen,
      
      AntibioticName ==
        antibiotic
    ) %>%
    droplevels()
  
  
  model_label <- paste(
    pathogen,
    antibiotic,
    sep = " | "
  )
  
  
  cat(
    "\n",
    i,
    "/",
    nrow(
      primary_combinations
    ),
    " — ",
    model_label,
    "\n",
    sep = ""
  )
  
  
  beta_model <- tryCatch(
    
    glmmTMB(
      
      cbind(
        Resistant,
        NonResistant
      ) ~
        
        YearCentered +
        
        (
          1 |
            Iso3
        ),
      
      data =
        combination_data,
      
      family =
        betabinomial(
          link = "logit"
        )
    ),
    
    error = function(e) {
      
      betabinomial_failures[[model_label]] <<-
        e$message
      
      NULL
    }
  )
  
  
  if (
    is.null(
      beta_model
    )
  ) {
    
    cat(
      "MODEL FAILED\n"
    )
    
    next
  }
  
  
  result <-
    extract_betabinomial_year(
      
      model =
        beta_model,
      
      model_data =
        combination_data,
      
      model_label =
        model_label
    ) %>%
    mutate(
      
      PathogenName =
        pathogen,
      
      AntibioticName =
        antibiotic,
      
      .before =
        Model
    )
  
  
  betabinomial_results_list[[i]] <-
    result
  
  
  cat(
    "OR/year:",
    round(
      result$OddsRatioPerYear,
      3
    ),
    "| p:",
    signif(
      result$P_value,
      3
    ),
    "| Hessian OK:",
    result$PositiveDefiniteHessian,
    "| convergence code:",
    result$ConvergenceCode,
    "\n"
  )
}


# ----------------------------------------------------------
# 15. Inspect beta-binomial failures
# ----------------------------------------------------------

if (
  length(
    betabinomial_failures
  ) > 0
) {
  
  beta_failure_table <- tibble(
    
    Model =
      names(
        betabinomial_failures
      ),
    
    Error =
      unlist(
        betabinomial_failures
      )
  )
  
  
  write_csv(
    beta_failure_table,
    here(
      "outputs",
      "betabinomial_model_failures.csv"
    )
  )
  
  
  print(
    beta_failure_table
  )
}


# ----------------------------------------------------------
# 16. Combine beta-binomial results
# ----------------------------------------------------------

betabinomial_results <-
  bind_rows(
    betabinomial_results_list
  ) %>%
  mutate(
    
    FDR_P_value =
      p.adjust(
        P_value,
        method = "BH"
      ),
    
    FDR_below_0_05 =
      FDR_P_value <
      0.05
  ) %>%
  arrange(
    FDR_P_value
  )


write_csv(
  betabinomial_results,
  here(
    "outputs",
    "sensitivity_betabinomial_models.csv"
  )
)


cat(
  "\nBeta-binomial models fitted:",
  nrow(
    betabinomial_results
  ),
  "\n"
)


cat(
  "BH-adjusted p < 0.05:",
  sum(
    betabinomial_results$FDR_below_0_05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Models without positive-definite Hessian:",
  sum(
    !betabinomial_results$PositiveDefiniteHessian,
    na.rm = TRUE
  ),
  "\n"
)


############################################################
# PART E — SENSITIVITY 2:
# STABLE-COUNTRY ANALYSIS
############################################################


# ----------------------------------------------------------
# 17. Identify countries contributing each combination
#     in all four years
# ----------------------------------------------------------
#
# This directly reduces bias from countries entering or
# leaving the surveillance dataset over time.
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "SENSITIVITY 2: STABLE-COUNTRY ANALYSIS\n",
  "========================================\n"
)


stable_country_results_list <-
  list()


stable_country_coverage_list <-
  list()


for (
  i in seq_len(
    nrow(
      primary_combinations
    )
  )
) {
  
  
  pathogen <-
    as.character(
      primary_combinations$PathogenName[i]
    )
  
  
  antibiotic <-
    as.character(
      primary_combinations$AntibioticName[i]
    )
  
  
  combination_data <- primary_analysis_data %>%
    filter(
      
      PathogenName ==
        pathogen,
      
      AntibioticName ==
        antibiotic
    )
  
  
  stable_countries <- combination_data %>%
    distinct(
      Iso3,
      Year
    ) %>%
    count(
      Iso3,
      name = "YearsReported"
    ) %>%
    filter(
      YearsReported == 4
    ) %>%
    pull(
      Iso3
    )
  
  
  stable_data <- combination_data %>%
    filter(
      Iso3 %in%
        stable_countries
    ) %>%
    droplevels()
  
  
  number_stable_countries <-
    n_distinct(
      stable_data$Iso3
    )
  
  
  model_label <- paste(
    pathogen,
    antibiotic,
    sep = " | "
  )
  
  
  stable_country_coverage_list[[i]] <-
    tibble(
      
      PathogenName =
        pathogen,
      
      AntibioticName =
        antibiotic,
      
      StableCountries =
        number_stable_countries,
      
      StableCountriesAtLeast20 =
        number_stable_countries >= 20
    )
  
  
  cat(
    "\n",
    i,
    "/",
    nrow(
      primary_combinations
    ),
    " — ",
    model_label,
    " | stable countries: ",
    number_stable_countries,
    "\n",
    sep = ""
  )
  
  
  if (
    number_stable_countries < 10
  ) {
    
    cat(
      "Skipped: fewer than 10 stable countries.\n"
    )
    
    next
  }
  
  
  stable_model <- glm(
    
    cbind(
      Resistant,
      NonResistant
    ) ~
      
      YearCentered +
      
      factor(
        Iso3
      ),
    
    data =
      stable_data,
    
    family =
      binomial
  )
  
  
  stable_result <-
    extract_cluster_robust_year(
      
      model =
        stable_model,
      
      model_data =
        stable_data,
      
      model_label =
        model_label
    ) %>%
    mutate(
      
      PathogenName =
        pathogen,
      
      AntibioticName =
        antibiotic,
      
      StableCountries =
        number_stable_countries,
      
      StableCountriesAtLeast20 =
        number_stable_countries >= 20,
      
      .before =
        Model
    )
  
  
  stable_country_results_list[[i]] <-
    stable_result
}


# ----------------------------------------------------------
# 18. Combine stable-country results
# ----------------------------------------------------------

stable_country_coverage <-
  bind_rows(
    stable_country_coverage_list
  )


stable_country_results <-
  bind_rows(
    stable_country_results_list
  ) %>%
  mutate(
    
    FDR_P_value =
      p.adjust(
        P_value,
        method = "BH"
      ),
    
    FDR_below_0_05 =
      FDR_P_value <
      0.05
  ) %>%
  arrange(
    FDR_P_value
  )


write_csv(
  stable_country_coverage,
  here(
    "outputs",
    "stable_country_coverage.csv"
  )
)


write_csv(
  stable_country_results,
  here(
    "outputs",
    "sensitivity_stable_country_models.csv"
  )
)


cat(
  "\nStable-country coverage summary:\n"
)


print(
  summary(
    stable_country_coverage$StableCountries
  )
)


cat(
  "\nStable-country models fitted:",
  nrow(
    stable_country_results
  ),
  "\n"
)


cat(
  "Stable-country models with >=20 countries:",
  sum(
    stable_country_results$StableCountriesAtLeast20,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "BH-adjusted p < 0.05:",
  sum(
    stable_country_results$FDR_below_0_05,
    na.rm = TRUE
  ),
  "\n"
)


############################################################
# PART F — SENSITIVITY 3:
# CATEGORICAL YEAR
############################################################


# ----------------------------------------------------------
# 19. Fit categorical-year models
# ----------------------------------------------------------
#
# Primary models treat year as a continuous variable:
#
#   0, 1, 2, 3
#
# This assumes a linear change in log odds.
#
# Here year is categorical instead.
#
# The main coefficient extracted is:
#
#   2023 versus 2020
#
# so 2021 and 2022 are not forced to lie on a straight trend.
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "SENSITIVITY 3: CATEGORICAL YEAR\n",
  "========================================\n"
)


categorical_results_list <-
  list()


for (
  i in seq_len(
    nrow(
      primary_combinations
    )
  )
) {
  
  
  pathogen <-
    as.character(
      primary_combinations$PathogenName[i]
    )
  
  
  antibiotic <-
    as.character(
      primary_combinations$AntibioticName[i]
    )
  
  
  combination_data <- primary_analysis_data %>%
    filter(
      
      PathogenName ==
        pathogen,
      
      AntibioticName ==
        antibiotic
    ) %>%
    droplevels()
  
  
  model_label <- paste(
    pathogen,
    antibiotic,
    sep = " | "
  )
  
  
  categorical_model <- glm(
    
    cbind(
      Resistant,
      NonResistant
    ) ~
      
      YearFactor +
      
      factor(
        Iso3
      ),
    
    data =
      combination_data,
    
    family =
      binomial
  )
  
  
  result <-
    extract_categorical_2023(
      
      model =
        categorical_model,
      
      model_data =
        combination_data,
      
      model_label =
        model_label
    ) %>%
    mutate(
      
      PathogenName =
        pathogen,
      
      AntibioticName =
        antibiotic,
      
      .before =
        Model
    )
  
  
  categorical_results_list[[i]] <-
    result
}


# ----------------------------------------------------------
# 20. Combine categorical-year results
# ----------------------------------------------------------

categorical_year_results <-
  bind_rows(
    categorical_results_list
  ) %>%
  mutate(
    
    FDR_P_value =
      p.adjust(
        P_value,
        method = "BH"
      ),
    
    FDR_below_0_05 =
      FDR_P_value <
      0.05
  ) %>%
  arrange(
    FDR_P_value
  )


write_csv(
  categorical_year_results,
  here(
    "outputs",
    "sensitivity_categorical_year_2023_vs_2020.csv"
  )
)


cat(
  "\nCategorical-year models fitted:",
  nrow(
    categorical_year_results
  ),
  "\n"
)


cat(
  "BH-adjusted 2023-vs-2020 p < 0.05:",
  sum(
    categorical_year_results$FDR_below_0_05,
    na.rm = TRUE
  ),
  "\n"
)


############################################################
# PART G — SENSITIVITY 4:
# EQUAL COUNTRY WEIGHTING
############################################################


# ----------------------------------------------------------
# 21. Calculate one resistance percentage per country,
#     combination and year
# ----------------------------------------------------------
#
# Pooled AST analyses give greater weight to surveillance
# systems contributing more AST results.
#
# Here each country contributes one percentage to the
# annual mean, regardless of AST volume.
#
# This answers a different descriptive question:
#
#   What was the average reported country-level resistance
#   percentage among contributing countries?
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "SENSITIVITY 4: EQUAL COUNTRY WEIGHTING\n",
  "========================================\n"
)


country_level_resistance <- primary_analysis_data %>%
  group_by(
    PathogenName,
    AntibioticName,
    Iso3,
    Year
  ) %>%
  summarise(
    
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
    
    CountryResistancePercentage =
      Resistant /
      InterpretableAST *
      100,
    
    .groups = "drop"
  )


# ----------------------------------------------------------
# 22. Equal-weight annual country mean
# ----------------------------------------------------------

country_weighted_yearly <- country_level_resistance %>%
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
    
    MeanCountryResistancePercentage =
      mean(
        CountryResistancePercentage,
        na.rm = TRUE
      ),
    
    MedianCountryResistancePercentage =
      median(
        CountryResistancePercentage,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


write_csv(
  country_weighted_yearly,
  here(
    "outputs",
    "sensitivity_equal_country_weighted_yearly.csv"
  )
)


# ----------------------------------------------------------
# 23. Equal-country-weighted 2020-to-2023 change
# ----------------------------------------------------------

country_weighted_change <-
  country_weighted_yearly %>%
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
    MeanCountryResistancePercentage,
    MedianCountryResistancePercentage
  ) %>%
  pivot_wider(
    
    names_from =
      Year,
    
    values_from = c(
      MeanCountryResistancePercentage,
      MedianCountryResistancePercentage
    ),
    
    names_glue =
      "{.value}_{Year}"
  ) %>%
  mutate(
    
    EqualWeightMeanChangePP =
      MeanCountryResistancePercentage_2023 -
      MeanCountryResistancePercentage_2020,
    
    EqualWeightMedianChangePP =
      MedianCountryResistancePercentage_2023 -
      MedianCountryResistancePercentage_2020
  )


write_csv(
  country_weighted_change,
  here(
    "outputs",
    "sensitivity_equal_country_weighted_change_2020_2023.csv"
  )
)


# ----------------------------------------------------------
# 24. Equal-country weighting among stable countries only
# ----------------------------------------------------------

stable_country_keys <- primary_analysis_data %>%
  distinct(
    PathogenName,
    AntibioticName,
    Iso3,
    Year
  ) %>%
  group_by(
    PathogenName,
    AntibioticName,
    Iso3
  ) %>%
  summarise(
    
    YearsReported =
      n_distinct(
        Year
      ),
    
    .groups = "drop"
  ) %>%
  filter(
    YearsReported == 4
  ) %>%
  select(
    PathogenName,
    AntibioticName,
    Iso3
  )


stable_country_level_resistance <-
  country_level_resistance %>%
  semi_join(
    stable_country_keys,
    by = c(
      "PathogenName",
      "AntibioticName",
      "Iso3"
    )
  )


stable_country_weighted_yearly <-
  stable_country_level_resistance %>%
  group_by(
    PathogenName,
    AntibioticName,
    Year
  ) %>%
  summarise(
    
    StableCountries =
      n_distinct(
        Iso3
      ),
    
    MeanStableCountryResistancePercentage =
      mean(
        CountryResistancePercentage,
        na.rm = TRUE
      ),
    
    MedianStableCountryResistancePercentage =
      median(
        CountryResistancePercentage,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


stable_country_weighted_change <-
  stable_country_weighted_yearly %>%
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
    StableCountries,
    MeanStableCountryResistancePercentage,
    MedianStableCountryResistancePercentage
  ) %>%
  pivot_wider(
    
    names_from =
      Year,
    
    values_from = c(
      StableCountries,
      MeanStableCountryResistancePercentage,
      MedianStableCountryResistancePercentage
    ),
    
    names_glue =
      "{.value}_{Year}"
  ) %>%
  mutate(
    
    StableEqualWeightMeanChangePP =
      MeanStableCountryResistancePercentage_2023 -
      MeanStableCountryResistancePercentage_2020,
    
    StableEqualWeightMedianChangePP =
      MedianStableCountryResistancePercentage_2023 -
      MedianStableCountryResistancePercentage_2020
  )


write_csv(
  stable_country_weighted_change,
  here(
    "outputs",
    "sensitivity_stable_equal_country_weighted_change.csv"
  )
)


############################################################
# PART H — BUILD CROSS-MODEL ROBUSTNESS TABLE
############################################################


# ----------------------------------------------------------
# 25. Load Script 04 and Script 05 results
# ----------------------------------------------------------

script04_path <- here(
  "outputs",
  "final_temporal_trend_results_primary.csv"
)


script05_path <- here(
  "outputs",
  "pathogen_antibiotic_primary_temporal_glmm.csv"
)


if (
  !file.exists(
    script04_path
  )
) {
  
  stop(
    "Script 04 result file not found. Run Script 04 first."
  )
}


if (
  !file.exists(
    script05_path
  )
) {
  
  stop(
    "Script 05 result file not found. Run Script 05 first."
  )
}


script04_results <- read_csv(
  script04_path,
  show_col_types = FALSE
)


script05_results <- read_csv(
  script05_path,
  show_col_types = FALSE
)


# ----------------------------------------------------------
# 26. Select and rename results for comparison
# ----------------------------------------------------------

primary_for_comparison <- script04_results %>%
  select(
    
    PathogenName,
    
    AntibioticName,
    
    Primary_OR_per_year =
      OddsRatioPerYear,
    
    Primary_CI_low =
      CI_low,
    
    Primary_CI_high =
      CI_high,
    
    Primary_FDR =
      FDR_P_value,
    
    Primary_AbsoluteChangePP =
      AbsoluteChangePercentagePoints
  )


binomial_glmm_for_comparison <-
  script05_results %>%
  select(
    
    PathogenName,
    
    AntibioticName,
    
    BinomialGLMM_OR_per_year =
      OddsRatioPerYear,
    
    BinomialGLMM_FDR =
      FDR_P_value,
    
    BinomialGLMM_Dispersion =
      PearsonDispersion
  )


beta_for_comparison <-
  betabinomial_results %>%
  select(
    
    PathogenName,
    
    AntibioticName,
    
    BetaBinomial_OR_per_year =
      OddsRatioPerYear,
    
    BetaBinomial_CI_low =
      CI_low,
    
    BetaBinomial_CI_high =
      CI_high,
    
    BetaBinomial_FDR =
      FDR_P_value,
    
    BetaBinomial_HessianOK =
      PositiveDefiniteHessian
  )


stable_for_comparison <-
  stable_country_results %>%
  select(
    
    PathogenName,
    
    AntibioticName,
    
    StableCountries,
    
    StableCountry_OR_per_year =
      OddsRatioPerYear,
    
    StableCountry_CI_low =
      CI_low,
    
    StableCountry_CI_high =
      CI_high,
    
    StableCountry_FDR =
      FDR_P_value
  )


categorical_for_comparison <-
  categorical_year_results %>%
  select(
    
    PathogenName,
    
    AntibioticName,
    
    Categorical_OR_2023_vs_2020 =
      OR_2023_vs_2020,
    
    Categorical_CI_low =
      CI_low,
    
    Categorical_CI_high =
      CI_high,
    
    Categorical_FDR =
      FDR_P_value
  )


country_weight_for_comparison <-
  country_weighted_change %>%
  select(
    
    PathogenName,
    
    AntibioticName,
    
    EqualWeightMeanChangePP,
    
    EqualWeightMedianChangePP
  )


stable_weight_for_comparison <-
  stable_country_weighted_change %>%
  select(
    
    PathogenName,
    
    AntibioticName,
    
    StableEqualWeightMeanChangePP,
    
    StableEqualWeightMedianChangePP
  )


# ----------------------------------------------------------
# 27. Join all sensitivity results
# ----------------------------------------------------------

robustness_summary <-
  primary_for_comparison %>%
  
  left_join(
    binomial_glmm_for_comparison,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  
  left_join(
    beta_for_comparison,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  
  left_join(
    stable_for_comparison,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  
  left_join(
    categorical_for_comparison,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  
  left_join(
    country_weight_for_comparison,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  ) %>%
  
  left_join(
    stable_weight_for_comparison,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  )


# ----------------------------------------------------------
# 28. Add direction and robustness indicators
# ----------------------------------------------------------
#
# These indicators are descriptive aids only.
#
# They are NOT a formal scoring system for deciding whether
# a result is "true".
# ----------------------------------------------------------

robustness_summary <- robustness_summary %>%
  mutate(
    
    PrimaryDirection =
      case_when(
        
        Primary_OR_per_year > 1 ~
          "Increase",
        
        Primary_OR_per_year < 1 ~
          "Decrease",
        
        TRUE ~
          "No change"
      ),
    
    BetaDirection =
      case_when(
        
        BetaBinomial_OR_per_year > 1 ~
          "Increase",
        
        BetaBinomial_OR_per_year < 1 ~
          "Decrease",
        
        TRUE ~
          "No change"
      ),
    
    StableDirection =
      case_when(
        
        StableCountry_OR_per_year > 1 ~
          "Increase",
        
        StableCountry_OR_per_year < 1 ~
          "Decrease",
        
        TRUE ~
          "No change"
      ),
    
    CategoricalDirection =
      case_when(
        
        Categorical_OR_2023_vs_2020 > 1 ~
          "Increase",
        
        Categorical_OR_2023_vs_2020 < 1 ~
          "Decrease",
        
        TRUE ~
          "No change"
      ),
    
    EqualWeightDirection =
      case_when(
        
        EqualWeightMeanChangePP > 0 ~
          "Increase",
        
        EqualWeightMeanChangePP < 0 ~
          "Decrease",
        
        TRUE ~
          "No change"
      ),
    
    DirectionAgreementAcrossSensitivities =
      PrimaryDirection ==
      BetaDirection &
      
      PrimaryDirection ==
      StableDirection &
      
      PrimaryDirection ==
      CategoricalDirection &
      
      PrimaryDirection ==
      EqualWeightDirection,
    
    PrimaryFDRBelow05 =
      Primary_FDR <
      0.05,
    
    BetaFDRBelow05 =
      BetaBinomial_FDR <
      0.05,
    
    StableFDRBelow05 =
      StableCountry_FDR <
      0.05,
    
    CategoricalFDRBelow05 =
      Categorical_FDR <
      0.05
  ) %>%
  arrange(
    Primary_FDR
  )


write_csv(
  robustness_summary,
  here(
    "outputs",
    "final_sensitivity_robustness_summary.csv"
  )
)


############################################################
# PART I — INSPECT KEY SCRIPT 04 SIGNALS
############################################################


# ----------------------------------------------------------
# 29. Focus on combinations significant in primary Script 04
# ----------------------------------------------------------

primary_signals <- robustness_summary %>%
  filter(
    PrimaryFDRBelow05
  )


cat(
  "\n========================================\n",
  "PRIMARY SCRIPT 04 SIGNALS ACROSS SENSITIVITIES\n",
  "========================================\n"
)


print(
  primary_signals %>%
    select(
      
      PathogenName,
      
      AntibioticName,
      
      Primary_OR_per_year,
      
      Primary_CI_low,
      
      Primary_CI_high,
      
      Primary_FDR,
      
      BetaBinomial_OR_per_year,
      
      BetaBinomial_CI_low,
      
      BetaBinomial_CI_high,
      
      BetaBinomial_FDR,
      
      StableCountries,
      
      StableCountry_OR_per_year,
      
      StableCountry_CI_low,
      
      StableCountry_CI_high,
      
      StableCountry_FDR,
      
      Categorical_OR_2023_vs_2020,
      
      Categorical_CI_low,
      
      Categorical_CI_high,
      
      Categorical_FDR,
      
      EqualWeightMeanChangePP,
      
      StableEqualWeightMeanChangePP,
      
      DirectionAgreementAcrossSensitivities
    ),
  
  n = Inf,
  
  width = Inf
)


write_csv(
  primary_signals,
  here(
    "outputs",
    "primary_signals_across_sensitivities.csv"
  )
)


############################################################
# PART J — OVERALL SENSITIVITY SUMMARY
############################################################


# ----------------------------------------------------------
# 30. Print key counts
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "SENSITIVITY ANALYSIS SUMMARY\n",
  "========================================\n"
)


cat(
  "\nPrimary combinations:",
  nrow(
    robustness_summary
  ),
  "\n"
)


cat(
  "Script 04 FDR < 0.05:",
  sum(
    robustness_summary$PrimaryFDRBelow05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Beta-binomial FDR < 0.05:",
  sum(
    robustness_summary$BetaFDRBelow05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Stable-country FDR < 0.05:",
  sum(
    robustness_summary$StableFDRBelow05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Categorical-year FDR < 0.05:",
  sum(
    robustness_summary$CategoricalFDRBelow05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Direction agreement across all sensitivity approaches:",
  sum(
    robustness_summary$DirectionAgreementAcrossSensitivities,
    na.rm = TRUE
  ),
  "/",
  nrow(
    robustness_summary
  ),
  "\n"
)


############################################################
# PART K — FINAL MESSAGE
############################################################


cat(
  "\n========================================\n",
  "SENSITIVITY ANALYSES COMPLETE\n",
  "========================================\n"
)


cat(
  "\nMain outputs:\n",
  "- sensitivity_betabinomial_models.csv\n",
  "- stable_country_coverage.csv\n",
  "- sensitivity_stable_country_models.csv\n",
  "- sensitivity_categorical_year_2023_vs_2020.csv\n",
  "- sensitivity_equal_country_weighted_yearly.csv\n",
  "- sensitivity_equal_country_weighted_change_2020_2023.csv\n",
  "- sensitivity_stable_equal_country_weighted_change.csv\n",
  "- final_sensitivity_robustness_summary.csv\n",
  "- primary_signals_across_sensitivities.csv\n"
)


cat(
  "\nIMPORTANT:\n",
  paste(
    "Sensitivity analyses should be interpreted together.",
    "A finding is more convincing when its direction and",
    "approximate magnitude remain similar across reasonable",
    "analytical choices. Statistical significance alone is",
    "not used as the criterion for scientific interpretation."
  ),
  "\n"
)