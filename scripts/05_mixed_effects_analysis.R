############################################################
# Project:
# Temporal and Geographic Patterns in Antimicrobial
# Resistance Reported to WHO GLASS, 2020–2023
#
# Script:
# 05_mixed_effects_analysis.R
#
# Purpose:
#   Assess the robustness of temporal AMR associations using
#   binomial mixed-effects models with country-level random
#   intercepts.
#
# Script 04 used:
#
#   - country fixed effects
#   - country-clustered robust standard errors
#
# Script 05 uses:
#
#   - country random intercepts
#
# Comparing these approaches helps determine whether the
# main conclusions are robust to reasonable differences in
# modelling country-level heterogeneity.
#
# Important interpretation:
#   Models describe temporal associations in reported
#   AST-level resistance among contributing WHO GLASS
#   surveillance systems.
#
#   They do NOT estimate population-representative global
#   AMR prevalence or causal effects.
#
# Input:
#   data/processed/AMR_GLASS_clean.csv
#
# Script 04 comparison input:
#   outputs/final_temporal_trend_results_primary.csv
#
# Outputs:
#   Mixed-effects model results and diagnostics in outputs/
############################################################


# ----------------------------------------------------------
# 1. Load packages
# ----------------------------------------------------------

library(tidyverse)
library(here)
library(lme4)


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
# NonResistant is required for grouped binomial modelling.
#
# Year is centred at 2020:
#
#   2020 = 0
#   2021 = 1
#   2022 = 2
#   2023 = 3
#
# This does NOT change the temporal slope.
#
# It improves numerical stability and makes the intercept
# correspond to the modelled baseline in 2020.
# ----------------------------------------------------------

analysis_data <- amr_data %>%
  mutate(
    
    NonResistant =
      InterpretableAST - Resistant,
    
    YearCentered =
      Year - 2020,
    
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
      ),
    
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
  "MIXED-EFFECTS ANALYSIS\n",
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
  length(
    observed_years
  ) !=
  length(
    expected_years
  ) ||
  any(
    observed_years !=
    expected_years
  )
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
# 8. Validate resistant/non-resistant counts
# ----------------------------------------------------------

invalid_count_rows <- analysis_data %>%
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
  "\nInvalid count rows:",
  nrow(
    invalid_count_rows
  ),
  "\n"
)


if (
  nrow(
    invalid_count_rows
  ) > 0
) {
  
  stop(
    "Invalid resistant/non-resistant counts detected."
  )
}


############################################################
# PART B — REPRODUCE PRIMARY ELIGIBILITY DEFINITION
############################################################


# ----------------------------------------------------------
# 9. Coverage by pathogen-antibiotic combination and year
# ----------------------------------------------------------
#
# We deliberately reproduce the eligibility definition from
# Script 04 so both modelling approaches analyse the same
# pathogen-antibiotic combinations.
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
# 10. Summarise coverage
# ----------------------------------------------------------

overall_country_coverage <- analysis_data %>%
  group_by(
    PathogenName,
    AntibioticName
  ) %>%
  summarise(
    
    OverallCountries =
      n_distinct(
        Iso3
      ),
    
    .groups = "drop"
  )


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
    
    MedianCountriesPerYear =
      median(
        Countries
      ),
    
    MaxCountriesPerYear =
      max(
        Countries
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
  ) %>%
  left_join(
    overall_country_coverage,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  )


# ----------------------------------------------------------
# 11. Apply primary eligibility definition
# ----------------------------------------------------------
#
# Same primary definition as Script 04:
#
#   - represented in all four study years
#   - 2020 represented
#   - 2023 represented
#   - at least 20 contributing countries in EVERY year
#
# This produces the same 37-combination analysis family.
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
  "\n========================================\n",
  "PRIMARY ANALYSIS ELIGIBILITY\n",
  "========================================\n"
)


cat(
  "\nTotal pathogen-antibiotic combinations:",
  nrow(
    combination_coverage
  ),
  "\n"
)


cat(
  "Eligible pathogen-antibiotic combinations:",
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
  
  warning(
    "Expected 37 primary combinations based on Script 04, ",
    "but found ",
    nrow(
      primary_combinations
    ),
    ". Check whether the dataset or eligibility definition changed."
  )
}


# ----------------------------------------------------------
# 12. Construct primary inferential dataset
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
  droplevels()


cat(
  "\nPrimary analysis rows:",
  nrow(
    primary_analysis_data
  ),
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
# PART C — GLMM HELPER FUNCTIONS
############################################################


# ----------------------------------------------------------
# 13. Extract convergence message
# ----------------------------------------------------------
#
# lme4 stores convergence warnings inside the fitted model.
#
# No message means that lme4 did not record a convergence
# warning.
# ----------------------------------------------------------

get_convergence_message <- function(
    model
) {
  
  messages <-
    model@optinfo$conv$lme4$messages
  
  
  if (
    is.null(
      messages
    )
  ) {
    
    return(
      NA_character_
    )
  }
  
  
  paste(
    messages,
    collapse = " | "
  )
}


# ----------------------------------------------------------
# 14. Extract optimizer convergence code
# ----------------------------------------------------------
#
# An optimizer code of 0 generally indicates successful
# optimizer termination.
# ----------------------------------------------------------

get_optimizer_code <- function(
    model
) {
  
  optimizer_code <-
    model@optinfo$conv$opt
  
  
  if (
    is.null(
      optimizer_code
    )
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  as.numeric(
    optimizer_code
  )
}


# ----------------------------------------------------------
# 15. Calculate Pearson dispersion
# ----------------------------------------------------------
#
# Approximate diagnostic:
#
#   sum(Pearson residuals^2) / residual degrees of freedom
#
# Values substantially above 1 suggest extra-binomial
# variation.
#
# This is a diagnostic, not an automatic pass/fail rule.
# ----------------------------------------------------------

calculate_glmm_dispersion <- function(
    model
) {
  
  pearson_residuals <- residuals(
    model,
    type = "pearson"
  )
  
  
  dispersion_ratio <-
    sum(
      pearson_residuals^2,
      na.rm = TRUE
    ) /
    df.residual(
      model
    )
  
  
  dispersion_ratio
}


# ----------------------------------------------------------
# 16. Extract temporal effect from a fitted GLMM
# ----------------------------------------------------------
#
# This function:
#
#   - extracts the YearCentered coefficient
#   - calculates the Wald standard error and 95% CI
#   - converts log odds to odds ratios
#   - calculates the implied 2020-to-2023 OR
#   - records country random-intercept variation
#   - checks singularity
#   - records convergence information
#   - calculates Pearson dispersion
# ----------------------------------------------------------

extract_glmm_year_effect <- function(
    model,
    model_data,
    model_label
) {
  
  fixed_effects <- fixef(
    model
  )
  
  
  if (
    !"YearCentered" %in%
    names(
      fixed_effects
    )
  ) {
    
    stop(
      "YearCentered coefficient missing in model: ",
      model_label
    )
  }
  
  
  estimate <-
    unname(
      fixed_effects[
        "YearCentered"
      ]
    )
  
  
  standard_errors <-
    sqrt(
      diag(
        vcov(
          model
        )
      )
    )
  
  
  standard_error <-
    unname(
      standard_errors[
        "YearCentered"
      ]
    )
  
  
  if (
    is.na(
      standard_error
    ) ||
    !is.finite(
      standard_error
    )
  ) {
    
    stop(
      "Invalid YearCentered standard error in model: ",
      model_label
    )
  }
  
  
  z_value <-
    estimate /
    standard_error
  
  
  p_value <-
    2 *
    pnorm(
      abs(
        z_value
      ),
      lower.tail = FALSE
    )
  
  
  ci_log_low <-
    estimate -
    1.96 *
    standard_error
  
  
  ci_log_high <-
    estimate +
    1.96 *
    standard_error
  
  
  random_effect_sd <-
    as.numeric(
      attr(
        VarCorr(
          model
        )$Iso3,
        "stddev"
      )[1]
    )
  
  
  singular <-
    lme4::isSingular(
      model,
      tol = 1e-4
    )
  
  
  convergence_message <-
    get_convergence_message(
      model
    )
  
  
  optimizer_code <-
    get_optimizer_code(
      model
    )
  
  
  dispersion_ratio <-
    calculate_glmm_dispersion(
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
    
    PearsonDispersion =
      dispersion_ratio,
    
    Singular =
      singular,
    
    OptimizerCode =
      optimizer_code,
    
    ConvergenceMessage =
      convergence_message
  )
}


############################################################
# PART D — OVERALL ADJUSTED GLMM
############################################################


# ----------------------------------------------------------
# 17. Overall mixed-effects model
# ----------------------------------------------------------
#
# Model:
#
#   resistance ~
#     year
#     + pathogen-antibiotic combination
#     + country random intercept
#
# This estimates a common temporal association while
# adjusting for large baseline differences between specific
# pathogen-antibiotic combinations.
#
# Country baseline heterogeneity is represented by a random
# intercept.
#
# The model remains AST-volume weighted.
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "OVERALL ADJUSTED GLMM\n",
  "========================================\n"
)


overall_glmm <- glmer(
  
  cbind(
    Resistant,
    NonResistant
  ) ~
    
    YearCentered +
    
    PathogenAntibiotic +
    
    (
      1 |
        Iso3
    ),
  
  data =
    primary_analysis_data,
  
  family =
    binomial,
  
  control =
    glmerControl(
      
      optimizer =
        "bobyqa",
      
      optCtrl =
        list(
          maxfun = 200000
        )
    )
)


overall_glmm_result <-
  extract_glmm_year_effect(
    
    model =
      overall_glmm,
    
    model_data =
      primary_analysis_data,
    
    model_label =
      "Overall adjusted GLMM"
  )


print(
  overall_glmm_result,
  width = Inf
)


write_csv(
  overall_glmm_result,
  here(
    "outputs",
    "overall_adjusted_temporal_glmm.csv"
  )
)


############################################################
# PART E — PATHOGEN-SPECIFIC ADJUSTED GLMMs
############################################################


# ----------------------------------------------------------
# 18. Pathogen-specific models
# ----------------------------------------------------------
#
# Within each pathogen:
#
#   resistance ~
#     year
#     + antibiotic identity
#     + country random intercept
#
# Adjusting for antibiotic identity reduces the risk that
# changes in antibiotic-testing composition are interpreted
# as temporal resistance changes.
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "PATHOGEN-SPECIFIC ADJUSTED GLMMs\n",
  "========================================\n"
)


pathogens <- sort(
  unique(
    as.character(
      primary_analysis_data$PathogenName
    )
  )
)


pathogen_glmm_results_list <- list()


for (
  pathogen in pathogens
) {
  
  
  pathogen_data <- primary_analysis_data %>%
    filter(
      PathogenName ==
        pathogen
    ) %>%
    droplevels()
  
  
  cat(
    "\n----------------------------------------\n"
  )
  
  
  cat(
    "Pathogen:",
    pathogen,
    "\n"
  )
  
  
  pathogen_model <- glmer(
    
    cbind(
      Resistant,
      NonResistant
    ) ~
      
      YearCentered +
      
      AntibioticName +
      
      (
        1 |
          Iso3
      ),
    
    data =
      pathogen_data,
    
    family =
      binomial,
    
    control =
      glmerControl(
        
        optimizer =
          "bobyqa",
        
        optCtrl =
          list(
            maxfun = 200000
          )
      )
  )
  
  
  result <-
    extract_glmm_year_effect(
      
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
      
      .before =
        Model
    )
  
  
  print(
    result %>%
      select(
        PathogenName,
        Countries,
        Antibiotics,
        OddsRatioPerYear,
        CI_low,
        CI_high,
        P_value,
        PearsonDispersion,
        Singular,
        OptimizerCode,
        ConvergenceMessage
      ),
    width = Inf
  )
  
  
  pathogen_glmm_results_list[[as.character(pathogen)]] <-
    result
}


# ----------------------------------------------------------
# 19. Combine pathogen results and apply BH correction
# ----------------------------------------------------------

pathogen_glmm_results <- bind_rows(
  pathogen_glmm_results_list
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


write_csv(
  pathogen_glmm_results,
  here(
    "outputs",
    "pathogen_adjusted_temporal_glmm.csv"
  )
)


############################################################
# PART F — PRIMARY PATHOGEN-ANTIBIOTIC GLMMs
############################################################


# ----------------------------------------------------------
# 20. Fit one GLMM per primary combination
# ----------------------------------------------------------
#
# For each pathogen-antibiotic combination:
#
#   resistance ~
#     year
#     + country random intercept
#
# These are the mixed-effects counterparts of the
# combination-specific Script 04 models.
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "PATHOGEN–ANTIBIOTIC GLMMs\n",
  "========================================\n"
)


combination_glmm_results_list <- list()


model_failures <- list()


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
  
  
  model <- tryCatch(
    
    glmer(
      
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
        binomial,
      
      control =
        glmerControl(
          
          optimizer =
            "bobyqa",
          
          optCtrl =
            list(
              maxfun = 200000
            )
        )
    ),
    
    error = function(e) {
      
      model_failures[[model_label]] <<-
        e$message
      
      NULL
    }
  )
  
  
  if (
    is.null(
      model
    )
  ) {
    
    cat(
      "MODEL FAILED\n"
    )
    
    next
  }
  
  
  result <-
    extract_glmm_year_effect(
      
      model =
        model,
      
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
  
  
  combination_glmm_results_list[[i]] <-
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
    "| dispersion:",
    round(
      result$PearsonDispersion,
      2
    ),
    "| singular:",
    result$Singular,
    "| optimizer code:",
    result$OptimizerCode,
    "\n"
  )
  
  
  if (
    !is.na(
      result$ConvergenceMessage
    )
  ) {
    
    cat(
      "Convergence message:",
      result$ConvergenceMessage,
      "\n"
    )
  }
}


# ----------------------------------------------------------
# 21. Check for failed primary models
# ----------------------------------------------------------

if (
  length(
    model_failures
  ) > 0
) {
  
  failure_table <- tibble(
    
    Model =
      names(
        model_failures
      ),
    
    Error =
      unlist(
        model_failures
      )
  )
  
  
  write_csv(
    failure_table,
    here(
      "outputs",
      "glmm_model_failures.csv"
    )
  )
  
  
  print(
    failure_table
  )
  
  
  stop(
    "One or more primary GLMMs failed. ",
    "FDR correction has not been finalised."
  )
}


# ----------------------------------------------------------
# 22. Combine primary GLMM results
# ----------------------------------------------------------

combination_glmm_results <- bind_rows(
  combination_glmm_results_list
)


if (
  nrow(
    combination_glmm_results
  ) !=
  nrow(
    primary_combinations
  )
) {
  
  stop(
    "Not all primary combinations produced GLMM results."
  )
}


if (
  any(
    is.na(
      combination_glmm_results$P_value
    )
  )
) {
  
  stop(
    "At least one primary GLMM produced a missing p-value. ",
    "Investigate diagnostics before FDR correction."
  )
}


############################################################
# PART G — MULTIPLE TESTING
############################################################


# ----------------------------------------------------------
# 23. Benjamini-Hochberg FDR correction
# ----------------------------------------------------------
#
# The correction is applied across the same family of
# 37 primary pathogen-antibiotic combinations used in
# Script 04.
# ----------------------------------------------------------

combination_glmm_results <- combination_glmm_results %>%
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
  combination_glmm_results,
  here(
    "outputs",
    "pathogen_antibiotic_primary_temporal_glmm.csv"
  )
)


############################################################
# PART H — GLMM DIAGNOSTICS
############################################################


# ----------------------------------------------------------
# 24. Diagnostic summary
# ----------------------------------------------------------

glmm_diagnostics <- combination_glmm_results %>%
  select(
    
    PathogenName,
    
    AntibioticName,
    
    Countries,
    
    Observations,
    
    CountryRandomInterceptSD,
    
    PearsonDispersion,
    
    Singular,
    
    OptimizerCode,
    
    ConvergenceMessage
  )


write_csv(
  glmm_diagnostics,
  here(
    "outputs",
    "pathogen_antibiotic_glmm_diagnostics.csv"
  )
)


cat(
  "\n========================================\n",
  "GLMM DIAGNOSTIC SUMMARY\n",
  "========================================\n"
)


cat(
  "\nSingular models:",
  sum(
    glmm_diagnostics$Singular,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Models with non-zero optimizer code:",
  sum(
    glmm_diagnostics$OptimizerCode != 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Models with convergence messages:",
  sum(
    !is.na(
      glmm_diagnostics$ConvergenceMessage
    )
  ),
  "\n"
)


cat(
  "\nPearson dispersion summary:\n"
)


print(
  summary(
    glmm_diagnostics$PearsonDispersion
  )
)


# ----------------------------------------------------------
# 25. Identify models with potential diagnostic concerns
# ----------------------------------------------------------

diagnostic_concerns <- glmm_diagnostics %>%
  filter(
    
    Singular |
      
      OptimizerCode != 0 |
      
      !is.na(
        ConvergenceMessage
      ) |
      
      PearsonDispersion > 2
  ) %>%
  arrange(
    desc(
      PearsonDispersion
    )
  )


cat(
  "\nModels requiring diagnostic attention:",
  nrow(
    diagnostic_concerns
  ),
  "\n"
)


write_csv(
  diagnostic_concerns,
  here(
    "outputs",
    "glmm_diagnostic_concerns.csv"
  )
)


############################################################
# PART I — FDR-SIGNIFICANT GLMM RESULTS
############################################################


# ----------------------------------------------------------
# 26. Detailed results with BH-adjusted p < 0.05
# ----------------------------------------------------------
#
# This is a convenience table for inspection.
#
# FDR < 0.05 alone does NOT establish that a result is:
#
#   - robust,
#   - clinically important,
#   - causal,
#   - or suitable for final emphasis.
#
# Diagnostics and sensitivity analyses must also be
# considered.
# ----------------------------------------------------------

fdr_glmm_results <- combination_glmm_results %>%
  filter(
    FDR_P_value <
      0.05
  )


cat(
  "\n========================================\n",
  "GLMM RESULTS WITH BH-ADJUSTED P < 0.05\n",
  "========================================\n"
)


print(
  fdr_glmm_results %>%
    select(
      
      PathogenName,
      
      AntibioticName,
      
      Countries,
      
      OddsRatioPerYear,
      
      CI_low,
      
      CI_high,
      
      OR_2023_vs_2020,
      
      OR_2023_vs_2020_CI_low,
      
      OR_2023_vs_2020_CI_high,
      
      P_value,
      
      FDR_P_value,
      
      PearsonDispersion,
      
      Singular,
      
      OptimizerCode,
      
      ConvergenceMessage
    ),
  
  n = Inf,
  
  width = Inf
)


write_csv(
  fdr_glmm_results,
  here(
    "outputs",
    "fdr_below_0_05_primary_glmm.csv"
  )
)


############################################################
# PART J — COMPARE SCRIPT 04 AND SCRIPT 05
############################################################


# ----------------------------------------------------------
# 27. Load Script 04 robust-model results
# ----------------------------------------------------------
#
# Script 04 and Script 05 address closely related questions
# using different treatments of country-level heterogeneity.
#
# Script 04:
#   country fixed effects + cluster-robust SE
#
# Script 05:
#   country random intercept
#
# Agreement is evidence of robustness.
#
# Disagreement is also scientifically important and must be
# investigated rather than hidden.
# ----------------------------------------------------------

robust_results_path <- here(
  "outputs",
  "final_temporal_trend_results_primary.csv"
)


if (
  file.exists(
    robust_results_path
  )
) {
  
  
  robust_results <- read_csv(
    robust_results_path,
    show_col_types = FALSE
  )
  
  
  robust_for_comparison <- robust_results %>%
    select(
      
      PathogenName,
      
      AntibioticName,
      
      Robust_OR_per_year =
        OddsRatioPerYear,
      
      Robust_CI_low =
        CI_low,
      
      Robust_CI_high =
        CI_high,
      
      Robust_P =
        P_value,
      
      Robust_FDR =
        FDR_P_value
    )
  
  
  glmm_for_comparison <- combination_glmm_results %>%
    select(
      
      PathogenName,
      
      AntibioticName,
      
      GLMM_OR_per_year =
        OddsRatioPerYear,
      
      GLMM_CI_low =
        CI_low,
      
      GLMM_CI_high =
        CI_high,
      
      GLMM_P =
        P_value,
      
      GLMM_FDR =
        FDR_P_value,
      
      GLMM_Dispersion =
        PearsonDispersion,
      
      GLMM_Singular =
        Singular,
      
      GLMM_OptimizerCode =
        OptimizerCode,
      
      GLMM_ConvergenceMessage =
        ConvergenceMessage
    )
  
  
  model_comparison <- robust_for_comparison %>%
    inner_join(
      
      glmm_for_comparison,
      
      by = c(
        "PathogenName",
        "AntibioticName"
      )
    ) %>%
    mutate(
      
      RobustDirection =
        case_when(
          
          Robust_OR_per_year > 1 ~
            "Increase",
          
          Robust_OR_per_year < 1 ~
            "Decrease",
          
          TRUE ~
            "No change"
        ),
      
      GLMMDirection =
        case_when(
          
          GLMM_OR_per_year > 1 ~
            "Increase",
          
          GLMM_OR_per_year < 1 ~
            "Decrease",
          
          TRUE ~
            "No change"
        ),
      
      DirectionAgreement =
        RobustDirection ==
        GLMMDirection,
      
      RobustFDRBelow05 =
        Robust_FDR <
        0.05,
      
      GLMMFDRBelow05 =
        GLMM_FDR <
        0.05,
      
      SignificantInBoth =
        RobustFDRBelow05 &
        GLMMFDRBelow05
    ) %>%
    arrange(
      Robust_FDR
    )
  
  
  if (
    nrow(
      model_comparison
    ) != 37
  ) {
    
    warning(
      "Script 04 vs Script 05 comparison contains ",
      nrow(
        model_comparison
      ),
      " combinations rather than the expected 37."
    )
  }
  
  
  write_csv(
    model_comparison,
    here(
      "outputs",
      "model_comparison_robust_vs_glmm.csv"
    )
  )
  
  
  cat(
    "\n========================================\n",
    "SCRIPT 04 vs SCRIPT 05 COMPARISON\n",
    "========================================\n"
  )
  
  
  cat(
    "\nCombinations compared:",
    nrow(
      model_comparison
    ),
    "\n"
  )
  
  
  cat(
    "Direction agreement:",
    sum(
      model_comparison$DirectionAgreement,
      na.rm = TRUE
    ),
    "/",
    nrow(
      model_comparison
    ),
    "\n"
  )
  
  
  cat(
    "FDR < 0.05 in Script 04:",
    sum(
      model_comparison$RobustFDRBelow05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "FDR < 0.05 in GLMM:",
    sum(
      model_comparison$GLMMFDRBelow05,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "FDR < 0.05 in BOTH approaches:",
    sum(
      model_comparison$SignificantInBoth,
      na.rm = TRUE
    ),
    "\n"
  )
  
  
  cat(
    "\nCombinations significant in either approach:\n"
  )
  
  
  print(
    model_comparison %>%
      filter(
        RobustFDRBelow05 |
          GLMMFDRBelow05
      ) %>%
      select(
        
        PathogenName,
        
        AntibioticName,
        
        Robust_OR_per_year,
        
        Robust_CI_low,
        
        Robust_CI_high,
        
        Robust_FDR,
        
        GLMM_OR_per_year,
        
        GLMM_CI_low,
        
        GLMM_CI_high,
        
        GLMM_FDR,
        
        SignificantInBoth,
        
        GLMM_Dispersion,
        
        GLMM_Singular,
        
        GLMM_OptimizerCode,
        
        GLMM_ConvergenceMessage
      ),
    
    n = Inf,
    
    width = Inf
  )
  
  
} else {
  
  
  warning(
    paste(
      "Script 04 results were not found.",
      "The GLMM analysis completed, but the",
      "cross-model comparison was skipped."
    )
  )
}


############################################################
# PART K — FINAL MODEL SUMMARIES
############################################################


# ----------------------------------------------------------
# 28. Pathogen-specific GLMM summary
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "PATHOGEN-SPECIFIC GLMM SUMMARY\n",
  "========================================\n"
)


print(
  pathogen_glmm_results %>%
    select(
      
      PathogenName,
      
      Countries,
      
      Antibiotics,
      
      OddsRatioPerYear,
      
      CI_low,
      
      CI_high,
      
      P_value,
      
      FDR_P_value,
      
      OR_2023_vs_2020,
      
      PearsonDispersion,
      
      CountryRandomInterceptSD,
      
      Singular,
      
      OptimizerCode,
      
      ConvergenceMessage
    ),
  
  n = Inf,
  
  width = Inf
)


# ----------------------------------------------------------
# 29. Overall adjusted GLMM summary
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "OVERALL GLMM SUMMARY\n",
  "========================================\n"
)


print(
  overall_glmm_result %>%
    select(
      
      Countries,
      
      OddsRatioPerYear,
      
      CI_low,
      
      CI_high,
      
      P_value,
      
      OR_2023_vs_2020,
      
      OR_2023_vs_2020_CI_low,
      
      OR_2023_vs_2020_CI_high,
      
      CountryRandomInterceptSD,
      
      PearsonDispersion,
      
      Singular,
      
      OptimizerCode,
      
      ConvergenceMessage
    ),
  
  width = Inf
)


############################################################
# PART L — FINAL MESSAGE
############################################################


# ----------------------------------------------------------
# 30. Completion summary
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "MIXED-EFFECTS ANALYSIS COMPLETE\n",
  "========================================\n"
)


cat(
  "\nPrimary combinations modelled:",
  nrow(
    combination_glmm_results
  ),
  "\n"
)


cat(
  "BH-adjusted p < 0.05:",
  sum(
    combination_glmm_results$FDR_below_0_05,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Singular primary models:",
  sum(
    combination_glmm_results$Singular,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Primary models with convergence messages:",
  sum(
    !is.na(
      combination_glmm_results$ConvergenceMessage
    )
  ),
  "\n"
)


cat(
  "\nResults saved to:\n",
  here(
    "outputs"
  ),
  "\n"
)


cat(
  "\nMain files:\n",
  "- overall_adjusted_temporal_glmm.csv\n",
  "- pathogen_adjusted_temporal_glmm.csv\n",
  "- pathogen_antibiotic_primary_temporal_glmm.csv\n",
  "- pathogen_antibiotic_glmm_diagnostics.csv\n",
  "- glmm_diagnostic_concerns.csv\n",
  "- fdr_below_0_05_primary_glmm.csv\n",
  "- model_comparison_robust_vs_glmm.csv\n"
)


cat(
  "\nIMPORTANT:\n",
  paste(
    "The GLMM results are a robustness analysis.",
    "Final scientific interpretation should consider",
    "agreement between modelling approaches, effect sizes,",
    "confidence intervals, model diagnostics, surveillance",
    "coverage, and additional sensitivity analyses."
  ),
  "\n"
)