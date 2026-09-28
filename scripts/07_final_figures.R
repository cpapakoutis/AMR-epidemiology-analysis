############################################################
# Project:
# Temporal and Geographic Patterns in Antimicrobial
# Resistance Reported to WHO GLASS, 2020–2023
#
# Script:
# 07_final_figures.R
#
# Purpose:
#   Create final candidate figures based on the locked
#   descriptive, primary inferential, and sensitivity
#   analyses.
#
# This script DOES NOT fit new statistical models.
#
# Figures:
#
#   Figure 6:
#     Temporal trajectories for pathogen-antibiotic
#     combinations identified in the primary analysis
#
#   Figure 7:
#     Sensitivity comparison for primary signals
#
#   Supplementary Figure S1:
#     Forest plot of all 37 primary temporal estimates
#
# Inputs:
#   data/processed/AMR_GLASS_clean.csv
#   outputs/final_temporal_trend_results_primary.csv
#   outputs/final_sensitivity_robustness_summary.csv
#
# Outputs:
#   figures/06_selected_temporal_trajectories.png
#   figures/07_primary_signal_sensitivity_comparison.png
#   figures/S1_primary_temporal_forest_plot.png
#
# Supporting figure datasets are saved to outputs/.
############################################################


# ----------------------------------------------------------
# 1. Load packages
# ----------------------------------------------------------

library(tidyverse)
library(here)


# ----------------------------------------------------------
# 2. Create directories if needed
# ----------------------------------------------------------

dir.create(
  here("figures"),
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  here("outputs"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ----------------------------------------------------------
# 3. Define consistent plotting theme
# ----------------------------------------------------------
#
# White background and restrained formatting are used to
# match the earlier project figures.
#
# Slightly larger plot margins are used to prevent titles,
# subtitles and captions from being clipped in exported
# images.
# ----------------------------------------------------------

theme_amr <- theme_minimal(
  base_size = 12
) +
  theme(
    
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major.x =
      element_line(
        linewidth = 0.3
      ),
    
    panel.grid.major.y =
      element_line(
        linewidth = 0.3
      ),
    
    plot.background =
      element_rect(
        fill = "white",
        colour = NA
      ),
    
    panel.background =
      element_rect(
        fill = "white",
        colour = NA
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 16
      ),
    
    plot.subtitle =
      element_text(
        size = 11,
        margin = margin(
          b = 10
        )
      ),
    
    plot.caption =
      element_text(
        size = 9,
        hjust = 0,
        margin = margin(
          t = 12
        )
      ),
    
    plot.title.position =
      "plot",
    
    plot.caption.position =
      "plot",
    
    plot.margin =
      margin(
        t = 15,
        r = 20,
        b = 20,
        l = 15
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 11
      ),
    
    axis.title =
      element_text(
        face = "bold"
      ),
    
    legend.position =
      "bottom",
    
    legend.title =
      element_blank()
  )


############################################################
# PART A — LOAD AND VALIDATE RESULTS
############################################################


# ----------------------------------------------------------
# 4. Define required input paths
# ----------------------------------------------------------

clean_data_path <- here(
  "data",
  "processed",
  "AMR_GLASS_clean.csv"
)


primary_results_path <- here(
  "outputs",
  "final_temporal_trend_results_primary.csv"
)


sensitivity_results_path <- here(
  "outputs",
  "final_sensitivity_robustness_summary.csv"
)


# ----------------------------------------------------------
# 5. Check required files exist
# ----------------------------------------------------------

required_files <- c(
  clean_data_path,
  primary_results_path,
  sensitivity_results_path
)


missing_files <- required_files[
  !file.exists(
    required_files
  )
]


if (
  length(
    missing_files
  ) > 0
) {
  
  stop(
    "Required files are missing:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ----------------------------------------------------------
# 6. Load files
# ----------------------------------------------------------

analysis_data <- read_csv(
  clean_data_path,
  show_col_types = FALSE
)


primary_results <- read_csv(
  primary_results_path,
  show_col_types = FALSE
)


sensitivity_results <- read_csv(
  sensitivity_results_path,
  show_col_types = FALSE
)


# ----------------------------------------------------------
# 7. Basic validation
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "FINAL FIGURES\n",
  "========================================\n"
)


cat(
  "\nClean data rows:",
  nrow(
    analysis_data
  ),
  "\n"
)


cat(
  "Primary model combinations:",
  nrow(
    primary_results
  ),
  "\n"
)


cat(
  "Sensitivity summary combinations:",
  nrow(
    sensitivity_results
  ),
  "\n"
)


if (
  nrow(
    primary_results
  ) != 37
) {
  
  warning(
    "Expected 37 primary model results but found ",
    nrow(
      primary_results
    ),
    "."
  )
}


if (
  nrow(
    sensitivity_results
  ) != 37
) {
  
  warning(
    "Expected 37 sensitivity-summary rows but found ",
    nrow(
      sensitivity_results
    ),
    "."
  )
}


############################################################
# PART B — IDENTIFY PRIMARY SIGNALS
############################################################


# ----------------------------------------------------------
# 8. Identify combinations meeting primary BH-FDR < 0.05
# ----------------------------------------------------------
#
# These are selected from the predefined primary analysis.
#
# They are NOT manually selected after examining the
# sensitivity analyses.
# ----------------------------------------------------------

primary_signal_keys <- primary_results %>%
  filter(
    FDR_P_value <
      0.05
  ) %>%
  select(
    PathogenName,
    AntibioticName
  )


cat(
  "\nPrimary combinations with BH-FDR < 0.05:",
  nrow(
    primary_signal_keys
  ),
  "\n"
)


print(
  primary_signal_keys
)


if (
  nrow(
    primary_signal_keys
  ) == 0
) {
  
  stop(
    "No primary FDR < 0.05 combinations found. ",
    "Selected-signal figures cannot be created."
  )
}


############################################################
# PART C — FIGURE 6:
# SELECTED TEMPORAL TRAJECTORIES
############################################################


# ----------------------------------------------------------
# 9. Restrict clean data to primary signals
# ----------------------------------------------------------

selected_data <- analysis_data %>%
  semi_join(
    primary_signal_keys,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  )


# ----------------------------------------------------------
# 10. Calculate pooled AST-weighted percentages
# ----------------------------------------------------------
#
# This is the same general estimand used in the descriptive
# pooled summaries:
#
#   total resistant AST /
#   total interpretable AST
#
# Countries contributing more AST results therefore receive
# greater weight.
# ----------------------------------------------------------

pooled_selected_yearly <- selected_data %>%
  group_by(
    PathogenName,
    AntibioticName,
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
    
    ResistancePercentage =
      Resistant /
      InterpretableAST *
      100,
    
    .groups = "drop"
  ) %>%
  mutate(
    Approach =
      "Pooled AST-weighted"
  )


# ----------------------------------------------------------
# 11. Identify stable countries for each selected combination
# ----------------------------------------------------------
#
# Stable countries contribute that specific
# pathogen-antibiotic combination in all four study years.
# ----------------------------------------------------------

stable_selected_country_keys <- selected_data %>%
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


# ----------------------------------------------------------
# 12. Calculate stable-country equal-weight percentages
# ----------------------------------------------------------
#
# First calculate resistance within each country-year.
#
# Then take the simple mean across stable countries.
#
# Each country therefore receives equal weight regardless
# of its AST volume.
# ----------------------------------------------------------

stable_country_year <- selected_data %>%
  semi_join(
    stable_selected_country_keys,
    by = c(
      "PathogenName",
      "AntibioticName",
      "Iso3"
    )
  ) %>%
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


stable_equal_weight_yearly <- stable_country_year %>%
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
    
    ResistancePercentage =
      mean(
        CountryResistancePercentage,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Approach =
      "Stable countries, equal weight"
  )


# ----------------------------------------------------------
# 13. Combine trajectory datasets
# ----------------------------------------------------------

trajectory_plot_data <- bind_rows(
  
  pooled_selected_yearly %>%
    select(
      PathogenName,
      AntibioticName,
      Year,
      ResistancePercentage,
      Approach
    ),
  
  stable_equal_weight_yearly %>%
    select(
      PathogenName,
      AntibioticName,
      Year,
      ResistancePercentage,
      Approach
    )
) %>%
  mutate(
    
    Combination =
      paste(
        PathogenName,
        AntibioticName,
        sep = " — "
      ),
    
    Approach =
      factor(
        Approach,
        levels = c(
          "Pooled AST-weighted",
          "Stable countries, equal weight"
        )
      )
  )


write_csv(
  trajectory_plot_data,
  here(
    "outputs",
    "figure06_selected_trajectory_data.csv"
  )
)


# ----------------------------------------------------------
# 14. Create Figure 6
# ----------------------------------------------------------
#
# Panels deliberately use independent y-axis scales because
# the absolute resistance levels differ substantially between
# the two pathogen-antibiotic combinations.
# ----------------------------------------------------------

figure_06 <- ggplot(
  trajectory_plot_data,
  aes(
    x = Year,
    y = ResistancePercentage,
    group = Approach,
    linetype = Approach,
    shape = Approach
  )
) +
  
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 3
  ) +
  
  facet_wrap(
    vars(
      Combination
    ),
    scales = "free_y"
  ) +
  
  scale_x_continuous(
    breaks = 2020:2023
  ) +
  
  scale_y_continuous(
    labels = function(x) {
      paste0(
        round(
          x,
          1
        ),
        "%"
      )
    }
  ) +
  
  labs(
    
    title =
      "Temporal patterns in the two primary pathogen–antibiotic signals",
    
    subtitle =
      stringr::str_wrap(
        paste(
          "Pooled AST-weighted estimates are compared with",
          "equal-weight estimates among countries reporting",
          "the combination in all four years."
        ),
        width = 105
      ),
    
    x =
      "Year",
    
    y =
      "Reported resistance (%)",
    
    caption =
      stringr::str_wrap(
        paste(
          "Primary signals were selected using the predefined",
          "Script 04 BH-FDR < 0.05 criterion.",
          "Panels use independent y-axis scales.",
          "Percentages describe reported surveillance data and",
          "should not be interpreted as population-representative",
          "global prevalence."
        ),
        width = 130
      )
  ) +
  
  theme_amr


ggsave(
  filename =
    here(
      "figures",
      "06_selected_temporal_trajectories.png"
    ),
  plot =
    figure_06,
  width =
    12,
  height =
    8,
  dpi =
    320,
  bg =
    "white"
)


############################################################
# PART D — SUPPLEMENTARY FIGURE S1:
# PRIMARY FOREST PLOT
############################################################


# ----------------------------------------------------------
# 15. Prepare Supplementary Figure S1 forest-plot data
# ----------------------------------------------------------
#
# This figure shows ALL 37 primary estimates.
#
# Showing all predefined models avoids displaying only
# combinations that crossed a statistical threshold.
#
# This figure is retained as Supplementary Figure S1
# because it displays the complete set of 37 estimates.
# ----------------------------------------------------------

forest_plot_data <- primary_results %>%
  mutate(
    
    CombinationID =
      paste(
        PathogenName,
        AntibioticName,
        sep = "___"
      ),
    
    CombinationID =
      forcats::fct_reorder(
        CombinationID,
        OddsRatioPerYear
      ),
    
    FDRStatus =
      if_else(
        FDR_P_value < 0.05,
        "BH-FDR < 0.05",
        "BH-FDR ≥ 0.05"
      )
  )


write_csv(
  forest_plot_data,
  here(
    "outputs",
    "S1_primary_forest_plot_data.csv"
  )
)


# ----------------------------------------------------------
# 16. Create Supplementary Figure S1
# ----------------------------------------------------------

figure_s1 <- ggplot(
  forest_plot_data,
  aes(
    x = OddsRatioPerYear,
    y = CombinationID
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      y = CombinationID,
      yend = CombinationID
    ),
    linewidth = 0.7
  ) +
  
  geom_point(
    aes(
      shape = FDRStatus
    ),
    size = 2.8
  ) +
  
  facet_grid(
    rows = vars(
      PathogenName
    ),
    scales = "free_y",
    space = "free_y"
  ) +
  
  scale_x_log10(
    breaks = c(
      0.5,
      0.75,
      1,
      1.25,
      1.5,
      2,
      2.5
    ),
    labels =
      scales::label_number(
        accuracy = 0.01
      )
  ) +
  
  scale_y_discrete(
    labels = function(x) {
      stringr::str_replace(
        x,
        "^.*___",
        ""
      )
    }
  ) +
  
  scale_shape_manual(
    values = c(
      "BH-FDR < 0.05" = 16,
      "BH-FDR ≥ 0.05" = 1
    )
  ) +
  
  labs(
    
    title =
      "Primary temporal estimates across pathogen–antibiotic combinations",
    
    subtitle =
      stringr::str_wrap(
        paste(
          "Odds ratios per additional year, with country fixed effects",
          "and country-clustered robust standard errors."
        ),
        width = 100
      ),
    
    x =
      "Odds ratio per year (log scale)",
    
    y =
      NULL,
    
    caption =
      stringr::str_wrap(
        paste(
          "Vertical dashed line indicates OR = 1.",
          "BH-FDR correction was applied across the 37",
          "predefined primary pathogen–antibiotic models.",
          "Filled points indicate BH-FDR < 0.05."
        ),
        width = 120
      )
  ) +
  
  theme_amr +
  
  theme(
    
    strip.text.y =
      element_text(
        angle = 0,
        face = "bold"
      ),
    
    axis.text.y =
      element_text(
        size = 9
      )
  )


ggsave(
  filename =
    here(
      "figures",
      "S1_primary_temporal_forest_plot.png"
    ),
  plot =
    figure_s1,
  width =
    12,
  height =
    15.5,
  dpi =
    320,
  bg =
    "white"
)


############################################################
# PART E — FIGURE 7:
# SENSITIVITY COMPARISON
############################################################


# ----------------------------------------------------------
# 17. Restrict sensitivity summary to primary signals
# ----------------------------------------------------------

selected_sensitivity <- sensitivity_results %>%
  semi_join(
    primary_signal_keys,
    by = c(
      "PathogenName",
      "AntibioticName"
    )
  )


# ----------------------------------------------------------
# 18. Convert continuous-year models to a common
#     2023-versus-2020 scale
# ----------------------------------------------------------
#
# The following models estimate an OR per year:
#
#   - Primary robust
#   - Beta-binomial
#   - Stable-country
#
# Because 2023 is three years after 2020:
#
#   OR_2023_vs_2020 = OR_per_year ^ 3
#
# Their confidence limits are transformed in the same way.
#
# The categorical-year model already directly estimates
# 2023 versus 2020.
#
# This puts all four approaches on the same effect scale.
# ----------------------------------------------------------

sensitivity_plot_data <- bind_rows(
  
  selected_sensitivity %>%
    transmute(
      
      PathogenName,
      
      AntibioticName,
      
      Method =
        "Primary robust",
      
      OR_2023_vs_2020 =
        Primary_OR_per_year^3,
      
      CI_low =
        Primary_CI_low^3,
      
      CI_high =
        Primary_CI_high^3
    ),
  
  
  selected_sensitivity %>%
    transmute(
      
      PathogenName,
      
      AntibioticName,
      
      Method =
        "Beta-binomial",
      
      OR_2023_vs_2020 =
        BetaBinomial_OR_per_year^3,
      
      CI_low =
        BetaBinomial_CI_low^3,
      
      CI_high =
        BetaBinomial_CI_high^3
    ),
  
  
  selected_sensitivity %>%
    transmute(
      
      PathogenName,
      
      AntibioticName,
      
      Method =
        "Stable-country",
      
      OR_2023_vs_2020 =
        StableCountry_OR_per_year^3,
      
      CI_low =
        StableCountry_CI_low^3,
      
      CI_high =
        StableCountry_CI_high^3
    ),
  
  
  selected_sensitivity %>%
    transmute(
      
      PathogenName,
      
      AntibioticName,
      
      Method =
        "Categorical year",
      
      OR_2023_vs_2020 =
        Categorical_OR_2023_vs_2020,
      
      CI_low =
        Categorical_CI_low,
      
      CI_high =
        Categorical_CI_high
    )
) %>%
  mutate(
    
    Combination =
      paste(
        PathogenName,
        AntibioticName,
        sep = " — "
      ),
    
    Method =
      factor(
        Method,
        levels = rev(
          c(
            "Primary robust",
            "Beta-binomial",
            "Stable-country",
            "Categorical year"
          )
        )
      )
  )


write_csv(
  sensitivity_plot_data,
  here(
    "outputs",
    "figure07_sensitivity_comparison_data.csv"
  )
)


# ----------------------------------------------------------
# 19. Create Figure 7
# ----------------------------------------------------------

figure_07 <- ggplot(
  sensitivity_plot_data,
  aes(
    x = OR_2023_vs_2020,
    y = Method
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      y = Method,
      yend = Method
    ),
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 3
  ) +
  
  facet_wrap(
    vars(
      Combination
    ),
    ncol = 1
  ) +
  
  scale_x_log10(
    breaks = c(
      0.75,
      1,
      1.25,
      1.5,
      2
    ),
    labels =
      scales::label_number(
        accuracy = 0.01
      )
  ) +
  
  labs(
    
    title =
      "Sensitivity analyses of the two primary temporal signals",
    
    subtitle =
      stringr::str_wrap(
        paste(
          "Estimates from alternative model specifications",
          "are shown on a common 2023-versus-2020",
          "odds-ratio scale."
        ),
        width = 100
      ),
    
    x =
      "Odds ratio: 2023 versus 2020 (log scale)",
    
    y =
      NULL,
    
    caption =
      stringr::str_wrap(
        paste(
          "Primary robust: country fixed effects with",
          "country-clustered robust SEs.",
          "Beta-binomial: country random intercept with",
          "extra-binomial variation.",
          "Stable-country: restricted to countries reporting",
          "the combination in all four years.",
          "For these three continuous-year models, per-year",
          "odds ratios were converted to 2023-versus-2020",
          "effects over three years.",
          "Categorical year directly compares 2023 with 2020."
        ),
        width = 125
      )
  ) +
  
  theme_amr


ggsave(
  filename =
    here(
      "figures",
      "07_primary_signal_sensitivity_comparison.png"
    ),
  plot =
    figure_07,
  width =
    11,
  height =
    9,
  dpi =
    320,
  bg =
    "white"
)


############################################################
# PART F — FINAL VALIDATION
############################################################


# ----------------------------------------------------------
# 20. Confirm expected files were created
# ----------------------------------------------------------

expected_figure_files <- c(
  
  here(
    "figures",
    "06_selected_temporal_trajectories.png"
  ),
  
  here(
    "figures",
    "S1_primary_temporal_forest_plot.png"
  ),
  
  here(
    "figures",
    "07_primary_signal_sensitivity_comparison.png"
  )
)


figure_created <- file.exists(
  expected_figure_files
)


figure_validation <- tibble(
  
  Figure =
    basename(
      expected_figure_files
    ),
  
  Created =
    figure_created
)


cat(
  "\n========================================\n",
  "FIGURE CREATION CHECK\n",
  "========================================\n"
)


print(
  figure_validation
)


if (
  !all(
    figure_created
  )
) {
  
  stop(
    "One or more final figure files were not created."
  )
}


############################################################
# PART G — FINAL MESSAGE
############################################################


cat(
  "\n========================================\n",
  "FINAL FIGURES CREATED SUCCESSFULLY\n",
  "========================================\n"
)


cat(
  "\nCreated:\n",
  "- figures/06_selected_temporal_trajectories.png\n",
  "- figures/07_primary_signal_sensitivity_comparison.png\n",
  "- figures/S1_primary_temporal_forest_plot.png\n"
)


cat(
  "\nSupporting data saved to:\n",
  "- outputs/figure06_selected_trajectory_data.csv\n",
  "- outputs/figure07_sensitivity_comparison_data.csv\n",
  "- outputs/S1_primary_forest_plot_data.csv\n"
)


cat(
  "\nCURRENT FIGURE PLAN:\n",
  paste(
    "Figures 6 and 7 are retained for the main report.",
    "Supplementary Figure S1 provides the complete set",
    "of 37 primary temporal estimates."
  ),
  "\n"
)


cat(
  "\nIMPORTANT:\n",
  paste(
    "This script only visualises results generated in",
    "earlier scripts. No statistical models are fitted",
    "or re-estimated here."
  ),
  "\n"
)