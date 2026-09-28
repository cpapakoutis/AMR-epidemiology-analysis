############################################################
# Project: WHO GLASS Antimicrobial Resistance Analysis
# Script: 03_exploratory_analysis.R
#
# Purpose:
#   Describe the structure of the cleaned WHO GLASS
#   surveillance dataset and explore resistance patterns
#   before inferential modelling.
#
# Important interpretation:
#   InterpretableAST counts antibiotic susceptibility tests
#   (AST results), not necessarily unique isolates.
#
#   Therefore, summaries that pool across antibiotics represent
#   the proportion of reported AST results classified as
#   resistant.
#
#   They should NOT be interpreted as the proportion of
#   infections or unique isolates that are resistant.
#
# Input:
#   data/processed/AMR_GLASS_clean.csv
#
# Outputs:
#   Exploratory summary tables in outputs/
#   Exploratory figures in figures/
############################################################


# ----------------------------------------------------------
# 1. Load packages
# ----------------------------------------------------------

library(tidyverse)
library(here)
library(viridis)


# ----------------------------------------------------------
# 2. Load cleaned analytical dataset
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
# 3. Confirm required variables
# ----------------------------------------------------------

required_variables <- c(
  "Specimen",
  "PathogenName",
  "AntibioticName",
  "Iso3",
  "CountryTerritoryArea",
  "WHORegionName",
  "InterpretableAST",
  "Resistant",
  "ResistancePercentage",
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
# 4. Create output directories
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
# 5. Create common figure theme
# ----------------------------------------------------------
#
# Using one shared theme keeps all project figures visually
# consistent.
# ----------------------------------------------------------

theme_amr <- theme_minimal(
  base_size = 12
) +
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 11,
      margin = margin(
        b = 10
      )
    ),
    
    plot.caption = element_text(
      size = 9,
      hjust = 0,
      margin = margin(
        t = 10
      )
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom",
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.margin = margin(
      15,
      20,
      15,
      15
    )
  )


############################################################
# PART A — SURVEILLANCE STRUCTURE
############################################################


# ----------------------------------------------------------
# 6. Basic dataset overview
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "EXPLORATORY ANALYSIS\n",
  "========================================\n"
)


cat(
  "\nRows:",
  nrow(amr_data),
  "\n"
)


cat(
  "Countries/territories:",
  n_distinct(amr_data$Iso3),
  "\n"
)


cat(
  "Years:",
  paste(
    sort(
      unique(
        amr_data$Year
      )
    ),
    collapse = ", "
  ),
  "\n"
)


cat(
  "Pathogens:",
  n_distinct(
    amr_data$PathogenName
  ),
  "\n"
)


cat(
  "Antibiotics:",
  n_distinct(
    amr_data$AntibioticName
  ),
  "\n"
)


# ----------------------------------------------------------
# 7. Surveillance coverage by year
# ----------------------------------------------------------

coverage_by_year <- amr_data %>%
  group_by(
    Year
  ) %>%
  summarise(
    
    Records = n(),
    
    Countries = n_distinct(
      Iso3
    ),
    
    Pathogens = n_distinct(
      PathogenName
    ),
    
    Antibiotics = n_distinct(
      AntibioticName
    ),
    
    InterpretableAST = sum(
      InterpretableAST,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


print(
  coverage_by_year
)


write_csv(
  coverage_by_year,
  here(
    "outputs",
    "surveillance_coverage_by_year.csv"
  )
)


# ----------------------------------------------------------
# 8. Surveillance coverage by pathogen and year
# ----------------------------------------------------------

coverage_by_pathogen_year <- amr_data %>%
  group_by(
    PathogenName,
    Year
  ) %>%
  summarise(
    
    Countries = n_distinct(
      Iso3
    ),
    
    Antibiotics = n_distinct(
      AntibioticName
    ),
    
    Records = n(),
    
    InterpretableAST = sum(
      InterpretableAST,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


print(
  coverage_by_pathogen_year
)


write_csv(
  coverage_by_pathogen_year,
  here(
    "outputs",
    "surveillance_coverage_by_pathogen_year.csv"
  )
)


# ----------------------------------------------------------
# 9. Surveillance coverage by WHO region and year
# ----------------------------------------------------------

coverage_by_region_year <- amr_data %>%
  group_by(
    WHORegionName,
    Year
  ) %>%
  summarise(
    
    Countries = n_distinct(
      Iso3
    ),
    
    InterpretableAST = sum(
      InterpretableAST,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


write_csv(
  coverage_by_region_year,
  here(
    "outputs",
    "surveillance_coverage_by_region_year.csv"
  )
)


# ----------------------------------------------------------
# 10. Country participation consistency
# ----------------------------------------------------------
#
# This identifies how many of the four study years each
# country/territory contributes data to.
#
# It will later help motivate and construct the complete-case
# sensitivity analysis restricted to countries contributing
# data across all four years.
# ----------------------------------------------------------

country_year_presence <- amr_data %>%
  distinct(
    Iso3,
    CountryTerritoryArea,
    Year
  )


country_participation <- country_year_presence %>%
  group_by(
    Iso3,
    CountryTerritoryArea
  ) %>%
  summarise(
    YearsReported = n_distinct(
      Year
    ),
    .groups = "drop"
  )


country_participation_consistency <- country_participation %>%
  count(
    YearsReported,
    name = "Countries"
  ) %>%
  arrange(
    YearsReported
  )


cat(
  "\n--- Number of study years contributed by countries ---\n"
)


print(
  country_participation_consistency
)


write_csv(
  country_participation_consistency,
  here(
    "outputs",
    "country_participation_consistency.csv"
  )
)


# ----------------------------------------------------------
# 11. Figure 1:
#     Participating countries/territories by year
# ----------------------------------------------------------

plot_country_coverage <- ggplot(
  coverage_by_year,
  aes(
    x = Year,
    y = Countries
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = Countries
    ),
    vjust = -0.9,
    size = 4
  ) +
  scale_x_continuous(
    breaks = 2020:2023
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.05,
        0.15
      )
    )
  ) +
  labs(
    title = "Country participation in the analysed WHO GLASS data",
    subtitle = "Bloodstream-infection surveillance, 2020–2023",
    x = "Year",
    y = "Countries/territories reporting relevant data",
    caption = paste(
      "Increasing country participation means that the surveillance",
      "population is not identical across study years."
    )
  ) +
  theme_amr


print(
  plot_country_coverage
)


ggsave(
  here(
    "figures",
    "01_country_participation_over_time.png"
  ),
  plot = plot_country_coverage,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ----------------------------------------------------------
# 12. Figure 2:
#     Reported AST surveillance volume by year
# ----------------------------------------------------------
#
# Country participation alone does not describe surveillance
# volume.
#
# This figure shows the total number of interpretable AST
# results represented in the analytical dataset each year.
# ----------------------------------------------------------

coverage_by_year <- coverage_by_year %>%
  mutate(
    InterpretableASTMillions =
      InterpretableAST / 1000000
  )


plot_ast_volume <- ggplot(
  coverage_by_year,
  aes(
    x = Year,
    y = InterpretableASTMillions
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2fM",
        InterpretableASTMillions
      )
    ),
    vjust = -0.9,
    size = 4
  ) +
  scale_x_continuous(
    breaks = 2020:2023
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.16
      )
    )
  ) +
  labs(
    title = "Volume of reported antimicrobial susceptibility testing",
    subtitle = "Interpretable AST results represented in the analysed WHO GLASS data",
    x = "Year",
    y = "Interpretable AST results (millions)",
    caption = paste(
      "Testing volume and country participation do not change identically over time;",
      "both influence the composition of the surveillance dataset."
    )
  ) +
  theme_amr


print(
  plot_ast_volume
)


ggsave(
  here(
    "figures",
    "02_reported_AST_volume_over_time.png"
  ),
  plot = plot_ast_volume,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)


############################################################
# PART B — POOLED AST DESCRIPTIVE SUMMARIES
############################################################


# ----------------------------------------------------------
# 13. Pooled resistant AST results by year
# ----------------------------------------------------------
#
# IMPORTANT:
#
# This summary pools AST results across different pathogens,
# antibiotics, countries, and surveillance systems.
#
# Therefore it describes the proportion of reported AST
# results classified as resistant.
#
# It is NOT a population-representative estimate of:
#
#   - global AMR prevalence,
#   - the percentage of infections resistant,
#   - or the percentage of unique isolates resistant.
#
# It is sensitive to changes in surveillance volume and
# pathogen/antibiotic composition.
# ----------------------------------------------------------

pooled_ast_resistance_by_year <- amr_data %>%
  group_by(
    Year
  ) %>%
  summarise(
    
    ResistantAST = sum(
      Resistant,
      na.rm = TRUE
    ),
    
    InterpretableAST = sum(
      InterpretableAST,
      na.rm = TRUE
    ),
    
    PooledResistancePercentage =
      ResistantAST /
      InterpretableAST *
      100,
    
    .groups = "drop"
  )


print(
  pooled_ast_resistance_by_year
)


write_csv(
  pooled_ast_resistance_by_year,
  here(
    "outputs",
    "pooled_ast_resistance_by_year.csv"
  )
)


# ----------------------------------------------------------
# 14. Figure 3:
#     Pooled proportion of resistant AST results over time
# ----------------------------------------------------------

plot_pooled_ast_year <- ggplot(
  pooled_ast_resistance_by_year,
  aes(
    x = Year,
    y = PooledResistancePercentage
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        PooledResistancePercentage
      )
    ),
    vjust = -0.9,
    size = 4
  ) +
  scale_x_continuous(
    breaks = 2020:2023
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.08,
        0.18
      )
    )
  ) +
  labs(
    title = "Pooled proportion of reported AST results classified as resistant",
    subtitle = "Descriptive summary across included pathogens, antibiotics and countries",
    x = "Year",
    y = "Resistant AST results (%)",
    caption = paste(
      "Composition-dependent surveillance summary.",
      "Not a population-representative estimate of global AMR prevalence."
    )
  ) +
  theme_amr


print(
  plot_pooled_ast_year
)


ggsave(
  here(
    "figures",
    "03_pooled_resistant_AST_over_time.png"
  ),
  plot = plot_pooled_ast_year,
  width = 9,
  height = 6,
  dpi = 300,
  bg = "white"
)


# ----------------------------------------------------------
# 15. Pooled AST resistance by pathogen and year
# ----------------------------------------------------------
#
# These summaries still combine different antibiotics
# within each pathogen.
#
# They are therefore useful descriptively but remain
# sensitive to antibiotic-testing composition.
# ----------------------------------------------------------

pooled_ast_by_pathogen_year <- amr_data %>%
  group_by(
    PathogenName,
    Year
  ) %>%
  summarise(
    
    ResistantAST = sum(
      Resistant,
      na.rm = TRUE
    ),
    
    InterpretableAST = sum(
      InterpretableAST,
      na.rm = TRUE
    ),
    
    PooledResistancePercentage =
      ResistantAST /
      InterpretableAST *
      100,
    
    .groups = "drop"
  )


write_csv(
  pooled_ast_by_pathogen_year,
  here(
    "outputs",
    "pooled_ast_resistance_by_pathogen_year.csv"
  )
)


# ----------------------------------------------------------
# 16. Figure 4:
#     Pooled AST resistance by pathogen
# ----------------------------------------------------------

plot_pathogen_trends <- ggplot(
  pooled_ast_by_pathogen_year,
  aes(
    x = Year,
    y = PooledResistancePercentage,
    group = PathogenName,
    color = PathogenName
  )
) +
  geom_line(
    linewidth = 1
  ) +
  geom_point(
    size = 2.7
  ) +
  scale_x_continuous(
    breaks = 2020:2023
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.05,
        0.08
      )
    )
  ) +
  labs(
    title = "Pooled resistant AST results by pathogen",
    subtitle = "Descriptive summaries combining antibiotics within each pathogen",
    x = "Year",
    y = "Resistant AST results (%)",
    color = "Pathogen",
    caption = paste(
      "Interpret cautiously: temporal differences may partly reflect",
      "changes in antibiotic-testing composition."
    )
  ) +
  theme_amr


print(
  plot_pathogen_trends
)


ggsave(
  here(
    "figures",
    "04_pooled_AST_resistance_by_pathogen.png"
  ),
  plot = plot_pathogen_trends,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)


############################################################
# PART C — PATHOGEN × ANTIBIOTIC ANALYSIS
############################################################


# ----------------------------------------------------------
# 17. Pathogen-antibiotic coverage and resistance
# ----------------------------------------------------------

pathogen_antibiotic_summary <- amr_data %>%
  group_by(
    PathogenName,
    AntibioticName
  ) %>%
  summarise(
    
    Countries = n_distinct(
      Iso3
    ),
    
    Years = n_distinct(
      Year
    ),
    
    Records = n(),
    
    ResistantAST = sum(
      Resistant,
      na.rm = TRUE
    ),
    
    InterpretableAST = sum(
      InterpretableAST,
      na.rm = TRUE
    ),
    
    ResistancePercentage =
      ResistantAST /
      InterpretableAST *
      100,
    
    .groups = "drop"
  )


print(
  pathogen_antibiotic_summary
)


write_csv(
  pathogen_antibiotic_summary,
  here(
    "outputs",
    "pathogen_antibiotic_summary.csv"
  )
)


# ----------------------------------------------------------
# 18. Pathogen-antibiotic resistance in 2023
# ----------------------------------------------------------
#
# A single-year heatmap avoids blending temporal changes
# across the full four-year study period.
#
# Values are still pooled across contributing surveillance
# systems and therefore remain surveillance-volume weighted.
# ----------------------------------------------------------

pathogen_antibiotic_2023 <- amr_data %>%
  filter(
    Year == 2023
  ) %>%
  group_by(
    PathogenName,
    AntibioticName
  ) %>%
  summarise(
    
    Countries = n_distinct(
      Iso3
    ),
    
    ResistantAST = sum(
      Resistant,
      na.rm = TRUE
    ),
    
    InterpretableAST = sum(
      InterpretableAST,
      na.rm = TRUE
    ),
    
    ResistancePercentage =
      ResistantAST /
      InterpretableAST *
      100,
    
    .groups = "drop"
  )


# ----------------------------------------------------------
# 19. Construct complete heatmap grid
# ----------------------------------------------------------
#
# Explicitly including all pathogen-antibiotic combinations
# means combinations without observations appear as grey
# cells rather than simply disappearing from the plot.
# ----------------------------------------------------------

pathogen_order <- c(
  "Acinetobacter spp.",
  "Escherichia coli",
  "Klebsiella pneumoniae",
  "Salmonella spp.",
  "Streptococcus pneumoniae"
)


antibiotic_order <- sort(
  unique(
    amr_data$AntibioticName
  )
)


heatmap_2023 <- pathogen_antibiotic_2023 %>%
  complete(
    PathogenName = pathogen_order,
    AntibioticName = antibiotic_order
  ) %>%
  mutate(
    
    PathogenName = factor(
      PathogenName,
      levels = pathogen_order
    ),
    
    AntibioticName = factor(
      AntibioticName,
      levels = antibiotic_order
    )
  )


# ----------------------------------------------------------
# 20. Figure 5:
#     Pathogen × antibiotic heatmap, 2023
# ----------------------------------------------------------

plot_heatmap_2023 <- ggplot(
  heatmap_2023,
  aes(
    x = AntibioticName,
    y = PathogenName,
    fill = ResistancePercentage
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_viridis_c(
    limits = c(
      0,
      100
    ),
    na.value = "grey90"
  ) +
  scale_x_discrete(
    labels = function(x) {
      
      case_when(
        
        str_detect(
          x,
          regex(
            "third.?generation.*cephalospor",
            ignore_case = TRUE
          )
        ) ~ "3rd-gen\ncephalosporins",
        
        TRUE ~ str_wrap(
          x,
          width = 14
        )
      )
    }
  ) +
  labs(
    title = "Reported resistance by pathogen and antibiotic, 2023",
    subtitle = "Pooled across contributing WHO GLASS surveillance systems",
    x = "Antibiotic",
    y = "Pathogen",
    fill = "Resistance (%)",
    caption = paste(
      "Grey cells indicate pathogen–antibiotic combinations without",
      "observations in the analysed 2023 data."
    )
  ) +
  theme_amr +
  theme(
    
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    legend.position = "right"
  )


print(
  plot_heatmap_2023
)


ggsave(
  here(
    "figures",
    "05_pathogen_antibiotic_heatmap_2023.png"
  ),
  plot = plot_heatmap_2023,
  width = 14,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

# ----------------------------------------------------------
# 21. Descriptive change by pathogen-antibiotic combination
# ----------------------------------------------------------
#
# Comparing within the same pathogen-antibiotic combination
# avoids the major problem of pooling different antibiotics
# together.
#
# Country composition can still change between years.
#
# Therefore these values remain descriptive and will later
# be evaluated using country-aware statistical models and
# sensitivity analyses.
# ----------------------------------------------------------

pathogen_antibiotic_year <- amr_data %>%
  group_by(
    PathogenName,
    AntibioticName,
    Year
  ) %>%
  summarise(
    
    Countries = n_distinct(
      Iso3
    ),
    
    ResistantAST = sum(
      Resistant,
      na.rm = TRUE
    ),
    
    InterpretableAST = sum(
      InterpretableAST,
      na.rm = TRUE
    ),
    
    ResistancePercentage =
      ResistantAST /
      InterpretableAST *
      100,
    
    .groups = "drop"
  )


change_2020_2023 <- pathogen_antibiotic_year %>%
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
    ResistancePercentage
  ) %>%
  pivot_wider(
    names_from = Year,
    values_from = ResistancePercentage,
    names_prefix = "Resistance_"
  ) %>%
  mutate(
    AbsoluteChange =
      Resistance_2023 -
      Resistance_2020
  ) %>%
  arrange(
    desc(
      AbsoluteChange
    )
  )


write_csv(
  change_2020_2023,
  here(
    "outputs",
    "pathogen_antibiotic_change_2020_2023_descriptive.csv"
  )
)


############################################################
# PART D — FINAL SUMMARY
############################################################


cat(
  "\n========================================\n",
  "EXPLORATORY ANALYSIS COMPLETE\n",
  "========================================\n"
)


cat(
  "\nFigures saved to:\n",
  here("figures"),
  "\n"
)


cat(
  "\nSummary tables saved to:\n",
  here("outputs"),
  "\n"
)


cat(
  "\nIMPORTANT:\n",
  paste(
    "Cross-antibiotic pooled resistance summaries are",
    "descriptive AST-level surveillance summaries and",
    "must not be interpreted as population-representative",
    "global AMR prevalence."
  ),
  "\n"
)