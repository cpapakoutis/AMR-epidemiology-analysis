Temporal and Pathogen-Antibiotic Patterns in Antimicrobial Resistance Reported to WHO GLASS, 2020-2023
An epidemiological analysis of WHO GLASS bloodstream-infection surveillance data using R.
Overview
Antimicrobial resistance (AMR) is often summarized using broad global percentages, but surveillance data are shaped by which countries report, how much testing they contribute, which pathogens and antibiotics are represented, and how those patterns change over time.
This project examines reported AMR patterns in WHO GLASS bloodstream-infection surveillance data from 2020-2023, with particular attention to:
- changing surveillance participation;
- pathogen- and antibiotic-specific resistance patterns;
- repeated observations within countries;
- unequal AST reporting volumes;
- overdispersion;
- multiple testing;
- robustness to alternative modelling assumptions.
The goal is not to estimate a single population-representative measure of global AMR prevalence. Instead, the analysis asks:
How did reported antimicrobial resistance patterns change between 2020 and 2023 among countries and territories contributing bloodstream-infection surveillance data to WHO GLASS?

Dataset
The final cleaned analytical dataset contains:
- 9,098 country-level pathogen-antibiotic surveillance observations;
- 103 countries and territories;
- 4 years: 2020-2023;
- 5 pathogen groups:
  - Acinetobacter spp.
  - Escherichia coli
  - Klebsiella pneumoniae
  - Salmonella spp.
  - Streptococcus pneumoniae
Key variables include year, country/territory, WHO region, pathogen, antibiotic, number of interpretable antimicrobial susceptibility tests (AST), and number of resistant AST results.
InterpretableAST is treated as a count of reported susceptibility-test results and is not assumed to represent unique patients or unique bacterial isolates.
Data Source and Redistribution
The analysis uses bloodstream-infection surveillance exports from the World Health Organization Global Antimicrobial Resistance and Use Surveillance System (WHO GLASS).
The original WHO source files and the processed analytical dataset are intentionally excluded from version control:
data/raw/
data/processed/
The repository therefore contains the analytical code, documentation, model workflow, and final figures rather than redistributing the source data.
The exact WHO source page, date of access, export/version identifier, and data-use terms that applied to the files used in this project should be recorded before public release or publication. These details are flagged in [`docs/dataset_notes.md`](docs/dataset_notes.md) rather than reconstructed from memory.
The repository LICENSE applies to the project code and documentation and should not be interpreted as relicensing the underlying WHO data.
Analytical Workflow
The project is organised as a reproducible R pipeline:
scripts/
├── 00_setup_packages.R
├── 01_data_understanding.R
├── 02_data_import_cleaning.R
├── 03_exploratory_analysis.R
├── 04_statistical_analysis.R
├── 05_mixed_effects_analysis.R
├── 06_sensitivity_analysis.R
└── 07_final_figures.R
The workflow covers:
1. data inspection and parsing;
2. cleaning and validation;
3. descriptive epidemiology;
4. primary country-adjusted statistical modelling;
5. mixed-effects robustness analysis;
6. sensitivity analyses;
7. final figure generation.
Descriptive Findings
Country participation increased from 81 countries/territories in 2020 to 95 in 2023, while annual AST volume did not increase monotonically.
Across all included observations, the pooled proportion of reported AST results classified as resistant was:
- 18.4% in 2020
- 21.0% in 2021
- 20.5% in 2022
- 23.2% in 2023
These pooled percentages are composition-dependent AST-level surveillance summaries, not estimates of global AMR prevalence.
Pathogen-specific and pathogen-antibiotic analyses showed substantial heterogeneity, reinforcing the limitations of interpreting a single pooled resistance percentage.
Primary Inferential Analysis
Forty-eight pathogen-antibiotic combinations were present in the cleaned data.
For the primary analysis, combinations were required to:
- appear in all four study years;
- include both 2020 and 2023;
- have at least 20 contributing countries in every year.
This produced 37 predefined pathogen-antibiotic combinations.
Primary models used:
- grouped-binomial logistic regression;
- country fixed effects;
- country-clustered robust standard errors;
- continuous year centred at 2020;
- Benjamini-Hochberg false discovery rate correction across the 37 primary tests.
Main Findings
Two pathogen-antibiotic combinations met the predefined BH-FDR < 0.05 criterion in the primary analysis.
Escherichia coli - Meropenem
- OR per year: 1.13
- 95% CI: 1.07-1.19
- BH-FDR: 0.0022
- pooled resistance: approximately 0.47% in 2020 → 1.77% in 2023
This was the most consistently supported temporal increase across the primary and sensitivity analyses.
Klebsiella pneumoniae - Levofloxacin
- OR per year: 1.09
- 95% CI: 1.03-1.15
- BH-FDR: 0.0395
- pooled resistance: approximately 12.3% in 2020 → 17.5% in 2023
The increasing direction and approximate effect magnitude were consistent across sensitivity analyses, although multiplicity-adjusted evidence was weaker in the stable-country analysis.
Robustness and Sensitivity Analyses
The project deliberately compared multiple reasonable analytical approaches rather than selecting whichever model produced the smallest p-value.
Sensitivity analyses included:
- ordinary binomial random-intercept GLMMs;
- beta-binomial mixed models;
- restriction to countries reporting a given combination in all four years;
- categorical-year models comparing 2023 directly with 2020;
- equal-country-weighted descriptive summaries.
The ordinary binomial GLMMs produced similar point estimates but much narrower standard errors and widespread statistical significance in the presence of substantial residual overdispersion. Their nominal p-values were therefore not treated as the primary inferential evidence.
Across the main sensitivity analyses, the two primary signals retained the same increasing direction.
Selected Figures
Surveillance Participation
 
Pooled AST-Level Resistance
 
Pathogen-Antibiotic Resistance in 2023
 
Temporal Patterns in the Two Primary Signals
 
Sensitivity Analyses of the Two Primary Signals
 
The complete forest plot of all 37 primary pathogen-antibiotic estimates is retained as Supplementary Figure S1:
figures/S1_primary_temporal_forest_plot.png
Interpretation
The project distinguishes between:
Descriptive evidence
Patterns visible in the reported surveillance data.
Inferential evidence
Statistical support for temporal associations after accounting for repeated country observations, uncertainty, multiplicity, and alternative analytical assumptions.
Observed associations are not interpreted as causal effects.
The contrast between crude pooled resistance and adjusted model results is itself an important finding: broad pooled AMR summaries can be strongly influenced by surveillance composition, testing volume, pathogen mix, and antibiotic mix.
Limitations
Important limitations include:
- changing country participation over time;
- non-random surveillance coverage;
- unequal AST reporting volumes;
- variation in laboratory and surveillance capacity;
- differences in pathogen and antibiotic composition;
- repeated observations within countries;
- residual heterogeneity and overdispersion;
- inability to identify unique patients or isolates from aggregate AST counts;
- possible repeated contribution of the same isolate across antibiotics;
- short four-year follow-up;
- observational surveillance design.
Country adjustment and sensitivity analyses reduce some of these concerns but cannot eliminate them.
Repository Structure
AMR-epidemiology-analysis/
├── data/
│   ├── raw/                  # original WHO files; not tracked
│   └── processed/            # cleaned analytical data; not tracked
├── docs/
│   ├── analysis_plan.md
│   ├── dataset_notes.md
│   └── figure_descriptions.md
├── figures/                  # final and supplementary figures
├── outputs/                  # intermediate tables/model outputs; not tracked
├── reports/
│   ├── methods.md
│   ├── results.md
│   └── discussion.md
├── scripts/
│   ├── 00_setup_packages.R
│   ├── 01_data_understanding.R
│   ├── 02_data_import_cleaning.R
│   ├── 03_exploratory_analysis.R
│   ├── 04_statistical_analysis.R
│   ├── 05_mixed_effects_analysis.R
│   ├── 06_sensitivity_analysis.R
│   └── 07_final_figures.R
├── .gitignore
├── LICENSE
├── AMR-epidemiology-analysis.Rproj
└── README.md
Reproducibility
Scripts are designed to be run in numerical order.
The cleaning script expects the 20 WHO GLASS source CSVs documented in [`docs/dataset_notes.md`](docs/dataset_notes.md). After those files are placed under data/raw/, the workflow can recreate the cleaned analytical dataset and subsequent outputs.
The project uses relative paths through the here package to avoid machine-specific file paths.
A final local 00 → 07 reproducibility run should be completed after repository cleanup. Package versions can then be frozen with renv so that the tested software environment is recorded.
Documentation
Detailed project records are available in:
- [`docs/analysis_plan.md`](docs/analysis_plan.md) - final analytical strategy and decision record
- [`docs/dataset_notes.md`](docs/dataset_notes.md) - source structure, parsing decisions, variables, validation, and reproducibility requirements
- [`docs/figure_descriptions.md`](docs/figure_descriptions.md) - calculation, interpretation, and limitations for each figure
- [`reports/methods.md`](reports/methods.md) - manuscript-style methods
- [`reports/results.md`](reports/results.md) - manuscript-style results
- [`reports/discussion.md`](reports/discussion.md) - interpretation, strengths, limitations, implications, and conclusion
Tools
Analysis was conducted in R.
Key packages include:
- tidyverse
- here
- sandwich
- lmtest
- lme4
- glmmTMB
Project Status
Core data cleaning, descriptive analysis, primary modelling, sensitivity analyses, final figure generation, and manuscript-style Methods, Results, and Discussion documentation are complete.
Remaining repository tasks are the final local reproducibility run, package-version locking with renv, and GitHub/portfolio release preparation.