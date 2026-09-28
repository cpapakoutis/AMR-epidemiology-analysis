Dataset Notes
Project
WHO GLASS Antimicrobial Resistance Analysis
This document records the source-data structure, analytical variables, parsing decisions, validation checks, and reproducibility requirements for the WHO GLASS bloodstream-infection dataset used in this project.
Data Source
Data were obtained from the World Health Organization Global Antimicrobial Resistance and Use Surveillance System (WHO GLASS).
Official GLASS initiative page:
https://www.who.int/initiatives/glass
The WHO GLASS initiative page provides access to the GLASS Data Visualization Dashboard, which presents downloadable antimicrobial resistance (AMR) surveillance data for countries, territories, and areas.
The current WHO dashboard states that it was last updated on 25 September 2025 and contains data for 2016–2023 submitted by the end of 2024.
This project uses country/territory-level bloodstream-infection surveillance exports for 2020–2023.
The original date on which the project files were downloaded was not preserved and is therefore not reconstructed from memory. The official WHO GLASS source and current dashboard were re-verified in September 2026.
Source Links
- WHO GLASS initiative: https://www.who.int/initiatives/glass
- WHO GLASS Data Visualization Dashboard: https://worldhealthorg.shinyapps.io/glass-dashboard/
Redistribution
The original WHO data files are not tracked in this Git repository.
Both of the following directories are excluded from version control:
data/raw/
data/processed/
This repository therefore contains the analytical code and documentation but does not redistribute the original WHO files.
Before redistributing any WHO-derived dataset, the terms that applied to the exact source/export used for this project should be checked and documented.
Study Scope
Years
2020
2021
2022
2023
Specimen / Infection Type
The project analyses bloodstream-infection surveillance data.
Included Pathogen Groups
- Acinetobacter spp.
- Escherichia coli
- Klebsiella pneumoniae
- Salmonella spp.
- Streptococcus pneumoniae
Pseudomonas aeruginosa was considered during project development but was excluded from the final analytical scope.
Required Raw Files
The cleaning pipeline in scripts/02_data_import_cleaning.R expects exactly 20 CSV files:
5 pathogen groups × 4 years = 20 files
The expected project structure is:
data/raw/
├── Acinetobacter/
│   ├── Acinetobacter_2020.csv
│   ├── Acinetobacter_2021.csv
│   ├── Acinetobacter_2022.csv
│   └── Acinetobacter_2023.csv
├── Escherichia_coli/
│   ├── Ecoli_2020.csv
│   ├── Ecoli_2021.csv
│   ├── Ecoli_2022.csv
│   └── Ecoli_2023.csv
├── Klebsiella_pneumoniae/
│   ├── Klebsiella_2020.csv
│   ├── Klebsiella_2021.csv
│   ├── Klebsiella_2022.csv
│   └── Klebsiella_2023.csv
├── Salmonella/
│   ├── Salmonella_2020.csv
│   ├── Salmonella_2021.csv
│   ├── Salmonella_2022.csv
│   └── Salmonella_2023.csv
└── Streptococcus_pneumoniae/
    ├── Spneumoniae_2020.csv
    ├── Spneumoniae_2021.csv
    ├── Spneumoniae_2022.csv
    └── Spneumoniae_2023.csv
The cleaning script searches recursively within data/raw/, so the pathogen subdirectories are useful for organisation but the filenames and total file count are the important reproducibility requirements.
Raw WHO File Structure
The source exports are not conventional single-table CSV files.
Each file contains multiple sections. The country-level observations used for this project occur after a line containing:
Data for boxplots
The following line is expected to contain this header:
Specimen,PathogenName,AntibioticName,Iso3,CountryTerritoryArea,WHORegionName,InterpretableAST,Resistant,ResistancePercentage
The analysis uses the rows following this header.
The cleaning script checks that the "Data for boxplots" marker occurs exactly once in every file and that the expected header structure is present. The script stops with an error if either condition fails.
Why a Custom Parser Is Required
A standard comma-based import cannot safely parse all records because some country/territory names themselves contain commas.
For example, a value such as:
occupied Palestinian territory, including east Jerusalem
contains a comma that is part of the country/territory name rather than a delimiter between analytical variables.
A naive comma split would therefore create too many fields and shift the remaining variables into the wrong columns.
The custom parser uses the known structure of each row:
- the first four fields are fixed:
  - Specimen
  - PathogenName
  - AntibioticName
  - Iso3
- the final four fields are fixed:
  - WHORegionName
  - InterpretableAST
  - Resistant
  - ResistancePercentage
- all fields between Iso3 and WHORegionName are reconstructed as CountryTerritoryArea
This preserves country and territory names containing embedded commas.
Analytical Variables
The cleaned dataset produced by the current cleaning script contains the following core variables.
Specimen
Specimen/infection category represented by the WHO export.
The project is restricted to bloodstream-infection surveillance observations.
PathogenName
Reported bacterial pathogen group.
The five pathogen groups included in the final project are listed above.
AntibioticName
Antibiotic for which susceptibility-test results were reported.
Analyses involving resistance should generally be interpreted at the pathogen–antibiotic level because resistance to different antibiotics is biologically and epidemiologically distinct.
Iso3
Three-character country/territory identifier supplied in the WHO data.
This variable is used as the country identifier in statistical models and clustering.
CountryTerritoryArea
Country or territory name.
This field may contain commas and is one of the reasons custom parsing is required.
WHORegionName
WHO regional classification supplied in the source data.
This variable is used for descriptive surveillance-coverage summaries.
InterpretableAST
Number of reported antimicrobial susceptibility tests with interpretable results for the corresponding country–year–pathogen–antibiotic observation.
Important: this variable should not automatically be interpreted as a count of unique patients or unique bacterial isolates.
The same isolate may be tested against several antibiotics. Therefore, summing InterpretableAST across antibiotics can count susceptibility-test results from the same isolate more than once.
Resistant
Number of reported AST results classified as resistant within the observation.
For statistical modelling, non-resistant results are derived as:
NonResistant = InterpretableAST - Resistant
ResistancePercentage
Resistance percentage reported in the WHO source file.
This is validated during cleaning against:
Resistant / InterpretableAST × 100
Year
Study year extracted from the source filename.
Allowed values are restricted to:
2020:2023
SourceFile
Name of the WHO source file from which each observation was parsed.
This preserves provenance within the combined analytical dataset.
Unit of Analysis
The analytical unit is a:
country/territory × year × pathogen × antibiotic surveillance observation

The dataset is aggregate surveillance data rather than patient-level or isolate-level data.
Accordingly, the analysis does not assume that rows correspond to independent individual patients.
Cleaning and Validation Checks
scripts/02_data_import_cleaning.R performs the following checks before saving the processed dataset.
File-Count Validation
Exactly 20 CSV files must be present under data/raw/.
Structural Validation
For every file:
- the "Data for boxplots" marker must occur exactly once;
- the expected nine-field header must match exactly;
- every data row must contain enough fields to reconstruct the expected variables.
Missing-Data Validation
The script checks for missing values in the key analytical variables:
- Specimen
- PathogenName
- AntibioticName
- Iso3
- CountryTerritoryArea
- WHORegionName
- InterpretableAST
- Resistant
- ResistancePercentage
- Year
The script stops if missing values occur in these critical variables.
Study-Year Validation
The combined dataset must contain exactly the four expected study years:
2020
2021
2022
2023
Logical Count Validation
Rows are considered invalid if:
InterpretableAST <= 0
Resistant < 0
Resistant > InterpretableAST
ResistancePercentage < 0
ResistancePercentage > 100
The cleaning script stops if any such observations are detected.
Resistance-Percentage Validation
For every row, the WHO-provided resistance percentage is compared with the value recalculated from the counts:
CalculatedResistance = Resistant / InterpretableAST × 100
The maximum absolute discrepancy must be no greater than 1e-6.
During the completed analysis, the agreement was effectively exact apart from floating-point precision.
Country-Name Validation
The cleaning workflow explicitly inspects country/territory names containing commas to confirm that they survived reconstruction correctly.
Duplicate-Key Inspection
Possible duplicate analytical keys are inspected using:
Year
Iso3
Specimen
PathogenName
AntibioticName
This check is diagnostic: potential duplicates are displayed for inspection rather than silently removed.
Final Analytical Dataset
The validated dataset used for the completed analysis contained:
- 9,098 observations
- 103 countries and territories
- 4 years
- 5 pathogen groups
Observations by year were:
Year	Observations
2020	2,097
2021	2,135
2022	2,343
2023	2,523


Observations by pathogen were:
Pathogen	Observations
Escherichia coli	3,294
Klebsiella pneumoniae	2,830
Acinetobacter spp.	1,389
Salmonella spp.	907
Streptococcus pneumoniae	678


The processed dataset is written to:
data/processed/AMR_GLASS_clean.csv
Because this directory is excluded by .gitignore, another user reproducing the project must obtain the source files and run the cleaning script locally.
Important Interpretation Rules
AST Counts Are Not Automatically Unique Isolate Counts
InterpretableAST records susceptibility-test results. Without isolate-level identifiers, the analysis cannot determine how many unique isolates or patients generated those tests.
Cross-Antibiotic Pooling Is Composition-Dependent
When AST counts are pooled across antibiotics, the resulting percentage depends on which antibiotics were tested and how frequently they were represented.
Therefore, pooled values across antibiotics are described as:
the proportion of reported AST results classified as resistant

They are not described as the percentage of infections or isolates that were resistant.
Country Contributions Are Unequal
Countries contribute very different numbers of AST results.
AST-weighted pooled estimates therefore give greater influence to high-volume surveillance systems.
For this reason, the project also includes equal-country-weighted descriptive sensitivity analyses.
Surveillance Participation Changes Over Time
The countries and territories contributing GLASS observations changed between 2020 and 2023.
Temporal comparisons therefore require attention to changing surveillance composition, which is addressed through country adjustment and stable-country sensitivity analyses.
The Data Are Not a Representative Global Sample
WHO GLASS surveillance coverage is not equivalent to probability sampling of the global population.
Results are therefore interpreted as patterns among contributing surveillance data, not as population-representative estimates of worldwide AMR prevalence.
Relationship to Project Scripts
The main scripts relevant to this document are:
- scripts/01_data_understanding.R
  Inspects the raw WHO file structure and identifies the parsing problem caused by comma-containing country names.
- scripts/02_data_import_cleaning.R
  Parses, combines, validates, and saves the analytical dataset.
- scripts/03_exploratory_analysis.R
  Uses the cleaned dataset for descriptive surveillance summaries.
- scripts/04_statistical_analysis.R onward
  Uses the cleaned dataset for inferential and sensitivity analyses.
Reproducibility Note
The raw and processed datasets are intentionally excluded from Git version control.
A reproducible run therefore requires:
1. obtaining the same WHO GLASS source exports;
2. arranging the 20 files under data/raw/;
3. running scripts/00_setup_packages.R;
4. running scripts/01_data_understanding.R as an inspection step;
5. running scripts/02_data_import_cleaning.R to regenerate data/processed/AMR_GLASS_clean.csv;
6. continuing with scripts 03 through 07.
All scripts use project-relative paths through the here package.
The repository also includes an renv.lock file to record the R package environment used during the final reproducibility run.
Record-Status Note
This document describes the final analytical workflow reflected in the current R scripts.
The exact original download date of the source CSVs was not preserved and is therefore not reconstructed from memory.
The official WHO GLASS initiative page and current dashboard were re-verified in September 2026. If the source data are re-downloaded or the analysis is prepared for formal publication, the exact source URL, access date, version/export identifier, and applicable data-use terms for the newly obtained files should be recorded before treating them as equivalent to the original project files.
