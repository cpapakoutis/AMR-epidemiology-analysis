Dataset Notes
Project
WHO GLASS Antimicrobial Resistance Analysis
This document records the source-data structure, analytical variables, parsing decisions, validation checks, and reproducibility requirements for the WHO GLASS bloodstream-infection dataset used in this project.
Data source
Data were obtained from the World Health Organization Global Antimicrobial Resistance and Use Surveillance System (WHO GLASS).
The analysis uses country/territory-level bloodstream-infection surveillance exports for 2020-2023.
The exact WHO download page, date of access, and dataset/version identifier were not preserved in the current project files. These details should be recorded here before public release or publication rather than reconstructed from memory.
Redistribution
The original WHO data files are not tracked in this Git repository.
Both of the following directories are excluded from version control:
data/raw/
data/processed/
This repository therefore contains the analytical code and documentation but does not redistribute the original WHO files.
Before redistributing any WHO-derived dataset, the terms that applied to the exact source/export used for this project should be checked and documented.
Study scope
Years
2020
2021
2022
2023
Specimen / infection type
The project analyses bloodstream-infection surveillance data.
Included pathogen groups
- Acinetobacter spp.
- Escherichia coli
- Klebsiella pneumoniae
- Salmonella spp.
- Streptococcus pneumoniae
Pseudomonas aeruginosa was considered during project development but was excluded from the final analytical scope.
Required raw files
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
Raw WHO file structure
The source exports are not conventional single-table CSV files.
Each file contains multiple sections. The country-level observations used for this project occur after a line containing:
Data for boxplots
The following line is expected to contain this header:
Specimen,PathogenName,AntibioticName,Iso3,CountryTerritoryArea,WHORegionName,InterpretableAST,Resistant,ResistancePercentage
The analysis uses the rows following this header.
The cleaning script checks that the "Data for boxplots" marker occurs exactly once in every file and that the expected header structure is present. The script stops with an error if either condition fails.
Why a custom parser is required
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
Analytical variables
The cleaned dataset produced by the current cleaning script contains the following core variables.
Specimen
Specimen/infection category represented by the WHO export.
The project is restricted to bloodstream-infection surveillance observations.
PathogenName
Reported bacterial pathogen group.
The five pathogen groups included in the final project are listed above.
AntibioticName
Antibiotic for which susceptibility-test results were reported.
Analyses involving resistance should generally be interpreted at the pathogen-antibiotic level because resistance to different antibiotics is biologically and epidemiologically distinct.
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
Number of reported antimicrobial susceptibility tests with interpretable results for the corresponding country-year-pathogen-antibiotic observation.
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
Unit of analysis
The analytical unit is a:
country/territory × year × pathogen × antibiotic surveillance observation

The dataset is aggregate surveillance data rather than patient-level or isolate-level data.
Accordingly, the analysis does not assume that rows correspond to independent individual patients.
Cleaning and validation checks
scripts/02_data_import_cleaning.R performs the following checks before saving the processed dataset.
File-count validation
Exactly 20 CSV files must be present under data/raw/.
Structural validation
For every file:
- the "Data for boxplots" marker must occur exactly once;
- the expected nine-field header must match exactly;
- every data row must contain enough fields to reconstruct the expected variables.
Missing-data validation
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
Study-year validation
The combined dataset must contain exactly the four expected study years:
2020
2021
2022
2023
Logical count validation
Rows are considered invalid if:
InterpretableAST <= 0
Resistant < 0
Resistant > InterpretableAST
ResistancePercentage < 0
ResistancePercentage > 100
The cleaning script stops if any such observations are detected.
Resistance-percentage validation
For every row, the WHO-provided resistance percentage is compared with the value recalculated from the counts:
CalculatedResistance = Resistant / InterpretableAST × 100
The maximum absolute discrepancy must be no greater than 1e-6.
During the completed analysis, the agreement was effectively exact apart from floating-point precision.
Country-name validation
The cleaning workflow explicitly inspects country/territory names containing commas to confirm that they survived reconstruction correctly.
Duplicate-key inspection
Possible duplicate analytical keys are inspected using:
Year
Iso3
Specimen
PathogenName
AntibioticName
This check is diagnostic: potential duplicates are displayed for inspection rather than silently removed.
Final analytical dataset
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
Important interpretation rules
AST counts are not automatically unique isolate counts
InterpretableAST records susceptibility-test results. Without isolate-level identifiers, the analysis cannot determine how many unique isolates or patients generated those tests.
Cross-antibiotic pooling is composition-dependent
When AST counts are pooled across antibiotics, the resulting percentage depends on which antibiotics were tested and how frequently they were represented.
Therefore, pooled values across antibiotics are described as:
the proportion of reported AST results classified as resistant

They are not described as the percentage of infections or isolates that were resistant.
Country contributions are unequal
Countries contribute very different numbers of AST results.
AST-weighted pooled estimates therefore give greater influence to high-volume surveillance systems.
For this reason, the project also includes equal-country-weighted descriptive sensitivity analyses.
Surveillance participation changes over time
The countries and territories contributing GLASS observations changed between 2020 and 2023.
Temporal comparisons therefore require attention to changing surveillance composition, which is addressed through country adjustment and stable-country sensitivity analyses.
The data are not a representative global sample
WHO GLASS surveillance coverage is not equivalent to probability sampling of the global population.
Results are therefore interpreted as patterns among contributing surveillance data, not as population-representative estimates of worldwide AMR prevalence.
Relationship to project scripts
The main scripts relevant to this document are:
- scripts/01_data_understanding.R
  Inspects the raw WHO file structure and identifies the parsing problem caused by comma-containing country names.
- scripts/02_data_import_cleaning.R
  Parses, combines, validates, and saves the analytical dataset.
- scripts/03_exploratory_analysis.R
  Uses the cleaned dataset for descriptive surveillance summaries.
- scripts/04_statistical_analysis.R onward
  Uses the cleaned dataset for inferential and sensitivity analyses.
Reproducibility note
The raw and processed datasets are intentionally excluded from Git version control.
A reproducible run therefore requires:
1. obtaining the same WHO GLASS source exports;
2. arranging the 20 files under data/raw/;
3. running scripts/00_setup_packages.R;
4. running scripts/01_data_understanding.R as an inspection step;
5. running scripts/02_data_import_cleaning.R to regenerate data/processed/AMR_GLASS_clean.csv;
6. continuing with scripts 03 through 07.
All scripts use project-relative paths through the here package.
Record-status note
This document describes the final analytical workflow reflected in the current R scripts.
If the WHO source files are re-downloaded in the future, the exact source URL, access date, version/export identifier, and applicable data-use terms should be added to this document before treating the new files as equivalent to the originals.