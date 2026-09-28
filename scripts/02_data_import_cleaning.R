############################################################
# Project: WHO GLASS Antimicrobial Resistance Analysis
# Script: 02_data_import_cleaning.R
#
# Purpose:
#   Import, parse, clean, combine, and validate the raw
#   WHO GLASS bloodstream-infection surveillance files.
#
# Input:
#   Raw WHO GLASS CSV exports stored under:
#   data/raw/
#
# Output:
#   data/processed/AMR_GLASS_clean.csv
#
# Main steps:
#   1. Locate all raw CSV files.
#   2. Locate the "Data for boxplots" section in each file.
#   3. Parse country-level pathogen-antibiotic observations.
#   4. Correctly handle country names containing commas.
#   5. Add year and source-file metadata.
#   6. Combine all files.
#   7. Run data-quality and consistency checks.
#   8. Save the validated analytical dataset.
############################################################


# ----------------------------------------------------------
# 1. Load packages
# ----------------------------------------------------------

library(tidyverse)
library(here)


# ----------------------------------------------------------
# 2. Function: parse one WHO GLASS data row
# ----------------------------------------------------------

parse_amr_row <- function(x) {
  
  # Split the raw row at commas.
  #
  # Most rows contain 9 fields.
  # Some contain additional commas inside the country/
  # territory name, so they produce more than 9 pieces.
  
  pieces <- str_split(
    x,
    fixed(",")
  )[[1]]
  
  # Remove quotation marks surrounding fields and
  # standardise leading/trailing whitespace.
  
  pieces <- pieces %>%
    str_remove('^"') %>%
    str_remove('"$') %>%
    str_trim()
  
  
  # A valid row must contain at least the expected
  # nine fields.
  
  if (length(pieces) < 9) {
    
    stop(
      "Malformed WHO GLASS row: expected at least 9 fields.\n",
      "Row:\n",
      x
    )
  }
  
  
  # --------------------------------------------------------
  # Fixed fields at the beginning of the row
  # --------------------------------------------------------
  
  Specimen <- pieces[1]
  
  PathogenName <- pieces[2]
  
  AntibioticName <- pieces[3]
  
  Iso3 <- pieces[4]
  
  
  # --------------------------------------------------------
  # Country / territory field
  # --------------------------------------------------------
  
  # The final four fields have fixed meanings:
  #
  # WHORegionName
  # InterpretableAST
  # Resistant
  # ResistancePercentage
  #
  # Therefore everything between Iso3 and those final
  # four fields belongs to CountryTerritoryArea.
  #
  # This allows country names containing commas to be
  # reconstructed correctly.
  
  country_end <- length(pieces) - 4
  
  CountryTerritoryArea <- str_c(
    pieces[5:country_end],
    collapse = ", "
  )
  
  
  # --------------------------------------------------------
  # Fixed fields at the end of the row
  # --------------------------------------------------------
  
  WHORegionName <- pieces[length(pieces) - 3]
  
  InterpretableAST <- parse_double(
    pieces[length(pieces) - 2]
  )
  
  Resistant <- parse_double(
    pieces[length(pieces) - 1]
  )
  
  ResistancePercentage <- parse_double(
    pieces[length(pieces)]
  )
  
  
  # --------------------------------------------------------
  # Return one parsed observation
  # --------------------------------------------------------
  
  tibble(
    Specimen,
    PathogenName,
    AntibioticName,
    Iso3,
    CountryTerritoryArea,
    WHORegionName,
    InterpretableAST,
    Resistant,
    ResistancePercentage
  )
}


# ----------------------------------------------------------
# 3. Function: clean one WHO GLASS file
# ----------------------------------------------------------

clean_glass_file <- function(file_path, year) {
  
  # Read the WHO export line by line rather than assuming
  # that the entire file forms one conventional CSV table.
  
  raw_lines <- read_lines(file_path)
  
  
  # Create a normalised copy for finding structural markers.
  # The original lines are retained for actual parsing.
  
  normalized_lines <- raw_lines %>%
    str_remove('^"') %>%
    str_remove('"$') %>%
    str_trim()
  
  
  # --------------------------------------------------------
  # Locate the country-level section
  # --------------------------------------------------------
  
  marker_position <- which(
    normalized_lines == "Data for boxplots"
  )
  
  
  # The marker should occur exactly once.
  
  if (length(marker_position) != 1) {
    
    stop(
      "Expected exactly one 'Data for boxplots' marker in:\n",
      file_path,
      "\nFound: ",
      length(marker_position)
    )
  }
  
  
  # --------------------------------------------------------
  # Validate the expected header
  # --------------------------------------------------------
  
  header_line <- raw_lines[
    marker_position + 1
  ]
  
  header_fields <- str_split(
    header_line,
    fixed(",")
  )[[1]] %>%
    str_remove('^"') %>%
    str_remove('"$') %>%
    str_trim()
  
  
  expected_header <- c(
    "Specimen",
    "PathogenName",
    "AntibioticName",
    "Iso3",
    "CountryTerritoryArea",
    "WHORegionName",
    "InterpretableAST",
    "Resistant",
    "ResistancePercentage"
  )
  
  
  if (!identical(
    header_fields,
    expected_header
  )) {
    
    stop(
      "Unexpected WHO GLASS header structure in:\n",
      file_path
    )
  }
  
  
  # --------------------------------------------------------
  # Extract actual observations
  # --------------------------------------------------------
  
  data_start <- marker_position + 2
  
  boxplot_raw <- raw_lines[
    data_start:length(raw_lines)
  ]
  
  
  # Remove completely empty lines if present.
  
  boxplot_raw <- boxplot_raw[
    str_trim(boxplot_raw) != ""
  ]
  
  
  # --------------------------------------------------------
  # Parse all observations
  # --------------------------------------------------------
  
  clean_data <- map_dfr(
    boxplot_raw,
    parse_amr_row
  )
  
  
  # --------------------------------------------------------
  # Add metadata
  # --------------------------------------------------------
  
  clean_data <- clean_data %>%
    mutate(
      Year = as.integer(year),
      SourceFile = basename(file_path)
    )
  
  
  clean_data
}


# ----------------------------------------------------------
# 4. Locate all raw WHO GLASS files
# ----------------------------------------------------------

files <- list.files(
  path = here("data", "raw"),
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
) %>%
  sort()


cat(
  "\n--- Raw WHO GLASS files ---\n"
)

cat(
  "Number of CSV files found:",
  length(files),
  "\n"
)


# This project should contain:
#
# 5 pathogens × 4 years = 20 files.

if (length(files) != 20) {
  
  stop(
    "Expected 20 raw CSV files but found ",
    length(files),
    "."
  )
}


# ----------------------------------------------------------
# 5. Parse and combine all files
# ----------------------------------------------------------

all_data <- map_dfr(
  files,
  function(file) {
    
    # Extract the four-digit year from the filename.
    
    year <- str_extract(
      basename(file),
      "20[0-9]{2}"
    )
    
    
    # Fail if no valid year was found.
    
    if (
      is.na(year) ||
      !as.integer(year) %in% 2020:2023
    ) {
      
      stop(
        "Could not identify a valid study year from file:\n",
        file
      )
    }
    
    
    clean_glass_file(
      file_path = file,
      year = year
    )
  }
)


# Sort rows into a stable and interpretable order.

all_data <- all_data %>%
  arrange(
    Year,
    PathogenName,
    Iso3,
    AntibioticName
  )


############################################################
# QUALITY CONTROL
############################################################


# ----------------------------------------------------------
# 6. Dataset dimensions
# ----------------------------------------------------------

cat(
  "\n--- Dataset dimensions ---\n"
)

cat(
  "Rows:",
  nrow(all_data),
  "\n"
)

cat(
  "Columns:",
  ncol(all_data),
  "\n"
)


# ----------------------------------------------------------
# 7. Dataset structure
# ----------------------------------------------------------

cat(
  "\n--- Dataset structure ---\n"
)

str(all_data)


# ----------------------------------------------------------
# 8. Missing values
# ----------------------------------------------------------

cat(
  "\n--- Missing values ---\n"
)

missing_values <- all_data %>%
  summarise(
    across(
      everything(),
      ~ sum(is.na(.x))
    )
  ) %>%
  pivot_longer(
    everything(),
    names_to = "Variable",
    values_to = "Missing"
  )

print(missing_values)


# Key analytical variables should not be missing.

critical_variables <- c(
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

critical_missing <- missing_values %>%
  filter(
    Variable %in% critical_variables,
    Missing > 0
  )


if (nrow(critical_missing) > 0) {
  
  stop(
    "Missing values detected in critical analytical variables."
  )
}


# ----------------------------------------------------------
# 9. Study years
# ----------------------------------------------------------

cat(
  "\n--- Years ---\n"
)

print(
  sort(
    unique(all_data$Year)
  )
)


if (!identical(
  sort(unique(all_data$Year)),
  2020:2023
)) {
  
  stop(
    "Dataset does not contain exactly the expected years 2020-2023."
  )
}


# ----------------------------------------------------------
# 10. Pathogens
# ----------------------------------------------------------

cat(
  "\n--- Pathogens ---\n"
)

print(
  sort(
    unique(all_data$PathogenName)
  )
)


# ----------------------------------------------------------
# 11. Infection types
# ----------------------------------------------------------

cat(
  "\n--- Infection types ---\n"
)

print(
  unique(all_data$Specimen)
)


# ----------------------------------------------------------
# 12. Observations by year
# ----------------------------------------------------------

cat(
  "\n--- Observations per year ---\n"
)

print(
  table(all_data$Year)
)


# ----------------------------------------------------------
# 13. Observations by pathogen
# ----------------------------------------------------------

cat(
  "\n--- Observations per pathogen ---\n"
)

print(
  table(all_data$PathogenName)
)


# ----------------------------------------------------------
# 14. Validate biologically/logically possible values
# ----------------------------------------------------------

invalid_values <- all_data %>%
  filter(
    InterpretableAST <= 0 |
      Resistant < 0 |
      Resistant > InterpretableAST |
      ResistancePercentage < 0 |
      ResistancePercentage > 100
  )


cat(
  "\n--- Invalid resistance/count values ---\n"
)

cat(
  "Rows with impossible values:",
  nrow(invalid_values),
  "\n"
)


if (nrow(invalid_values) > 0) {
  
  print(invalid_values)
  
  stop(
    "Impossible resistance/count values detected."
  )
}


# ----------------------------------------------------------
# 15. Validate WHO resistance percentages
# ----------------------------------------------------------

percentage_check <- all_data %>%
  mutate(
    CalculatedResistance =
      Resistant / InterpretableAST * 100,
    
    PercentageDifference =
      abs(
        ResistancePercentage -
          CalculatedResistance
      )
  )


max_percentage_difference <- max(
  percentage_check$PercentageDifference,
  na.rm = TRUE
)


cat(
  "\n--- WHO percentage validation ---\n"
)

cat(
  "Maximum absolute difference:",
  max_percentage_difference,
  "\n"
)


# Small numerical differences can occur because of
# floating-point precision.

if (max_percentage_difference > 1e-6) {
  
  stop(
    "WHO resistance percentages do not match values calculated from counts."
  )
}


# ----------------------------------------------------------
# 16. Confirm that comma-containing country names survived
# ----------------------------------------------------------

countries_with_commas <- all_data %>%
  filter(
    str_detect(
      CountryTerritoryArea,
      fixed(",")
    )
  ) %>%
  distinct(
    CountryTerritoryArea
  )


cat(
  "\n--- Country/territory names containing commas ---\n"
)

print(
  countries_with_commas
)


# ----------------------------------------------------------
# 17. Examine possible duplicate analytical keys
# ----------------------------------------------------------

duplicate_keys <- all_data %>%
  count(
    Year,
    Iso3,
    Specimen,
    PathogenName,
    AntibioticName,
    name = "n"
  ) %>%
  filter(
    n > 1
  )


cat(
  "\n--- Potential duplicate analytical keys ---\n"
)

cat(
  "Number of duplicated keys:",
  nrow(duplicate_keys),
  "\n"
)

if (nrow(duplicate_keys) > 0) {
  
  print(
    duplicate_keys,
    n = 20
  )
}


# ----------------------------------------------------------
# 18. Inspect first observations
# ----------------------------------------------------------

cat(
  "\n--- First observations ---\n"
)

print(
  all_data %>%
    select(
      Specimen,
      PathogenName,
      AntibioticName,
      Iso3,
      CountryTerritoryArea,
      WHORegionName,
      InterpretableAST,
      Resistant,
      ResistancePercentage,
      Year,
      SourceFile
    ) %>%
    head()
)


# ----------------------------------------------------------
# 19. Save processed analytical dataset
# ----------------------------------------------------------

processed_dir <- here(
  "data",
  "processed"
)

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


write_csv(
  all_data,
  here(
    "data",
    "processed",
    "AMR_GLASS_clean.csv"
  )
)


# ----------------------------------------------------------
# 20. Final message
# ----------------------------------------------------------

cat(
  "\n========================================\n",
  "Cleaning and validation complete\n",
  "Files processed: ", length(files), "\n",
  "Rows: ", nrow(all_data), "\n",
  "Columns: ", ncol(all_data), "\n",
  "========================================\n"
)