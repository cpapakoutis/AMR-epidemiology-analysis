############################################################
# Project: WHO GLASS Antimicrobial Resistance Analysis
# Script: 01_data_understanding.R
#
# Purpose:
#   Inspect the structure of a representative raw WHO GLASS
#   file before building the data-cleaning pipeline.
#
# Questions addressed:
#   1. How is the raw WHO export structured?
#   2. Where does the country-level "Data for boxplots"
#      section begin?
#   3. What variables are present in that section?
#   4. Can the rows be safely parsed by simply splitting
#      every comma?
#
# Input:
#   One representative raw WHO GLASS CSV file:
#   Acinetobacter spp., 2020.
#
# Output:
#   No files are written.
#   This script is exploratory and informs the cleaning
#   strategy implemented in 02_data_import_cleaning.R.
############################################################


# ----------------------------------------------------------
# 1. Load packages
# ----------------------------------------------------------

library(here)
library(readr)
library(dplyr)
library(stringr)


# ----------------------------------------------------------
# 2. Define representative raw file
# ----------------------------------------------------------

file_path <- here(
  "data",
  "raw",
  "Acinetobacter",
  "Acinetobacter_2020.csv"
)

# Stop immediately if the expected file is not available.
stopifnot(file.exists(file_path))


# ----------------------------------------------------------
# 3. Inspect the raw file as text
# ----------------------------------------------------------

# The WHO export contains several sections rather than a
# single conventional rectangular CSV table.
#
# Reading the file line-by-line allows us to inspect its
# original structure before making assumptions about columns.

raw_lines <- read_lines(file_path)

# Number of lines in the file
length(raw_lines)

# Inspect the first 20 rows/lines
head(raw_lines, 20)


# ----------------------------------------------------------
# 4. Locate the country-level data section
# ----------------------------------------------------------

# The country-level observations used in this project occur
# after the marker "Data for boxplots".

normalized_lines <- raw_lines %>%
  str_remove_all('"') %>%
  str_trim()

boxplot_start <- which(
  normalized_lines == "Data for boxplots"
)

boxplot_start


# Check that the marker occurs exactly once.
if (length(boxplot_start) != 1) {
  stop(
    "Expected exactly one 'Data for boxplots' marker, but found ",
    length(boxplot_start),
    "."
  )
}


# ----------------------------------------------------------
# 5. Inspect the header and first observations
# ----------------------------------------------------------

boxplot_header <- raw_lines[boxplot_start + 1]

boxplot_header


# Extract the observations following the header.
boxplot_lines <- raw_lines[
  (boxplot_start + 2):length(raw_lines)
]

# Remove completely empty lines, if present.
boxplot_lines <- boxplot_lines[
  str_trim(boxplot_lines) != ""
]

# Inspect the first 10 data rows.
head(boxplot_lines, 10)


# ----------------------------------------------------------
# 6. Examine the expected column structure
# ----------------------------------------------------------

# The expected variables are:
#
# 1. Specimen
# 2. PathogenName
# 3. AntibioticName
# 4. Iso3
# 5. CountryTerritoryArea
# 6. WHORegionName
# 7. InterpretableAST
# 8. Resistant
# 9. ResistancePercentage
#
# A conventional row containing nine fields should therefore
# contain eight comma delimiters.

expected_commas <- 8

comma_counts <- str_count(
  boxplot_lines,
  fixed(",")
)

table(comma_counts)


# ----------------------------------------------------------
# 7. Identify rows that do not have the expected number
#    of comma delimiters
# ----------------------------------------------------------

unexpected_rows <- tibble(
  line_number = seq_along(boxplot_lines),
  raw_line = boxplot_lines,
  comma_count = comma_counts
) %>%
  filter(comma_count != expected_commas)

unexpected_rows


# ----------------------------------------------------------
# 8. Investigate why these rows differ
# ----------------------------------------------------------

# Some country/territory names themselves contain commas.
#
# Therefore a naive operation such as:
#
# separate(..., sep = ",")
#
# cannot reliably distinguish between:
#
#   commas separating variables
#
# and
#
#   commas that are part of the country/territory name.
#
# This means that the raw data require a custom parsing
# strategy rather than a simple comma split.

unexpected_rows %>%
  select(
    line_number,
    comma_count,
    raw_line
  ) %>%
  print(n = Inf)


# ----------------------------------------------------------
# 9. Data-understanding conclusion
# ----------------------------------------------------------

message(
  paste0(
    "Raw-file inspection complete. ",
    length(boxplot_lines),
    " country-level rows were identified in the representative file. ",
    nrow(unexpected_rows),
    " rows did not contain the expected ",
    expected_commas,
    " comma delimiters and require special handling during parsing."
  )
)