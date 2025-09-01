# NOTE FOR REVIEWERS / APPENDIX USE:
# This script is provided for methods transparency. It is NOT runnable as-is:
# - All file paths and filenames are redacted.
# - Input data are protected (CPRD) and cannot be shared.


###################################################################################
###################################################################################
# Antibiotic extraction

library(data.table)
library(tictoc)


# Load the antibiotic code list
abx_codes <- fread("APPENDIX_8_ANTIBIOTIC_CODES_FINAL.txt")
abx_prodcodes <- abx_codes[Include == 1, as.character(ProdCodeId)]

# Check length
cat("Number of antibiotic product codes included: ", length(abx_prodcodes), "\n")

# Combine drug list files
drug_list <- c(drug_list_1, drug_list_2)
selected_drug_list <- list()

tic("Loading and filtering antibiotic prescriptions")
for (f in drug_list) {
  tmp <- fread(f, colClasses = 'character')
  tmp_abx <- tmp[prodcodeid %in% abx_prodcodes, 
                 .(patid, issuedate = as.Date(issuedate), prodcodeid)]
  selected_drug_list <- append(selected_drug_list, list(tmp_abx))
}
drug_abx <- rbindlist(selected_drug_list)
toc()

# Filter out any abnormal dates before subset
drug_abx <- drug_abx[issuedate >= as.IDate("2000-01-01") & issuedate <= as.IDate("2025-12-31")]


# Check that the dates are plausible
summary(gi_obs$eventdate)

# Load and filter GI events
tic("Loading and filtering GI observation data")
gi_obs <- readRDS("study_observations.rds")
gi_obs[, eventdate := as.Date(as.character(obsdate), format = "%d/%m/%Y")]
gi_obs <- gi_obs[medcodeid %in% GI3_medcodes_vec, .(patid, medcodeid, eventdate)]
toc()

# Subset for testing 
subset_patid <- unique(drug_abx$patid)[1:5000]
drug_abx_sub <- drug_abx[patid %in% subset_patid]
gi_obs_sub   <- gi_obs[patid %in% subset_patid]
str(gi_obs_sub)


# Create windows: 0–7 days before abx prescribing
drug_abx_sub[, start := issuedate - 7L]
drug_abx_sub[, end := issuedate]
setkey(drug_abx_sub, patid, start, end)

gi_obs_sub[, start := eventdate]
gi_obs_sub[, end := eventdate]
setkey(gi_obs_sub, patid, start, end)

# Perform join using foverlaps
#prepare by removing NA 
# Drop rows with NA in range columns (start or end)
gi_obs_sub <- gi_obs_sub[!is.na(start) & !is.na(end)]
drug_abx_sub <- drug_abx_sub[!is.na(start) & !is.na(end)]
cat("Dropped GI rows with missing dates: ", sum(is.na(gi_obs_sub$start) | is.na(gi_obs_sub$end)), "\n")
cat("Dropped ABX rows with missing dates: ", sum(is.na(drug_abx_sub$start) | is.na(drug_abx_sub$end)), "\n")


tic("Running foverlaps join")
abx_with_gi <- foverlaps(gi_obs_sub, drug_abx_sub,
                         by.x = c("patid", "start", "end"),
                         by.y = c("patid", "start", "end"),
                         type = "within",
                         nomatch = 0L)
toc()

str(abx_with_gi)


#double check extraction with date format 
# Combine drug list files
drug_list <- c(drug_list_1, drug_list_2)
selected_drug_list <- list()

tic("Loading and filtering antibiotic prescriptions")
for (f in drug_list) {
  tmp <- fread(f, colClasses = 'character')
  tmp_abx <- tmp[prodcodeid %in% abx_prodcodes, 
                 .(patid, issuedate = as.Date(issuedate, format = "%d/%m/%Y"), prodcodeid)]
  selected_drug_list <- append(selected_drug_list, list(tmp_abx))
}
drug_abx <- rbindlist(selected_drug_list)
toc()

# Sense checks on extracted antibiotic prescription data
cat("\n--- Antibiotic Extraction Summary ---\n")
cat("Total prescriptions extracted: ", nrow(drug_abx), "\n")
cat("Unique patients: ", length(unique(drug_abx$patid)), "\n")
cat("Date range: ", format(min(drug_abx$issuedate, na.rm = TRUE)), "to", format(max(drug_abx$issuedate, na.rm = TRUE)), "\n")
cat("Missing issuedate entries: ", sum(is.na(drug_abx$issuedate)), "\n")
cat("Missing prodcodeid entries: ", sum(is.na(drug_abx$prodcodeid)), "\n")
cat("Sample of prodcodeid values:\n")


#final step cleaning 
#keep only the first GI event in the 7 day window, if multiple on the same day randomly select one

set.seed(123)  # for reproducibility 

abx_with_gi_dedup <- abx_with_gi %>%
  mutate(
    days_diff = as.integer(as.Date(issuedate) - as.Date(eventdate))
  ) %>%
  group_by(patid, issuedate) %>%
  filter(days_diff == min(days_diff)) %>%
  slice_sample(n = 1) %>%
  ungroup()

#recheck for duplicates 
abx_gi_event_counts_dedup <- abx_with_gi_dedup %>%
  group_by(patid, issuedate) %>%
  summarise(
    n_gi_events = n_distinct(eventdate),  # or medcodeid if preferred
    .groups = "drop"
  )

multi_event_abx_dedup <- abx_gi_event_counts_dedup %>%
  filter(n_gi_events > 1)

# Count how many remain
nrow(multi_event_abx_dedup)

#filter dates 2010-2020
abx_with_gi_filtered <- abx_with_gi_dedup %>%
  filter(issuedate >= as.Date("2010-01-01") & issuedate <= as.Date("2020-1-1"))

range(abx_with_gi_filtered$issuedate)  # Confirm date range

#remove columns not needed for join to anderson gill and rename
abx_with_gi_clean <- abx_with_gi_filtered %>%
  rename(abx_gi_event_date = issuedate) %>%
  select(patid, abx_gi_event_date)

rm(abx_with_gi_dedup, abx_with_gi_filtered)



