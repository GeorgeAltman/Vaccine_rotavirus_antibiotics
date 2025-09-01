# NOTE FOR REVIEWERS / APPENDIX USE:
# This script is provided for methods transparency. It is NOT runnable as-is:
# - All file paths and filenames are redacted.
# - Input data are protected (CPRD) and cannot be shared.


### Step 1: Study observation and drug issue extraction

## Load packages

library(data.table)
library(lubridate)
library(dplyr)
library(stringr)
library(tictoc)      



### Load medcode lookup
tic()
medcode_look <- fread("load medcode_lookup")
toc()

### Load the required medcodes (observations) and prodcodes (drug issue)
source('ref_medcodes_prodcodes_IVAR.R')
GI3_medcodes <- readRDS("GI3_medcodes location")
str(GI3_medcodes)
GI3_medcodes_vec <- as.character(GI3_medcodes$med_code_id)
length(GI3_medcodes_vec)  


# Combine list of desired medcodes
study_medcodes <- c(rotavacc_medcodes, PCV_medcodes, GI3_medcodes_vec)
length(study_medcodes)


### Load observation data to test
tic()
obs_path <- file.path("file location")
obs01 <- fread(obs_path, nrow = 100000)
toc()

# convert medcodeid to character for matching
obs01$medcodeid <- as.character(obs01$medcodeid)

sel_obs <- obs01 %>% filter(medcodeid %in% study_medcodes)
########################################################################

obs_list_1 <- Sys.glob("location")
obs_list_2 <- Sys.glob("location")

# Combine the file lists from both folders
obs_list <- c(obs_list_1, obs_list_2)
length(obs_list)

# Create an empty list to store selected observations
selected_obs_list <- list()

# Start the overall timing
tic("Processing files")

# Loop over each file and select only specific medcodes
for (f in obs_list) {
  tic(paste("Processing file", f))
  
  # Read the file as a character data.table to preserve leading zeros, etc.
  tmp <- fread(f, colClasses = 'character')
  
  # Select only columns 1, 2, 3, 4, 5, 6, 9 by column number
  tmp_selected <- tmp[, c(1, 2, 3, 4, 5, 6, 9), with = FALSE]
  
  # Filter by medcodeid (assuming 'study_medcodes' is predefined)
  sel_obs <- tmp_selected %>% filter(medcodeid %in% study_medcodes)
  
  # Append the selected rows to the list
  selected_obs_list <- append(selected_obs_list, list(sel_obs))
  
  toc()  # End timing for the current file
}

# Combine all the data.frames into one
df <- bind_rows(selected_obs_list)


# End overall timing
toc()


gc()
########################################################################

######## Step 2: Drug issue

### Load product medcode lookup
tic()
prod_path <- "location"
prodcode_look <- fread(prod_path)
toc()

# Combine list of desired medcodes
study_prodcodes <- c(rota_prodcodes, PCV_prodcodes)
length(study_prodcodes)
length(rota_prodcodes)

# Get list of observation files from two directories
drug_list_1 <- Sys.glob("location")
drug_list_2 <- Sys.glob("location")


# Combine the file lists from both folders
drug_list <- c(drug_list_1, drug_list_2)

# Create a list to store selected observations
selected_obs_list <- list()

# Start the overall timing
tic("Processing drug files")

# Loop over each file and select only specific prodcodes
for (f in drug_list) {
  tic(paste("Processing file", f))
  
  # Read the file
  tmp <- fread(f)
  
  # Only convert prodcodeid to character if it's not already
  if (!is.character(tmp$prodcodeid)) {
    tmp$prodcodeid <- as.character(tmp$prodcodeid)
  }
  
  # Filter by prodcodeid 
  sel_obs <- tmp %>% filter(prodcodeid %in% study_prodcodes)
  
  # Append to the list of selected observations
  selected_obs_list <- append(selected_obs_list, list(sel_obs))
  
  toc()  
}

# Combine all the data.frames into one at the end
df <- bind_rows(selected_obs_list)

# End overall timing
toc()

###################################################################################
###################################################################################
# Step 3: Antibiotic extraction

library(data.table)
library(tictoc)

# Set working directory


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

# Optional flag in drug_abx if GI event found
drug_abx_sub[, gi_flag := patid %in% abx_with_gi$patid & 
               issuedate %in% abx_with_gi$issuedate]
