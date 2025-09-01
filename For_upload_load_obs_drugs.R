#Step 1 - Libraries and code lists
# Load required packages
library(data.table)
library(dplyr)
library(lubridate)

# Load code lists (: GI_medcodes, rotavacc_medcodes, pcv_medcodes)
getwd()
source("location")
source("ref_medcodes_prodcodes_IVAR.R")

# Remove prodcodes for non-PCV pneumococcal vaccines (given to children at risk only > age 2)
non_pcv_prodcodes <- c("12706141000033114", "2182541000033114", "1094441000033111")
PCV_prodcodes <- setdiff(PCV_prodcodes, non_pcv_prodcodes)


#Step 2: load and filter study observations
#goals: 
#1. Load the file
#2. Keep only rows where medcodeid is in one of the relevant lists
#3. Keep all rows, no summarising
#4. Add a column to identify what each code is (GI, rota vacc, or PCV)
#5. Convert obsdate to proper Date format



# Load observations
obs_full <- fread("study_observations.csv")

# Convert medcode to character for matching
obs_full[, medcodeid := as.character(medcodeid)]

# Create a tagging table to label exposures
medcode_tags <- rbindlist(list(
  data.table(medcodeid = rotavacc_medcodes, exposure = "rotavirus vaccination"),
  data.table(medcodeid = PCV_medcodes, exposure = "PCV vaccination"),
  data.table(medcodeid = GI3_medcodes, exposure = "GI")
), use.names = TRUE, fill = TRUE)

# Join tags to observations
obs_full <- merge(obs_full, medcode_tags, by = "medcodeid", all.x = FALSE)

# Convert date
obs_full[, obsdate := dmy(obsdate)]
str(obs_full)

# Keep relevant variables (retain medcodeid and enterdate for later use)
obs_full <- obs_full[, .(patid, obsdate, medcodeid, exposure, enterdate)]

# Rename obsdate to standard name
setnames(obs_full, "obsdate", "event_date")

#Step 3: load and clean
#goals 1. Loads drug data.
# 2.Links prodcodeids to exposures.
# 3.Converts date formats.
# 4.Keeps only essential variables.
# 5.Renames date to event_date for consistency.
### Load selected drug issues

tic("Loading drug issue data")
drugs_full <- fread("location")
toc()

# Convert prodcodeid to character (for safe joining)
drugs_full$prodcodeid <- as.character(drugs_full$prodcodeid)

### Sense-check: Number of patients with PCV or rotavirus codes
print(rota_prodcodes)
drugs_full %>% filter(prodcodeid %in% PCV_prodcodes) %>% distinct(patid) %>% tally()
drugs_full %>% filter(prodcodeid %in% rota_prodcodes) %>% distinct(patid) %>% tally()

#(delete later) check all drug issues 
nrow(drugs_full)
length(unique(drugs_full$patid))
summary(drugs_full)

#(delete later) check recorded as medcodes
obs_full %>% 
  filter(medcodeid %in% rotavacc_medcodes) %>% 
  distinct(patid) %>% 
  tally()


### Create exposure labels
rota_prodcodes_df <- data.frame(prodcode = rota_prodcodes, exposure = "rotavirus vaccination")
pcv_prodcodes_df  <- data.frame(prodcode = PCV_prodcodes,  exposure = "PCV vaccination")

# Convert prodcodeid to character (for safe joining)
drugs_full$prodcodeid <- as.character(drugs_full$prodcodeid)

# Drop event_date if it exists from previous steps
if ("event_date" %in% colnames(drugs_full)) {
  drugs_full <- drugs_full[, !("event_date"), with = FALSE]
}

#  safely join
study_prodcodes <- bind_rows(rota_prodcodes_df, pcv_prodcodes_df)
drugs_full <- left_join(drugs_full, study_prodcodes, by = c("prodcodeid" = "prodcode"))


# Remove non-PCV codes (exposure already assigned)
drugs_full <- drugs_full %>%
  filter(!(prodcodeid %in% non_pcv_prodcodes))

# Convert issuedate to Date
drugs_full$issuedate <- as.Date(dmy(drugs_full$issuedate))

# Keep needed columns
drugs_full <- drugs_full[, .(patid, event_date = issuedate, prodcodeid, exposure)]


### Step 4 - Combine exposures from observations and drugs

exp_full <- bind_rows(
  obs_full[, .(patid, event_date, exposure)],
  drugs_full[, .(patid, event_date, exposure)]
)

# Check structure
str(exp_full)


### Step 6 - Generate binary vaccine flags per person

binary_vax_status <- exp_full %>%
  filter(exposure %in% c("rotavirus vaccination", "PCV vaccination")) %>%
  distinct(patid, exposure) %>%
  mutate(status = 1) %>%
  pivot_wider(names_from = exposure, values_from = status, values_fill = 0)



### Study observation and drug issue reshaping for analysis
## Updated for rotavirus + PCV, Andersen-Gill, multiple events

#sense checks 
str(binary_vax_status)
# Rotavirus vaccination coverage
mean(binary_vax_status$`rotavirus vaccination`) * 100

# PCV vaccination coverage
mean(binary_vax_status$`PCV vaccination`) * 100

#check no None-PCV codes remaining
drugs_full %>%
  filter(prodcodeid %in% non_pcv_prodcodes)

