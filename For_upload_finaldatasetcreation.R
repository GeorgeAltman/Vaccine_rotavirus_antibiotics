# NOTE FOR REVIEWERS / APPENDIX USE:
# This script is provided for methods transparency. It is NOT runnable as-is:
# - All file paths and filenames are redacted.
# - Input data are protected (CPRD) and cannot be shared.

### Study population extraction

## Load packages

library(data.table)
library(lubridate)
library(dplyr)
library(tictoc)       
library(ggplot2)
library(tidyr)


### Step 1: Load patient and practice data

# 1.1 Load practice data
tic()
prac1 <- fread("")
prac2 <- fread("")
#bind
prac <- rbind(prac1, prac2)
toc()

# 1.2 change dates to date format
prac$lcd <- dmy(prac$lcd)

# 1.3 load patient data
tic()
# Load both sets
pat1 <- fread("")
pat2 <- fread("")
# Bind together
pat <- rbind(pat1, pat2)
toc()
gc()

# 1.4 change dates to date format
pat$emis_ddate <- dmy(pat$emis_ddate)
pat$regstartdate <- dmy(pat$regstartdate)
pat$regenddate <- dmy(pat$regenddate)
pat$cprd_ddate <- dmy(pat$cprd_ddate)

# 1.5 remove partial objects to free memory
rm(pat1, pat2, prac1, prac2)

# 1.6 check for duplicated in prac
prac %>%
  group_by(pracid) %>%
  summarise(n = n()) %>%
  filter(n > 1)

#1.7 deduplicate prac
prac <- prac %>%
  arrange(pracid) %>%
  distinct(pracid, .keep_all = TRUE)

#1.8 join practice ID to patient data
patprac <- left_join(pat, prac, by = "pracid")

#1.9check number of missing pracid in patprac
sum(is.na(patprac$pracid))


# 1.11 convert patient id to character variable for matching
patprac$patid <- as.character(patprac$patid)

# 1.12 number of gp surgeries
length(unique(patprac$pracid))

# 1.13 number of unique individuals
length(unique(patprac$patid))

# 1.14 remove individual patient and practice data.frames
rm("pat", "prac")

#####################################################

## Step 2: Join exposure data

# 2.1 filter outside 2010-2020 early 
patprac <- patprac %>%
  filter(yob >= 2010 & yob <= 2020)

# 2.2 recheck number of unique individuals
length(unique(patprac$patid))


## 2.4 NEW - Read in exposure datapatient data
str(exp_full)

# 2.5 Ensure IDs are character and dates are parsed
exp_full$patid <- as.character(exp_full$patid)
exp_full$event_date <- ymd(exp_full$event_date)


# 2.6: Filter relevant exposures and get first date per exposure
relevant_exposures <- c("PCV vaccination", "rotavirus vaccination")

#2.7 remove all but first event 
exp_first <- exp_full %>%
  filter(exposure %in% relevant_exposures) %>%
  mutate(patid = as.character(patid)) %>%
  group_by(patid, exposure) %>%
  summarise(first_date = min(event_date, na.rm = FALSE), .groups = "drop")  # allow NAs


# Step 2.8: Pivot to wide format
exp_wide <- exp_first %>%
  pivot_wider(
    names_from = exposure,
    values_from = first_date
  ) %>%
  rename_with(~ gsub(" ", "_", tolower(.x))) %>%
  rename(
    rota_vacc_date = rotavirus_vaccination,
    pcv_vacc_date  = pcv_vaccination
  )

# Step 2.9: Join to all patients in `patprac` (left join keeps unvaccinated)
patprac <- patprac %>%
  mutate(patid = as.character(patid)) %>%
  left_join(exp_wide, by = "patid")

# Step 2.10: Create binary flags (still NA = unvaccinated, no events)
patprac <- patprac %>%
  mutate(
    rota_vacc = !is.na(rota_vacc_date),
    pcv_vacc  = !is.na(pcv_vacc_date)
  )

#####################################################

## Step 3: Join healthcare contact information

## 3.1 Read in number of consultations per patient
tic()
all_cons <- fread("")
toc()


# 3.2 change patient id to character for matching
all_cons$patid <- as.character(all_cons$patid)

# 3.3 make sure patid is unique
all_cons <- all_cons %>% 
  group_by(patid) %>% 
  summarise(num_consults = sum(consults))

# 3.4 join to main data
patprac <- left_join(patprac, all_cons, by = c("patid" = "patid"))

# 3.5 replace NA with 0
patprac$num_consults[is.na(patprac$num_consults)] <- 0


#####################################################

## Step 4 Join IMD data to patient data

## 4.1 Read in IMD data
tic()
imd <- fread("")
toc()
gc()

# 4.2 change patient id to character for matching
imd$patid <- as.character(imd$patid)

# 4.3 remove pracid variable
imd$pracid <- NULL

# 4.4 join to main data
patprac <- left_join(patprac, imd, by = c("patid" = "patid"))

# 4.5 check how many don't match
sum(is.na(patprac$e2019_imd_5))

###################################################
## Step 5 - Join Urban-rural data to patient data

## 5.1 Read in urban-rural data
tic()
ur <- fread(")
toc()
gc()

# 5.2 change patient id to character for matching
ur$patid <- as.character(ur$patid)

# 5.3 remove pracid variable
ur$pracid <- NULL

# 5.4 join to main data
patprac <- left_join(patprac, ur, by = c("patid" = "patid"))

# 5.5 check how many don't match
sum(is.na(patprac$e2011_urban_rural))

#5.6 remove for memory 
rm("ur","imd","all_cons")

#######


## Step 6 - Create DOB , then work out rota group, then start date

# 6.1. create proxy date of birth
patprac$dob <- paste0(patprac$yob, "/", patprac$mob, "/", "01")

# 6.2 parse using lubridate
patprac$dob <- ymd(patprac$dob)

# 6.3 create rotagroup if born after April 2013 will have been eligible for vaccine
patprac <- patprac %>%
  mutate(rota_group = if_else(dob >= as.Date("2013-05-01"), 1, 0))

# 6.4 Create start date as dob + 24 weeks (24*7=168) !!! changed from 12 weeks (12*7)
patprac$start_date <- as.Date(patprac$dob + days(168))

# 6.5 Difference between date of birth and registration start date
patprac$regdelay <- as.numeric(patprac$regstartdate - patprac$dob)

# 6.6 Show registration delay
ggplot(patprac, aes(x = regdelay)) +
  geom_histogram() +
  scale_x_log10()

# 6.7 Number without a registration end date
sum(is.na(patprac$regenddate))

# 6.8 Create new end date variable based on registration end date
patprac$end_date <- patprac$regenddate

# 6.9 If NA (no end date), replace with practice last collection date
patprac$end_date <- dplyr::coalesce(patprac$end_date, patprac$lcd)
gc()


#####################################################
#Step 7 - getting correct end date censoring for end of study 2019, last collection date , last reg date and death


# 7.1  Recalculate start_date: 6 months from date of birth
patprac <- patprac %>%
  mutate(start_date = dob + months(6))

# 7.2 Create censoring date at 7th birthday
patprac <- patprac %>%
  mutate(date_age_7 = start_date + years(7))

# 7.3 Derive end_date as the earliest of regenddate, cprd_ddate, lcd, and age 7, and end of 2019
#old patprac <- patprac %>%
#old  mutate(end_date = pmin(regenddate, cprd_ddate, lcd, date_age_7, na.rm = TRUE))
patprac <- patprac %>%
  mutate(
    date_age_7 = dob + years(7),
    end_date = pmin(regenddate, cprd_ddate, lcd, date_age_7, as.Date("2019-12-30"), na.rm = TRUE)
  )

# 7.4 calculate study consultation rate, uses DOB as assumed study consults is all from birth - to be checked 
patprac <- patprac %>%
  mutate(
    rate_consults = round((num_consults / as.numeric((end_date - dob) / 7)) * 52, 1),
    rate_consults = pmin(pmax(rate_consults, 0), 52)  # cap between 0 and 52
  )

# 7.4 Calculate follow-up time
patprac <- patprac %>%
  mutate(follow_days = as.numeric(end_date - start_date))

# 7.5 Follow-up duration in days (total for person)
patprac <- patprac %>%
  mutate(
    followup_total_days = as.integer(end_date - start_date)
  )

# Difference in days between end_date and regstartdate
# old - patprac$follow_days <- as.numeric(patprac$end_date - patprac$start_date)

# 
# 7.6 Examine the characteristics of follow-up days
#

patprac <- patprac %>% 
  mutate(regend_valid = if_else(is.na(regenddate), 0,1))

ggplot(patprac, aes(x = follow_days)) +
  geom_histogram() +
  facet_wrap(~regend_valid)

ggplot(patprac, aes(x = follow_days)) +
  geom_histogram() +
  facet_wrap(~rota_group)

ggplot(patprac, aes(x = end_date)) +
  geom_histogram() +
  facet_wrap(~regend_valid)

# 7.7. Count individuals with negative follow-up time
sum(patprac$follow_days < 0, na.rm = TRUE)

# 7.8 View rows with negative follow-up
patprac %>%
  filter(follow_days < 0)

# 7.9 Summary stats of follow-up time
summary(patprac$follow_days)


###### 
#Step 8 - check vaccine uptake

table(patprac$rota_vacc, useNA = "ifany")
table(patprac$pcv_vacc, useNA = "ifany")

table(patprac$rota_group, patprac$rota_vacc)
table(patprac$rota_group, patprac$pcv_vacc)

patprac %>%
  summarise(
    total = n(),
    rota_uptake = mean(rota_vacc, na.rm = TRUE),
    pcv_uptake  = mean(pcv_vacc, na.rm = TRUE)
  )


vaccine_uptake_by_yob <- patprac %>%
  group_by(yob) %>%
  summarise(
    n = n(),
    rota_vaccinated = sum(rota_vacc, na.rm = TRUE),
    rota_uptake = mean(rota_vacc, na.rm = TRUE),
    pcv_vaccinated  = sum(pcv_vacc, na.rm = TRUE),
    pcv_uptake = mean(pcv_vacc, na.rm = TRUE)
  ) %>%
  arrange(yob)

print(vaccine_uptake_by_yob, n = Inf)

################################
# step 9 - Prepare clean datasets for joining to anderson gill  

# 9.1 filter out unnecessary columns from pat prac
patprac_ag <- patprac %>%
  select(-patienttypeid, -acceptable, -lcd, -uts)

# 9.2 remove those with negative follow up time 
patprac_ag <- patprac_ag %>%
  filter(follow_days >= 0)

#9.3 log how many removed 
n_removed <- sum(patprac$follow_days < 0)
cat("Removed", n_removed, "patients with negative follow-up time.\n")


# 9.4 Combine into final anderson gill - Join abx/GI events to cleaned patprac
ag_base <- patprac_ag[, .(patid, start_date, end_date, rota_group, rota_vacc, rota_vacc_date, pcv_vacc, pcv_vacc_date, yob, dob, gender,
                          e2019_imd_5, e2011_urban_rural, num_consults, rate_consults,follow_days)]


abx_events <- abx_with_gi_clean[, c("patid", "abx_gi_event_date"), with = FALSE]


# 9.5 Join to get start/end/fixed covariates
ag_joined <- merge(abx_events, ag_base, by = "patid", all.x = TRUE)

# 9.6 Remove events that occur before start or end date
ag_joined <- ag_joined[abx_gi_event_date >= start_date & abx_gi_event_date <= end_date]

# 9.7 Sort events per person and assign row numbers
ag_joined <- ag_joined[order(patid, abx_gi_event_date)]
ag_joined[, event_number := 1:.N, by = patid]

# 9.8 Create start/stop intervals
ag_long <- copy(ag_joined)

# 9.9 Ensure consistent date format
ag_long[, abx_gi_event_date := as.Date(abx_gi_event_date)]
ag_long[, start_date := as.Date(start_date)]

# 9.10 Create lagged start_time
ag_long[, start_time := shift(abx_gi_event_date, n = 1, type = "lag"), by = patid]

# 9.11 Fill first row per patient with start_date
ag_long[, start_time := fifelse(is.na(start_time), start_date, start_time)]

# 9.12 Set stop_time and event flag
ag_long[, stop_time := abx_gi_event_date]
ag_long[, event := 1]


# 9.13 add censoring row for each patient(if event doesnt end follow up)
# Get last event per patient
last_event <- ag_long[, .SD[.N], by = patid]

# 9.14 Only add censoring if stop_time < end_date
censored_rows <- last_event[stop_time < end_date, .(
  patid,
  start_time = stop_time,
  stop_time = end_date,
  event = 0
)]

# 9.15 Copy over covariates
censored_rows <- merge(censored_rows, ag_base, by = "patid", all.x = TRUE)


# 9.16 Identify columns in ag_long but missing in censored_rows
missing_cols <- setdiff(names(ag_long), names(censored_rows))

# 9.17 Add them as NA in censored_rows
for (col in missing_cols) {
  censored_rows[, (col) := NA]
}

# 9.18 combine safely
ag_ready <- rbindlist(list(ag_long, censored_rows), use.names = TRUE)
setorder(ag_ready, patid, start_time)


#Step 10 - build long dataset
# 10.1 define core covariates to remain
keep_cols <- c(
  "patid", "start_date", "end_date", "dob", "yob", "gender", "mob", "pracid", "usualgpstaffid",
  "region", "rota_vacc", "rota_vacc_date", "rota_group", "pcv_vacc", "pcv_vacc_date",
  "e2019_imd_5", "e2011_urban_rural", "num_consults", "rate_consults"
)

# 10.2 build ag_base from patprac
ag_base <- patprac_ag[, ..keep_cols]

#10.3 join clean abx with gi to ag base
names(abx_with_gi_clean)
setDT(abx_with_gi_clean)
abx_events <- abx_with_gi_clean[, .(patid, abx_gi_event_date)]
ag_joined <- merge(abx_events, ag_base, by = "patid", all.x = TRUE)
ag_joined <- ag_joined[abx_gi_event_date >= start_date & abx_gi_event_date <= end_date]

#10.4 create recurrent events structure
setorder(ag_joined, patid, abx_gi_event_date)
ag_joined[, start_date := as.Date(start_date)]
ag_joined[, start_time := as.Date(start_time)]
ag_joined[, event_number := 1:.N, by = patid]
ag_joined[, start_time := shift(abx_gi_event_date, type = "lag"), by = patid]
ag_joined[, start_time := fifelse(is.na(start_time), start_date, start_time)]
ag_joined[, stop_time := abx_gi_event_date]
ag_joined[, event := 1]


#10.5 add censoring row for patients with events 
last_event <- ag_joined[, .SD[.N], by = patid]
censored_rows <- last_event[stop_time < end_date, .(
  patid,
  start_time = stop_time,
  stop_time = end_date,
  event = 0
)]
censored_rows <- merge(censored_rows, ag_base, by = "patid", all.x = TRUE)


#10.6 add patients with no events 
patients_with_events <- unique(ag_joined$patid)
patients_without_events <- setdiff(ag_base$patid, patients_with_events)

no_event_rows <- ag_base[patid %in% patients_without_events]
no_event_rows[, `:=` (
  abx_gi_event_date = as.Date(NA),
  event_number = NA_integer_,
  start_time = start_date,
  stop_time = end_date,
  event = 0
)]


#10.7 combine all 3 datasets 
ag_ready <- rbindlist(list(ag_joined, censored_rows, no_event_rows), use.names = TRUE, fill = TRUE)
setorder(ag_ready, patid, start_time)


#10.8 check final numbers 
cat("Final row count: ", nrow(ag_ready), "\n")
cat("Unique patients: ", length(unique(ag_ready$patid)), "\n")
str(ag_ready[1:5])  # Look at first few entries

#10.9 change date format
ag_ready[, `:=`(
  start_date = as.Date(start_date),
  start_time = as.Date(start_time),
  stop_time  = as.Date(stop_time)
)]

#10.10 fill missing start_time safely
ag_ready[is.na(start_time), start_time := start_date]

# 10.11  intervals relative to start_date
ag_ready[, start_time_days := as.integer(start_time - start_date)]
ag_ready[, stop_time_days  := as.integer(stop_time  - start_date)]
ag_ready[, followup_days   := stop_time_days - start_time_days]









