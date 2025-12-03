# Harmonising and cleaning the displacement data 

# Load packages
library(dplyr)
library(readr)
library(lubridate)

# Load data 
idps <- read_delim("01_data_processing/displacement/irq_idps.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Example GID_1 codes
gid_1_codes <- c("IRQ.15_1", "IRQ.2_1", "IRQ.11_1", "IRQ.14_1")

# Generate full sequence of weeks
start_date <- as.Date("2022-01-03") # first Monday of 2022 (epiweek 1)
end_date <- as.Date("2025-12-28")   # last Sunday of 2025

# create a weekly sequence
all_weeks <- seq.Date(start_date, end_date, by = "week")

# extract year and epiweek
week_grid <- data.frame(
  week_start = all_weeks
) %>%
  mutate(
    year = isoyear(week_start),
    epiweek = isoweek(week_start)
  )

# create full grid with GID_1
full_grid <- expand.grid(
  GID_1 = gid_1_codes,
  year = unique(week_grid$year),
  epiweek = 1:53  # max epiweek (some years have 52)
) %>%
  left_join(week_grid, by = c("year", "epiweek")) %>%
  arrange(GID_1, year, epiweek)

# Add months
full_grid <- full_grid %>%
  mutate(month = month(week_start)) %>% select(-week_start)
full_grid$month[is.na(full_grid$month)]<-12

# Join
idps <- idps %>% rename("year" = Year, "month" = Month)
full_weekly_data <- full_grid %>%
  left_join(idps, by = c("GID_1", "year", "month")) %>%
  arrange(GID_1, year, month) %>%  # sort by actual week start date
  group_by(GID_1) %>%
  fill(idp_no, .direction = "downup") %>%
  ungroup()

# Remove un-nessessary columns 
full_weekly_data <- full_weekly_data %>% select(-NAME_1, -month)

# Save data 
write_rds(full_weekly_data, "01_data_processing/displacement/irq_idps.rds")



















