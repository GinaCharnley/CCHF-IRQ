# Generating a religious holidays dataframe for Iraq 

# Load required libraries
library(dplyr)
library(lubridate)

# Define the start and end dates of your study period 
start_date <- as.Date("2022-01-01")
end_date <- as.Date("2025-12-31")

# Create a weekly sequence (epidemiological weeks) 
weeks_df <- data.frame(
  date = seq(from = floor_date(start_date, "week", week_start = 1),
                       to = floor_date(end_date, "week", week_start = 1),
                       by = "1 week")
) %>%
  mutate(
    year = year(date),
    epiweek = epiweek(date)
  )

# Manually define major Islamic holiday dates for Iraq 
# (approximate based on official announcements; actual local sighting may vary by ±1 day)
eid_dates <- as.Date(c(
  # Eid al-Fitr
  "2022-05-02", "2023-04-21", "2024-04-10", "2025-03-31",
  # Eid al-Adha
  "2022-07-09", "2023-06-28", "2024-06-16", "2025-06-06"
))

# Create a helper function to check if a date falls in a given epi week 
weeks_df <- weeks_df %>%
  rowwise() %>%
  mutate(
    holiday_this_week = any(eid_dates >= date & eid_dates < date + 7),
    holiday_previous_week = any(eid_dates >= (date - 7) & eid_dates < date)
  ) %>%
  ungroup() %>%
  mutate(
    holiday_this_week = as.integer(holiday_this_week),
    holiday_previous_week = as.integer(holiday_previous_week)
  ) %>% select(-date)

# Save 
write_csv(weeks_df, "01_data_processing/holidays/irq_relig_hol.csv")

