# Harmonise the CCHF data to the common spatial/temporal scale 

# Load packages 
library(dplyr)
library(readr)
library(sf)
library(lubridate)

# Load the pre-cleaned data file 
cchf <- read_excel("01_data_processing/cchf/cchf_clean.xlsx", sheet = "Sheet1")

# Load the Iraq shapefile 
irq <- st_read("02_data_harmonisation/gadm41_IRQ_shp/gadm41_IRQ_1.shp", quiet = TRUE) %>% st_transform(4326)

# Load place name key 
districts <- read_delim("02_data_harmonisation/irq_places_key/district-Table 1.csv", 
                        delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Merge key to CCHF data 
cchf <- left_join(cchf, districts)

# Separate the date column and convert to epiweek 
cchf <- cchf %>%
  mutate(
    year = year(Dateofonset),
    month = month(Dateofonset),
    day = day(Dateofonset),
    epiweek = epiweek(Dateofonset)  # epidemiological week
  )

# Add an age category 
cchf <- cchf %>% mutate(age_group = case_when(Ageyear < 16 ~ "Children",
                                              Ageyear > 15 & Ageyear < 60 ~ "Adults",
                                              Ageyear > 59 ~ "Older Adults")) 

# Tally for the different groups of cases 
cchf2 <- cchf %>% 
  summarise(
    total_cases = n(),
    male_cases = sum(Gender == "Male", na.rm = TRUE),
    female_cases = sum(Gender == "Female", na.rm = TRUE),
    child_cases = sum(age_group == "Children", na.rm = TRUE),
    adult_cases = sum(age_group == "Adults", na.rm = TRUE),
    elderly_cases = sum(age_group == "Older Adults", na.rm = TRUE),
    .by = c(GID_1, year, epiweek)
  )

# Define full epiweek/year range
epiweeks <- expand.grid(
  GID_1 = unique(cchf2$GID_1),
  year = 2022:2025,
  epiweek = 1:53
)

# Join and fill missing combinations with zeros ---
cchf2 <- epiweeks %>%
  left_join(cchf2, by = c("GID_1", "year", "epiweek")) %>%
  mutate(across(c(total_cases, male_cases, female_cases,
                  child_cases, adult_cases, elderly_cases),
                ~ tidyr::replace_na(., 0)))

# Save the cleaned file 
write_csv(cchf2, "01_data_processing/cchf/irq_cchf.csv")














