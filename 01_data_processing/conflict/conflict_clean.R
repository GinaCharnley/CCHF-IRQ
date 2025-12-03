# Process the ACLED conflict data 

# Load packages 
library(dplyr)
library(sf)
library(lubridate)

# Load data 
conflict_data_irq <- read_csv("01_data_processing/conflict/conflict_data_irq.csv")

# Select the columns of interest 
conf <- conflict_data_irq %>% select(date_start, latitude, longitude)
conf <- conf[-c(1),]

# Read shapefile
irq <- st_read("02_data_harmonisation/gadm41_IRQ_shp/gadm41_IRQ_1.shp", quiet = TRUE) %>% st_transform(4326)

# Convert grid to sf object (points)
grid_sf <- st_as_sf(conf, coords = c("longitude", "latitude"), crs = 4326)

# Make sure both have the same CRS
irq <- st_transform(irq, crs = st_crs(grid_sf))

# Spatial join: adds polygon attributes to points
grid_with_gid <- st_join(grid_sf, irq, join = st_within)

# Tidy dates 
conf <- grid_with_gid %>% select(date_start, GID_1) %>% st_drop_geometry()
conf <- conf %>%
  # First, convert the character column to POSIXct
  mutate(date_start = ymd_hms(date_start)) %>%
  # Extract year, month, day, and epiweek
  mutate(
    year = year(date_start),
    month = month(date_start),
    day = day(date_start),
    epiweek = epiweek(date_start)
  )

# Tally number of events for GID and epiweek 
conf_sum <- conf %>% group_by(GID_1, year, epiweek) %>% tally()

# Save 
saveRDS(conf_sum, "01_data_processing/conflict/conflict_clean.rds")







