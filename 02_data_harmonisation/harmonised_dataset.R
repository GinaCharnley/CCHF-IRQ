# Combine all the cleaned data to one harmonised dataset 

# Load packages 
library(dplyr)
library(readr)

# Load all clean data 
land <- read_rds("01_data_processing/land_cover/irq_land.rds")
elev <- read_csv("01_data_processing/elevation/irq_elev.csv")
conf <- read_rds("01_data_processing/conflict/conflict_clean.rds")
clim <- read_rds("01_data_processing/climate/irq_climate.rds")
hols <- read_csv("01_data_processing/holidays/irq_relig_hol.csv")
cchf <- read_csv("01_data_processing/cchf/irq_cchf.csv")
idps <- read_delim("01_data_processing/displacement/irq_idps.rds")
pov <- read_delim("01_data_processing/poverty/irq_poverty.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Merge into one dataframe 
harmonised <- left_join(cchf, land)
harmonised <- left_join(harmonised, elev)
harmonised <- left_join(harmonised, conf)
harmonised <- left_join(harmonised, clim)
harmonised <- left_join(harmonised, hols)
harmonised <- left_join(harmonised, pov)
harmonised <- left_join(harmonised, idps)

# Tidy 
harmonised <- harmonised %>% rename("conf_events" = n)
harmonised$conf_events[is.na(harmonised$conf_events)]<-0

# Save
write_rds(harmonised, "02_data_harmonisation/harmonised_set.rds")










