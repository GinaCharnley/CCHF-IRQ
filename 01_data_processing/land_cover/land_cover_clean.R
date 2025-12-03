# Format the processed global Copernicus Land Monitoring Service file 

# Load packages
library(dplyr)
library(sf)

# Load the file 
# This is already processed from a different project 
land_copern <- readRDS("01_data_processing/land_cover/copernicus_global_land_use_2019_0.25.rds")

# Crop to southern Iraq 
land_copern <- land_copern %>% filter(x > 40 & x < 50 & y > 28 & y < 35)

# Read shapefile
irq <- st_read("02_data_harmonisation/gadm41_IRQ_shp/gadm41_IRQ_1.shp", quiet = TRUE) %>% st_transform(4326)

# Convert grid to sf object (points)
grid_sf <- st_as_sf(land_copern, coords = c("x", "y"), crs = 4326)

# Make sure both have the same CRS
irq <- st_transform(irq, crs = st_crs(grid_sf))

# Spatial join: adds polygon attributes to points
grid_with_gid <- st_join(grid_sf, irq, join = st_within)

# Save land cover classes
grid_with_gid <- grid_with_gid %>% st_drop_geometry()
classes <- colnames(grid_with_gid[1:11])

# Calculate spatial coverage by GID_1 
grid_summary <- grid_with_gid %>% 
  group_by(GID_1) %>%
  summarise(across(all_of(classes), mean, na.rm = TRUE))

# Tidy and save 
grid_summary$Forest[is.nan(grid_summary$Forest)]<-0
grid_summary <- na.omit(grid_summary)

saveRDS(grid_summary, "01_data_processing/land_cover/irq_land.rds")




