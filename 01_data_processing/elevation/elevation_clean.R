# Process NASA Shuttle Radar Topography Mission Global 
# https://lpdaac.usgs.gov/products/srtmgl3v003/

# Load packages 
library(dplyr)
library(terra)
library(sf)

# Load file
rast <- rast("01_data_processing/elevation/wc2.1_10m_elev.tif")

# Set the desired resolution 
# 0.1° x 0.1° resolution
new_res <- c(0.25, 0.25)  

# Convert the raster to the correct spatial scale
template <- rast(ext(rast), res = new_res, crs = crs(rast))
rast2 <- resample(rast, template, method = "bilinear")

# Convert to a data frame 
rast2 <- as.data.frame(rast2, xy = TRUE)

# Crop to southern Iraq 
rast2 <- rast2 %>% filter(x > 40 & x < 50 & y > 28 & y < 35)

# Read shapefile
irq <- st_read("02_data_harmonisation/gadm41_IRQ_shp/gadm41_IRQ_1.shp", quiet = TRUE) %>% st_transform(4326)

# Convert grid to sf object (points)
grid_sf <- st_as_sf(rast2, coords = c("x", "y"), crs = 4326)

# Make sure both have the same CRS
irq <- st_transform(irq, crs = st_crs(grid_sf))

# Spatial join: adds polygon attributes to points
grid_with_gid <- st_join(grid_sf, irq, join = st_within)

# Spatially aggregate 
grid_with_gid <- grid_with_gid %>%
  group_by(GID_1) %>%
  summarise(elevation = mean(elevation, na.rm = TRUE)) %>%
  ungroup() %>% st_drop_geometry()

# Save 
readr::write_csv(grid_with_gid, "01_data_processing/elevation/irq_elev.csv")




