# Function for transforming the daily soil moisture to weekly

process_soil_moisture <- function(
    input_dir = "01_data_processing/climate/ncs/soil_moisture",
    shapefile_path = "02_data_harmonisation/gadm41_IRQ_shp/gadm41_IRQ_1.shp"
) {
  library(terra)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(lubridate)
  library(epiR)
  
  message("📂 Reading all daily soil moisture files...")
  files <- list.files(input_dir, full.names = TRUE, pattern = "\\.nc$")
  if (length(files) == 0) stop("No NetCDF files found in input directory.")
  
  message(paste("Found", length(files), "files"))
  
  #--- Read shapefile and ensure correct CRS
  irq <- read_sf(shapefile_path) %>% st_transform(4326)
  
  #--- Define crop extent (Iraq region)
  region_extent <- terra::ext(41, 49, 29, 34)
  
  #--- Function to read, crop, and extract date
  process_file <- function(f) {
    # Extract date (8 digits after "DAILY-")
    date_str <- sub(".*DAILY-([0-9]{8}).*", "\\1", basename(f))
    date_val <- as.Date(date_str, format = "%Y%m%d")
    if (is.na(date_val)) stop(paste("Could not parse date from:", f))
    
    # Read raster
    r <- terra::rast(f)
    terra::crs(r) <- "EPSG:4326"
    
    # Crop to Iraq region
    r <- terra::crop(r, region_extent)
    
    # Extract first layer (soil moisture)
    r <- r[[1]]
    
    # Name it by date
    names(r) <- as.character(date_val)
    
    return(r)
  }
  
  #--- Read all rasters
  message("🧩 Reading and cropping rasters...")
  rlist <- lapply(files, process_file)
  
  #--- Combine into a single SpatRaster
  rstack <- terra::rast(rlist)
  
  #--- Convert to data.frame (x, y, date, soil_moisture)
  df <- as.data.frame(rstack, xy = TRUE)
  
  # Pivot longer: columns = x, y, date, soil_moisture
  df_long <- df %>%
    pivot_longer(cols = starts_with("20"), names_to = "date", values_to = "soil_moisture") %>%
    mutate(date = as.Date(date))
  
  #--- Compute epiweek and year
  df_long <- df_long %>%
    mutate(
      year = year(date),
      epiweek = epiweek(date)
    )
  
  #--- Average daily → weekly
  df_weekly <- df_long %>%
    group_by(x, y, year, epiweek) %>%
    summarise(soil_moisture = mean(soil_moisture, na.rm = TRUE), .groups = "drop")
  
  #--- Snap coordinates to grid
  snap_to_grid <- function(x, res = 0.25) round(x / res) * res
  df_weekly$lon_grid <- snap_to_grid(df_weekly$x)
  df_weekly$lat_grid <- snap_to_grid(df_weekly$y)
  
  #--- Create polygons for each grid cell
  make_cell <- function(lon, lat, res = 0.25) {
    st_polygon(list(rbind(
      c(lon - res/2, lat - res/2),
      c(lon + res/2, lat - res/2),
      c(lon + res/2, lat + res/2),
      c(lon - res/2, lat + res/2),
      c(lon - res/2, lat - res/2)
    )))
  }
  
  era5_cells <- st_sfc(
    mapply(make_cell, df_weekly$lon_grid, df_weekly$lat_grid, SIMPLIFY = FALSE),
    crs = 4326
  )
  
  #--- Spatial join
  overlaps <- st_intersects(era5_cells, irq)
  df_weekly$GID_1 <- sapply(overlaps, function(idx) {
    if (length(idx) == 0) return(NA)
    paste(irq$GID_1[idx], collapse = ",")
  })
  
  df_weekly <- df_weekly %>%
    separate_rows(GID_1, sep = ",") %>%
    filter(!is.na(GID_1))
  
  #--- Aggregate by shapefile region
  data_final <- df_weekly %>%
    group_by(GID_1, year, epiweek) %>%
    summarise(soil_moisture = mean(soil_moisture, na.rm = TRUE), .groups = "drop")
  
  message("✅ Soil moisture processing complete.")
  return(data_final)
}

sm_all <- process_soil_moisture("01_data_processing/climate/ncs/soil_moisture")

irq_climate <- readRDS("01_data_processing/climate/irq_climate.rds")

irq_climate <- left_join(irq_climate, sm_all)

saveRDS(irq_climate, "01_data_processing/climate/irq_climate.rds")



