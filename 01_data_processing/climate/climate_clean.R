# Function for transforming the daily ERA5 to weekly

# The files are pulled from Copernicus 
# As they are daily, they only allow you to download a year at once 
# So I have each year saved in separate .nc files in the ncs folder
# The function reads in each year, transforms them to weekly and saves it as a combined file 
# Aggregated to GADM codes 

process_variable_year <- function(year, 
                                  input_dir = "01_data_processing/climate/ncs",
                                  shapefile_path = "02_data_harmonisation/gadm41_IRQ_shp/gadm41_IRQ_1.shp") {
  library(raster)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(sf)
  library(epiR)
  library(terra)
  
  message(paste("Processing variable data for year", year, "..."))
  
  #--- Load raster files for the specified year
  files <- list.files(input_dir, full.names = TRUE, pattern = as.character(year))
  if (length(files) == 0) stop(paste("No raster files found for year", year))
  
  #--- Read raster(s)
  r <- terra::rast(files)
  
  #--- Always assign known Iraq region extent (since Copernicus clipping sometimes omits it)
  message("🗺️  Assigning explicit extent for Iraq region: 41–49E, 29–34N")
  terra::ext(r) <- terra::ext(41, 49, 29, 34)
  
  #--- Handle 0–360° longitude if needed
  ext_r <- terra::ext(r)
  if (ext_r$xmin >= 0 && ext_r$xmax > 180) {
    message("🧭 Detected 0–360° longitude range — shifting to –180–180...")
    r <- terra::shift(r, dx = -180)
  }
  
  #--- Ensure CRS
  terra::crs(r) <- "EPSG:4326"
  
  #--- Convert to raster brick
  brick <- raster::brick(r)
  n_layers <- raster::nlayers(brick)
  
  #--- Convert to dataframe
  data <- as.data.frame(brick, xy = TRUE)
  
  #--- Generate dynamic date names
  dates <- c("x", "y", as.character(seq.Date(
    from = as.Date(paste0(year, "-01-01")),
    by   = "day",
    length.out = n_layers
  )))
  colnames(data) <- dates
  
  #--- Reshape to long format
  data <- data %>%
    tidyr::gather(dates, precipitation, 3:ncol(data)) %>% # change of variable of interest
    tidyr::separate(dates, into = c("year", "month", "day"), sep = "-")
  
  #--- Convert to epiweeks
  data$date <- as.Date(with(data, paste(year, month, day, sep = "-")))
  data$epiweek <- epiweek(data$date)
  
  #--- Aggregate daily → weekly
  data <- data %>%
    group_by(x, y, year, epiweek) %>%
    summarise(
      prec_me = mean(precipitation, na.rm = TRUE),  # change of variable of interest
      prec_mn = min(precipitation, na.rm = TRUE),
      prec_mx = max(precipitation, na.rm = TRUE),
      .groups = "drop"
    )
  
  #--- Load shapefile
  irq <- read_sf(shapefile_path, stringsAsFactors = FALSE) %>%
    st_transform(4326)
  
  #--- Snap coordinates to 0.25° ERA5 grid
  snap_to_grid <- function(x, res = 0.25) round(x / res) * res
  data$lon_grid <- snap_to_grid(data$x)
  data$lat_grid <- snap_to_grid(data$y)
  
  #--- Create grid polygons
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
    mapply(make_cell, data$lon_grid, data$lat_grid, SIMPLIFY = FALSE),
    crs = 4326
  )
  
  #--- Find intersections with Iraq polygons
  overlaps <- st_intersects(era5_cells, irq)
  data$GID_1 <- sapply(overlaps, function(idx) {
    if (length(idx) == 0) return(NA)
    paste(irq$GID_1[idx], collapse = ",")
  })
  
  #--- Expand and clean
  data <- data %>%
    tidyr::separate_rows(GID_1, sep = ",") %>%
    filter(!is.na(GID_1))
  
  #--- Aggregate by region/week
  data_processed <- data %>%
    group_by(year, epiweek, GID_1) %>%
    summarise(
      prec_me = mean(prec_me, na.rm = TRUE),  # change of variable of interest
      prec_mn = mean(prec_mn, na.rm = TRUE),
      prec_mx = mean(prec_mx, na.rm = TRUE),
      .groups = "drop"
    )
  
  message(paste("✅ Finished processing variable data for", year))
  return(data_processed)
}

# Process years 
data_2022 <- process_variable_year(2022)
data_2023 <- process_variable_year(2023)
data_2024 <- process_variable_year(2024)
data_2025 <- process_variable_year(2025)

# Combine all
irq_climate <- readRDS("01_data_processing/climate/irq_climate.rds")
irq_climate2 <- bind_rows(data_2022, data_2023, data_2024, data_2025) %>% distinct()
irq_climate <- left_join(irq_climate, irq_climate2)

# Save once at the end
saveRDS(irq_climate, "01_data_processing/climate/irq_climate.rds")

# Clear the netCDF folder 
file.remove(list.files("01_data_processing/climate/ncs", full.names = TRUE))


