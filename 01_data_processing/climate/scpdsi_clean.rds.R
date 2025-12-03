# Function for transforming the daily PDSI to weekly

# --- Load packages 
library(terra)
library(dplyr)
library(tidyr)
library(sf)
library(lubridate)
library(epiR)

# --- Load raster
r <- terra::rast("01_data_processing/climate/ncs/scPDSI_cru_ts.nc")
terra::crs(r) <- "EPSG:4326"
  
# --- Crop to Iraq region
region_extent <- terra::ext(41, 49, 29, 34)
r <- terra::crop(r, region_extent)

# --- Extract or construct time dimension
time_vals <- terra::time(r)
if (is.null(time_vals)) {
    # fallback: use layer names
    layer_names <- names(r)
    # try to extract YYYY-MM-DD from layer names if possible
    parsed_dates <- as.Date(sub(".*?(\\d{4}-\\d{2}-\\d{2}).*", "\\1", layer_names), format = "%Y-%m-%d")
    if (all(is.na(parsed_dates))) {
      # if parsing fails, assume first layer is 1901-01-01 and monthly sequence
      start_date <- as.Date("1901-01-01")
      parsed_dates <- seq.Date(start_date, by = "month", length.out = terra::nlyr(r))
    }
    time_vals <- parsed_dates
}
  
# --- Convert raster to dataframe
# Get number of layers
n_layers <- terra::nlyr(r)
# Construct monthly dates starting from known first month (1901-01-01)
time_vals <- seq.Date(from = as.Date("1901-01-01"), by = "month", length.out = n_layers)
# Assign as column names
df <- as.data.frame(r, xy = TRUE)
colnames(df) <- c("x", "y", as.character(time_vals))

# Pivot to long format
df_long <- df %>%
    pivot_longer(cols = -c(x, y), names_to = "date", values_to = "scpdsi") %>%
    mutate(date = as.Date(date)) %>%
    mutate(year = lubridate::year(date),
           month = lubridate::month(date)) %>%
    select(x, y, date, year, month, scpdsi)
  
# --- Precompute month -> epiweek mapping 
months_seq <- seq(min(df_long$date), max(df_long$date), by = "month")
month_epi <- lapply(months_seq, function(m) {
    days <- seq(m, m %m+% months(1) - days(1), by = "day")
    tibble(
      year = year(m),
      month = month(m),
      epiweek = unique(epiweek(days))
    )
  }) %>% bind_rows()
  
# --- Join months to raster values to get one row per (x, y, year, epiweek)
df_epi <- df_long %>%
left_join(month_epi, by = c("year", "month"))
  
# --- Read shapefile
irq <- st_read("02_data_harmonisation/gadm41_IRQ_shp/gadm41_IRQ_1.shp", quiet = TRUE) %>% st_transform(4326)

# --- Snap coordinates to 0.5° grid
snap_to_grid <- function(x, res = 0.5) round(x / res) * res
df_epi$lon_grid <- snap_to_grid(df_epi$x)
df_epi$lat_grid <- snap_to_grid(df_epi$y)

# --- Create grid polygons
make_cell <- function(lon, lat, res = 0.5) {
  st_polygon(list(rbind(
    c(lon - res/2, lat - res/2),
    c(lon + res/2, lat - res/2),
    c(lon + res/2, lat + res/2),
    c(lon - res/2, lat + res/2),
    c(lon - res/2, lat - res/2)
  )))
}
grid_cells <- st_sfc(
  mapply(make_cell, df_epi$lon_grid, df_epi$lat_grid, SIMPLIFY = FALSE),
  crs = 4326
)

# --- Spatial join to assign GID_1
overlaps <- st_intersects(grid_cells, irq)
df_epi$GID_1 <- sapply(overlaps, function(idx) {
  if (length(idx) == 0) return(NA)
  paste(irq$GID_1[idx], collapse = ",")
})
df_epi <- df_epi %>%
  separate_rows(GID_1, sep = ",") %>%
  filter(!is.na(GID_1))

# --- Aggregate by GID_1, year, epiweek
df_final <- df_epi %>%
  group_by(GID_1, year, epiweek) %>%
  summarise(scpdsi = mean(scpdsi, na.rm = TRUE), .groups = "drop")

irq_climate <- readRDS("01_data_processing/climate/irq_climate.rds")

irq_climate <- left_join(irq_climate, df_final)

saveRDS(irq_climate, "01_data_processing/climate/irq_climate.rds")

