library(terra)

# Path to your soil moisture folder
input_dir <- "01_data_processing/climate/ncs"

# Get the first NetCDF file
files <- list.files(input_dir, full.names = TRUE, pattern = "\\.nc$")
first_file <- files[1]

cat("📄 Checking file:", basename(first_file), "\n")

# Read it in
r <- terra::rast(first_file)

# Print metadata
cat("\n--- Basic Info ---\n")
print(r)

cat("\n--- Layer Names ---\n")
print(names(r))

cat("\n--- Number of layers ---\n")
print(terra::nlyr(r))

cat("\n--- Extent ---\n")
print(terra::ext(r))

cat("\n--- CRS ---\n")
print(terra::crs(r))

cat("\n--- Variable summary (first few cells) ---\n")
print(terra::minmax(r))

# Optional: check attributes directly
cat("\n--- File attributes ---\n")
print(terra::describe(r))
