# Reviewer comment 1 sensitivity analysis
# Refit the FINAL main model only at monthly resolution

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ISOweek)
  library(mgcv)
  library(readr)
})

out_dir <- "06_review_sensitivity/monthly_sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# Load weekly harmonised data
# ----------------------
df <- readRDS("02_data_harmonisation/harmonised_set.rds") %>%
  na.omit() %>%
  mutate(
    GID_1 = factor(GID_1),
    week_date = ISOweek::ISOweek2date(sprintf("%s-W%02d-1", year, epiweek)),
    month_date = floor_date(week_date, unit = "month"),
    year_month = format(month_date, "%Y-%m")
  )

# ----------------------
# Monthly aggregation
# Adjust names here if needed to match your harmonised dataset exactly
# ----------------------
monthly_df <- df %>%
  group_by(GID_1, month_date, year_month) %>%
  summarise(
    total_cases = sum(total_cases, na.rm = TRUE),
    
    # final selected predictors
    scPDSI = mean(scpdsi, na.rm = TRUE),   # monthly variable
    IDPs = mean(idp_no, na.rm = TRUE),       # quarterly variable repeated within quarter
    
    .groups = "drop"
  ) %>%
  mutate(
    month = lubridate::month(month_date),
    GID_1 = factor(GID_1)
  )

write_csv(monthly_df, file.path(out_dir, "monthly_refit_dataset.csv"))

# ----------------------
# Refit final model at monthly resolution
# ----------------------
monthly_model <- gam(
  total_cases ~ s(month, bs = "cc", k = 12) + s(GID_1, bs = "re") + scPDSI + IDPs,
  data = monthly_df,
  family = nb(),
  method = "REML"
)

# ----------------------
# Save outputs
# ----------------------
capture.output(
  summary(monthly_model),
  file = file.path(out_dir, "monthly_final_model_summary.txt")
)

model_metrics <- tibble(
  model = "monthly_refit_final_model",
  AIC = AIC(monthly_model),
  deviance_explained = summary(monthly_model)$dev.expl
)

write_csv(model_metrics, file.path(out_dir, "monthly_final_model_metrics.csv"))

coef_table <- as.data.frame(summary(monthly_model)$p.table) %>%
  tibble::rownames_to_column("term")

write_csv(coef_table, file.path(out_dir, "monthly_final_model_coefficients.csv"))

message("Done. Outputs saved to: ", out_dir)
