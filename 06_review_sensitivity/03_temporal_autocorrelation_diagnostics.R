# Reviewer comment 3 diagnostics
# Examine residual temporal dependence beyond AR(1) for the MAIN MODEL.
#
# This does NOT force higher-order AR terms into the model.
# Instead, it shows whether there is empirical residual dependence at lags > 1
# after fitting the selected GAM.

suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
  library(readr)
  library(ggplot2)
})

out_dir <- "06_review_sensitivity/temporal_autocorrelation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# Load data and selected covariates
# ----------------------
df <- readRDS("02_data_harmonisation/harmonised_set.rds") %>%
  na.omit() %>%
  arrange(GID_1, year, epiweek) %>%
  mutate(GID_1 = factor(GID_1))

top_covars <- read_csv("03_prelim_analysis/top_covars.csv", show_col_types = FALSE)
main_covs <- top_covars %>%
  filter(Outcome == "total_cases") %>%
  pull(Feature) %>%
  unique()

# ----------------------
# Refit main total-cases model using the selected weekly covariates
# ----------------------
stepwise_gam <- function(outcome, covs, data) {
  base_formula <- as.formula(
    paste0(outcome, " ~ s(epiweek, bs='cc', k=12) + s(GID_1, bs='re')")
  )
  best_model <- bam(base_formula, data = data, family = nb(), discrete = TRUE)
  best_aic <- AIC(best_model)
  selected <- character(0)
  remaining <- covs
  improved <- TRUE

  while (improved && length(remaining) > 0) {
    improved <- FALSE
    best_candidate <- NULL
    best_candidate_model <- NULL
    best_candidate_aic <- best_aic

    for (cov in remaining) {
      test_formula <- as.formula(
        paste0(
          outcome,
          " ~ s(epiweek, bs='cc', k=12) + s(GID_1, bs='re') + ",
          paste(c(selected, cov), collapse = " + ")
        )
      )
      fit <- try(bam(test_formula, data = data, family = nb(), discrete = TRUE), silent = TRUE)
      if (inherits(fit, "try-error")) next
      fit_aic <- AIC(fit)
      if (fit_aic < best_candidate_aic - 2) {
        best_candidate <- cov
        best_candidate_model <- fit
        best_candidate_aic <- fit_aic
      }
    }

    if (!is.null(best_candidate)) {
      selected <- c(selected, best_candidate)
      remaining <- setdiff(remaining, best_candidate)
      best_model <- best_candidate_model
      best_aic <- best_candidate_aic
      improved <- TRUE
    }
  }

  list(model = best_model, selected_covars = selected, final_aic = best_aic)
}

main_fit <- stepwise_gam("total_cases", main_covs, df)
mod <- main_fit$model
resid_p <- residuals(mod, type = "pearson")

capture.output(summary(mod), file = file.path(out_dir, "main_model_summary.txt"))

# ----------------------
# Median within-region ACF by lag
# ----------------------
compute_group_acf_lags <- function(residuals, groups, lag_max = 6) {
  tmp <- tibble(res = residuals, id = groups)

  bind_rows(lapply(1:lag_max, function(L) {
    acf_vals <- tmp %>%
      group_by(id) %>%
      summarise(
        acf_lag = {
          r <- res[!is.na(res)]
          if (length(r) <= (L + 8)) return(NA_real_)
          acf(r, plot = FALSE, lag.max = lag_max)$acf[L + 1]
        },
        .groups = "drop"
      ) %>%
      pull(acf_lag)

    tibble(
      lag = L,
      median_acf = median(acf_vals, na.rm = TRUE),
      mean_acf = mean(acf_vals, na.rm = TRUE)
    )
  }))
}

acf_tbl <- compute_group_acf_lags(resid_p, df$GID_1, lag_max = 6)
write_csv(acf_tbl, file.path(out_dir, "median_within_region_acf_lags1to6.csv"))

p_acf <- ggplot(acf_tbl, aes(x = lag, y = median_acf)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Median within-region Pearson residual ACF by lag",
    x = "Lag (weeks)",
    y = "Median ACF"
  )

ggsave(file.path(out_dir, "median_within_region_acf_lags1to6.png"), p_acf, width = 7, height = 4.5)

# ----------------------
# Pooled residual ACF/PACF as a simple visual supplement
# ----------------------
png(file.path(out_dir, "pooled_residual_acf_pacf.png"), width = 1200, height = 500)
par(mfrow = c(1, 2))
acf(resid_p, lag.max = 12, main = "Pooled Pearson residual ACF")
pacf(resid_p, lag.max = 12, main = "Pooled Pearson residual PACF")
dev.off()

# ----------------------
# Text summary for rebuttal drafting
# ----------------------
summary_lines <- c(
  paste0("Selected covariates in refit: ", paste(main_fit$selected_covars, collapse = ", ")),
  paste0("Final AIC: ", round(main_fit$final_aic, 3)),
  paste0("Deviance explained: ", round(summary(mod)$dev.expl, 4)),
  paste0("Median lag-1 ACF: ", round(acf_tbl$median_acf[acf_tbl$lag == 1], 4)),
  paste0("Median lag-2 ACF: ", round(acf_tbl$median_acf[acf_tbl$lag == 2], 4)),
  paste0("Median lag-3 ACF: ", round(acf_tbl$median_acf[acf_tbl$lag == 3], 4)),
  paste0("Median lag-4 ACF: ", round(acf_tbl$median_acf[acf_tbl$lag == 4], 4))
)
writeLines(summary_lines, file.path(out_dir, "autocorrelation_summary_for_rebuttal.txt"))

message("Temporal autocorrelation diagnostics complete. Outputs saved to: ", out_dir)
