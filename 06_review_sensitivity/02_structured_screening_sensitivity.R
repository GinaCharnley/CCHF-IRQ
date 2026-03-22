# Reviewer comment 2 sensitivity analysis
# Screening that respects seasonality and region effects for the MAIN OUTCOME.
#
# Strategy:
#   1. Fit the same base GAM used in the main analysis:
#      total_cases ~ s(epiweek, bs='cc') + s(GID_1, bs='re')
#   2. Add each candidate covariate individually and rank by delta AIC and delta deviance explained.
#   3. Run the same forward AIC selection, but starting from the full candidate covariate set
#      instead of the glmnet-screened subset.
#
# This directly answers whether scPDSI / IDPs stay competitive once the time-space structure
# is already in the baseline model.

suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
  library(readr)
  library(purrr)
})

out_dir <- "06_review_sensitivity/structured_screening"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# Load data
# ----------------------
df <- readRDS("02_data_harmonisation/harmonised_set.rds") %>%
  na.omit() %>%
  mutate(GID_1 = factor(GID_1))

outcomes <- c("total_cases", "male_cases", "female_cases", "child_cases", "adult_cases", "elderly_cases")

# Keep candidate covariates broad.
# This mirrors your pipeline but avoids the outcome columns and identifiers.
exclude <- c("GID_1", "year", "epiweek", outcomes)
all_covariates <- setdiff(names(df), exclude)

# ----------------------
# Base model for the main outcome
# ----------------------
base_formula <- total_cases ~ s(epiweek, bs = "cc", k = 12) + s(GID_1, bs = "re")
base_model <- bam(base_formula, data = df, family = nb(), discrete = TRUE)
base_aic <- AIC(base_model)
base_dev <- summary(base_model)$dev.expl

capture.output(summary(base_model), file = file.path(out_dir, "base_gam_summary.txt"))

# ----------------------
# Single-covariate additions on top of the structured base model
# ----------------------
rank_one_by_one <- map_dfr(all_covariates, function(cov) {
  f <- as.formula(
    paste0("total_cases ~ s(epiweek, bs = 'cc', k = 12) + s(GID_1, bs = 're') + ", cov)
  )
  fit <- try(bam(f, data = df, family = nb(), discrete = TRUE), silent = TRUE)

  if (inherits(fit, "try-error")) {
    return(tibble(
      covariate = cov,
      AIC = NA_real_,
      delta_AIC = NA_real_,
      deviance_explained = NA_real_,
      delta_deviance = NA_real_,
      converged = FALSE
    ))
  }

  tibble(
    covariate = cov,
    AIC = AIC(fit),
    delta_AIC = AIC(fit) - base_aic,
    deviance_explained = summary(fit)$dev.expl,
    delta_deviance = summary(fit)$dev.expl - base_dev,
    converged = TRUE
  )
}) %>%
  arrange(delta_AIC)

write_csv(rank_one_by_one, file.path(out_dir, "one_by_one_structured_screening.csv"))

# ----------------------
# Forward AIC selection from FULL candidate list
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
        best_candidate_aic <- fit_aic
        best_candidate_model <- fit
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

full_structured <- stepwise_gam("total_cases", all_covariates, df)

capture.output(summary(full_structured$model),
               file = file.path(out_dir, "structured_forward_selection_summary.txt"))

write_csv(
  tibble(
    analysis = c("original_weekly_prescreened", "structured_forward_from_full_set"),
    selected_covariates = c(
      paste(read_csv("03_prelim_analysis/top_covars.csv", show_col_types = FALSE) %>%
              filter(Outcome == "total_cases") %>% pull(Feature) %>% unique(), collapse = ", "),
      paste(full_structured$selected_covars, collapse = ", ")
    )
  ),
  file.path(out_dir, "selected_covariate_comparison.csv")
)

# Optional: check original reviewer-targeted covariates directly
focus_covs <- c("scpdsi", "idp_no")
focus_tbl <- rank_one_by_one %>% filter(covariate %in% focus_covs)
write_csv(focus_tbl, file.path(out_dir, "structured_screening_focus_covariates.csv"))

message("Structured screening sensitivity complete. Outputs saved to: ", out_dir)
