# Reviewer comment 4 analysis
# Leave-one-year-out (LOYO) robustness analysis for the MAIN OUTCOME.
#
# For each held-out year:
#   1. train on all remaining years
#   2. rerun forward AIC selection for total_cases
#   3. predict the held-out year
#   4. save selected covariates, effect estimates, train AIC/deviance,
#      and held-out predictive metrics

suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
  library(readr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

out_dir <- "06_review_sensitivity/leave_one_year_out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------
# Load data and candidate covariates
# ----------------------
df <- readRDS("02_data_harmonisation/harmonised_set.rds") %>%
  na.omit() %>%
  mutate(GID_1 = factor(GID_1))

top_covars <- read_csv("03_prelim_analysis/top_covars.csv", show_col_types = FALSE)
main_covs <- top_covars %>%
  filter(Outcome == "total_cases") %>%
  pull(Feature) %>%
  unique()

years <- sort(unique(df$year))

# ----------------------
# Forward AIC function
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

# ----------------------
# Metrics helper
# ----------------------
calc_metrics <- function(obs, pred) {
  tibble(
    rmse = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    mae = mean(abs(obs - pred), na.rm = TRUE),
    spearman = suppressWarnings(cor(obs, pred, method = "spearman", use = "complete.obs"))
  )
}

# ----------------------
# Run LOYO
# ----------------------
loyo_results <- map(years, function(yr_out) {
  train <- df %>% filter(year != yr_out)
  test  <- df %>% filter(year == yr_out)

  fit_obj <- stepwise_gam("total_cases", main_covs, train)
  mod <- fit_obj$model

  # Predictions on held-out year
  pred <- predict(mod, newdata = test, type = "response")
  metrics <- calc_metrics(test$total_cases, pred)

  # Parametric coefficients for selected covariates only
  ptab <- as.data.frame(summary(mod)$p.table)
  ptab$term <- rownames(ptab)
  rownames(ptab) <- NULL
  
  # Detect correct column names (t vs z)
  stat_col <- intersect(c("z value", "t value"), colnames(ptab))[1]
  p_col <- intersect(c("Pr(>|z|)", "Pr(>|t|)"), colnames(ptab))[1]
  
  coef_tbl <- ptab %>%
    rename(
      estimate = Estimate,
      std_error = `Std. Error`
    ) %>%
    mutate(
      stat = .data[[stat_col]],
      p_value = .data[[p_col]]
    ) %>%
    select(term, estimate, std_error, stat, p_value) %>%
    filter(term %in% fit_obj$selected_covars) %>%
    mutate(held_out_year = yr_out)
  
  list(
    overview = tibble(
      held_out_year = yr_out,
      selected_covariates = paste(fit_obj$selected_covars, collapse = ", "),
      train_AIC = fit_obj$final_aic,
      train_deviance_explained = summary(mod)$dev.expl,
      rmse = metrics$rmse,
      mae = metrics$mae,
      spearman = metrics$spearman
    ),
    coefs = coef_tbl,
    predictions = tibble(
      held_out_year = yr_out,
      GID_1 = test$GID_1,
      year = test$year,
      epiweek = test$epiweek,
      observed = test$total_cases,
      predicted = pred
    )
  )
})

overview_tbl <- bind_rows(map(loyo_results, "overview"))
coef_tbl     <- bind_rows(map(loyo_results, "coefs"))
pred_tbl     <- bind_rows(map(loyo_results, "predictions"))

write_csv(overview_tbl, file.path(out_dir, "loyo_overview.csv"))
write_csv(coef_tbl, file.path(out_dir, "loyo_coefficients.csv"))
write_csv(pred_tbl, file.path(out_dir, "loyo_predictions.csv"))

# Stability summary for the two headline covariates
focus_covs <- c("scpdsi", "idp_no")
focus_tbl <- coef_tbl %>%
  filter(term %in% focus_covs) %>%
  arrange(term, held_out_year)
write_csv(focus_tbl, file.path(out_dir, "loyo_focus_covariates.csv"))

# Plot coefficient stability across held-out years
if (nrow(focus_tbl) > 0) {
  p_coef <- ggplot(focus_tbl, aes(x = factor(held_out_year), y = estimate, group = term)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point() +
    geom_line() +
    facet_wrap(~ term, scales = "free_y") +
    theme_minimal() +
    labs(x = "Held-out year", y = "Coefficient estimate", title = "LOYO coefficient stability")

  ggsave(file.path(out_dir, "loyo_focus_covariate_stability.png"), p_coef, width = 7, height = 4.5)
}

# Plot observed vs predicted for held-out years
p_pred <- ggplot(pred_tbl, aes(x = observed, y = predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap(~ held_out_year, scales = "free") +
  theme_minimal() +
  labs(title = "LOYO observed vs predicted: total cases", x = "Observed weekly cases", y = "Predicted weekly cases")

ggsave(file.path(out_dir, "loyo_observed_vs_predicted.png"), p_pred, width = 8, height = 6)

message("Leave-one-year-out analysis complete. Outputs saved to: ", out_dir)
