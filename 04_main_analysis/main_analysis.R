# Running the main analysis on the harmonised dataset 

# ================================================
# Main analysis: Stepwise GAM with detailed summary
# ================================================

library(mgcv)
library(dplyr)
library(glmmTMB)
library(purrr)
library(tidyr)

# ----------------------
# Load data
# ----------------------
df <- readRDS("02_data_harmonisation/harmonised_set.rds")
df$GID_1 <- factor(df$GID_1)
df <- na.omit(df)

top_covars <- read.csv("03_prelim_analysis/top_covars.csv", stringsAsFactors = FALSE)

outcomes <- c("total_cases","male_cases","female_cases",
              "child_cases","adult_cases","elderly_cases")

# ----------------------
# Helper: compute median within-group AR1
# ----------------------
compute_group_acf1 <- function(residuals, groups){
  df_tmp <- data.frame(res=residuals, id=groups)
  df_tmp %>%
    group_by(id) %>%
    summarise(acf1 = {
      r <- res
      r <- r[!is.na(r)]
      if(length(r) < 10) return(NA_real_)
      acf(r, plot=FALSE)$acf[2]
    }) %>%
    pull(acf1) %>%
    median(na.rm=TRUE)
}

# ----------------------
# STEPWISE AIC GAM selection function
# ----------------------
stepwise_gam <- function(outcome, covs, df){
  
  message("  Starting stepwise AIC for ", outcome)
  
  # Base formula (no covariates)
  base_formula <- as.formula(
    paste0(outcome, " ~ s(epiweek, bs='cc', k=12) + s(GID_1, bs='re')")
  )
  
  # Initialize step failure flag
  theta_warning <- FALSE
  best_model <- try(bam(base_formula, data=df, family=nb(), discrete=TRUE), silent=TRUE)
  if(inherits(best_model, "try-error")) {
    warning("Base model failed for ", outcome)
    return(list(model=NULL, selected_covars=NULL, final_aic=NA, theta_warning=TRUE))
  }
  
  best_aic <- AIC(best_model)
  
  remaining <- covs
  selected <- c()
  improved <- TRUE
  
  while(improved && length(remaining) > 0){
    improved <- FALSE
    best_candidate <- NULL
    best_candidate_aic <- best_aic
    best_candidate_model <- NULL
    
    for(cov in remaining){
      test_formula <- as.formula(
        paste0(
          outcome, " ~ s(epiweek, bs='cc', k=12) + s(GID_1, bs='re') + ",
          paste(c(selected, cov), collapse=" + ")
        )
      )
      
      test_model <- try(bam(test_formula, data=df, family=nb(), discrete=TRUE), silent=TRUE)
      test_aic <- if(!inherits(test_model, "try-error")) AIC(test_model) else Inf
      
      if(inherits(test_model, "try-error")) {
        theta_warning <- TRUE
        next
      }
      
      if(test_aic < best_candidate_aic - 2){
        best_candidate <- cov
        best_candidate_aic <- test_aic
        best_candidate_model <- test_model
      }
    }
    
    if(!is.null(best_candidate)){
      selected <- c(selected, best_candidate)
      remaining <- setdiff(remaining, best_candidate)
      best_aic <- best_candidate_aic
      best_model <- best_candidate_model
      improved <- TRUE
      message("    + Added: ", best_candidate,
              "  (AIC=", round(best_aic,1), ")")
    }
  }
  
  return(list(
    model = best_model,
    selected_covars = selected,
    final_aic = best_aic,
    theta_warning = theta_warning
  ))
}

# ----------------------
# Initialize results and summary table
# ----------------------
results_list <- list()
summary_table <- data.frame()

# ----------------------
# MAIN LOOP
# ----------------------
for(outcome in outcomes){
  
  message("\n\n=======================================")
  message("          OUTCOME: ", outcome)
  message("=======================================\n")
  
  covs <- top_covars$Feature[top_covars$Outcome == outcome]
  
  if(length(covs) < 1){
    warning("Outcome ", outcome, " has <1 covariate. Skipping.")
    next
  }
  
  # STEPWISE model selection
  step_res <- stepwise_gam(outcome, covs, df)
  final_model <- step_res$model
  theta_warning <- step_res$theta_warning
  
  if(is.null(final_model)){
    results_list[[outcome]] <- NULL
    next
  }
  
  # Random-effect variance
  re_var <- final_model$sp[grep("GID_1", names(final_model$sp))]
  random_effect_needed <- ifelse(length(re_var)==0, FALSE, re_var > 1e-4)
  
  # AR1 diagnostic
  pearson_res <- residuals(final_model, type="pearson")
  median_acf1 <- compute_group_acf1(pearson_res, df$GID_1)
  ar1_needed <- median_acf1 > 0.25
  
  # Optional AR1 GLMM
  glmm_fit <- NULL
  glmm_failed <- FALSE
  if(ar1_needed){
    message("  -> AR1 needed: fitting glmmTMB...")
    
    if(length(step_res$selected_covars) == 0){
      fixed_part <- "1"
    } else {
      fixed_part <- paste(step_res$selected_covars, collapse=" + ")
    }
    
    glmm_formula <- as.formula(
      paste0(
        outcome, " ~ ", fixed_part,
        " + (1|GID_1) + ar1(formula = ~ epiweek | GID_1)"
      )
    )
    
    glmm_fit <- try(glmmTMB(glmm_formula, family=nbinom2, data=df),
                    silent=TRUE)
    if(inherits(glmm_fit, "try-error")){
      warning("glmmTMB failed for AR1.")
      glmm_fit <- NULL
      glmm_failed <- TRUE
    }
  }
  
  # Deviance explained
  dev_expl <- summary(final_model)$dev.expl
  
  # ----------------------
  # Parametric coefficients
  # ----------------------
  param_coef <- summary(final_model)$p.table
  if(!is.null(param_coef)){
    param_df <- data.frame(
      Outcome = outcome,
      Term = rownames(param_coef),
      Type = "Parametric",
      Estimate = param_coef[,1],
      Std_Error = param_coef[,2],
      p_value = param_coef[,4],
      stringsAsFactors = FALSE
    )
  } else param_df <- NULL
  
  # ----------------------
  # Smooth terms
  # ----------------------
  smooth_summary <- summary(final_model)$s.table
  if(!is.null(smooth_summary)){
    smooth_df <- data.frame(
      Outcome = outcome,
      Term = rownames(smooth_summary),
      Type = "Smooth",
      edf = smooth_summary[,1],
      F = smooth_summary[,2],
      p_value = smooth_summary[,4],
      stringsAsFactors = FALSE
    )
  } else smooth_df <- NULL
  
  # Combine parametric and smooth into one table
  outcome_summary <- bind_rows(param_df, smooth_df)
  if(!is.null(outcome_summary)){
    outcome_summary <- outcome_summary %>%
      mutate(
        Selected_Covariates = paste(step_res$selected_covars, collapse=", "),
        AIC = step_res$final_aic,
        Deviance_Explained = dev_expl,
        Random_Effect = random_effect_needed,
        RE_Variance = ifelse(length(re_var)==0, NA, re_var),
        AR1_Recommended = ar1_needed,
        Median_ACF1 = median_acf1,
        Theta_Warning = theta_warning,
        GLMM_Failed = glmm_failed
      )
    summary_table <- bind_rows(summary_table, outcome_summary)
  }
  
  # Save full results
  results_list[[outcome]] <- list(
    selected_covariates = step_res$selected_covars,
    best_gam = final_model,
    AIC = step_res$final_aic,
    deviance_explained = dev_expl,
    random_effect_recommended = random_effect_needed,
    ranef_variance = re_var,
    ar1_recommended = ar1_needed,
    median_acf1 = median_acf1,
    theta_warning = theta_warning,
    glmm_model = glmm_fit,
    glmm_failed = glmm_failed,
    parametric_coeffs = param_df,
    smooth_terms = smooth_df
  )
}

message("\n\nSTEPWISE MODELS COMPLETE AND FULL SUMMARY TABLE CREATED.")

# Save summary table as CSV if needed
write.csv(summary_table, "04_main_analysis/stepwise_model_summary.csv", row.names=FALSE)

# Diagnostics and plots for total_cases ---
summary(results_list$total_cases$best_gam)
gam.check(results_list$total_cases$best_gam)
plot(results_list$total_cases$best_gam, pages=1, scheme=1)






