library(ggplot2)
library(dplyr)
library(mgcv)

# ---------------------------------------------
# Marginal effects
# ---------------------------------------------

plot_marginal_fixed_effect <- function(outcome_name, results_list, df, covariate, n_points = 100){
  
  gam_model <- results_list[[outcome_name]]$best_gam
  if(is.null(gam_model)) stop("No GAM model available for this outcome")
  
  if(!(covariate %in% names(df))) stop("Covariate not found in data")
  
  # Sequence of values for covariate
  x_seq <- seq(min(df[[covariate]], na.rm=TRUE),
               max(df[[covariate]], na.rm=TRUE),
               length.out = n_points)
  
  # Build newdata: hold all other covariates constant
  newdata <- data.frame(matrix(ncol = length(names(gam_model$model)), nrow = n_points))
  colnames(newdata) <- names(gam_model$model)
  
  for(col in names(newdata)){
    if(col == covariate) {
      newdata[[col]] <- x_seq
    } else if(is.factor(df[[col]])){
      newdata[[col]] <- factor(levels(df[[col]])[1], levels=levels(df[[col]]))
    } else {
      newdata[[col]] <- mean(df[[col]], na.rm=TRUE)
    }
  }
  
  # Predict
  pred <- predict(gam_model, newdata=newdata, type="response", se.fit=TRUE)
  
  df_plot <- data.frame(
    covariate_value = x_seq,
    predicted = pred$fit,
    lower = pred$fit - 2*pred$se.fit,
    upper = pred$fit + 2*pred$se.fit
  )
  
  ggplot(df_plot, aes(x=covariate_value, y=predicted)) +
    geom_line(color="steelblue", linewidth=1.2) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill="steelblue") +
    labs(title=paste0(outcome_name, ": Marginal effect of ", covariate),
         x=covariate, y="Predicted counts") +
    theme_minimal(base_size=14)
}

# ------------------------
# Example usage
# ------------------------
plot_marginal_fixed_effect("total_cases", results_list, df, "scpdsi")
plot_marginal_fixed_effect("total_cases", results_list, df, "idp_no")

ggpubr::ggarrange(
plot_marginal_fixed_effect("total_cases", results_list, df, "scpdsi"),
plot_marginal_fixed_effect("total_cases", results_list, df, "idp_no")
)



