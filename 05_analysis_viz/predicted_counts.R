library(ggplot2)
library(dplyr)
library(mgcv)

# ---------------------------------------------
# Predicted case counts 
# ---------------------------------------------

plot_predicted_counts <- function(outcome_name, results_list, df, n_points = 100){
  
  gam_model <- results_list[[outcome_name]]$best_gam
  if(is.null(gam_model)) stop("No GAM model available for this outcome")
  
  # Create sequence over epiweek
  x_seq <- seq(min(df$epiweek, na.rm=TRUE), max(df$epiweek, na.rm=TRUE), length.out = n_points)
  
  # Build newdata: hold all other covariates constant
  newdata <- data.frame(matrix(ncol=length(names(gam_model$model)), nrow=n_points))
  colnames(newdata) <- names(gam_model$model)
  
  for(col in names(newdata)){
    if(col == "epiweek") newdata[[col]] <- x_seq
    else if(is.factor(df[[col]])) newdata[[col]] <- factor(levels(df[[col]])[1], levels=levels(df[[col]]))
    else newdata[[col]] <- mean(df[[col]], na.rm=TRUE)
  }
  
  # Predict
  pred_counts <- predict(gam_model, newdata=newdata, type="response", se.fit=TRUE)
  df_plot <- data.frame(
    epiweek = x_seq,
    predicted = pred_counts$fit,
    lower = pred_counts$fit - 2*pred_counts$se.fit,
    upper = pred_counts$fit + 2*pred_counts$se.fit
  )
  
  # Plot
  ggplot(df_plot, aes(x=epiweek, y=predicted)) +
    geom_line(color="steelblue", linewidth=1.2) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, fill="steelblue") +
    labs(title=paste0(outcome_name, ": Predicted counts over time"),
         x="Epiweek", y="Predicted counts") +
    theme_minimal(base_size=14)
}

# ---- Call it outside the function ----
plot_predicted_counts("total_cases", results_list, df)


