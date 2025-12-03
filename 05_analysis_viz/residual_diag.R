library(ggplot2)
library(dplyr)
library(mgcv)

# ---------------------------------------------
# Residual diagnostics
# ---------------------------------------------
plot_residual_diagnostics <- function(outcome_name, results_list){
  gam_model <- results_list[[outcome_name]]$best_gam
  if(is.null(gam_model)) stop("No GAM model available for this outcome")
  
  res <- residuals(gam_model, type="pearson")
  fitted <- fitted(gam_model)
  df_res <- data.frame(fitted=fitted, residual=res)
  
  p1 <- ggplot(df_res, aes(x=fitted, y=residual)) +
    geom_point(alpha=0.5) +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(title=paste0(outcome_name, ": Residuals vs Fitted"), x="Fitted", y="Pearson residuals") +
    theme_minimal(base_size=14)
  
  # QQ plot
  qq <- qqnorm(res, plot.it=FALSE)
  df_qq <- data.frame(theoretical=qq$x, sample=qq$y)
  p2 <- ggplot(df_qq, aes(x=theoretical, y=sample)) +
    geom_point(alpha=0.5) +
    geom_abline(intercept=0, slope=1, linetype="dashed") +
    labs(title=paste0(outcome_name, ": QQ plot of residuals"), x="Theoretical quantiles", y="Sample quantiles") +
    theme_minimal(base_size=14)
  
  return(list(residuals_vs_fitted = p1, qq_plot = p2))
}

res_plots <- plot_residual_diagnostics("total_cases", results_list)
res_plots$residuals_vs_fitted
res_plots$qq_plot
ggpubr::ggarrange(res_plots$residuals_vs_fitted,
                  res_plots$qq_plot)
