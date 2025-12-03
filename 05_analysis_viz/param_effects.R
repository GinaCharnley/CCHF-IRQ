library(ggplot2)
library(dplyr)

# ---------------------------------------------
# Parameter effects 
# ---------------------------------------------

plot_gam_effects <- function(outcome_name, results_list, n_points = 100) {
  
  if(!(outcome_name %in% names(results_list))){
    stop("Outcome not found in results_list")
  }
  
  model_info <- results_list[[outcome_name]]
  gam_model <- model_info$best_gam
  if(is.null(gam_model)){
    stop("No GAM model available for this outcome")
  }
  
  smooth_plots <- list()
  
  for(smooth_obj in gam_model$smooth){
    st <- smooth_obj$label
    term_var <- smooth_obj$term
    # Check if numeric
    if(is.numeric(df[[term_var]])){
      # numeric smooth: line + ribbon
      x_seq <- seq(min(df[[term_var]], na.rm=TRUE), max(df[[term_var]], na.rm=TRUE), length.out = n_points)
      # minimal newdata
      newdata <- data.frame(matrix(ncol = length(names(gam_model$model)), nrow = n_points))
      colnames(newdata) <- names(gam_model$model)
      for(col in names(newdata)){
        if(col == term_var){
          newdata[[col]] <- x_seq
        } else if(is.factor(df[[col]])){
          newdata[[col]] <- factor(levels(df[[col]])[1], levels=levels(df[[col]]))
        } else {
          newdata[[col]] <- mean(df[[col]], na.rm=TRUE)
        }
      }
      pred <- predict(gam_model, newdata=newdata, type="terms", se.fit=TRUE)
      df_plot <- data.frame(
        x = x_seq,
        fit = pred$fit[, st],
        se = pred$se.fit[, st]
      )
      p <- ggplot(df_plot, aes(x=x, y=fit)) +
        geom_line(color="steelblue", linewidth=1.2) +
        geom_ribbon(aes(ymin=fit-2*se, ymax=fit+2*se), alpha=0.2, fill="steelblue") +
        theme_minimal() + 
        labs(title=paste0(outcome_name, ": Smooth effect of ", st),
             x=term_var, y="Effect on log(counts)") +
        theme(text = element_text(face = "bold"))
      smooth_plots[[st]] <- p
    } else if(is.factor(df[[term_var]])){
      # factor/re random-effect: plot as points
      re_coef <- coef(gam_model)[grep(term_var, names(coef(gam_model)))]
      df_plot <- data.frame(
        level = names(re_coef),
        estimate = as.numeric(re_coef)
      )
      p <- ggplot(df_plot, aes(x=reorder(level, estimate), y=estimate)) +
        geom_point(color="steelblue", size=4) +
        coord_flip() + theme_minimal() + 
        labs(title=paste0(outcome_name, ": Random effect ", st),
             x=term_var, y="Estimated effect") +
        theme(text = element_text(face = "bold"))
      smooth_plots[[st]] <- p
    }
  }
  
  # Parametric coefficients
  param_df <- model_info$parametric_coeffs
  param_plot <- NULL
  if(!is.null(param_df) && nrow(param_df) > 0){
    param_df <- param_df %>%
      filter(Term != "(Intercept)") %>%
      mutate(lower = Estimate - 1.96 * Std_Error,
             upper = Estimate + 1.96 * Std_Error)
    if(nrow(param_df) > 0){
      param_plot <- ggplot(param_df, aes(x=reorder(Term, Estimate), y=Estimate)) +
        geom_point(color="darkred", size=4) +
        geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, size = 2, color="darkred", alpha = .4) +
        coord_flip() + theme_minimal() + 
        labs(title=paste0(outcome_name, ": Parametric effect sizes"),
             x="Covariate", y="Effect (log scale)") +
        theme(text = element_text(face = "bold"))
    }
  }
  
  return(list(
    smooth_plots = smooth_plots,
    param_plot = param_plot
  ))
}


# -------------------------------
# Example usage:
# -------------------------------
# To plot for "total_cases":
plots <- plot_gam_effects("total_cases", results_list)
plots$smooth_plots[[1]]  # first smooth term
plots$smooth_plots[[2]]  # first smooth term
plots$param_plot          # parametric coefficients

ggpubr::ggarrange(plots$smooth_plots[[1]],  # first smooth term
                  plots$smooth_plots[[2]],  # first smooth term
                  plots$param_plot)

