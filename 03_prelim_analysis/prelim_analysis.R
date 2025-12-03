# Preliminary analysis on the harmonised dataset 

# Load the data
harmonised_set <- readRDS("02_data_harmonisation/harmonised_set.rds")
harmonised_set <- na.omit(harmonised_set)

# Libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(lubridate)
library(tidyr)
library(glmnet)

# Plot a time-series of the outcome 

# Create a proper weekly date
df <- harmonised_set %>%
  mutate(week_date = ISOweek::ISOweek2date(paste0(year, "-W", sprintf("%02d", epiweek), "-1")))
# "-1" is Monday of the ISO week

# Plot time series faceted by GID_1
ggplot(df, aes(x = week_date, y = total_cases, fill = GID_1)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x="Week",
       y="Total Cases") +
  scale_x_date(date_labels = "%Y-%W", date_breaks = "4 weeks") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(face = "bold"))

# Identify correlations among the covariates 

# Select covariates
covariate_data <- harmonised_set[, 10:41]
# Select only the continuous 
covariate_data <- covariate_data %>% select(-elevation, -holiday_this_week, -holiday_previous_week,
                                            -conf_events)

# Compute correlation matrix
cor_matrix <- cor(covariate_data, use="pairwise.complete.obs", method="pearson")

# Convert to long format for ggplot
cor_long <- melt(cor_matrix)
colnames(cor_long) <- c("Var1", "Var2", "Correlation")

cor_matrix[upper.tri(cor_matrix)] <- NA
cor_long <- melt(cor_matrix, na.rm=TRUE)
colnames(cor_long) <- c("Var1", "Var2", "Correlation")

pretty_names <- c(
  temp_mn = "TemperatureMin",
  temp_mx = "TemperatureMax",
  temp_me = "TemperatureMean",
  hurs_me = "HumidityMean",
  hurs_mn = "HumidityMin",
  hurs_mx = "HumidityMax",
  prec_me = "PrecipMean",
  prec_mn = "PrecipMin",
  prec_mx = "PrecipMax",
  soil_moisture = "SoilMoisture",
  scpdsi = "scPDSI",
  mpi = "MPI",
  pov_headcount = "Poverty Headcount",
  pov_intensity = "Poverty Intensity",
  pov_vulnerability = "Poverty Vulnerability",
  pov_severe = "Poverty Severe",
  idp_no = "IDPs"
)

# Plot heatmap
ggplot(cor_long, aes(x=Var1, y=Var2, fill=Correlation)) +
  geom_tile(color="white") +
  geom_text(aes(label=round(Correlation, 2)), size=3) +
  scale_fill_gradient2(
    low="blue", mid="white", high="red",
    midpoint=0, limits=c(-1, 1),
    name="Pearson\nCorrelation"
  ) +
  scale_x_discrete(labels = function(x) ifelse(x %in% names(pretty_names), pretty_names[x], x)) +
  scale_y_discrete(labels = function(x) ifelse(x %in% names(pretty_names), pretty_names[x], x)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        text = element_text(face="bold")) +
  labs(x = "", y = "")


# Calculate variable importance 

# Define outcomes
outcomes <- c("total_cases", "male_cases", "female_cases", "child_cases", "adult_cases", "elderly_cases")

# Covariates in your dataset (columns 10:35)
covariates <- colnames(df)[10:41]

# Ensure correct data format 
df$GID_1 <- as.factor(df$GID_1)
df <- na.omit(df)
x <- model.matrix(~ ., df[, covariates])
y <- df$total_cases

# Function to extract importance using LASSO and Elastic Net 

get_importance_glmnet <- function(df, outcome, covariates,
                                  alpha_enet = 0.5,
                                  family = "poisson") {
  # Load required package
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required but not installed.")
  }
  
  # Build x matrix and outcome y vector
  x <- model.matrix(~ ., df[, covariates])
  y <- df[[outcome]]
  
  # -------------------------
  # 1. Fit LASSO (alpha = 1)
  # -------------------------
  lasso_fit <- glmnet::cv.glmnet(
    x, y,
    family = family,
    alpha = 1
  )
  
  lasso_coef <- coef(lasso_fit, s = "lambda.min")
  lasso_selected <- rownames(lasso_coef)[lasso_coef[,1] != 0]
  lasso_selected <- lasso_selected[lasso_selected != "(Intercept)"]
  
  # Rank by coefficient magnitude
  lasso_ranking <- sort(abs(lasso_coef[lasso_selected,1]), decreasing = TRUE)
  
  # -------------------------------
  # 2. Fit Elastic Net (alpha = X)
  # -------------------------------
  enet_fit <- glmnet::cv.glmnet(
    x, y,
    family = family,
    alpha = alpha_enet
  )
  
  enet_coef <- coef(enet_fit, s = "lambda.min")
  enet_selected <- rownames(enet_coef)[enet_coef[,1] != 0]
  enet_selected <- enet_selected[enet_selected != "(Intercept)"]
  
  # Rank by coefficient magnitude
  enet_ranking <- sort(abs(enet_coef[enet_selected,1]), decreasing = TRUE)
  
  # -------------------------
  # 3. Compare selections
  # -------------------------
  both_selected <- intersect(lasso_selected, enet_selected)
  only_lasso   <- setdiff(lasso_selected, enet_selected)
  only_enet    <- setdiff(enet_selected, lasso_selected)
  union_all    <- union(lasso_selected, enet_selected)
  
  # -------------------------
  # 4. Return results
  # -------------------------
  return(list(
    outcome = outcome,
    
    lasso = list(
      model = lasso_fit,
      selected = lasso_selected,
      ranking = lasso_ranking,
      lambda_min = lasso_fit$lambda.min
    ),
    
    elastic_net = list(
      model = enet_fit,
      alpha = alpha_enet,
      selected = enet_selected,
      ranking = enet_ranking,
      lambda_min = enet_fit$lambda.min
    ),
    
    comparison = list(
      both_selected = both_selected,
      only_lasso = only_lasso,
      only_enet = only_enet,
      union_all = union_all
    )
  ))
}

run_importance_all_outcomes <- function(df, outcomes, covariates,
                                        alpha_enet = 0.5,
                                        save_plots = FALSE,
                                        plot_dir = "03_prelim_analysis/importance_plots") {
  
  # ---- Pretty labels for features ----
  pretty_labels <- c(
    # climate
    hurs_me       = "HumidityMean",
    hurs_mn       = "HumidityMin",
    hurs_mx       = "HumidityMax",
    temp_me       = "TemperatureMean",
    temp_mn       = "TemperatureMin",
    temp_mx       = "TemperatureMax",
    prec_me       = "PrecipitationMean",
    prec_mn       = "PrecipitationMin",
    prec_mx       = "PrecipitationMax",
    soil_moisture = "SoilMoisture",
    scpdsi        = "scPDSI",
    
    # holidays / events / geography
    elevation             = "Elevation",
    conf_events           = "ConflictEvents",
    holiday_this_week     = "HolidayThisWeek",
    holiday_previous_week = "HolidayPreviousWeek",
    
    # socioeconomics
    mpi               = "MPI",
    pov_headcount     = "PovertyHeadcount",
    pov_intensity     = "PovertyIntensity",
    pov_vulnerability = "PovertyVulnerability",
    pov_severe        = "PovertySevere",
    idp_no            = "IDPs"
  )
  
  if (save_plots && !dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  results_all <- list()
  
  for (outcome in outcomes) {
    message("Running variable importance for outcome: ", outcome)
    
    res <- get_importance_glmnet(
      df         = df,
      outcome    = outcome,
      covariates = covariates,
      alpha_enet = alpha_enet
    )
    
    results_all[[outcome]] <- res
    
    # Intersection set = what you use for top_covariates_df
    stable_vars <- res$comparison$both_selected
    
    if (length(stable_vars) == 0) {
      message("  No variables selected by BOTH LASSO and Elastic Net for outcome: ", outcome)
      next
    }
    
    # -----------------------------
    # Build long-format ranking df
    # -----------------------------
    df_plot <- dplyr::bind_rows(
      data.frame(
        Feature    = names(res$lasso$ranking),
        Importance = as.numeric(res$lasso$ranking),
        Method     = "LASSO",
        stringsAsFactors = FALSE
      ),
      data.frame(
        Feature    = names(res$elastic_net$ranking),
        Importance = as.numeric(res$elastic_net$ranking),
        Method     = paste0("ElasticNet α=", alpha_enet),
        stringsAsFactors = FALSE
      )
    ) %>%
      # 🔍 keep only variables that are in both_selected
      dplyr::filter(Feature %in% stable_vars) %>%
      dplyr::mutate(
        Importance_plot = log1p(Importance),
        Feature_pretty  = dplyr::recode(Feature, !!!pretty_labels)
      )
    
    # If for some reason nothing left after filtering, skip
    if (nrow(df_plot) == 0) {
      message("  No overlapping ranked features for outcome: ", outcome)
      next
    }
    
    # -----------------------------
    # Make the plot
    # -----------------------------
    p <- ggplot(
      df_plot,
      aes(
        x    = reorder(Feature_pretty, Importance_plot),
        y    = Importance_plot,
        fill = Method
      )
    ) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      theme_minimal() +
      labs(
        title = paste("Variable Importance (stable predictors) for", outcome),
        x     = "Feature",
        y     = "Importance (log scale)"
      ) +
      theme(text = element_text(face = "bold")) +
      scale_fill_manual(values = c("dodgerblue", "#FF6F61"))
    
    print(p)
    
    if (save_plots) {
      ggsave(
        filename = file.path(plot_dir, paste0(outcome, "_importance_stable.png")),
        plot     = p,
        width    = 8,
        height   = 6
      )
    }
  }
  
  return(results_all)
}


results_all <- run_importance_all_outcomes(
  df = df,
  outcomes = outcomes,
  covariates = covariates,
  alpha_enet = 0.5,
  save_plots = TRUE,         # set to FALSE if you don’t want files
  plot_dir = "03_prelim_analysis/importance_plots"
)

# Export the top set 
top_covariates_df <- data.frame(
  Outcome = character(),
  Feature = character(),
  Method = character(),
  stringsAsFactors = FALSE
)

for (outcome in names(results_all)) {
  res <- results_all[[outcome]]
  
  # You can choose which set you want: both_selected, union_all, etc.
  top_vars <- res$comparison$both_selected  # most stable predictors
  
  if(length(top_vars) > 0){
    tmp <- data.frame(
      Outcome = outcome,
      Feature = top_vars,
      Method = "Both LASSO + Elastic Net",
      stringsAsFactors = FALSE
    )
    top_covariates_df <- rbind(top_covariates_df, tmp)
  }
}
top_covariates_df <- top_covariates_df %>%
  rowwise() %>%
  mutate(
    LASSO_coef = ifelse(Feature %in% names(results_all[[Outcome]]$lasso$ranking),
                        results_all[[Outcome]]$lasso$ranking[Feature], NA_real_),
    ENET_coef = ifelse(Feature %in% names(results_all[[Outcome]]$elastic_net$ranking),
                       results_all[[Outcome]]$elastic_net$ranking[Feature], NA_real_)
  ) %>%
  ungroup()

readr::write_csv(top_covariates_df, "03_prelim_analysis/top_covars.csv")




