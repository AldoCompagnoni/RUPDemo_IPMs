# Function to get survival probabilities and confidence intervals
fun_sur_prob_ci <- function(age_class) {
  # Subset the data for the given age class
  data_subset <- subset(df_suv, age == age_class)
  
  # Check if the subset has more than one observation
  if (nrow(data_subset) < 2) {
    # If there are fewer than 2 rows, skip this age class (or return NA)
    return(data.frame(
      age = age_class,
      mean_survival = NA,
      lower_ci      = NA,
      upper_ci      = NA,
      n_obs         = 0  # Add number of observations
    ))
  }
  
  # Fit the logistic regression model
  model <- glm(survives ~ 1, family = 'binomial', data = data_subset)
  
  # Get confidence intervals for the model's intercept (as a vector)
  ci <- confint(model)
  
  # Extract lower and upper bounds for the intercept
  lower_ci <- boot::inv.logit(ci[1])  # 2.5% percentile
  upper_ci <- boot::inv.logit(ci[2])  # 97.5% percentile
  
  # Mean survival probability (logit to probability transformation)
  mean_survival <- boot::inv.logit(coef(model))
  
  # Get the number of observations in this age class
  n_obs <- nrow(data_subset)
  
  # Return a data frame with age class, mean survival, 
  #  confidence intervals, and number of observations
  return(data.frame(
    age           = age_class,
    mean_survival = mean_survival,
    lower_ci      = lower_ci,
    upper_ci      = upper_ci,
    n_obs         = n_obs
  ))
}