# function to plot your survival data "binned" per year (instead of "jittered")

# Arguments: data frame, number of bins, 
# UNQUOTED name of size variable, 
# UNQUOTED name of response variable

df_binned_prop_year <- function(ii, df_in, n_bins, siz_var, rsp_var, years) {
  
  # Subset data for the specific year
  df <- subset(df_in, year == years$year[ii])
  
  if (nrow(df) == 0) return(NULL)  # Check if there is any data for the year
  
  size_var <- deparse(substitute(siz_var))
  resp_var <- deparse(substitute(rsp_var))
  
  # Binning the size variable
  h    <- (max(df[,size_var], na.rm = TRUE) - min(df[,size_var], na.rm = TRUE)) / n_bins
  lwr  <- min(df[,size_var], na.rm = TRUE) + (h * c(0:(n_bins - 1)))
  upr  <- lwr + h
  mid  <- lwr + (1/2 * h)
  
  # Standard error of a Bernoulli process with handling for small sample sizes
  se_bern <- function(x, lwr_upr) {
    if (length(x) == 0) return(c(NA, NA))  # Return NA if no observations in the bin
    surv_n <- sum(x, na.rm = TRUE)
    tot_n  <- length(x)
    if (tot_n == 0) return(c(NA, NA))  # Return NA if there are no observations in the bin
    binom.confint(surv_n, tot_n, methods = c("wilson"))[, lwr_upr]
  }
  
  # Binned proportion calculation
  binned_prop <- function(lwr_x, upr_x, response, lwr_upr) {
    id  <- which(df[, size_var] > lwr_x & df[, size_var] < upr_x) 
    tmp <- df[id,]
    
    if (response == 'prob') {   
      return(sum(tmp[, resp_var], na.rm = TRUE) / nrow(tmp)) 
    }
    if (response == 'n_size') { 
      return(nrow(tmp)) 
    }
    if (response == 'se') {    
      return(se_bern(tmp[, resp_var], lwr_upr)) 
    }
  }
  
  # Get binned proportions and standard errors for the year
  y_binned   <- Map(binned_prop, lwr, upr, 'prob', 'mean') %>% unlist
  x_binned   <- mid
  y_n_size   <- Map(binned_prop, lwr, upr, 'n_size') %>% unlist
  y_se_lwr   <- Map(binned_prop, lwr, upr, 'se', 'lower') %>% unlist
  y_se_upr   <- Map(binned_prop, lwr, upr, 'se', 'upper') %>% unlist
  
  # Ensure that all vectors have the same length (fill NA where necessary)
  max_len <- max(length(y_binned), length(y_n_size), length(y_se_lwr), length(y_se_upr))
  
  # Fill the shorter vectors with NA to ensure equal length
  y_binned   <- c(y_binned, rep(NA, max_len - length(y_binned)))
  y_n_size   <- c(y_n_size, rep(NA, max_len - length(y_n_size)))
  y_se_lwr   <- c(y_se_lwr, rep(NA, max_len - length(y_se_lwr)))
  y_se_upr   <- c(y_se_upr, rep(NA, max_len - length(y_se_upr)))
  x_binned   <- c(x_binned, rep(NA, max_len - length(x_binned)))
  
  # Create and return the data frame for this year
  data.frame(
    xx   = x_binned,
    yy   = y_binned,
    nn   = y_n_size,
    lwr  = y_se_lwr,
    upr  = y_se_upr
  ) %>% 
    mutate(year = years$year[ii]) %>% 
    setNames(c(size_var, resp_var, 'n_size', 'lwr', 'upr', 'year')) %>%
    drop_na()
}
