# function to plot your survival data "binned" per year (instead of "jittered")

# Arguments: data frame, number of bins, 
# UNQUOTED name of size variable, 
# UNQUOTED name of response variable

df_binned_prop_year <- function(ii, df_in, n_bins, siz_var, rsp_var, years){
  
  # make sub-selection of data
  df   <- subset(df_in, year == years$year[ii])
  
  if(nrow(df) == 0) return(NULL)
  
  size_var <- deparse(substitute(siz_var))
  resp_var <- deparse(substitute(rsp_var))
  
  # binned survival probabilities
  h    <- (max(df[,size_var], na.rm = T) - min(df[,size_var], na.rm = T)) / n_bins
  lwr  <-  min(df[,size_var], na.rm = T) + (h * c(0:(n_bins - 1)))
  upr  <- lwr + h
  mid  <- lwr + (1/2 * h)
  
  binned_prop <- function(lwr_x, upr_x, response){
    
    id  <- which(df[,size_var] > lwr_x & df[,size_var] < upr_x) 
    tmp <- df[id,]
    
    if(response == 'prob'){   return(sum (tmp[,resp_var], na.rm = T)/nrow(tmp))}
    if(response == 'n_size'){ return(nrow(tmp))}
    
  }
  
  y_binned <- Map(binned_prop, lwr, upr, 'prob') %>% unlist
  x_binned <- mid
  y_n_size <- Map(binned_prop, lwr, upr, 'n_size') %>% unlist
  
  # output data frame
  data.frame(xx  = x_binned, 
             yy  = y_binned,
             nn  = y_n_size) %>% 
    setNames(c( size_var, resp_var, 'n_size')) %>% 
    mutate(year  = years$year[ii]) %>% 
    drop_na()
  
}