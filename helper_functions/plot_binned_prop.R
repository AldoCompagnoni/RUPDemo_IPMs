# function to plot your survival data "binned" (instead of "jittered")

# Arguments: data frame, number of bins, 
# UNQUOTED name of size variable, 
# UNQUOTED name of response variable
plot_binned_prop <- function(df, n_bins, siz_var, rsp_var){
  
  size_var <- deparse( substitute(siz_var) )
  resp_var <- deparse( substitute(rsp_var) )
  
  # remove all NAs
  na_ids   <- c( which( is.na(df[,size_var]) ),
                 which( is.na(df[,resp_var]) )
  ) %>% unique
  
  # remove NAs only if there is at least one NA
  if( length(na_ids) > 0 ) df <- df[-na_ids,]
  
  # binned survival probabilities
  h    <- (max(df[,size_var],na.rm=T) - min(df[,size_var],na.rm=T)) / n_bins
  lwr  <-  min(df[,size_var],na.rm=T) + (h*c(0:(n_bins-1)))
  upr  <- lwr + h
  mid  <- lwr + (1/2*h)
  
  # standard error of a bernoulli process
  # https://stats.stackexchange.com/questions/82720/confidence-interval-around-binomial-estimate-of-0-or-1
  se_bern <- function( x, lwr_upr ){
    
    # do not suppress potential convergence warnings
    surv_n <- sum( x )
    tot_n  <- length( x )
    binom.confint( surv_n, tot_n, methods=c("wilson") )[,lwr_upr]
    
  }
  
  binned_prop <- function(lwr_x, upr_x, response, lwr_upr ){
    
    id  <- which(df[,size_var] >= lwr_x & df[,size_var] < upr_x) 
    tmp <- as.data.frame(df[id,])
    
    if( response == 'prob' ){   return( sum(tmp[,resp_var],na.rm=T) / nrow(tmp) ) }
    if( response == 'n_size' ){ return( nrow(tmp) ) }
    if( response == 'se' ){     return( se_bern(tmp[,resp_var], lwr_upr) ) }
    
  }
  
  y_binned <- Map(binned_prop, lwr, upr, 'prob')        %>% unlist
  x_binned <- mid
  y_n_size <- Map(binned_prop, lwr, upr, 'n_size')      %>% unlist
  y_se_lwr <- Map(binned_prop, lwr, upr, 'se', 'lower') %>% unlist
  y_se_upr <- Map(binned_prop, lwr, upr, 'se', 'upper') %>% unlist
  
  data.frame(x_binned, 
             y_binned,
             n_s  = y_n_size,
             lwr  = y_se_lwr,
             upr  = y_se_upr ) %>% 
    mutate( n_prob = y_n_size/sum(y_n_size) ) %>% 
    setNames( c(size_var, resp_var, 'n_s', 'lwr', 'upr', 'n_prob') )
  
}
