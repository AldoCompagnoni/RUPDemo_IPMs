# function to plot your survival data "binned" (instead of "jittered")
# https://stats.stackexchange.com/questions/82720/confidence-interval-around-binomial-estimate-of-0-or-1

# Arguments: data frame, number of bins, 
# UNQUOTED name of discrete variable (usually age)
# UNQUOTED name of response variable (usually survival)
plot_binned_prop <- function(df, age_var, rsp_var){
  
  # set up the variables
  size_var <- deparse( substitute(age_var) )
  resp_var <- deparse( substitute(rsp_var) )
  
  # remove all NAs
  na_ids   <- c( which( is.na(df[,size_var]) ),
                 which( is.na(df[,resp_var]) )
  ) %>% unique
  
  # remove NAs only if there is at least one NA
  if( length(na_ids) > 0 ) df <- df[-na_ids,]
  
  # final data frame
  df %>% 
    # group by age
    group_by( age ) %>% 
    # compute numbers of surviving individuals, and total individuals 
    summarise( surv_n = sum(survives),
               tot_n  = n() ) %>% 
    ungroup %>% 
    # compute proportions of survivors, and upper/lower 95% C.I.
    mutate( surv_p = surv_n / tot_n,
            lwr    = binom.confint( surv_n, tot_n, methods=c("wilson") )$lower,
            upr    = binom.confint( surv_n, tot_n, methods=c("wilson") )$upper
            ) %>% 
    dplyr::select( -surv_n, -tot_n ) %>% 
    setNames( c(size_var,resp_var,'lwr','upr') ) %>% 
    as.data.frame
  
}
