# MPM year specific Adler 2007 Kansas Ratibida columnifera

# Author: Niklas Neisse
# Email: neisse.n@protonmail.com &
#  aldo.compagnoni@gmail.com
# Date: 2024.10.30

# Comments ---------------------------------------------------------------------
#  find all the results in the result folder of the respective species 
#  check the global environment for all the objects
#   if this is the format we are opting for in the longrun 
#    maybe it would be great to list and describe important objects


# Define the species variable --------------------------------------------------
species <- 'Ratibida columnifera'


# MPM pipeline -- with model comparison -- -------------------------------------
source('helper_functions/adler07_mpm_forbs.R')


# Published vs our lambdas -----------------------------------------------------
rmse( compute_df$obs_pgr,
      compute_df$pub_lam )
rmse( compute_df$obs_pgr,
      compute_df$lam )
rmse( compute_df$obs_pgr,
      compute_df$pub_proj_lam )
rmse( compute_df$obs_pgr,
      compute_df$proj_lam )