# MPM year specific Adler 2007 Kansas Hedyotis nigricans

# Author: Niklas Neisse
# Email: neisse.n@protonmail.com &
#  aldo.compagnoni@gmail.com
# Date: 2024.10.30

# Comments ---------------------------------------------------------------------
# 1. the pipeline runs plant tracker only if the data does not exist already
# 2. find all the graphics in the result folder of the respective species
# 2.1 and the survival and recruitment data in its respective folder
# 3. check the global environment for all the other objects
# 3.1 recr_pred_df: view(recr_pred_df)
#      year specific recruitment predictions
# 3.2 pop_counts:   view(pop_counts)
#      year specific lambdas and underlying survival and recruitment predictions
#      alongside observed population count and growth rate
# 3.3 model survival: summary(mod_surv)
# 3.4 model recruit1: summary(yesno_mod)
#      recruitment being there or not
# 3.4 model recruit2: summary(pcr_mod)
#      per-capita recruitment (pcr) *conditional* on recruitment happening
# 3.5 model recruit3: summary(pcr_uc_mod)
#      per capita recruitment *unconditinoal*


# Define the species variable --------------------------------------------------
species <- 'Hedyotis nigricans'


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