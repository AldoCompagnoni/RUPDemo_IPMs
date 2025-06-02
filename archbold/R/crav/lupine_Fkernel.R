# Building an F kernel
# 28.05.2025 AW


# Notes on this specific model:

  # The code below refers only to a mean IPM - not year- or site-specific IPMs

  # Life cycle of this species includes two seed banks and a continuous size class
    # of adult individuals. I have tried to simplify this where I can.

  # We have data on the number of fruits are produced by each reproductive stem
    # (called a "raceme" - "fruit_rac"), and the number of seeds in each fruit 
    # ("seed_fruit"), on average. These data come from another source and are
    # constants in this model.

  # "flow_b0" and "flow_b1" refer to the probability of flowering for an
    # individual of a given size, which we have modeled as
    # glm( flow_t0 ~ log_area_t0, family = "binomial" ) where flow_t0 is a
    # binary variable referring to if an individual flowered in a given year (1)
    # or not (0)

  # "fert_b0" and "fert_b1" refer to the number of racemes produced by a
    # reproductive individual of a given size, which we have modeled as
    # glm( numrac_t0 ~ log_area_t0, family = "poisson" ) where numrac_t0 is a 
    # count variable referring to the number of racemes produced by an individual
    # in the given year

  # "recr_sz" is the mean size of new recruits calculated from the raw data
    # aggregated across all recruits (not year- or site-specific)

  # "recr_sd" is the standard deviation of new recruit size calculated from the
    # raw data aggregated across all recruits (not year- or site-specific)


# Fecundity function

fx <- function( x, pars ){
  
  # Calculate total racemes produced
    # Product of probability of flowering of an individual of size x (1st term) 
    # and number of racemes produced by an individual of size x (2nd term)
  tot_rac  <- inv_logit( pars$flow_b0 + pars$flow_b1*x ) * 
              exp(       pars$fert_b0 + pars$fert_b1*x )
  
  # Calculate total viable seeds produced from these racemes
  viab_sd  <- tot_rac * pars$fruit_rac * pars$seed_fruit
  
  return( viab_sd )
}

# Size distribution of new recruits
  # Probably similar to other models - this just gives each new recruit a realistic
  # size as it is added to the population

recs <- function( y, pars ){
  dnorm( y, mean = pars$recr_sz, sd = pars$recr_sd )
}

# Full fecundity function
  # Multiply the output from the fecundity function (number of seeds) and the 
  # size distribution of recruits to add the appropriate amount of realistically-
  # sized recruits into the population

fxy <- function( y, x, pars ){
  fx( x, pars ) * recs( y, pars )
}



# Implementation into IPM kernel function
  # Ignore pretty much everything with "##" as it shouldn't be much different
  # from things you've already seen

kernel <- function( pars ){
  
  ## set up IPM domain ---------------------------------------------------------
  
  n <- pars$mat_siz
  L <- pars$L
  U <- pars$U
  h <- ( U - L ) / n
  b <- L + c( 0:n ) * h
  y <- 0.5 * ( b[1:n] + b[2:( n + 1 )] )
  
  # populate kernel ------------------------------------------------------------
  
  # seeds mini matrix (2x2 to account for both seed bank stages)
  s_mat <- matrix( 0, 2, 2 )
  
  # seeds that enter 1 yr-old seed bank (fecundity function multiplied by a 
    # germination factor)
  plant_s1 <- fx( y, pars ) * ( 1 - pars$g0 )
  
  # no seeds go directly to 2 yr-old seed bank
  plant_s2 <- numeric( n )
  
  # seeds that go directly to seedlings germinate right away (full fecuntity 
    # function multiplied by a germination factor)
  Fmat <- ( outer( y, y, fxy, pars ) * pars$g0 * h )
  
  # recruits from the 1 yr-old seedbank (recruit size distribution multiplied
    # by a germination factor)
  s1_rec <- h * recs( y, pars ) * pars$g1
  
  # seeds that enter the 2 yr-old seed bank (just a constant germination factor)
  s_mat[2,1] <- ( 1 - pars$g1 )
  
  # recruits from the 2 yr-old seedbank (recruit size distribution multiplied by
    # a different germination factor)
  s2_rec <- h * recs( y, pars ) * pars$g2
  
  ## survival and growth of adult plants
  Tmat <- ( outer( y, y, pxy_ns, pars ) * h )
  
  small_K <- Tmat + Fmat
  
  # assemble the kernel --------------------------------------------------------
  
  # top two vectors (seeds entering the seed banks) and survival/growth of
    # continuous size class
  from_plant <- rbind( rbind( plant_s1, plant_s2 ),
                       small_K )
  
  # leftmost vectors (seeds exiting seed banks into continuous size class)
  from_seed <- rbind( s_mat,
                      cbind( s1_rec, s2_rec ) )
  
  # stick 'em all together
  k_yx <- cbind( from_seed, from_plant )
  
  return( k_yx )
  
}




