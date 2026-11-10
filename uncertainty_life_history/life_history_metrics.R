# Functions to calculate life history metrics from ipmr objects

# Need to use the development version of ipmr
# Run this whole section to ensure the correct version is attached

ipmr_required <- "0.0.8"
ipmr_repo     <- "anon-lupineuncertainty/ipmr"
ipmr_ref      <- "v0.0.8"

# Version currently installed, or NA if the package isn't there at all
ipmr_current <- tryCatch(
  as.character(utils::packageVersion("ipmr")),
  error = function(e) NA_character_
)

if (is.na(ipmr_current) || ipmr_current != ipmr_required) {
  
  message(
    "ipmr ", ifelse(is.na(ipmr_current), "not installed", paste0(ipmr_current, " found")),
    " — installing ", ipmr_required, " from ", ipmr_repo, "@", ipmr_ref
  )
  
  # remotes is needed to do the install
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  
  # Can't overwrite a package whose namespace is loaded
  if ("ipmr" %in% loadedNamespaces()) unloadNamespace("ipmr")
  
  if (!is.na(ipmr_current)) {
    remove.packages("ipmr", lib = dirname(find.package("ipmr")))
  }
  
  remotes::install_github(ipmr_repo, ref = ipmr_ref, upgrade = "never")
}

library(ipmr)

stopifnot(as.character(packageVersion("ipmr")) == ipmr_required)

# Now you can run the functions below to test them with your models


# Helper functions -------------------------------------------------------------

# Dominant eigenvalue
  # Not exported from ipmr, oops

dominant_eigen <- function(mat) {
  eig <- eigen(mat)
  Re(eig$values[which.max(Mod(eig$values))])
}

# Extract P and F kernels from any IPM
  # Users must specify which of the sub-kernels belong to the P and F kernels
  # Users must also specify the mega matrix structure, as in the uncertainty
    # tools
  # Default behavior should work for most simple IPMs as long as the P and F
    # kernels are called P and F

get_P_F_kernels <- function( ipm,
                             P_kernels = "P",
                             F_kernels = "F",
                             mega_mat = NULL ) {
  
  # Validate kernel specifications
  
  if( !is.character(P_kernels) || length(P_kernels) == 0 ) {
    stop("'P_kernels' must be a non-empty character vector of subkernel names.")
  }
  
  if( !is.character(F_kernels) || length(F_kernels) == 0 ) {
    stop("'F_kernels' must be a non-empty character vector of subkernel names.")
  }
  
  overlap <- intersect(P_kernels, F_kernels)
  
  if( length(overlap) > 0 ) {
    stop(
      "The following subkernels are included in both 'P_kernels' ",
      "and 'F_kernels': ",
      paste(overlap, collapse = ", "),
      ". Each subkernel must be assigned to either P or F."
    )
  }
  
  available <- names(ipm$sub_kernels)
  
  missing_P <- setdiff(P_kernels, available)
  missing_F <- setdiff(F_kernels, available)
  
  if( length(missing_P) > 0 ) {
    stop(
      "The following names in 'P_kernels' are not present in ",
      "'ipm$sub_kernels': ",
      paste(missing_P, collapse = ", ")
    )
  }
  
  if( length(missing_F) > 0 ) {
    stop(
      "The following names in 'F_kernels' are not present in ",
      "'ipm$sub_kernels': ",
      paste(missing_F, collapse = ", ")
    )
  }
  
  
  # Simple IPM
  
  # All subkernels occupy the same matrix space, so no mega_mat is required.
  
  if( inherits(ipm, "simple_ipm") ||
      grepl("simple", class(ipm)[1]) ) {
    
    P <- Reduce(
      `+`,
      lapply(P_kernels, function(x) ipm$sub_kernels[[x]])
    )
    
    F <- Reduce(
      `+`,
      lapply(F_kernels, function(x) ipm$sub_kernels[[x]])
    )
    
    return(
      list(
        P = P,
        F = F
      )
    )
  }
  
  
  # General IPM 
  
  if( is.null(mega_mat) ) {
    stop(
      "'mega_mat' must be supplied for general IPMs because the ",
      "subkernels do not necessarily occupy the same matrix space."
    )
  }
  
  
  # Normalize mega_mat
  
  mega_mat <- rlang::enquo(mega_mat)
  
  if( rlang::quo_is_call(mega_mat) ) {
    
    mega_expr <- rlang::get_expr(mega_mat)
    
  } else {
    
    text <- rlang::eval_tidy(mega_mat)
    
    if( length(text) > 1 ) {
      
      mega_expr <- rlang::call2(
        "c",
        !!!rlang::syms(text)
      )
      
    } else {
      
      mega_expr <- rlang::parse_expr(text)
      
    }
  }
  
  
  # Check mega_mat structure
  
  if( !rlang::is_call(mega_expr, "c") ) {
    stop(
      "'mega_mat' must be specified using 'c()'. For example: ",
      "'c(P, F, A, B)' or c(\"P\", \"F\", \"A\", \"B\")."
    )
  }
  
  blocks <- rlang::call_args(mega_expr)
  
  
  # Identify kernels occurring in mega_mat
  
  block_kernels <- lapply(
    blocks,
    function(x) {
      .args_from_txt(rlang::expr_text(x))
    }
  )
  
  all_mega_kernels <- unique(unlist(block_kernels))
  
  structural <- c("0", "`0`", "I")
  
  mega_kernels <- setdiff(
    all_mega_kernels,
    structural
  )
  
  
  # Check that every kernel is classified
  
  unclassified <- setdiff(
    mega_kernels,
    c(P_kernels, F_kernels)
  )
  
  if( length(unclassified) > 0 ) {
    stop(
      "The following subkernels occur in 'mega_mat' but are not included ",
      "in either 'P_kernels' or 'F_kernels': ",
      paste(unclassified, collapse = ", "),
      ". Add them to the appropriate kernel list."
    )
  }
  
  
  # Recursively classify a mega-matrix block
  
  classify_block <- function( expr, keep, component ) {
    
    # Structural elements
    
    if( rlang::is_symbol(expr) ) {
      
      nm <- rlang::as_string(expr)
      
      # Identity represents persistence and therefore belongs to P.
      if( nm == "I" ) {
        if( component == "P" ) {
          return(expr)
        } else {
          return(rlang::expr(0))
        }
      }
      
      # A subkernel is retained if it belongs to the requested component.
      if( nm %in% keep ) {
        return(expr)
      }
      
      # Any other symbol is removed.
      return(rlang::expr(0))
    }
    
    
    # Explicit zero
    
    if( rlang::is_numeric(expr) &&
        length(expr) == 1 &&
        identical(expr, 0) ) {
      
      return(rlang::expr(0))
    }
    
    
    # Recursive expressions
    
    if( rlang::is_call(expr) ) {
      
      # Addition is the main expression we need to handle because multiple
      # subkernels can contribute to the same matrix block.
      
      if( rlang::call_name(expr) == "+" ) {
        
        args <- rlang::call_args(expr)
        
        classified <- lapply(
          args,
          classify_block,
          keep = keep,
          component = component
        )
        
        # Remove terms that contributed nothing to this component.
        keep_terms <- !vapply(
          classified,
          function(x) {
            identical(x, rlang::expr(0))
          },
          logical(1)
        )
        
        classified <- classified[keep_terms]
        
        # Nothing from this component occurs in the block.
        if( length(classified) == 0 ) {
          return(rlang::expr(0))
        }
        
        # Only one term remains.
        if( length(classified) == 1 ) {
          return(classified[[1]])
        }
        
        # Reconstruct the addition.
        return(
          rlang::call2(
            "+",
            !!!classified
          )
        )
      }
      
      # For expressions other than addition, retain the expression only if
      # all of its subexpressions belong to the requested component.
      #
      # This prevents us from silently changing the meaning of an arbitrary
      # expression that we do not know how to decompose.
      
      args <- rlang::call_args(expr)
      
      classified <- lapply(
        args,
        classify_block,
        keep = keep,
        component = component
      )
      
      if( all(vapply(
        classified,
        function(x) !identical(x, rlang::expr(0)),
        logical(1)
      )) ) {
        
        return(
          rlang::call2(
            rlang::call_name(expr),
            !!!classified
          )
        )
      }
      
      return(rlang::expr(0))
    }
    
    
    # Anything else is treated as an empty contribution.
    rlang::expr(0)
  }
  
  
  # Create P and F versions of every block
  
  P_blocks <- lapply( blocks, classify_block, 
                      keep = P_kernels, component = "P" )
  
  F_blocks <- lapply( blocks, classify_block,
                      keep = F_kernels, component = "F" )
  
  
  # Reconstruct complete mega_mat expressions
  
  P_expr <- rlang::call2(
    "c",
    !!!P_blocks
  )
  
  F_expr <- rlang::call2(
    "c",
    !!!F_blocks
  )
  
  
  # Assemble matrices using ipmr's matrix-construction machinery
  
  P <- .make_mega_mat( ipm = ipm, mega_mat = P_expr )
  
  F <- .make_mega_mat( ipm = ipm, mega_mat = F_expr )
  
  
  return( list( P = P,
                F = F ) )
}


# Calculate the stable stage distribution

get_stable_dist <- function( ipm ){
  
  out <- right_ev( ipm )[[1]]
  
  return( out )
}

# Calculate reproduction by state 

get_reproduction_by_state <- function( ipm,
                                       P_kernels = "P",
                                       F_kernels = "F",
                                       mega_mat = NULL ){
  
  kers <- get_P_F_kernels(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  Fm <- kers$F
  
  out <- colSums( Fm )
  
  return( out )
}

# Calculate newborn distribution

get_newborn_dist <- function( ipm,
                             P_kernels = "P",
                             F_kernels = "F",
                             mega_mat = NULL ){
  
  w <- get_stable_dist( ipm )
  
  kers <- get_P_F_kernels(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  Fm <- kers$F
  
  recruits <- as.vector( Fm %*% w )
  
  out <- recruits / sum( recruits )
  
  return( out )
}

# Calculate age-specific survivorship, fecundity, and reproduction

get_age_trajectory <- function( ipm,
                                P_kernels = "P",
                                F_kernels = "F",
                                mega_mat = NULL,
                                tolerance = 1e-8,
                                max_age = 100 ){
  
  kers <- get_P_F_kernels(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  Pm <- kers$P
  Fm <- kers$F
  
  c <- get_newborn_dist(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  n <- c
  
  age <- 0
  l   <- sum( n )
  b   <- sum( Fm %*% n )
  
  out <- data.frame(
    age   = age,
    l     = l,
    lx_fx = b
  )
  
  while( l > tolerance && age < max_age ) {
    
    n <- Pm %*% n
    age <- age + 1
    
    l <- sum( n )
    b <- sum( Fm %*% n )
    
    out <- rbind(
      out,
      data.frame(
        age   = age,
        l     = l,
        lx_fx = b
      )
    )
  }
  
  out$f <- out$lx_fx / out$l
  
  return( out )
}


# Calculate the survival-independent transition matrix

get_survival_independent_matrix <- function( ipm,
                                             P_kernels = "P",
                                             F_kernels = "F",
                                             mega_mat = NULL ){
  
  kers <- get_P_F_kernels(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  Pm <- kers$P
  
  sigma <- colSums( Pm )
  
  out <- sweep( Pm, 2, sigma, "/" )
  
  # Handle states with zero survival probability
  out[, sigma == 0] <- 0
  
  return( out )
}

# Calculate state-specific probability of reproducing

get_reproduction_probability <- function( ipm,
                                          P_kernels = "P",
                                          F_kernels = "F",
                                          mega_mat = NULL,
                                          repro_var = c("poisson",
                                                        "bernoulli",
                                                        "fixed") ){
  
  repro_var <- match.arg( repro_var )
  
  phi <- get_reproduction_by_state(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  if( repro_var == "poisson" ){
    
    p_repro <- 1 - exp(-phi)
    
  } else if( repro_var == "bernoulli" ){
    
    p_repro <- phi
    
  } else if( repro_var == "fixed" ){
    
    p_repro <- as.numeric(phi > 0)
  }
  
  p_repro <- pmin(pmax(p_repro, 0), 1)
  
  return( p_repro )
}

# Identify states that meet the reproductive maturity threshold

get_mature_states <- function( ipm,
                               P_kernels = "P",
                               F_kernels = "F",
                               mega_mat = NULL,
                               maturity_threshold = 0.5,
                               repro_var = "poisson" ){
  
  if( length(maturity_threshold) != 1 ||
      !is.numeric(maturity_threshold) ||
      maturity_threshold < 0 ||
      maturity_threshold > 1 ){
    
    stop("'maturity_threshold' must be a single number between 0 and 1.")
  }
  
  p_repro <- get_reproduction_probability(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat,
    repro_var = repro_var
  )
  
  mature <- p_repro >= maturity_threshold
  
  if( !any(mature) ){
    stop(
      "No states meet the specified maturity threshold of ",
      maturity_threshold, "."
    )
  }
  
  return( mature )
}


# Calculate first-passage information to reproductive maturity

get_maturity_info <- function( ipm,
                               P_kernels = "P",
                               F_kernels = "F",
                               mega_mat = NULL,
                               maturity_threshold = 0.5,
                               repro_var = "poisson" ){
  
  kers <- get_P_F_kernels(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  Pm <- kers$P
  
  mature <- get_mature_states(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat,
    maturity_threshold = maturity_threshold,
    repro_var = repro_var
  )
  
  c0 <- get_newborn_dist(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  transient <- !mature
  
  # Transition matrix among states that have not yet reached maturity.
  Q <- Pm[transient, transient, drop = FALSE]
  
  if( nrow(Q) == 0 ){
    # Everyone starts mature.
    return(
      list(
        mature = mature,
        maturity_probability = 1,
        mean_age = 0,
        first_mature_dist = c0,
        transient = transient
      )
    )
  }
  
  # Probability of entering a mature state in the next step.
  r <- colSums(
    Pm[mature, transient, drop = FALSE]
  )
  
  # Probability of eventually reaching maturity from each transient state.
  h <- solve(
    diag(nrow(Q)) - t(Q),
    r
  )
  
  # Starting probability of eventually reaching maturity.
  p_mature <- sum(c0[transient] * h) + sum(c0[mature])
  
  if( p_mature <= 0 ){
    stop(
      "No individuals in the derived newborn cohort have a ",
      "positive probability of reaching reproductive maturity."
    )
  }
  
  # Conditional transition matrix given eventual maturity.
  #
  # h gives the probability of eventual maturity from each transient state.
  # This is the standard Doob h-transform for conditioning on absorption.
  Q_cond <- Q * outer(h, 1 / h)
  Q_cond[!is.finite(Q_cond)] <- 0
  
  N_cond <- solve(
    diag(nrow(Q_cond)) - Q_cond
  )
  
  # Expected time to maturity, conditional on eventual maturity,
  # from each transient state.
  e <- rep(1, nrow(Q_cond))
  age_cond <- as.vector(e %*% N_cond)
  
  # Weight starting states by their probability of eventual maturity.
  start_weights <- c0[transient] * h
  
  if( any(c0[mature] > 0) ){
    start_weights <- c(
      start_weights,
      c0[mature]
    )
    
    age_values <- c(
      age_cond,
      rep(0, sum(mature))
    )
    
  } else {
    
    age_values <- age_cond
  }
  
  start_weights <- start_weights / sum(start_weights)
  
  mean_age <- sum(start_weights * age_values)
  
  # Distribution of the state in which maturity is first reached.
  #
  # First-maturity probabilities from each transient starting state:
  # N_cond is used to accumulate the conditional occupancy before absorption.
  first_mature_trans <- c0[transient] * h
  
  # This is easiest to obtain from the unconditional absorption
  # probabilities, followed by conditioning on eventual maturity.
  N <- solve(
    diag(nrow(Q)) - Q
  )
  
  first_mature <- rbind(
    Pm[mature, transient, drop = FALSE]
  ) %*% N %*% c0[transient]
  
  first_mature <- as.vector(first_mature)
  
  if( sum(first_mature) > 0 ){
    first_mature <- first_mature / sum(first_mature)
  }
  
  # Add individuals who start in a mature state.
  if( any(c0[mature] > 0) ){
    
    first_mature_full <- numeric(length(c0))
    first_mature_full[mature] <- first_mature
    
    first_mature_full[mature] <-
      first_mature_full[mature] * 
      (1 - sum(c0[mature]) / p_mature)
    
    first_mature_full[mature] <-
      first_mature_full[mature] +
      c0[mature] / p_mature
    
    first_mature <- first_mature_full
    
  } else {
    
    first_mature_full <- numeric(length(c0))
    first_mature_full[mature] <- first_mature
    first_mature <- first_mature_full
  }
  
  return(
    list(
      mature = mature,
      p_repro = get_reproduction_probability(
        ipm = ipm,
        P_kernels = P_kernels,
        F_kernels = F_kernels,
        mega_mat = mega_mat,
        repro_var = repro_var
      ),
      maturity_probability = p_mature,
      mean_age = mean_age,
      first_mature_dist = first_mature,
      transient = transient
    )
  )
}

# Generation time --------------------------------------------------------------

# T_a: Average time between recruitment events along an ancestral lineage.
#
# method = "bienvenu_legendre":
#   T_a = lambda * (v' w) / (v' F w)
#   This is the Bienvenu & Legendre (2015) measure, generalized to
#   continuous-state models by Ellner (2018).
#
# method = "R0":
#   T = log(R0) / log(lambda)
#   The generation-time measure used in Salguero-Gómez et al. (2015).
#
# timestep converts the result from projection intervals to the desired
# time unit. For example, timestep = 10 for a model censused every 10 years.

gen_time <- function( ipm,
                      P_kernels = "P",
                      F_kernels = "F",
                      mega_mat = NULL,
                      method = c("bienvenu_legendre", "R0"),
                      timestep = 1 ){
  
  method <- match.arg( method )
  
  if( length(timestep) != 1 || !is.numeric(timestep) || timestep <= 0 ){
    stop("'timestep' must be a single positive number.")
  }
  
  kers <- get_P_F_kernels(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  Pm <- kers$P
  Fm <- kers$F
  
  K <- Pm + Fm
  
  lam <- dominant_eigen( K )
  
  if( method == "R0" ){
    
    if( abs(lam - 1) < 1e-10 ){
      stop(
        "The R0-based generation time is undefined when ",
        "lambda is approximately 1."
      )
    }
    
    R0 <- R_0(
      ipm = ipm,
      P_kernels = P_kernels,
      F_kernels = F_kernels,
      mega_mat = mega_mat
    )
    
    out <- log(R0) / log(lam)
    
  } else {
    
    # Right eigenvector: stable distribution
    w <- Re(eigen(K)$vectors[, which.max(Mod(eigen(K)$values))])
    w <- w / sum(w)
    
    # Left eigenvector: reproductive value
    v <- Re(eigen(t(K))$vectors[, which.max(Mod(eigen(t(K))$values))])
    
    # Normalize so that v' w = 1
    v <- v / sum(v * w)
    
    denominator <- sum(v * (Fm %*% w))
    
    if( denominator <= 0 ){
      stop(
        "The Bienvenu-Legendre generation time cannot be calculated ",
        "because reproductive output is zero."
      )
    }
    
    out <- lam * sum(v * w) / denominator
  }
  
  return( out * timestep )
}


# Survivorship curve type ------------------------------------------------------

# H: shape of the age-specific survivorship curve as quantified by Keyfitz'
  # entropy, where H < 1 is type I, H = 1 is type II, and H > 1 is type III

surv_type <- function( ipm,
                       P_kernels = "P",
                       F_kernels = "F",
                       mega_mat = NULL,
                       tolerance = 1e-8,
                       max_age = 100 ){
  
  traj <- get_age_trajectory(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat,
    tolerance = tolerance,
    max_age = max_age
  )
  
  out <- -sum( traj$l * log( traj$l ) ) / sum( traj$l )
  
  return( out )
}



# Age at sexual maturity -------------------------------------------------------

# La: The number of years until the average individual reaches age of first
  # sexual reproduction among individuals that eventually mature
# More explicitly: the mean age at first entry into the mature state set,
  # conditional on eventually reaching maturity

age_maturity <- function( ipm,
                          P_kernels = "P",
                          F_kernels = "F",
                          mega_mat = NULL,
                          maturity_threshold = 0.5,
                          repro_var = "poisson",
                          timestep = 1 ){
  
  if( length(timestep) != 1 ||
      !is.numeric(timestep) ||
      timestep <= 0 ){
    stop("'timestep' must be a single positive number.")
  }
  
  info <- get_maturity_info(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat,
    maturity_threshold = maturity_threshold,
    repro_var = repro_var
  )
  
  return( info$mean_age * timestep )
}


# Progressive growth probability -----------------------------------------------

# Gamma: Mean probability of transitioning to a larger size, weighted by the
  # stable stage distribution

prog_growth <- function( ipm,
                         P_kernels = "P",
                         F_kernels = "F",
                         mega_mat = NULL ) {
  
  w <- get_stable_dist( ipm )
  
  Up <- get_survival_independent_matrix(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  gamma_j <- colSums( Up * lower.tri( Up ) )
  
  out <- sum( gamma_j * w )
  
  return( out )
}


# Retrogressive growth probability ---------------------------------------------

# Rho: Mean probability of transitioning to a smaller size, weighted by the
  # stable stage distribution

retro_growth <- function( ipm,
                          P_kernels = "P",
                          F_kernels = "F",
                          mega_mat = NULL ) {
  
  w <- get_stable_dist( ipm )
  
  Up <- get_survival_independent_matrix(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  rho_j <- colSums( Up * upper.tri( Up ) )
  
  out <- sum( rho_j * w )
  
  return( out )
}


# Mean sexual reproduction -----------------------------------------------------

# Mean per capita number of recruits produced across life stages, weighted by
  # the stable stage distribution

mean_rep <- function( ipm,
                      P_kernels = "P",
                      F_kernels = "F",
                      mega_mat = NULL ){
  
  w <- get_stable_dist( ipm )
  
  phi <- get_reproduction_by_state(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  out <- sum( phi * w )
  
  return( out )
}


# Degree of iteroparity --------------------------------------------------------

# S: Spread of reproduction across lifespan, as quantified by Demetrius' entropy

iteroparity <- function( ipm,
                         P_kernels = "P",
                         F_kernels = "F",
                         mega_mat = NULL,
                         tolerance = 1e-8,
                         max_age = 100 ){
  
  traj <- get_age_trajectory(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat,
    tolerance = tolerance,
    max_age = max_age
  )
  
  lam <- as.vector( lambda( ipm ) )
  
  q <- traj$lx_fx / lam^(traj$age + 1)
  q <- q[q > 0]
  
  
  # Remove negligible numerical noise.
  q[abs(q) < 1e-12] <- 0
  
  if (any(q < 0)) {
    stop("Negative reproductive contribution detected in iteroparity.")
  }
  
  q <- q[q > 0]
  
  if (length(q) == 0 || !all(is.finite(q))) {
    return(NA_real_)
  }
  
  # Normalize to protect against small numerical error.
  q <- q / sum(q)
  
  out <- -sum( q * log( q ) )
  
  return( out )
}


# Net reproductive rate --------------------------------------------------------

# Ro: Mean number of recruits produced during mean life expectancy

# Needs to be generalized to work with any ipmr IPM object

R_0 <- function( ipm,
                 P_kernels = "P",
                 F_kernels = "F",
                 mega_mat = NULL ){
  
  kers <- get_P_F_kernels(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  Pm <- kers$P
  Fm <- kers$F
  
  I <- diag( nrow( Pm ) )
  
  N <- solve( I - Pm )
  
  R <- Fm %*% N
  
  out <- dominant_eigen( R )
  
  return( out )
}



# Mature life expectancy -------------------------------------------------------

# Lw: Mean remaining lifespan after reaching reproductive maturity, conditional
  # on reaching maturity
# Could add an optional method later to implement the mature life expectancy as
  # the difference between mean life expectancy and age at maturity, if this is
  # desired

mature_expec <- function( ipm,
                          P_kernels = "P",
                          F_kernels = "F",
                          mega_mat = NULL,
                          maturity_threshold = 0.5,
                          repro_var = "poisson",
                          timestep = 1 ){
  
  if( length(timestep) != 1 ||
      !is.numeric(timestep) ||
      timestep <= 0 ){
    stop("'timestep' must be a single positive number.")
  }
  
  kers <- get_P_F_kernels(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat
  )
  
  Pm <- kers$P
  
  info <- get_maturity_info(
    ipm = ipm,
    P_kernels = P_kernels,
    F_kernels = F_kernels,
    mega_mat = mega_mat,
    maturity_threshold = maturity_threshold,
    repro_var = repro_var
  )
  
  # Mean remaining lifespan from every state.
  N <- solve(
    diag(nrow(Pm)) - Pm
  )
  
  life_expectancy <- as.vector(
    rep(1, nrow(Pm)) %*% N
  )
  
  # Average remaining lifespan from the state in which maturity
  # is first reached.
  out <- sum(
    life_expectancy[info$mature] *
      info$first_mature_dist[info$mature]
  )
  
  return( out * timestep )
}


