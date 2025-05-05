# Padrino ipm execution

# Author: Niklas Neisse
# Co    : Aspen Workman, Aldo Compagnoni
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.01.13


# Execution distributor --------------------------------------------------------
# Conditional execution based on ipm_type
if (v_ipm_type == 'mean') {
  source('pipeline/padrino_mean.R')
} else if (v_ipm_type == 'year_specific') {
  source('pipeline/padrino_year_specific.R')
} else {
  stop("Invalid ipm_type value. Choose either 'mean' or 'year_specific'.")
}