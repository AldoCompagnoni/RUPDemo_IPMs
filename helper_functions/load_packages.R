
# load or install packages 
load_packages <- function( ... ){
  
  # Define CRAN packages
  .cran_packages <- sapply( eval( substitute( alist(...) ) ),
                            deparse )
  
  # Check if CRAN packages are installed
  .inst <- .cran_packages %in% installed.packages() 
  
  if(any(!.inst)) {
    # Install missing CRAN packages
    install.packages(.cran_packages[!.inst]) 
  }
  
  # Load required packages
  sapply(.cran_packages, require, character.only = TRUE) 
  
}
