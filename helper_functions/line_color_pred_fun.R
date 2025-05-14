# Function to define the color depending on the number of parameters

# Author: Niklas Neisse
# Co    : Aspen Workman, Diāna Spurīte, Aldo Compagnoni*
# Email : neisse.n@protonmail.com
# Main  : aldo.compagnoni@idiv.de
# Web   : https://aldocompagnoni.weebly.com/
# Date  : 2025.02.21

line_color_pred_fun <- function(ranef) {
   # Determine the line color of a prediction 
    # based on the number of random effects
  
  
  # Check if the length of ranef is 2, and return 'red' if true
  if (length(ranef) == 2) {
    return('red')
  }
  # For 3 return 'green' 
  else if (length(ranef) == 3) {
    return('green')
  }
  # For 4 and return 'blue'
  else if (length(ranef) == 4) {
    return('blue')
  }
  # If none of the above conditions are met, return 'black'
  else {
    return('black')
  }
}