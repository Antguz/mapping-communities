# RGB color
library(stringr)

RGB <- function(X, Y, Z, min_range = -0.5, max_range = 0.5) {
  
  if(is.null(min_range) == TRUE & is.null(max_range) == TRUE) {
    
    x_new <-  ((X - min(X)) / ((max(X)) - (min(X)))) * (255 - 0) + 0
    y_new <-  ((Y - min(Y)) / ((max(Y)) - (min(Y)))) * (255 - 0) + 0
    z_new <-  ((Z - min(Z)) / ((max(Z)) - (min(Z)))) * (255 - 0) + 0
    
  } else {
    
    x_new <-  ((X - min_range) / (max_range - min_range)) * (255 - 0) + 0
    y_new <-  ((Y - min_range) / (max_range - min_range)) * (255 - 0) + 0
    z_new <-  ((Z - min_range) / (max_range - min_range)) * (255 - 0) + 0
    
  }
  
  # Convert colour channels to hex
  hex_value <- paste0(
    "#",
    # Convert red
    str_pad(
      string = as.hexmode(round(x_new, 0)),
      width = 2, pad = "0", side = "left"
    ),
    # Convert green
    str_pad(
      string = as.hexmode(round(y_new, 0)),
      width = 2, pad = "0", side = "left"
    ),
    # Convert blue
    str_pad(
      string = as.hexmode(round(z_new, 0)),
      width = 2, pad = "0", side = "left"
    )
  )
  
  return(hex_value)
  
}
