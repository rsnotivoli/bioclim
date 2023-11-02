#' Function to calculate Thornthwaite’s index
#'
#' @description This function calculates Thornthwaite’s index to refine the bioclimatic classification.
#' @param bh Water balance in data.frame format from watbal() function.
#' @return Character, describing the humid characteristics of the climate.
#' @examples
#' wb <- watbal(t = rnorm(12, 18, 6), p = rnorm(12, 50, 30), lat = 35, CC = 400)
#' ith(wb)
#' @export


ith <- function(bh){
  num <- (bh$TEAW[13]/bh$PET[13])*100
  if(num > 100){
    res <- 'Perhumid'
  } else if(num > 80 & num <= 100){
    res <- 'Very humid'
  } else if(num > 60 & num <= 80){
    res <- 'Highly humid'
  } else if(num > 40 & num <= 60){
    res <- 'Moderate humid'
  } else if(num > 20 & num <= 40){
    res <- 'Low humid'
  } else if(num > 0 & num <= 20){
    res <- 'Moist humid'
  } else if(num > -20 & num <= 0){
    res <- 'Dry humid'
  } else if(num > -40 & num <= -20){
    res <- 'Semiarid'
  } else if(num <= -40){
    res <- 'Arid'
  } 
  return(res)
}
