#' Function to calculate Thornthwaite’s index (raster format)
#'
#' @description This function calculates Thornthwaite’s index to refine the bioclimatic classification.
#' @param bh Water balance in SpatRaster format from watbalRaster() function.
#' @return Numeric, describing the humid characteristics of the climate. 1: 'HyperArid', 2: 'Arid', 3: 'Semiarid', 4: 'Dry humid', 5: 'Moist humid', 6 'Low humid', 7: 'Moderate humid', 8: 'Highly humid', 9: 'Very humid', 10: 'Perhumid'.
#' @examples
#' \donttest{
#' wb <- terra::rast(wbRast)
#' itr <- ithRaster(wb)
#' }
#' @export

ithRaster <- function(bh){
  num <- (sum(bh[[37:48]])/sum(bh[[25:36]]))*100

  m <- c(-99999999, -40, 1,
         -40, -20, 2,
         -20, 0, 3,
         0, 20, 4,
         20, 40, 5,
         40, 60, 6,
         60, 80, 7,
         80, 100, 8,
         100, 99999999, 9)
  m <- matrix(m, ncol=3, byrow=TRUE)
  clf <- terra::classify(num, m)

  # 1: 'Arid',
  # 2: 'Semiarid',
  # 3: 'Dry humid',
  # 4: 'Moist humid',
  # 5 'Low humid',
  # 6: 'Moderate humid',
  # 7: 'Highly humid',
  # 8: 'Very humid',
  # 9: 'Perhumid'

  return(clf)
}
