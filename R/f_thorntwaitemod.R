#' Function to calculate Thornthwaite potential evapotranspiration
#'
#' @description This function calculates Thornthwaite’s potential evapotranspiration.
#' @param Tave a numeric vector, matrix or time series of monthly mean temperatures, ºC.
#' @param lat a numeric vector with the latitude of the site or sites, in degrees.
#' @param na.rm optional, a logical value indicating whether NA values should be stripped from the computations.
#' @return A time series with the values of monthly potential or reference evapotranspiration, in mm. If the input is a matrix or a multivariate time series each column will be treated as independent data (e.g., diferent observatories), and the output will be a multivariate time series.
#' @details This function is a modified version of SPEI::thornthwaite() function.
#' @import stats
#' @author Santiago Begueria
#' @examples
#' thornthwaite(rnorm(12, 18, 6), 35)
#' @export


thornthwaite <- function (Tave, lat, na.rm = FALSE)
{
  if (sum(is.na(Tave)) != 0 & na.rm == FALSE) {
    stop("Data must not contain NAs")
  }
  if (is.ts(Tave) & frequency(Tave) != 12) {
    stop("Data should be a monthly time series (frequency = 12)")
  }
  if (!is.null(ncol(Tave))) {
    if (length(lat) != ncol(Tave)) {
      stop("Longitude of latitudes vector does not coincide with the number of columns in Tave")
    }
  }
  if (!is.ts(Tave)) {
    Tave <- ts(as.matrix(Tave), frequency = 12)
  }
  else {
    Tave <- ts(as.matrix(Tave), frequency = frequency(Tave),
               start = start(Tave))
  }
  PE <- Tave * NA
  n <- nrow(Tave)
  m <- ncol(Tave)
  c <- cycle(Tave)
  tanLat <- tan(lat/57.2957795)
  tanDelta <- c(-0.37012566, -0.23853358, -0.04679872, 0.16321764,
                0.32930908, 0.40677729, 0.3747741, 0.239063, 0.04044485,
                -0.16905776, -0.33306377, -0.40743608)
  tanLatM <- matrix(tanLat, nrow = 12, ncol = length(tanLat),
                    byrow = TRUE) * tanDelta
  tanLatM[tanLatM < {
    -1
  }] <- -1
  tanLatM[tanLatM > 1] <- 1
  omega <- acos(-tanLatM)
  N <- 24 * omega/pi
  days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  K <- N/12 * days/30
  T <- matrix(unlist(by(Tave, c, colMeans, na.rm = na.rm)),
              ncol = m, byrow = TRUE)
  T[T < 0] <- 0
  J <- colSums((T/5)^1.514, na.rm = na.rm)
  J <- matrix(J, n, m, TRUE)
  J2 <- J * J
  J3 <- J2 * J
  q <- 6.75e-07 * J3 - 7.71e-05 * J2 + 0.01792 * J + 0.49239
  Tave[Tave < 0] <- 0
  PE <- K[c, ] * 16 * (10 * Tave/J)^q
  PE[is.na(PE)] <- 0
  colnames(PE) <- rep("PET_tho", m)
  return(PE)
}
