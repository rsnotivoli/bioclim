#' Function to Compensated Thermal Index
#'
#' @description Computes Compensated Thermal Index from monthly temperature.
#' @param t Monthly average temperature data (12 nueric values).
#' @return Compensated Thermal Index
#' @examples
#' thermind(rnorm(12, 18, 6))
#' @export
#'
thermind <- function(t){
  It <- 10 * (mean(t) + 2 * min(t))
  Am <- diff(range(t))
  C1 <- 5 * (Am - 18)
  C2 <- 10 * (Am - 21)
  C3 <- 5 * (Am - 27)
  C4 <- 20 * (Am - 46)
  if(Am < 9) Itc <- It - 90 + 10 * Am
  if(Am >= 9 & Am <= 18) Itc <- It
  if(Am > 18) Itc <- It + C1 + C2 + C3 + C4
  Itc
}
