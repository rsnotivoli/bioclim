#' Function to Positive Temperature index
#'
#' @description Computes Positive Temperature index from monthly temperature.
#' @param t Monthly average temperature data (12 nueric values).
#' @return Positive Temperature index
#' @examples
#' postemp(rnorm(12, 18, 6))
#' @export
#'
postemp <- function(t){ sum(t[which(t > 0)])}
