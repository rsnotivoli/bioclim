#' Function to plot Walter and Lieth diagram
#'
#' @description Function to plot Walter and Lieth diagram.
#' @param t Monthly average temperature data (12 nueric values).
#' @param p Monthly average precipitation data (12 nueric values).
#' @return Plot of Walter and Lieth diagram
#' @importFrom berryFunctions climateGraph
#' @examples
#' plotWL(t = rnorm(12, 18, 6), p = rnorm(12, 50, 30))
#' @export
#'
plotWL <- function(t, p){
  climateGraph(t, p, main = '',
               mar=c(2,3,4,3), textprop = 0)
}
