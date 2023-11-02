#' Monthly precipitation
#'
#' A SpatRaster containing the monthly precipitation of the Alps.
#'
#' @format A PackedSpatRaster with 12 monthly values precipitation:
"preRast"

#' Monthly temperature
#'
#' A SpatRaster containing the monthly temperature of the Alps.
#'
#' @format A PackedSpatRaster with 12 monthly values temperature:
"tmpRast"

#' Water Balance
#'
#' A SpatRaster containing the water balance of the Alps.
#'
#' @format A PackedSpatRaster with 12 monthly values of 12 variables:
#'   \code{temp}, \code{prec}, \code{PET}, \code{P_PET}, \code{PPA},
#'   \code{ST}, \code{i_ST}, \code{RET}, \code{HD}, \code{HEX}, \code{r},
#'   \code{rP}.
"wbRast"

#' Bioclimatic Balance
#'
#' A SpatRaster containing the bioclimatic balance of the Alps.
#'
#' @format A PackedSpatRaster with 12 monthly values of 4 variables:
#' \code{B}, \code{b}, \code{bc} and \code{bl}.
"bbRast"
