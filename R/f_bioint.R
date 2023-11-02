#' Computation of Bioclimatic Intensities
#'
#' @description Computes bioclimatic intensities from bioclimatic balance.
#' @param bb Bioclimatic balance.
#' @return data.frame with 10 variables. See details.
#' @details The function yields 10 variables at monthly scale corresponding with the warm (w) and cold (c) variants of 5 bioclimatic intensities: PBI (Potential bioclimatic intensity), RBI (Real bioclimatic intensity), CBI (Conditioned bioclimatic intensity), FBI (Free bioclimatic intensity), and DBI (Dry bioclimatic intensity). 
#' @examples
#' wb <- watbal(t = rnorm(12, 18, 6), p = rnorm(12, 50, 30), lat = 35, CC = 400)
#' bb <- biobal(wb, 400)
#' bi <- bioint(bb)
#' @export
#'

bioint <- function(bb){
  intens <- matrix(NA, ncol = 14, nrow = 12)
  colnames(intens) <-  c('PBIw', 'PBIw_g', 'CBIw', 'FBIw', 'FBIw_g', 'RBIw', 'DBIw',
                         'PBIc', 'PBIc_g', 'CBIc', 'FBIc', 'FBIc_g', 'RBIc', 'DBIc')
  rownames(intens) <- month.name
  intens <- as.data.frame(intens)

  for(i in 1:12){
    if(bb$PBI[i] > 0) intens$PBIw[i] <- bb$PBI[i] else intens$PBIw[i] <- 0
    if(bb$CBI[i] > 0) intens$CBIw[i] <- bb$CBI[i] else intens$CBIw[i] <-0
    if(bb$PBI[i] > 0 & bb$RBI[i] > 0){if(bb$FBI[i] > 0) intens$FBIw[i] <- bb$FBI[i] else intens$FBIw[i] <- bb$RBI[i]} else intens$FBIw[i] <- 0
    # if(intens$FBIw[i] < intens$CBIw[i]) intens$FBIw_g[i] <- 0 else intens$FBIw_g[i] <- intens$FBIw[i] - intens$CBIw[i]
    # intens$PBIw_g[i] <- intens$PBIw[i] - (intens$CBIw[i] + intens$FBIw_g[i])
    intens$RBIw[i] <- intens$FBIw[i] + intens$CBIw[i]
    if(bb$PBI[i] > 0 & bb$RBI[i] < 0) intens$DBIw[i] <- bb$RBI[i] else intens$DBIw[i] <- 0
    if(bb$PBI[i] < 0) intens$PBIc[i] <- bb$PBI[i] else intens$PBIc[i] <- 0
    if(bb$CBI[3] < 0) intens$CBIc[i] <- bb$CBI[3] else intens$CBIc[i] <- 0
    if(bb$PBI[i] < 0 & bb$RBI[i] < 0) {if(bb$FBI[i] < 0) intens$FBIc[i] <- bb$FBI[i] else intens$FBIc[i] <- bb$RBI[i]} else intens$FBIc[i] <- 0
    # if(intens$CBIc[i] < intens$FBIc[i]) intens$FBIc_g[i] <- 0 else intens$FBIc_g[i] <- intens$FBIc[i] - intens$CBIc[i]
    # intens$PBIc_g[i] <- intens$PBIc[i] - (intens$CBIc[i] + intens$FBIc_g[i])
    intens$RBIc[i] <- intens$FBIc[i] + intens$CBIc[i]
    if(bb$PBI[i] < 0) {if(bb$RBI[i] > 0) intens$DBIc[i] <- bb$RBI[i] else intens$DBIc[i] <- 0} else intens$DBIc[i] <- 0
  }
  intens <- intens[, c(1,3,4,6,7,8,10,11,13,14)]
  intens
}
