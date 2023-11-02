#' Bioclimatic classification
#'
#' @description Calculates bioclimatic classification based on bioclimatic balance.
#' @param t Numeric. Monthly temperature required for water balance calculation.
#' @param p Numeric. Monthly precipitation required for water balance calculation.
#' @param lat Numeric. Latitude required for water balance calculation.
#' @param wb Water balance in data.frame format from watbal() function. If provided, 't' and 'p' are not used.
#' @param bb Bioclimatic balancein data.frame format from biobal() function. If provided, 't', 'p' and 'wb' are not used.
#' @param CC Field capacity. It depends on water retention capacity and depth of roots. 400 as default value.
#' @param mode Type of output: "TBR", "sub", or "zonal". See details.
#' @return character defining the type of climate.
#' @details Argument "mode" defines the type of return ("TBR": Types of Bioclimatic Regime; "zonal": zonal units; "sub": bioclimatic regime subtypes)
#' @examples
#' # calculation of water balance
#' wb <- watbal(t = rnorm(12, 18, 6), p = rnorm(12, 50, 30), lat = 35, CC = 400)
#' # calculation of bioclimatic balance
#' bb <- biobal(wb, 400)
#'
#' # bioclimatic classification at TBR levels
#' biotype(bb = bb, mode = 'TBR')
#'
#' # bioclimatic classification at zonal levels
#' biotype(bb = bb, mode = 'zonal')
#'
#' # bioclimatic classification at subtypes levels (requires water balance)
#' wb <- watbal(t = rnorm(12, 18, 6), p = rnorm(12, 50, 30), lat = 35, CC = 400)
#' biotype(wb = wb, CC = 400, mode = 'sub')
#'
#' @export
#'

biotype <- function(t=NULL, p=NULL, lat=NULL, wb=NULL, bb=NULL, CC=NULL, mode){
  if(!is.character(mode)) stop('mode must be set')

  bh <- NULL

  # bioclimatic balance
  if(is.null(bb)){
    # water balance
    if(!is.null(wb)){
      bh <- wb
      t <- bh[1:12,1]
      p <- bh[1:12,2]
    } else{
      bh <- watbal(t, p, lat, CC)
    }
    bb <- biobal(bh, CC)
  } else{
    t <- bb[1:12, 2]
    p <- bb[1:12, 1]
  }


  # Thornthwaite’s index
  if(mode=='sub'){
    if(is.null(bh)){
      stop('For mode = "sub", you need to provide wb or t,p,lat,CC values')
    }
    subtype <- ith(bh)
  }


  # bioclimatic intensities
  # intens <- bioint(bb)
  # other required variables
  # PVH <- length(which(intens$IBSc != 0))
  # PVT <- length(which(intens$IBRf != 0))
  PVH <- length(which((p-(0.2*bh$PET[1:12]))<0))
  PVT <- length(which(t<7.5))
  Tf <- min(t)
  P <- sum(p)


  # type <- NA
  if(Tf < 18){
    #no tropical
    if(Tf > 7.5){
      #subtropical
      zonal <- 'Subtropical'
      if(PVH == 0){
        if(P > 2000) type <- 'Euritermo-Ombrophyllo' else type <- 'Euritermo-Mesophyllo'
      } else{
        if(PVH > 0 & PVH <= 4) {type <- 'Euritermo-Tropophyllo'} else{
          if(PVH >= 5 & PVH < 8) type <- 'Euritermo-Xerophyllo' else type <- 'Euritermo-Hiperxerophyllo'
        }
      }
    } else{
      if(PVT >= 1 & PVT < 6){
        #templado
        zonal <- 'Mid latitudes'
        if(PVH == 0){
          if(P > 2000) type <- 'Cryo-Ombrophyllo' else type <- 'Cryo-Mesophyllo'
        } else{
          if(PVH > 0 & PVH <= 4) type <- 'Cryo-Tropophyllo' else {
            if(PVH >= 5 & PVH < 8) type <- 'Cryo-Xerophyllo' else type <- 'Cryo-Hiperxerophyllo'
          }
        }
      } else{
        if(PVT >= 6 & PVT < 11){
          #subpolar
          zonal <- 'Subpolar'
          if(PVH == 0){
            if(P > 1100) type <- 'Mesocryo-Ombrophyllo' else type <- 'Mesocryo-mesophyllo'
          } else if(PVH > 0 & PVH <= 4) {
            type <- 'Mesocryo-Tropophyllo'} else if(PVH >= 5 & PVH < 8) {
              type <- 'Mesocryo-Xerophyllo'} else type <- 'Mesocryo-Hiperxerophyllo'
        } else{
          #polar
          zonal <- 'Polar'
          if(PVH == 0){
            if(P > 1100) type <- 'Hypercryo-Ombrophyllo' else type <- 'Hypercryo-Mesophyllo'
          } else if(PVH > 0 & PVH <= 4) {
            type <- 'Hypercryo-Tropofilo'} else if(PVH >= 5 & PVH < 8) {
              type <- 'Hypercryo-Xerophyllo'} else type <- 'Hypercryo-Hypexerophyllo'
        }
      }
    }
  } else{
    #Tropical
    zonal <- 'Tropical'
    if(PVH == 0){
      if(P > 2000) type <- 'Ombrophyllo' else type <- 'Mesophyllo'
    } else if(PVH > 0 & PVH <= 4){
      if(P > 2000) type <- 'Ombro-Tropophyllo' else type <- 'Tropophyllo'
    } else if(PVH >= 4 & PVH < 8){
      if(P > 2000) type <- 'Ombro-Xerophyllo' else type <- 'Xerophyllo'
    } else if(PVH >= 8) type <- 'Hyperxerophyllo'
  }
  if(mode == 'zonal'){
    res <- zonal
  } else if(mode == 'TBR'){
    res <- type
  } else if(mode == 'sub'){
    res <- paste(subtype, type)
  }
  return(res)
}
