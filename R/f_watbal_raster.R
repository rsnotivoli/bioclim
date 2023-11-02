#' Water balance in raster format
#'
#' @description Computes water balance from temperature and precipitation in raster format.
#' @param temp SpatRaster containing 12 layers with monthly temperature from January to December.
#' @param prec SpatRaster containing 12 layers with monthlyprecipitation from January to December.
#' @param PET Optional. Potential evapotranspiration in raster format.
#' @param CC Field capacity. It depends on water retention capacity and depth of roots. 400 as default value. It can be a SpatRaster layer.
#' @param path Optional. Path (folder) where the output raster files and look-up-tables will be saved.
#' @param ncpu Number of cores used for the most demanding calculations.
#' @return SpatRaster with 144 layers corresponding to the 12 monthly values of 'temp', 'prec','PET','P_PET','PPA','ST','i_ST','RET','HD','HEX','r','rP'.
#' @import terra
#' @importFrom methods is
#' @examples
#' \donttest{
#' tmp <- terra::rast(tmpRast)
#' pre <- terra::rast(preRast)
#' wb <- watbalRaster(tmp, pre, PET = NULL, CC = 400, path=NULL, ncpu = 2)
#' }
#' @export



watbalRaster <- function(temp, prec, PET = NULL, CC, path=NULL, ncpu = 2){

  # extract latitudes
  e <- prec[[1]]
  e <- as.points(e)
  e$lat <- crds(e)[,2]
  lat <- terra::rasterize(e, prec, 'lat')

  # mask
  msk <- lat/lat

  if(is.numeric(CC)){
    if(length(CC)>1) stop("CC must be a single number") else{
      CC <- msk*CC
    }
  } else if(!is(CC, 'SpatRaster')) stop("CC must be numeric or a SpatRaster")


  if(!is.null(PET)){
    if(!is(PET, 'SpatRaster')) stop("PET must be a SpatRaster") else{
      pet <- PET
    }
  } else{
    # thorntwaite calculation
    lt <- c(temp,lat)
    ft <- function(lt){
      if(length(which(is.na(lt)))>0){
        rr <- rep(NA,12)} else{
          rr <- thornthwaite(lt[1:12],lt[13], na.rm = FALSE)
        }
      return(rr)
    }

    pet <- app(lt, ft, cores = ncpu)
  }
  #P-PET
  p_pet <- prec-pet

  #ppa (cumulated potential leakage)
  ppa_fun <- function(p_pet){
    ppa <- p_pet
    w1 <- which(p_pet > 0)
    w2 <- which(p_pet <= 0)
    if(length(w1) > 0){
      ppa[w1] <- 0
    }
    if(length(w2) > 0){
      for(h in w2){
        if(h == 1){
          ppa[h] <- p_pet[h] }else{
            ppa[h] <- ppa[(h - 1)] + p_pet[h]
          }
      }
    }
    return(ppa)
  }
  ppa <- app(p_pet, ppa_fun, cores = ncpu)

  #ST (soil pore water)
  stdata <- c(p_pet,ppa, CC) # join required data
  st_fun <- function(stdata){
    if(length(which(is.na(stdata))) > 0){
      st <- rep(NA, 12)
    } else{
      p_pet <- stdata[1:12]
      ppa <- stdata[13:24]
      CC <- stdata[25]
      st <- rep(NA, 12)
      for(i in 1:12){
        if(i == 1){ #if first month
          if(p_pet[i] > 0){
            if(p_pet[i] > CC) {st[i] <- CC} else {st[i] <- p_pet[i]}
          } else{
            st[i] <- CC * exp(ppa[i]/CC)
          }
        } else{
          if(p_pet[i] > 0){
            if(p_pet[i] + st[(i-1)] > CC) {st[i] <- CC} else {st[i] <- p_pet[i] + st[(i-1)]}
          } else{
            st[i] <- CC * exp(ppa[i]/CC)
          }
        }
      }
    }
    return(st)
  }
  st <- app(stdata, st_fun, cores = ncpu)

  #loop to update values until no changes
  dataloop <- c(p_pet, ppa, st, CC)

  floop <- function(dataloop){
    if(length(which(is.na(dataloop))) > 0){
      ppa <- st <- k <- rep(NA, 12)
    } else{
      p_pet <- dataloop[1:12]
      ppa <- dataloop[13:24]
      st <- dataloop[25:36]
      CC <- dataloop[37]
      seg <- TRUE
      it <- 0
      while(it < 20) { #should work with less than 20
        it <- it + 1
        # print(it)
        #ppa
        ppa2 <- ppa
        for(h in 1:12){
          if(h == 12){ #if last Dec, takes p_pet from Jan and Feb
            if(p_pet[1] >= 0) {
              if(p_pet[h] >= 0){ppa2[h] <- 0} else {
                ppa2[h] <- ppa2[(h-1)] + p_pet[h]
              }
            } else{
              if(p_pet[h] >= 0){if(st[h-1]+ p_pet[h] < CC) {ppa2[h] <- log(st[h]/CC) * CC} else 0} else{
                ppa2[h] <- ppa2[h-1] + p_pet[h]
              }
            }
          } else{
            if(p_pet[(h + 1)] >= 0) {
              if(p_pet[h] >= 0){ppa2[h] <- 0} else {
                if(h == 1){ppa2[h] <- ppa[length(ppa)] + p_pet[h]} else {ppa2[h] <- ppa2[(h-1)] + p_pet[h]}
              }
            } else{
              if(p_pet[h] >= 0){if(st[length(st)]+ p_pet[h] < CC) ppa2[h] <- log(st[h]/CC) * CC else 0} else{
                if(h == 1){ppa2[h] <- ppa[length(ppa)] + p_pet[h]} else {ppa2[h] <- ppa2[(h-1)] + p_pet[h]}
              }
            }
          }
        }


    #st
    st2 <- st
    for(i in 1:12){
      if(h == 12) ini_p_pet <- p_pet[1] else ini_p_pet <- p_pet[(i + 1)]
      if(ini_p_pet >= 0) {
        if(p_pet[i] >= 0){
          if(i == 1){if((st[length(st)] + p_pet[i]) > CC) {st2[i] <- CC} else {st2[i] <- (st[length(st)] + p_pet[i])}} else {
            if((st2[(i - 1)] + p_pet[i]) > CC) {st2[i] <- CC} else {st2[i] <- (st2[(i - 1)]  + p_pet[i])}
          }
        } else{st2[i] <- CC * exp(ppa2[i]/CC)}

      } else{
        if(p_pet[i] >= 0){
          if(i == 1){if((st[length(st)] + p_pet[i]) > CC) {st2[i] <- CC} else {st2[i] <- (st[length(st)] + p_pet[i])}} else {
            if((st2[(i - 1)] + p_pet[i]) > CC) {st2[i] <- CC} else {st2[i] <- (st2[(i - 1)]  + p_pet[i])}
          }
        } else{st2[i] <- CC * exp(ppa2[i]/CC)}
      }
    }

    k1 <- cbind(ppa2, st2)
    k <- k1
    ppa <- ppa2
    st <- st2
    rm(ppa2, st2)
      }
    }
    return(c(ppa,st))
  }

  ret <- app(dataloop, floop, cores = ncpu)
  ppa <- ret[[1:12]]
  st <- ret[[13:24]]
  rm(ret)

  #changes in soil humidity (i_st)
  i_st <-  c((st[[1]]-st[[12]]),(st[[2:12]]-st[[1:11]]))

  #real evapotranspiration
  w1 <- pet*(p_pet >= 0)
  w2 <- (prec+abs(i_st))*(p_pet < 0)
  etr <- w1 + w2

  #humidity deficit
  dh <- etr - pet

  #humidity exceedance
  w1 <- p_pet >=0 # & st >= CC
  S <- (prec - (etr + i_st)) * w1

  #runoff
  datarf <- c(p_pet, st, S, CC)
  frunoff <- function(datarf){
    if(length(which(is.na(datarf))) > 0){
      r <- rep(NA, 12)
    } else{
      p_pet <- datarf[1:12]
      st <- datarf[13:24]
      S <- datarf[25:36]
      CC <- datarf[37]
      w <- which(p_pet >= 0 & st >= CC)
      if(length(w) > 0){
        if(w[1] == 1) {rec <- 1:12} else {rec <- c(w[1]:12, 1:(w[1]-1))} #order of months
        rnf <- rep(0,12)
        seguir <- TRUE
          # initial
          for(i in rec){
            if(i == 1) {
              rnf[i] <- 0.5 * (S[i] + rnf[12])
            } else {
              rnf[i] <- 0.5 * (S[i] + rnf[i-1])
              }
          }
        while(seguir == TRUE){
          # check identical values
          rnf2 <- rnf
          for(i in rec){
            if(i == 1) {
              rnf2[i] <- 0.5 * (S[i] + rnf2[12])
            } else {
              rnf2[i] <- 0.5 * (S[i] + rnf2[i-1])
            }
          }
          if(identical(rnf, rnf2)){seguir <- FALSE} else{
            rnf <- rnf2
          }
        }
        r <- rnf
      } else{ r <- rep(0,12)}
    }
    return(r)
  }
  #
  r <- app(datarf, frunoff, cores = ncpu)

  # runoff in %
  rP <- (r * 100) / prec
  rP[is.infinite(rP)] <- 0
  w1 <- (is.na(rP) | is.nan(rP))
  w1 <- w1*msk
  w1[w1==0] <- NA
  w1[w1==1] <- 0
  rP <- cover(rP,w1)

  # return of final values
  reswatbal <- c(temp, prec, pet, p_pet, ppa, st, i_st, etr, dh, S, r, rP)
  names(reswatbal)[1:12] <- paste0('temp',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[13:24] <- paste0('prec',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[25:36] <- paste0('PET',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[37:48] <- paste0('TEAW',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[49:60] <- paste0('PALW',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[61:72] <- paste0('ST',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[73:84] <- paste0('i_ST',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[85:96] <- paste0('RET',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[97:108] <- paste0('MD',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[109:120] <- paste0('ME',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[121:132] <- paste0('r',formatC(1:12,width = 2, flag = '0'))
  names(reswatbal)[133:144] <- paste0('rP',formatC(1:12,width = 2, flag = '0'))

  if(!is.null(path)){
    for(i in 1:dim(reswatbal)[3]){
      terra::writeRaster(reswatbal[[i]], paste0(path,'/',names(reswatbal)[[i]],'.tif'), overwrite=TRUE)
    }
  }
  return(reswatbal)
}
