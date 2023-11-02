#' Computation of Bioclimatic Balance (raster mode)
#'
#' @description Computes bioclimatic balance from water balance in raster format.
#' @param bh Water balance in raster format.
#' @param CC Field capacity. It depends on water retention capacity and depth of roots. 400 as default value. It can be a SpatRaster layer.
#' @param path Optional. Path (folder) where the output raster files and look-up-tables will be saved.
#' @param ncpu Number of CPUs to use. By default, sequential mode (1 cpu) is used.
#' @return SpatRaster with 48 layers corresponding to the 12 monthly values of 'B', 'b','bc','bl'.
#' @import terra
#' @importFrom methods is
#' @examples
#' \donttest{
#' wb <- terra::rast(wbRast)
#' bb <- biobalRaster(wb, CC = 400, path=NULL, ncpu = 2)
#' }
#' @export
#'
biobalRaster <- function(bh, CC, path=NULL, ncpu = 1){

  # mask
  msk <- bh[[1]]/bh[[1]]
  # CC
  if(is.numeric(CC)){
    if(length(CC)>1) stop("CC must be a single number") else{
      CC <- msk*CC
    }
  } else if(!is(CC, 'SpatRaster')) stop("CC must be numeric or a SpatRaster")

  # precip util
  p <- (1-(bh[[133:144]]/100))*bh[[13:24]]
  names(p) <- paste0('p',formatC(1:12,width = 2, flag = '0'))

  # e
  e <- bh[[25:36]] / 5
  names(e) <- paste0('e',formatC(1:12,width = 2, flag = '0'))

  #D and S
  datavars <- c(p,e,bh[[25:36]], CC)
  fvars <- function(datavars){
    if(length(which(is.na(datavars))) > 0){
      D <- S <- e_D <- iS <- D_e <- sC <- Q <- rep(NA, 12)
    } else{
      p <- datavars[1:12]
      e <- datavars[13:24]
      pet <- datavars[25:36]
      CC <- datavars[37]
      #ini
      S <- D <- e_D <- iS <- rep(NA, 12)
      for(i in 1:12){
        if(i == 1) {
          D[i] <- p[i] #D
          if(e[i] > D[i]) e_D[i] <- e[i] - D[i] else e_D[i] <- 0 #e_D
          if(e_D[i] > 0) iS[i] <- e_D[i] else iS[i] <- 0 #iS
          if(D[i] - pet[i] < 0) S[i] <- 0 else S[i] <- D[i] - pet[i] #S
        } else{
          D[i] <- p[i] + S[i-1]
          if(e[i] > D[i]) e_D[i] <- e[i] - D[i] else e_D[i] <- 0 #e_D
          if(e_D[i] > 0) iS[i] <- e_D[i] + iS[i-1] else iS[i] <- 0 #iS
          if(D[i] - pet[i] < 0) S[i] <- 0 else S[i] <- D[i] - pet[i]
        }
      }

      #ini2
      iSold <-  iS
      D_e <- sC <- maxis <- Hm <- Q <- rep(NA, 12)
      e_Dold <- e_D
      Sold <- S
      Dold <- D
      for(i in 1:12){
        if(i == 1) wf <- 12 else wf <- i-1
        if(i == 1) D[i] <- p[i] else D[i] <- p[i] + S[i-1]#D
        if(e[i] > D[i]) e_D[i] <- e[i] - D[i] else e_D[i] <- 0 #e_D
        if(e_D[i] > 0) iS[i] <- e_D[i] + iS[wf] else iS[i] <- 0 #iS
        #D_e
        if(i == 1){
          if(e_D[i] <= 0){
            if(e_Dold[wf] <= 0) {if(e_Dold[wf] <= 0) D_e[i] <- 0 else D_e[i] <- D[i] - e[i]} else D_e[i] <- D[i] -e[i]
          } else D_e[i] <- 0
          if(D_e[i] > 0) sC[i] <- D_e[i] else sC[i] <- 0
        } else{
          if(e_D[i] <= 0){
            if(e_D[wf] <= 0) {if(sC[wf] >= maxis[wf]) D_e[i] <- 0 else if(sC[wf] <= 0) {D_e[i] <- 0} else {D_e[i] <- D[i] - e[i]}} else D_e[i] <- D[i] - e[i]
          } else D_e[i] <- 0
          if(D_e[i] > 0) sC[i] <- D_e[i] + sC[wf] else sC[i] <- 0
        }

        #maxis
        if(e_D[i] == 0) maxis[i] <- max(iS) else maxis[i] <- 0
        #Q
        if(sC[i] >= maxis[i]) Q[i] <- sC[i] - maxis[i] else Q[i] <- 0

        #S
        if(Q[i] > 0){
          cal <- (Q[i] * (D[i] - pet[i]) / D_e[i])
          if(length(which(is.infinite(cal)))>0) cal[which(is.infinite(cal))] <- 0
          if(cal < 0) S[i] <- 0 else{
            if(cal >= CC) S[i] <- CC else S[i] <- cal
          }
        } else{
          if((D[i] - pet[i]) > 0) {
            if((D[i] - pet[i]) > CC) S[i] <- CC else S[i] <- D[i] - pet[i]
          } else S[i] <- 0
        }
      }

      #tab.final
      iSold <-  iS
      e_Dold <- e_D
      Sold <- S
      Dold <- D
      maxisold <- maxis
      sCold <- sC
      D_e <- sC <- maxis <- Hm <- Q <- rep(NA, 12)
      for(i in 1:12){
        if(i == 1) wf <- 12 else wf <- i-1
        if(i == 1) D[i] <- p[i] + Sold[wf] else D[i] <- p[i] + S[i-1]#D
        if(e[i] > D[i]) e_D[i] <- e[i] - D[i] else e_D[i] <- 0 #e_D
        if(e_D[i] > 0) iS[i] <- e_D[i] + iS[wf] else iS[i] <- 0 #iS
        #maxis
        if(e_D[i] == 0) maxis[i] <- max(iS) else maxis[i] <- 0
        #D_e
        if(i == 1){
          if(e_D[i] <= 0){
            if(e_Dold[wf] <= 0) {if(sCold[wf] >= maxisold[wf]) D_e[i] <- 0 else {if(sCold[wf] <= 0) D_e[i] <- 0 else D_e[i] <- D[i+1] - e[i]}} else D_e[i] <- D[i] - e[i]
          } else D_e[i] <- 0
          if(D_e[i] > 0) sC[i] <- D_e[i] else sC[i] <- 0
        } else{
          if(e_D[i] <= 0){
            if(e_D[wf] <= 0) {if(sC[wf] >= maxis[wf]) D_e[i] <- 0 else if(sC[wf] <= 0) {D_e[i] <- 0} else {D_e[i] <- D[i] - e[i]}} else D_e[i] <- D[i] - e[i]
          } else D_e[i] <- 0
          if(D_e[i] > 0) sC[i] <- D_e[i] + sC[wf] else sC[i] <- 0
        }


        #Q
        if(sC[i] >= maxis[i]) Q[i] <- sC[i] - maxis[i] else Q[i] <- 0

        #S
        if(Q[i] > 0){
          cal <- (Q[i] * (D[i] - pet[i]) / D_e[i])
          if(length(which(is.infinite(cal)))>0) cal[which(is.infinite(cal))] <- 0
          if(cal < 0) S[i] <- 0 else{
            if(cal >= CC) S[i] <- CC else S[i] <- cal
          }
        } else{
          if((D[i] - pet[i]) > 0) {
            if((D[i] - pet[i]) > CC) S[i] <- CC else S[i] <- D[i] - pet[i]
          } else S[i] <- 0
        }
      }
    }
    return(c(D, S, e_D, iS, D_e, sC, Q))
  }

  vars <- app(datavars, fvars, cores = ncpu)
  names(vars)[1:12] <- paste0('D',formatC(1:12,width = 2, flag = '0'))
  names(vars)[13:24] <- paste0('S',formatC(1:12,width = 2, flag = '0'))
  names(vars)[25:36] <- paste0('e_D',formatC(1:12,width = 2, flag = '0'))
  names(vars)[37:48] <- paste0('iS',formatC(1:12,width = 2, flag = '0'))
  names(vars)[49:60] <- paste0('c_D_e',formatC(1:12,width = 2, flag = '0'))
  names(vars)[61:72] <- paste0('sC',formatC(1:12,width = 2, flag = '0'))
  names(vars)[73:84] <- paste0('Q',formatC(1:12,width = 2, flag = '0'))

  #x
  x <- vars[[73:84]] / (vars[[1:12]] - e)
  x[is.infinite(x)] <- 0
  names(x) <- paste0('x',formatC(1:12,width = 2, flag = '0'))

  #E_e
  E_e <- bh[[25:36]] - e
  names(E_e) <- paste0('E_e',formatC(1:12,width = 2, flag = '0'))

  #D_e
  D_e <- vars[[1:12]] - e
  names(D_e) <- paste0('D_e',formatC(1:12,width = 2, flag = '0'))

  #Cd
  w1 <- (bh[[25:36]] - e) == 0
  cd1 <- (vars[[1:12]] - e) * w1
  w2 <- (bh[[25:36]] - e) != 0
  interm <- ((vars[[1:12]] - e) / (bh[[25:36]] - e))
  interm[is.infinite(interm)] <- 0 # we need this to avoid 0/0 = inf
  cd2 <- interm * w2
  Cd <- cd1 + cd2
  names(Cd) <- paste0('CWA',formatC(1:12,width = 2, flag = '0'))

  #t75
  t75 <- bh[[1:12]] - 7.5
  names(t75) <- paste0('T75',formatC(1:12,width = 2, flag = '0'))

  #PBI
  B <- t75 / 5
  names(B) <- paste0('PBI',formatC(1:12,width = 2, flag = '0'))

  #RBI
  w1 <- Cd >1
  b1 <- B * w1
  w2 <- Cd <=1
  b2 <- (Cd * B) * w2
  b <- b1+b2
  names(b) <- paste0('RBI',formatC(1:12,width = 2, flag = '0'))

  #FBI
  databl <- c(vars[[61:72]], x, B, b)
  fbl <- function(databl){
    if(length(which(is.na(databl))) > 0){
      bl <- rep(NA, 12)
    } else{
      sC <- databl[1:12]
      x <- databl[13:24]
      B <- databl[25:36]
      b <- databl[37:48]
      bl <- rep(NA, 12)
      for(i in 1:12){
        if(sC[i] > 0){
          if(x[i] > 0) bl[i] <- x[i] * b[i] else bl[i] <- 0
        } else{
          if((B[i] > 0 & b[i] < 0) | (B[i] < 0 & b[i] > 0)) bl[i] <- 0 else bl[i] <- b[i]
        }
      }
    }
    return(bl)
  }

  bl <- app(databl, fbl, cores=ncpu)
  names(bl) <- paste0('FBI',formatC(1:12,width = 2, flag = '0'))

  #bCBI
  w1 <- vars[[61:72]] > 0
  w2 <- vars[[61:72]] <= 0
  bc <- ((b - bl) * w1) + ((b * 0) * w2)
  names(bc) <- paste0('CBI',formatC(1:12,width = 2, flag = '0'))

  # only bioclimatic intensities
  resbiobal <- c(B, b, bc, bl)


  if(!is.null(path)){
    for(i in 1:dim(resbiobal)[3]){
      terra::writeRaster(resbiobal[[i]], paste0(path,'/',names(resbiobal)[[i]],'.tif'), overwrite=TRUE)
    }
  }
  
  return(resbiobal)

}
