#' Function to calculate water balance
#'
#' @description Computes water balance from temperature and precipitation data.
#' @param t Monthly average temperature data (12 nueric values).
#' @param p Monthly average precipitation data (12 nueric values).
#' @param lat Latitude in degrees. For southern latitudes use negative values.
#' @param CC Field capacity. It depends on water retention capacity and depth of roots. Use 400 as default value.
#' @return data frame with 12 variables: 'Tmp', 'Pcp', 'PET', 'P_PET', 'ppa', 'ST', 'i_ST', 'ETR', 'Dh', 'S', 'r', 'rP'.
#' @examples
#' watbal(t = rnorm(12, 18, 6), p = rnorm(12, 50, 30),
#'        lat = 35, CC = 400)
#' @export
#'


watbal <- function(t, p, lat, CC){
  #
  balhid <- matrix(NA, ncol = 12, nrow = 12)
  colnames(balhid) <-  c('Tmp', 'Pcp', 'PET', 'P_PET', 'ppa', 'ST',
                         'i_ST', 'ETR', 'Dh', 'S', 'r', 'rP')
  rownames(balhid) <- month.name

  #temperature
  balhid[, 1] <- t
  #precipitation
  balhid[, 2] <- p
  #potential evapotranspiration (PET)
  balhid[, 3] <- thornthwaite(t, lat)
  #P-PET
  balhid[, 4] <- balhid[, 2] - balhid[, 3]

  #ppa (cumulated potential leakage)
  p_pet <- balhid[, 4]
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

  #ST (soil pore water)
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

  #loop to update values until no changes
  k <- cbind(ppa, st)
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
    # if(identical(k, k1)) {seg <- FALSE} else {
    k <- k1
    ppa <- ppa2
    st <- st2
    rm(ppa2, st2)
    # }

  }

  balhid[, 5] <- k[, 1]
  balhid[, 6] <- k[, 2]

  #changes in soil humidity
  balhid[, 7] <-  c((balhid[1, 6] - balhid[12, 6]), (balhid[2:12, 6] - balhid[1:11, 6]))

  #real evapotranspiration
  w1 <- which(p_pet >= 0)
  balhid[w1, 8] <- balhid[w1, 3]
  w2 <- which(p_pet < 0)
  balhid[w2, 8] <- balhid[w2, 2] + abs(balhid[w2, 7])

  #humidity deficit
  balhid[, 9] <- balhid[, 8] - balhid[, 3]

  #humidity excedance
  balhid[, 10] <- 0
  w1 <- which(p_pet >= 0) #& balhid[, 6] >= CC)
  balhid[w1, 10] <- balhid[w1, 2] - (balhid[w1, 8] + balhid[w1, 7])

  #runoff
  balhid[, 11] <- 0
  w <- which(p_pet >= 0 & balhid[, 6] >= CC)
  if(length(w) > 0){
    if(w[1] == 1) {rec <- 1:12} else {rec <- c(w[1]:12, 1:(w[1]-1))} #order of months
    rnf <- rep(0,12)
    seguir <- TRUE
    # initial
    for(i in rec){
      if(i == 1) {
        rnf[i] <- 0.5 * (balhid[i,10] + rnf[12])
      } else {
        rnf[i] <- 0.5 * (balhid[i,10] + rnf[i-1])
      }
    }
    while(seguir == TRUE){
      # check identical values
      rnf2 <- rnf
      for(i in rec){
        if(i == 1) {
          rnf2[i] <- 0.5 * (balhid[i,10] + rnf2[12])
        } else {
          rnf2[i] <- 0.5 * (balhid[i,10] + rnf2[i-1])
        }
      }
      if(identical(rnf, rnf2)){seguir <- FALSE} else{
        rnf <- rnf2
      }
    }
    r <- rnf
  } else{ r <- rep(0,12)}
    balhid[, 11] <- r


  #runoff in %
  balhid[, 12] <- (balhid[, 11] * 100) / balhid[, 2]
  w <- which(is.nan(balhid[,12]) | is.infinite(balhid[,12]))
  if(length(w)>0) balhid[w,12] <- 0

  tt <- data.frame(Tmp = mean(balhid[, 1]), Pcp = sum(balhid[, 2]), PET = sum(balhid[, 3]),
                   P_PET = sum(balhid[, 4]), ppa = NA, ST = NA, i_ST = NA, ETR = sum(balhid[, 8]),
                   Dh = sum(balhid[, 9]), S = sum(balhid[, 10]), r = sum(balhid[, 11]), rP = mean(balhid[, 12]))
  balhid <- rbind(balhid, tt)
  rownames(balhid)[13] <- 'TOTAL'
  
  colnames(balhid) <-  c('T', 'P', 'PET', 'TEAW', 'PALW', 'ST',
                         'i_ST', 'RET', 'MD', 'ME', 'r', 'rP')
  
  balhid
}

