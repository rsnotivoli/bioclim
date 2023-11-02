#' Computation of Bioclimatic Balance
#'
#' @description Computes bioclimatic balance from water balance.
#' @param balhid Water balance.
#' @param CC Field capacity. It depends on water retention capacity and depth of roots. 400 as default value.
#' @return data frame with 12 variables: 'p', 'Tm', 'PET', 'e', 'D', 'S', 'Cd', 'T_75', 'B', 'b', 'bl', 'bc'.
#' @examples
#' wb <- watbal(t = rnorm(12, 18, 6), p = rnorm(12, 50, 30), lat = 35, CC = 400)
#' biobal(wb, 400)
#' @export
#'

biobal <- function(balhid, CC){
  bh <- balhid[1:12, ]
  #
  bb <- matrix(NA, ncol = 20, nrow = 12)
  colnames(bb) <-  c('p', 'Tm', 'PET', 'e', 'D', 'S',
                     's_e_D', 'sum_s', 'c_D_e', 'sum_c',
                     'Q', 'x', 'E_e', 'D_e', 'Cd', 'T_75',
                     'B', 'b', 'bl', 'bc')
  rownames(bb) <- month.name

  #precip util
  Ce <- 1 - (bh$rP / 100)
  bb[, 1] <- Ce * bh$P

  #tmed
  bb[, 2] <- bh$T

  #pet
  bb[, 3] <- bh$PET

  #e
  bb[, 4] <- bh$PET / 5

  #D y S

  #ini
  S <- D <- e_D <- iS <- rep(NA, 12)
  for(i in 1:12){
    if(i == 1) {
      D[i] <- bb[i, 1] #D
      if(bb[i, 4] > D[i]) e_D[i] <- bb[i, 4] - D[i] else e_D[i] <- 0 #e_D
      if(e_D[i] > 0) iS[i] <- e_D[i] else iS[i] <- 0 #iS
      if(D[i] - bh$PET[i] < 0) S[i] <- 0 else S[i] <- D[i] - bh$PET[i] #S
    } else{
      D[i] <- bb[i, 1] + S[i-1]
      if(bb[i, 4] > D[i]) e_D[i] <- bb[i, 4] - D[i] else e_D[i] <- 0 #e_D
      if(e_D[i] > 0) iS[i] <- e_D[i] + iS[i-1] else iS[i] <- 0 #iS
      if(D[i] - bh$PET[i] < 0) S[i] <- 0 else S[i] <- D[i] - bh$PET[i]
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
    if(i == 1) D[i] <- bb[i, 1] else D[i] <- bb[i, 1] + S[i-1]#D
    if(bb[i, 4] > D[i]) e_D[i] <- bb[i, 4] - D[i] else e_D[i] <- 0 #e_D
    if(e_D[i] > 0) iS[i] <- e_D[i] + iS[wf] else iS[i] <- 0 #iS
    #D_e
    if(i == 1){
      if(e_D[i] <= 0){
        if(e_Dold[wf] <= 0) {if(e_Dold[wf] <= 0) D_e[i] <- 0 else D_e[i] <- D[i] - bb[i, 4]} else D_e[i] <- D[i] - bb[i,4]
      } else D_e[i] <- 0
      if(D_e[i] > 0) sC[i] <- D_e[i] else sC[i] <- 0
    } else{
      if(e_D[i] <= 0){
        if(e_D[wf] <= 0) {if(sC[wf] >= maxis[wf]) D_e[i] <- 0 else if(sC[wf] <= 0) {D_e[i] <- 0} else {D_e[i] <- D[i] - bb[i,4]}} else D_e[i] <- D[i] - bb[i,4]
      } else D_e[i] <- 0
      if(D_e[i] > 0) sC[i] <- D_e[i] + sC[wf] else sC[i] <- 0
    }

    #maxis
    if(e_D[i] == 0) maxis[i] <- max(iS) else maxis[i] <- 0
    #Q
    if(sC[i] >= maxis[i]) Q[i] <- sC[i] - maxis[i] else Q[i] <- 0

    #S
    if(Q[i] > 0){
      cal <- (Q[i] * (D[i] - bh$PET[i]) / D_e[i])
      if(length(which(is.infinite(cal)))>0) cal[which(is.infinite(cal))] <- 0
      if(cal < 0) S[i] <- 0 else{
        if(cal >= CC) S[i] <- CC else S[i] <- cal
      }
    } else{
      if((D[i] - bh$PET[i]) > 0) {
        if((D[i] - bh$PET[i]) > CC) S[i] <- CC else S[i] <- D[i] - bh$PET[i]
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
    if(i == 1) D[i] <- bb[i, 1] + Sold[wf] else D[i] <- bb[i, 1] + S[i-1]#D
    if(bb[i, 4] > D[i]) e_D[i] <- bb[i, 4] - D[i] else e_D[i] <- 0 #e_D
    if(e_D[i] > 0) iS[i] <- e_D[i] + iS[wf] else iS[i] <- 0 #iS
    #maxis
    if(e_D[i] == 0) maxis[i] <- max(iS) else maxis[i] <- 0
    #D_e
    if(i == 1){
      if(e_D[i] <= 0){
        if(e_Dold[wf] <= 0) {if(sCold[wf] >= maxisold[wf]) D_e[i] <- 0 else {if(sCold[wf] <= 0) D_e[i] <- 0 else D_e[i] <- D[i+1] - bb[i, 4]}} else D_e[i] <- D[i] - bb[i,4]
      } else D_e[i] <- 0
      if(D_e[i] > 0) sC[i] <- D_e[i] else sC[i] <- 0
    } else{
      if(e_D[i] <= 0){
        if(e_D[wf] <= 0) {if(sC[wf] >= maxis[wf]) D_e[i] <- 0 else if(sC[wf] <= 0) {D_e[i] <- 0} else {D_e[i] <- D[i] - bb[i,4]}} else D_e[i] <- D[i] - bb[i,4]
      } else D_e[i] <- 0
      if(D_e[i] > 0) sC[i] <- D_e[i] + sC[wf] else sC[i] <- 0
    }


    #Q
    if(sC[i] >= maxis[i]) Q[i] <- sC[i] - maxis[i] else Q[i] <- 0

    #S
    if(Q[i] > 0){
      cal <- (Q[i] * (D[i] - bh$PET[i]) / D_e[i])
      if(length(which(is.infinite(cal)))>0) cal[which(is.infinite(cal))] <- 0
      if(cal < 0) S[i] <- 0 else{
        if(cal >= CC) S[i] <- CC else S[i] <- cal
      }
    } else{
      if((D[i] - bh$PET[i]) > 0) {
        if((D[i] - bh$PET[i]) > CC) S[i] <- CC else S[i] <- D[i] - bh$PET[i]
      } else S[i] <- 0
    }
  }

  bb[, 5] <- D
  bb[, 6] <- S
  bb[, 7] <- e_D
  bb[, 8] <- iS
  bb[, 9] <- D_e
  bb[, 10] <- sC
  bb[, 11] <- Q

  #x
  bb[, 12] <- Q / (D - bb[,4])
  if(length(which(is.infinite(bh[,12])))>0) bb[which(is.infinite(bh[,12])), 12] <- 0

  #E_e
  bb[, 13] <- bb[, 3] - bb[, 4]

  #D_e
  bb[, 14] <- bb[, 5] - bb[, 4]

  #Cd
  cd <- rep(NA, 12)
  for(i in 1:12){
    if(bb[i, 3] - bb[i, 4] == 0){
      cd[i] <- bb[i, 5]- bb[i, 4]
    } else{
      interm <- (bb[i, 5]- bb[i, 4]) / (bb[i, 3] - bb[i, 4])
      if(length(which(is.infinite(interm)))>0) interm[which(is.infinite(interm))] <- 0
      cd[i] <- interm
      }
  }
  bb[, 15] <- cd

  #t75
  bb[, 16] <- bb[, 2] - 7.5

  #B
  bb[, 17] <- bb[, 16] / 5

  #b
  b <- rep(NA, 12)
  for(i in 1:12){ if(bb[i, 15] > 1) b[i] <- bb[i, 17] else b[i] <- bb[i, 15] * bb[i, 17]}
  bb[, 18] <- b

  #bl
  bl <- rep(NA, 12)
  for(i in 1:12){
    if(bb[i, 10] > 0){
      if(bb[i, 12] > 0) bl[i] <- bb[i, 12] * bb[i, 18] else bl[i] <- 0
    } else{
      if((bb[i, 17] > 0 & bb[i, 18] < 0) | (bb[i, 17] < 0 & bb[i, 18] > 0)) bl[i] <- 0 else bl[i] <- bb[i, 18]
    }
  }

  bb[, 19] <- bl

  #bc
  bc <- rep(NA, 12)
  for(i in 1:12){ if(bb[i, 10] > 0) bc[i] <- bb[i, 18] - bb[i, 19] else bc[i] <- 0}
  bb[, 20] <- bc
  
  # 
  bb <- bb[,c(1,2,3,4,5,6,15,16,17,18,19,20)]

  tt <- data.frame(p = sum(bb[,1]), Tm = mean(bb[,2]), 
                   PET = sum(bb[,3]), e = sum(bb[,4]),
                   D = NA, S = NA, Cd =  sum(bb[,7]),
                   T_75 =  sum(bb[,8]), B =  sum(bb[,9]), 
                   b =  sum(bb[,10]), bl =  sum(bb[,11]),
                   bc =  sum(bb[,12]))
  bb <- rbind(bb, tt)
  rownames(bb)[13] <- 'TOTAL'
  
  colnames(bb) <-  c('AIP', 'T', 'PET', 'RE', 'AW', 'S',
                     'CWA', 'T75', 'PBI', 'RBI', 'FBI', 'CBI')
  bb
}
