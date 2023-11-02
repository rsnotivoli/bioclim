#' Bioclimatic classification (raster mode)
#'
#' @description Calculates bioclimatic classification based on bioclimatic balance.
#' @param temp SpatRaster object with 12 layers representing temperature from January to December.
#' @param prec SpatRaster object with 12 layers representing precipitation from January to December.
#' @param CC Field capacity. It can be numeric (1 value) or a SpatRaster object.
#' @param PET Potential evapotranspiration. Optional. It must be a SpatRaster object.
#' @param bh Water balance. Optional. It must be a SpatRaster object.
#' @param path Optional. Path (folder) where the output raster files and look-up-tables will be saved.
#' @param ncpu number of cores to use in calculation. If not provided, sequential mode is used (1 core).
#' @return SpatRaster with 3 variables ("TBR": Types of Bioclimatic Regime; "zonal": zonal units; "sub": bioclimatic regime subtypes).
#' @examples
#' \donttest{
#' wb <- terra::rast(wbRast)
#' btr <- biotypeRaster(bh = wb)
#'}
#' @importFrom utils write.table
#' @importFrom methods is
#' @export
#'


biotypeRaster <- function(temp=NULL, prec=NULL, CC=NULL, path=NULL, ncpu = 1, PET = NULL, bh = NULL){

  # If water balance is not provided, CC is computed
  if(is.null(bh)){
    # CC
    if(is.numeric(CC)){
      if(length(CC)>1) stop("CC must be a single number") else{
        # mask
        msk <- prec[[1]]/prec[[1]]
        CC <- msk*CC
      }
    } else if(!is(CC, 'SpatRaster')) stop("CC must be numeric or a SpatRaster")
  } else{
    if(!is(bh, 'SpatRaster')) stop("bh must be a SpatRaster")
    prec <- bh[[13:24]]
    temp <- bh[[1:12]]
    PET <- bh[[25:36]]
  }

  # check if PET is provided, then water balance, then computes it
  if(!is.null(PET)){
    if(!is(PET, 'SpatRaster')) stop("PET must be a SpatRaster")
  } else{
    if(is.null(bh)){
      message(paste0("Computing water balance ",'[',Sys.time(),']'))
      bh <- watbalRaster(temp, prec, CC, ncpu)
      PET <- bh[[25:36]]
    }
  }

    # message("Computing Thornthwaite’s index")
    subtype <- ithRaster(bh)

  # other required variables
  PVT <- sum(temp < 7.5)
  PVH <- sum((prec-0.2*PET) < 0)

  Tf <- min(temp)
  P <- sum(prec)


  # classification ----
  message(paste0("Computing Classification ",'[',Sys.time(),']'))
  zonal <- type <- P * 0
  # non tropical
  wnt <- (Tf < 18)
    wst <- wnt * (Tf > 7.5)
    zonal <- zonal + (wst * 1) # subtropical
      wPVH <- (PVH == 0) * wst
      type <- type + (((P > 2000) * wPVH) * 1) # Euritermo-Ombrophyllo
      type <- type + (((P <= 2000) * wPVH) * 2) # Euritermo-Mesophyllo
      wPVH <- (PVH != 0) * wst
        type <- type + (((PVH > 0 & PVH <= 4) * wPVH) * 3) # Euritermo-Tropophyllo
        type <- type + (((PVH > 4 & PVH < 8) * wPVH) * 4) # Euritermo-Xerophyllo
        type <- type + (((PVH >= 8) * wPVH) * 5) # Euritermo-Hyperxerophyllo
    wnt <- wnt * (Tf <= 7.5)
      wtm <- (PVT >= 1 & PVT < 6) * wnt
        zonal <- zonal + (wtm * 2) # Mid-latitudes
        w2 <- wtm * (PVH == 0)
          type <- type + (((P > 2000) * w2) * 6) # Cryo-Ombrophyllo
          type <- type + (((P <= 2000) * w2) * 7) # Cryo-Mesophyllo
        w2 <- wtm * (PVH != 0)
          type <- type + (((PVH > 0 & PVH <= 4) * w2) * 8) # Cryo-Tropophyllo
          type <- type + (((PVH > 4 & PVH < 8) * w2) * 9) # Cryo-Xerophyllo
          type <- type + (((PVH >= 8) * w2) * 10) # Cryo-Hyperxerophyllo
      wsp <- (PVT >= 6 & PVT < 10) * wnt
        zonal <- zonal + (wsp * 3) # Subpolar
        w2 <- wsp * (PVH == 0)
          type <- type + (((P > 1100) * w2) * 11) # Mesocryo-Ombrophyllo
          type <- type + (((P <= 1100) * w2) * 12) # Mesocryo-mesophyllo
        w2 <- wsp * (PVH != 0)
          type <- type + (((PVH > 0 & PVH <= 4) * w2) * 13) # Mesocryo-Tropophyllo
          type <- type + (((PVH > 4 & PVH < 8) * w2) * 14) # Mesocryo-Xerophyllo
          type <- type + (((PVH >= 8) * w2) * 15) # Mesocryo-Hyperxerophyllo
      wp <- (PVT >= 10) * wnt
        zonal <- zonal + (wp * 4) # Polar
        w2 <- wp * (PVH == 0)
          type <- type + (((P > 1100) * w2) * 16) # Hypercryo-Ombrophyllo
          type <- type + (((P <= 1100) * w2) * 17) # Hypercryo-mesophyllo
        w2 <- wp * (PVH != 0)
          type <- type + (((PVH > 0 & PVH <= 4) * w2) * 18) # Hypercryo-Tropophyllo
          type <- type + (((PVH > 4 & PVH < 8) * w2) * 19) # Hypercryo-Xerophyllo
          type <- type + (((PVH >= 8) * w2) * 20) # Hypercryo-Hyperxerophyllo
  wt <- (Tf >= 18)
    zonal <- zonal + (wt * 5) # Tropical
    w2 <- wt * (PVH == 0)
      type <- type + (((P > 2000) * w2) * 21) # Ombrophyllo
      type <- type + (((P <= 2000) * w2) * 22) # Mesophyllo
    w2 <- wt * (PVH > 0 & PVH <= 4)
      type <- type + (((P > 2000) * w2) * 23) # Ombro-Tropophyllo !
      type <- type + (((P <= 2000) * w2) * 24) # Tropophyllo
    w2 <- wt * (PVH > 4 & PVH < 8)
      type <- type + (((P > 2000) * w2) * 25) # Ombro-Xerophyllo !
      type <- type + (((P <= 2000) * w2) * 26) # Xerophyllo
    w2 <- wt * (PVH >= 8)
      type <- type + (((PVH >= 8) * w2) * 27) # Hyperxerophyllo

  # remove zeros
  zonal[zonal == 0] <- NA
  type[type == 0] <- NA


  zonalnames <- c('Subtropical', 'Mid latitudes', 'Subpolar', 'Polar', 'Tropical')
  subtypenames <- c('Arid', 'Semiarid', 'Dry humid', 'Moist humid',
                    'Low humid', 'Moderate humid', 'Highly humid', 'Very humid','Perhumid')
  typenames <- c('Euritermo-Ombrophyllo', 'Euritermo-Mesophyllo', 'Euritermo-Tropophyllo',
                 'Euritermo-Xerophyllo', 'Euritermo-Hyperxerophyllo', 'Cryo-Ombrophyllo',
                 'Cryo-Mesophyllo', 'Cryo-Tropophyllo', 'Cryo-Xerophyllo', 'Cryo-Hyperxerophyllo',
                 'Mesocryo-Ombrophyllo', 'Mesocryo-Mesophyllo', 'Mesocryo-Tropophyllo', 'Mesocryo-Xerophyllo',
                 'Mesocryo-Hyperxerophyllo', 'Hypercryo-Ombrophyllo', 'Hypercryo-Mesophyllo',
                 'Hypercryo-Tropophyllo', 'Hypercryo-Xerophyllo', 'Hypercryo-Hyperxerophyllo',
                 'Ombrophyllo', 'Mesophyllo', 'Ombro-Tropophyllo', 'Tropophyllo', 'Ombro-Xerophyllo',
                 'Xerophyllo', 'Hyperxerophyllo')
  df <- data.frame(ID = as.numeric(sapply(101:109, paste0, 1:27)), int = 1:243,
                   category = as.character(sapply(subtypenames, paste, typenames)))

  levels(zonal) <- data.frame(ID = 1:5, category = zonalnames)
  # plot(zonal)

  levels(type) <- data.frame(ID = 1:27, category = typenames)
  # plot(type)

    levels(subtype) <- data.frame(ID = 101:109, category = subtypenames)

    funion <- function(uni){
      p <- as.numeric(paste0(uni[1]+100, uni[2]))
      return(p)
    }
    ct <- suppressWarnings(app(c(subtype, type), funion))
    d <- df[match(unique(ct)[[1]], df$ID), ]
    # m <- cbind(unique(ct)[[1]], 1:length(unique(ct)[[1]]))
    m <- cbind(unique(ct)[[1]], d[,2])
    ct <- terra::classify(ct, m)
    levels(ct) <- data.frame(ID = d[,2], category = d[,3])
    res <- c(zonal, type, ct)
    names(res) <- c('Zonal', 'TBR', 'Subtype')
    
  # LOOK-UP-TABLES
  # table zonal
    tabzon <- data.frame(ID = 1:5, category = zonalnames,
                         hex = c("#FFC000","#92D050","#00B0F0","#B4C6E7","#9966FF"))
    vals <- unique(values(zonal))
    if(length(is.na(vals))>0) vals <- vals[-which(is.na(vals))]
    tabzon <- tabzon[match(vals, tabzon$ID),]
    
  # table TBRs
    tabtbr <- data.frame(ID = 1:27, category = typenames,
                         hex = c("#C276B6","#FA8C6B","#FFA767","#FAC488","#FFF685",
                                 "#3A005D","#00508F","#006B71","#00864B","#529F2B",
                                 "#5C2799","#0070C0","#03A7B0","#5DC894","#AEDB88",
                                 "#E3CBE5","#ADC3EC","#BDE5E5","#BAE4D4","#E4F1D3",
                                 "#6F0056","#9A1B15","#E600A9","#9D5000","#FF73DF",
                                 "#FFA614","#FFFF00"))
    vals <- unique(values(type))
    if(length(is.na(vals))>0) vals <- vals[-which(is.na(vals))]
    tabtbr <- tabtbr[match(vals, tabtbr$ID),]
    
  # table subtypes
    tabsub <- data.frame(ID = 1:243, category = df$category,
                         hex = c("#F6DDCC","#FAE5D3","#FDEBD0","#FCF3CF","#F0F4C3",
                                 "#66FFCC","#DCEDC8","#C8E6C9","#D4EFDF","#D5F5E3",
                                 "#64FF57","#B2DFDB","#D0ECE7","#D1F2EB","#B2EBF2",
                                 "#C5CAE9","#D4E6F0","#D6EAF8","#BBDEFB","#B3E5FC",
                                 "#D1C4E9","#CDBEE7","#FFB9D0","#F2D7D5","#F3BFBF",
                                 "#FADBD8","#FFF9C4","#EDBB99","#F5CBA7","#FAD7A0",
                                 "#F9E79F","#E6EE9C","#66FF99","#C5E1A5","#A5D6A7",
                                 "#A9DFBF","#ABEBC6","#46F052","#80CBC4","#A2D9CE",
                                 "#A3E4D7","#80DEEA","#9FA8DA","#A9CCE3","#AED6F1",
                                 "#90CAF9","#81D4FA","#B39DDB","#CE93D8","#FF93B7",
                                 "#E6B0AA","#EEA4A4","#F5B7B1","#FFF59D","#E59866",
                                 "#F0B27A","#F8C471","#F7DC6F","#DCE775","#66FF66",
                                 "#AED581","#81C784","#7DCEA0","#82E0AA","#1ED24D",
                                 "#4DB6AC","#73C6B6","#76D7C4","#4DD0E1","#7986CB",
                                 "#7FB3D5","#85C1E9","#64B5F6","#4FC3F7","#9575CD",
                                 "#BA68C8","#FF6699","#D98880","#E57373","#F1948A",
                                 "#FFF176","#DC7633","#EB984E","#F5B041","#F4D03F",
                                 "#D4E157","#00FF00","#9CCC65","#66BB6A","#52BE80",
                                 "#58D68D","#00B448","#26A69A","#45B39D","#48C9B0",
                                 "#26C6DA","#5C6BC0","#5499C2","#5DADE2","#42A5F5",
                                 "#29B6F6","#7E57C2","#AB47BC","#EC407A","#CD6155",
                                 "#EF5350","#EC7063","#FCE112","#D35400","#E67E22",
                                 "#F39C12","#F1C40F","#CDDC39","#33CC33","#8BC34A",
                                 "#4CAF50","#27AE60","#2ECC71","#009643","#009688",
                                 "#16A085","#1ABC9C","#00BCD4","#3F51B5","#2980B9",
                                 "#3498DB","#2196F3","#03A9F4","#673AB7","#9C27B0",
                                 "#E91E63","#C0392B","#F44336","#E74C3C","#FACE12",
                                 "#BA4A00","#CA6F1E","#D68910","#D4AC0D","#C0CA33",
                                 "#00CC00","#7CB342","#43A047","#229954","#28B463",
                                 "#00783E","#00897B","#138D75","#17A589","#00ACC1",
                                 "#3949AB","#2471A3","#2E86C1","#1E88E5","#039BE5",
                                 "#5E35B1","#8E24AA","#D81B60","#A93226","#E53935",
                                 "#CB4335","#F9B709","#A04000","#AF601A","#B9770E",
                                 "#B7950B","#AFB42B","#008000","#689F38","#388E3C",
                                 "#1E8449","#239B56","#005D39","#00796B","#117A65",
                                 "#148F77","#0097A7","#303F9F","#1F618D","#2874A6",
                                 "#1976D2","#0288D1","#512DA8","#7B1FA2","#C2185B",
                                 "#922B21","#D32F2F","#B03A2E","#CD9403","#873600",
                                 "#935116","#9C640C","#9A7D0A","#9E9D24","#006600",
                                 "#558B2F","#2E7D32","#196F3D","#1D8348","#003E34",
                                 "#00695C","#0E6655","#117864","#00695C","#283593",
                                 "#1A5276","#21618C","#1565C0","#0277BD","#4527A0",
                                 "#6A1B9A","#AD1457","#7B241C","#C62828","#943126",
                                 "#AC7C02","#6E2C00","#784212","#7E5109","#7D6608",
                                 "#82777B","#154B19","#33691E","#1B5E20","#145A32",
                                 "#186A3B","#003E34","#97FF97","#0B5345","#0E6251",
                                 "#006064","#1A237E","#154360","#1B4F72","#0D47A1",
                                 "#01579B","#311B92","#4A148C","#880E4F","#641E16",
                                 "#B71C1C","#78281F","#AE5C02"))
    vals <- unique(values(ct))
    if(length(is.na(vals))>0) vals <- vals[-which(is.na(vals))]
    tabsub <- tabsub[match(vals, tabsub$ID),]
    
    
  if(!is.null(path)){
    for(i in 1:dim(res)[3]){
      terra::writeRaster(res[[i]], paste0(path,'/',names(res)[[i]],'.tif'), overwrite=TRUE)
    }
    write.table(tabzon, paste0(path,'/','lut_zonal.txt'), quote=F, sep='\t',row.names = F)
    write.table(tabtbr, paste0(path,'/','lut_TBR.txt'), quote=F, sep='\t',row.names = F)
    write.table(tabsub, paste0(path,'/','lut_subtypes.txt'), quote=F, sep='\t',row.names = F)
  }

  message(paste0("End ",'[',Sys.time(),']'))
  return(res)
}
