#' Computation of Bioclimatic Intensities (raster mode)
#'
#' @description Computes bioclimatic intensities from bioclimatic balance.
#' @param bb Bioclimatic balance in raster format.
#' @param path Optional. Path (folder) where the output raster files will be saved.
#' @return SpatRaster with 120 layers corresponding to the 12 monthly values of "IBPc","IBCc","IBLc","IBRc","IBSc","IBPf","IBCf","IBLf","IBRf","IBSf".
#' @examples
#' \donttest{
#' bb <- terra::rast(bbRast)
#' bi <- biointRaster(bb, path=NULL)
#' }
#' @export
#'

biointRaster <- function(bb, path=NULL){

  B <- bb[[1:12]]
  b <- bb[[13:24]]
  bc <- bb[[25:36]]
  bl <- bb[[37:48]]
  # PBIw
  PBIw <- B * (B > 0)
  names(PBIw) <- paste0('PBIw', formatC(1:12,width = 2, flag = '0'))
  # CBIw
  CBIw <- bc * (bc > 0)
  names(CBIw) <- paste0('CBIw', formatC(1:12,width = 2, flag = '0'))
  # FBIw
  w1 <- (B > 0 & b > 0)
  FBIw <- (bl * w1) * (bl > 0)
  names(FBIw) <- paste0('FBIw', formatC(1:12,width = 2, flag = '0'))

  # FBIw_g
  FBIw_g <- (FBIw > CBIw) * (FBIw - CBIw)
  names(FBIw_g) <- paste0('FBIw_g', formatC(1:12,width = 2, flag = '0'))

  # PBIw_g
  PBIw_g <- PBIw - (CBIw + FBIw_g)
  names(PBIw_g) <- paste0('PBIw_g', formatC(1:12,width = 2, flag = '0'))

  # RBIw
  RBIw <- FBIw - CBIw
  names(RBIw) <- paste0('RBIw', formatC(1:12,width = 2, flag = '0'))

  # DBIw
  DBIw <- b * (B > 0 & b < 0)
  names(DBIw) <- paste0('DBIw', formatC(1:12,width = 2, flag = '0'))

  # PBIc
  PBIc <- B * (B < 0)
  names(PBIc) <- paste0('PBIc', formatC(1:12,width = 2, flag = '0'))

  # CBIc
  CBIc <- bc * (bc[[3]] < 0)
  names(CBIc) <- paste0('CBIc', formatC(1:12,width = 2, flag = '0'))

  # FBIc
  w1 <- (B < 0 & b < 0)
  FBIc <- ((bl * (bl < 0)) * w1) + ((b * (bl > 0)) * w1)
  names(FBIc) <- paste0('FBIc', formatC(1:12,width = 2, flag = '0'))

  # FBIc_g
  FBIc_g <- (CBIc > FBIc) * (FBIc - CBIc)
  names(FBIc_g) <- paste0('FBIc_g', formatC(1:12,width = 2, flag = '0'))

  # PBIc_g
  PBIc_g <- PBIc + (CBIc - FBIc_g)
  names(PBIc_g) <- paste0('PBIc_g', formatC(1:12,width = 2, flag = '0'))

  # RBIc
  RBIc <- FBIc + CBIc
  names(RBIc) <- paste0('RBIc', formatC(1:12,width = 2, flag = '0'))

  # DBIc
  w1 <- (B < 0)
  DBIc <- (b * (b > 0)) * w1
  names(DBIc) <- paste0('DBIc', formatC(1:12,width = 2, flag = '0'))

  intens <- c(PBIw, CBIw, FBIw, RBIw, DBIw,
              PBIc, CBIc, FBIc, RBIc, DBIc)

  if(!is.null(path)){
    for(i in 1:dim(intens)[3]){
      terra::writeRaster(intens[[i]], paste0(path,'/',names(intens)[[i]],'.tif'), overwrite=TRUE)
    }
  }
  return(intens)
}
