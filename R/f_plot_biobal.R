#' Function to plot bioclimatic balance
#'
#' @description Function to plot bioclimatic balance.
#' @param intens bioclimatic intensities in data.frame format from bioint() function.
#' @return Plot of bioclimatic balance
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme_minimal theme guides xlab ylab scale_x_discrete ggtitle element_text element_blank guide_legend rel unit
#' @import reshape2
#' @examples
#' wb <- watbal(t = c(10, 11.5, 14, 16.5, 20, 24.5, 27.5, 28, 24.5, 19.5, 14.5, 11),
#'              p = c(55, 73, 84, 58, 33, 23, 2, 2, 28, 66, 94, 71), lat = 35, CC = 400)
#'bb <- biobal(wb, 400)
#'bi <- bioint(bb)
#'plotBiobal(bi)
#'
#' @export
#'

plotBiobal <- function(intens){
  #
  for(i in 1:12){
    if(intens$FBIw[i] < intens$CBIw[i]) intens$FBIw_g[i] <- 0 else intens$FBIw_g[i] <- intens$FBIw[i] - intens$CBIw[i]
    intens$PBIw_g[i] <- intens$PBIw[i] - (intens$CBIw[i] + intens$FBIw_g[i])
    if(intens$CBIc[i] < intens$FBIc[i]) intens$FBIc_g[i] <- 0 else intens$FBIc_g[i] <- intens$FBIc[i] - intens$CBIc[i]
    intens$PBIc_g[i] <- intens$PBIc[i] - (intens$CBIc[i] + intens$FBIc_g[i])
  }
  # "CBIw" "FBIw_g" "PBIw_g" "DBIw" "CBIc" "FBIc_g" "DBIc" "PBIc_g"
  dbb <- data.frame(PBIc_g = intens$PBIc_g,
                    FBIc_g = intens$FBIc_g,
                    CBIc = intens$CBIc, 
                    DBIc = intens$DBIc,
                    PBIw_g = intens$PBIw_g, 
                    FBIw_g = intens$FBIw_g, 
                    CBIw = intens$CBIw,
                    DBIw = intens$DBIw)
  # dbb <- dbb[,rev(c(1,2,3,4,5,6,7,8))]

  dbb$mon <- as.character(formatC(1:12, width = 2, flag = '0'))
  dbb <- melt(dbb, id.var="mon")

  ggplot(dbb[order(dbb$variable),], aes(x = mon, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("legend", values = c("PBIc_g" = '#35B2FF',
                                           "FBIc_g" = "#2B7FE0",
                                           'CBIc' = "#224CC2",
                                           "DBIc" = '#1919A4',
                                           'PBIw_g' = '#FFBC35',
                                           'FBIw_g' = '#F1842A',
                                           'CBIw' = '#E44C1F',
                                           'DBIw' = '#D71414'),
                      labels = c('PBIc', 'FBIc', 'CBIc', 'DBIc',
                                 'PBIw', 'FBIw', 'CBIw', 'DBIw')) +
    theme_minimal() +
    theme(legend.position="bottom",
          legend.title = element_blank(),
          axis.title.y = element_text(size = rel(0.8)),
          plot.title = element_text(size = rel(0.9), face = 'bold'),
          legend.spacing.x = unit(0.15, 'cm')) +
    guides(fill = guide_legend(nrow = 1)) +
    xlab('') + ylab('') +
    scale_x_discrete(labels = month.abb) +
    ggtitle("Bioclimatic balance")

}
