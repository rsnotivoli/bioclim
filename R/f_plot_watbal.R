#' Function to plot water balance
#'
#' @description Function to plot water balance.
#' @param bh Water balance in data.frame format from watbal() function.
#' @return Plot of water balance
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_manual theme_minimal theme guides xlab ylab scale_x_discrete ggtitle  element_text element_blank guide_legend  rel unit
#' @import reshape2
#' @importFrom grDevices rgb
#' @examples
#' wb <- watbal(t = c(10, 11.5, 14, 16.5, 20, 24.5, 27.5, 28, 24.5, 19.5, 14.5, 11),
#' p = c(55, 73, 84, 58, 33, 23, 2, 2, 28, 66, 94, 71), lat = 35, CC = 400)
#' plotWatbal(wb)
#' @export
#'

plotWatbal <- function(bh){
  #
  dbh <- matrix(NA, ncol = 5, nrow = 12)
  colnames(dbh) <-  c('white', 'water_exc', 'soil_moist', 'wat_def', 'soil_rech')
  rownames(dbh) <- month.name

  for(i in 1:12){
    dbh[i, 1] <- min(bh$PET[i], bh$P[i])
    if(bh$RET[i] == bh$PET[i]) {if(round(bh$ME[i],2) > 0) dbh[i, 2] <- bh$P[i] - dbh[i, 1] else dbh[i, 2] <- 0} else dbh[i, 2] <- 0
    if(bh$PET[i] > bh$P[i]) dbh[i, 3] <- bh$RET[i] - sum(dbh[i, 1:2]) else dbh[i, 3] <- 0
    if(bh$PET[i] > bh$RET[i]) dbh[i, 4] <- bh$PET[i] - sum(dbh[i, 1:3]) else dbh[i, 4] <- 0
    if(bh$PET[i] > bh$P[i]) dbh[i, 5] <- 0 else{
      if(bh$RET[i] == bh$PET[i]) {
        if(round(bh$ME[i],2) == 0) dbh[i, 5] <- bh$P[i] - sum(dbh[i, 1:4]) else dbh[i, 5] <- 0
      } else dbh[i, 5] <- 0
    }
  }

  dbh <- as.data.frame(dbh)
  dbh <- dbh[, c(4, 3, 5, 2, 1)]
  #names(dbh) <- paste(formatC(1:ncol(dbh), width = 2, flag = '0'), names(dbh))
  names(dbh) <- c('Water deficit', 'Soil water use', 'Soil water recharge', 'Water exceedance', ' ')
  dbh$mon <- as.character(formatC(1:12, width = 2, flag = '0'))
  dbh <- melt(dbh, id.var="mon")

  ggplot(dbh[order(dbh$variable),], aes(x = mon, y = value, fill = variable)) +
    geom_bar(stat = "identity") +
    scale_fill_manual("legend", values = c(" " = rgb(1, 1, 1, 0),
                                           "Water exceedance" = "darkblue",
                                           "Soil water use" = "orange",
                                           'Water deficit' = 'red',
                                           'Soil water recharge' = 'lightblue')) +
    theme_minimal() +
    theme(legend.position="bottom",
          legend.title = element_blank(),
          axis.title.y = element_text(size = rel(0.8)),
          plot.title = element_text(size = rel(0.9), face = 'bold'),
          legend.spacing.x = unit(0.15, 'cm')) +
    xlab('') + ylab('mm') +
    scale_x_discrete(labels = month.abb) +
    ggtitle("Water balance")

}
