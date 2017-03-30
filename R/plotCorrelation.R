#' Simple convenience function for plotting observed vs predicted log-fold-changes
#'
#' @param observed observed logFC
#' @param predicted predicted logFC
#' @author Lukasz Jan Kielpinski
#' @export plotCorrelation

plotCorrelation <- function(observed, predicted, col = "#00000010", ...){
  #Make the plot square:
  opar <- par(pty = "s")
  on.exit(par(opar))
  ###
  plot(x = predicted, y = observed, pch = 20, 
       xlab = expression(Predicted~log[2]~FC), 
       ylab = expression(Observed~log[2]~FC), 
       las = 1, col = col, ...)
  myCor <- round(cor(observed, predicted)**2,2)
  text(x = min(predicted), y=max(observed)*0.9, bquote(R^2 ~ "="~.(myCor)), pos = 4)
}
