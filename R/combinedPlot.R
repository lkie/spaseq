#' Plotting sequence preferences
#'
#' Function plots a heatmap of sequence preferences 
#' and optionally a barplot with information content
#'
#' @param fc_results Output of "enrichmentMatrix"
#' @param i_cont Output of "iContent"
#' @param main Title of the plot
#' @author Lukasz Jan Kielpinski
#' @seealso \code{\link{enrichmentMatrix}}
#' @seealso \code{\link{iContent}}
#' @import lattice
#' @export combinedPlot

combinedPlot <- function(fc_results, i_cont, main = ""){
  rgb.palette <- colorRampPalette(c("blue", "white", "red"), space = "Lab")
  r1 <- max(abs(fc_results))
  # rownames(fc_results)[4] <- "U"
  w1 = levelplot(t(fc_results)[,4:1], col.regions = rev(rgb.palette(100)), aspect = 4/ncol(fc_results), 
                 at = seq(-r1, r1, length.out = 100),
                 colorkey=c(), scales = list(tck = c(1,0), cex = 2, x = list(at = c())),
                 xlab = "", ylab = list("nt", cex = 2))
  
  if(!(missing(i_cont))){
    w2 = barchart((i_cont) ~ rep(1:ncol(fc_results), length(i_cont)/ncol(fc_results)), 
                groups = sort(rep(1:I(length(i_cont)/ncol(fc_results)), ncol(fc_results))), 
                scales = list(tck = c(1,0), x=list(draw=FALSE), cex = 2, 
                              tick.number = max(2, floor(max(i_cont, na.rm = T))/2) ), 
                origin = 0, xlab = "", ylab = list("Bit", cex = 2), horizontal = F, col = "black")
  }
  if(missing(i_cont)){
    print(w1)
  }else{
    print(w1, position=c(0, 0.2, 1, 1), more=TRUE)
    print(w2, position = c(0, 0, 1, .4))
    }
  grid::grid.text(label = main, x = 0.5, y = 1, just = "top")
}

