#' Plot correlation of logFC with target knockdown at different sites within antisense oligonucleotides
#'
#' Function plots correlation level with confidence intervals across set of equal-length oligos
#' It considers only bonds within the oligonucleotide.
#' 
#' @param cleavabilityMatrix Matrix with logFC scores calculated at different position of the oligonucleotides. 
#' Output of \code{\link{predictWithin}}
#' @param observed Experimental knockdown values
#' @param ... Further parameters to \code{\link{plot}}
#' @return Plot
#' @author Lukasz Jan Kielpinski
#' @seealso \code{\link{predictWithin}}
#' @export plotRunningCorrelation

plotRunningCorrelation <- function(cleavabilityMatrix, observed, conf.level = 0.99, ...){
  cor_test_mid <- c()
  cor_test_up <- c()
  cor_test_down <- c()
  for(i in 1:ncol(cleavabilityMatrix)){
    if(!any(is.na(cleavabilityMatrix[,i]))){
      cor_test_res <- cor.test(observed, cleavabilityMatrix[,i], conf.level = conf.level)
      cor_test_mid[i] <- cor_test_res$estimate
      cor_test_up[i] <- cor_test_res$conf.int[2]
      cor_test_down[i] <- cor_test_res$conf.int[1]
    }
  }
  #plot:
  local_x <- which(!is.na(cor_test_mid))
  plot(cor_test_mid[!is.na(cor_test_mid)] ~ local_x, type = "b", ylim = range(cor_test_up, cor_test_down, na.rm = T), xlab = "", ylab = "", ...)
  mtext(text = "Pearson's R", side = 2, line = 2)
  mtext(text = "Modeled RNase H cleavage site", side = 1, line = 2)
  segments(x0 = 1:length(cor_test_up), y0 = cor_test_up, y1 = cor_test_down)
  
  extremes <- c(cor_test_down, cor_test_up)
  
  text(x = min(local_x, na.rm = T), y = max(abs(extremes), na.rm = T)*0.9*sign(extremes[which.max(abs(extremes))]), pos = 4, labels = paste("N =", nrow(cleavabilityMatrix)))
  abline(h= 0, lty =  2)
  points(cor_test_mid[cor_test_up < 0 | cor_test_down > 0] ~ I((1:length(cor_test_mid))[cor_test_up < 0 | cor_test_down > 0]), pch = 20)
  points(cor_test_mid[which.max(abs(cor_test_mid))] ~ I((1:length(cor_test_mid))[which.max(abs(cor_test_mid))]), pch = 20, col = 2)

}
