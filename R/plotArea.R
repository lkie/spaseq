#' Plot HIV-1 RNase H cleavage scores together with a sequence
#'
#' Function plots cleavability score (based on a log2-fold-change of a heptamer in the probing experiment) 
#' of bonds within a provided sequence. It takes as input predicted cleavage scores or 
#' alternatively hbp list with which it runs the prediction internally.
#'                                                                                       
#'
#' @param rnaSequence Sequence of the RNA strand (T will be converted to U for display)
#' @param CleavageScores Cleavage scores to plot
#' @param hbp If CleavageScores are not provided, hbp and rnaSequence will be used to call \code{\link{predictCleavages}} and calcaulte the scores. Look there for more info. 
#' @param startPos From which position to plot
#' @param endPos To which position to plot
#' @param xlab X-axis label
#' @param show_sequence Should RNA sequence be displayed on the plot?
#' @param ... Further parameters to \code{\link{plot}}
#' @return Plot
#' @author Lukasz Jan Kielpinski
#' @seealso \code{\link{predictCleavages}}
#' @export plotArea

#Function for plotting RNase H cleavage scores in the sequence (HIV-style)
#CleavageScores OR hbp has to be provided:
plotArea <- function(rnaSequence, CleavageScores, hbp, startPos, endPos, xlab = "Bond number", show_sequence = TRUE, ...){
  if(missing(CleavageScores)){
    CleavageScores <- predictCleavages(rnaSequence, hbp)
    message("Calculating scores based on sequence and hbp provided")
  }
  if(missing(startPos)){startPos <- 1}
  if(missing(endPos)){endPos <- nchar(rnaSequence)}

  plot(CleavageScores[startPos:endPos] ~ I(startPos:endPos), 
       type = "h", ylab = expression(log[2]*(FC)), xlab = xlab, 
       las = 1, ylim = range(c(CleavageScores[startPos:endPos], 0), na.rm = T), ...)
  
  if(show_sequence){
      toRNA <- unlist(strsplit(as.character(rnaSequence[startPos:endPos]), ""))
  toRNA[toRNA == "T"] <- "U"
  text(x = seq(startPos - 0.5, endPos - 0.5,1), y = 0, 
       labels = toRNA, pos = 1)
  }
}
