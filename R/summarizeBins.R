#' Summarize Bins 
#'
#' Function counts occurences of numbers within intervals of interest.
#' It does so for a reference set and for a list of sampled sets, 
#' creating a list suitable for \code{\link{makeViolinPlot}}
#'
#' @param referenceSet Vector with numbers in the reference set
#' @param sampledList List of vectors with numbers in the background (e.g. sampled) sets
#' @param interval_starts Starts of the intervals [included]
#' @param interval_stops Stops of the intervals [included]
#' @return List with (a) counts in intervals for the reference and (b) for a background and containing names of bins as column names
#' @author Lukasz Jan Kielpinski
#' @export summarizeBins

summarizeBins <- function(referenceSet, sampledList, interval_starts, interval_stops){
  if(length(interval_starts) != length(interval_stops))stop("interval_starts and interval_stops must be of the same lengths")
  if(any(interval_starts > interval_stops))stop("interval_starts cannot be greater than interval_stops")
  
  #Summarize counts for the reference:
  countsInBinsGenome <- .summarizeBinsSingle(referenceSet, interval_starts = interval_starts, interval_stops = interval_stops)
  
  #Summarize counts for the background:
  countsInBinsSampling <- matrix(nrow = length(sampledList), ncol = length(countsInBinsGenome))
  for(j in 1:length(sampledList)){
    countsInBinsSampling[j,] <- .summarizeBinsSingle(sampledList[[j]], interval_starts = interval_starts, interval_stops = interval_stops)
  }
  
  #Make nice names of the bins:
  xaxis_names <- matrix(c(interval_starts, interval_stops), ncol = 2, byrow = F, nrow = length(interval_starts))
  colnames(countsInBinsSampling) <- apply(xaxis_names, 1, FUN = function(input){paste("[", input[1], ",", input[2], "]", sep = "")})
  #Return as a list:
  list(countsInBinsGenome, countsInBinsSampling)
}

.summarizeBinsSingle <- function(distances, interval_starts, interval_stops){
  # unlist(lapply(split(distances, cut(distances, c(1, interval_stops))), length)) #This line could be a standalone function but is actually slower
  binCount <- length(interval_starts)
  countsInBins <- rep(0, binCount)
  for(i in 1:binCount){
    countsInBins[i] <- sum(distances >= interval_starts[i] & distances <= interval_stops[i])
  }
  countsInBins
}
