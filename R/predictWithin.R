#' Predict cleavability of a set of antisene oligonucleotide target sites
#'
#' Function calculates the cleavability of oligo's target sites 
#' score using PWM of a set of oligonucleotides
#' 
#' @param oligoSequences Sequences of the antisense oligonucleotides
#' @param scores Matrix with scores; Subset of the output of the \code{\link{enrichmentMatrix}}
#' @param cleavageSiteInScores After which nucleotide of the scoring sequence the cleavage occurs. 
#' See details of \code{\link{scoreSequence}}
#' @return cleavabilityMatrix Matrix with cleavability scores at different position
#' @author Lukasz Jan Kielpinski
#' @import Biostrings
#' @export predictWithin

predictWithin <- function(oligoSequences, scores, cleavageSiteInScores){
  oligoLength <- nchar(oligoSequences[[1]])
  if(any(nchar(oligoSequences) != oligoLength))stop("All have to same the same length")
  
  targetSites <- as.character(reverseComplement(DNAStringSet(oligoSequences)))
  
  smallestUnit <- log(nrow(scores), 4)
  
  cleavabilityMatrix <- matrix(nrow = length(targetSites), ncol = oligoLength)
  for(oligo_no in 1:length(targetSites)){
    cleavabilityMatrix[oligo_no , ] <- rev(scoreSequence(targetSites[[oligo_no]], scores, cleavageSiteInScores)) # could be rewritten to directly use predictFC, probably would be faster
  }
  cleavabilityMatrix[,-1] #Do not include a cleavage outside of the oligo.
}




