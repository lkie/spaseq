#' Predict log2 fold changes based on a kmer-based-PWM
#' 
#' Function uses position-specific-kmer scoring matrix to score input sequences.
#' Input sequences must be of the same length and must be fully covered by the scoring matrix.
#' For single nucleotide model it is a sum of log2-fold-changes of those 
#' nucleotides from the scoring matrix that match the observed nucleotides. 
#' For longer models, kmers overlap, hence only a fraction of each kmers associated 
#' log2-fold-change is added to the prediction according to implemented weighing mechanism.
#' 
#'
#' @param charVector Sequences to be scored. All need to be of the same length. 
#' The length has to match the length of scores [number of columns + kmer length - 1]
#' @param scores Matrix with scores; Subset of the output of the \code{\link{enrichmentMatrix}}
#' @return predicted_scores
#' @author Lukasz Jan Kielpinski
#' @seealso For scoring longer sequence use: \code{\link{scoreSequence}}
#' @export predictFC

predictFC <- function(charVector, scores){
  
  if(length(unique(nchar(charVector))) != 1)
    stop("Characters have to be of the same length")

  smallestUnit <- log(nrow(scores), 4) #Is the model 1nt, 2nt, ...
  kmerLen <- nchar(charVector[1])
  
  if(kmerLen != ncol(scores) + smallestUnit - 1)
    stop("Scoring matrix has to fit in only one way")
  
  predicted_scores <- rep(0, length(charVector))
  
  weightsMatrix <- matrix(ncol = kmerLen, nrow = kmerLen - smallestUnit + 1)
  weightsMatrix[,] <- 0
  for(w in 1:nrow(weightsMatrix)){
    weightsMatrix[w,w:(w + smallestUnit - 1)] <- 1
  }
  
  for(i in 1:ncol(scores)){
    if(smallestUnit > 1 & smallestUnit < kmerLen){
      normalizer <- sum(1/colSums(weightsMatrix[, i:(i + smallestUnit - 1)])*1/smallestUnit)}else{
        normalizer <- 1}
    
    ntAtPos <- substr(charVector, i, i + smallestUnit - 1)
    j <- 0
    for(nt1 in rownames(scores)){
      j <- j+1
      predicted_scores[ntAtPos == nt1] <- predicted_scores[ntAtPos == nt1] + scores[j,i] * normalizer
    }
  }
  
  predicted_scores
}

