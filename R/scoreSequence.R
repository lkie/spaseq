#' Score sequence with PWM based scores
#'
#' Function to score a single sequence for cleavability
#' given kmer-based PWM of RNase H sequence preferences
#'
#' @param rnaSequence Sequence to be scored. Provide the RNA strand, but with T (not U)
#' @param scores Matrix with scores; Subset of the output of the \code{\link{enrichmentMatrix}}
#' @param cleavageSiteInScores After which nucleotide of the scoring sequence the cleavage occurs. See details
#' @details Assume that the scoring matrix (scores) has length i, and the atomic unit (kmer) length is j. 
#' Then the scoring sequence has the length i + j - 1. 
#' cleavageSiteInScores describes after which nucleotide within the scoring matrix the cleavage is expected to occur.
#' @return Vector with logFC scores for each bond of the input sequence (including after last nt)
#' @author Lukasz Jan Kielpinski
#' @export scoreSequence

scoreSequence <- function(rnaSequence, scores, cleavageSiteInScores){
  if(length(rnaSequence) != 1)stop("Input exactly one sequence")
  rnaSequence <- as.character(rnaSequence)
  if(any(!(unlist(strsplit(rnaSequence, "")) %in% c("A", "C", "G", "T"))))stop("Use only A, C, G, T")
  
  rnaSeqLen <- nchar(rnaSequence)
  smallestUnit <- log(nrow(scores), 4)
  pwm_length <- ncol(scores)
  kmerLen <- pwm_length + smallestUnit - 1
  localKmers <- substr(rep(rnaSequence, rnaSeqLen - kmerLen + 1), start = 1:(rnaSeqLen - kmerLen + 1), stop = kmerLen:rnaSeqLen)
  c(rep(NA, cleavageSiteInScores - 1), predictFC(localKmers, scores), rep(NA, kmerLen - cleavageSiteInScores))
}
