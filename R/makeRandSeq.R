#' Make random sequence
#'
#' Function plots a heatmap of sequence preferences 
#' and optionally a barplot with information content
#'
#' @param length length of the sequence
#' @param nts nucleotides to use
#' @param prob at what prob each nucleotide should be sampled
#' @author Lukasz Jan Kielpinski
#' @return Character with a randomized sequence
#' @export makeRandSeq

makeRandSeq <- function(length, nts = c("A", "C", "G", "T"), prob = rep(0.25, 4)){
  paste(sample(x = nts, size = length, prob = prob, replace = T), collapse = "")
}
