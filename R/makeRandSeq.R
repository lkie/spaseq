#' Make random sequence
#'
#'
#' @param length length of the sequence
#' @param nts nucleotides to use
#' @param prob at what prob each nucleotide should be sampled
#' @author Lukasz Jan Kielpinski
#' @return Character with a randomized sequence
#' @export makeRandSeq

makeRandSeq <- function(length, nts = c("A", "C", "G", "T"), prob = rep(0.25, 4), n = 1){
  generated_sequences <- rep(NA, i)
  for(i in 1:n){
    generated_sequences[i] <- paste(sample(x = nts, size = length, prob = prob, replace = T), collapse = "")
    }
}
