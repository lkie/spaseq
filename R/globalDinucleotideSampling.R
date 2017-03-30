#' Random sampling of the nucleotide sequence with statistical preservation of 
#' the global dinucleotide content
#'
#' Function generates a randomized version of the input DNA sequence, which is 
#' of the same length as the input sequence.
#' A new dinucleotide is selected according to:
#' (a) Narrow the set of dinucleotides to those which start from the second nucleotide of a previous dinucleotide.
#' (b) Randomly selected a new dinucleotide from the set defined in (a) with probabilities proportional to number 
#' of occurences of dinucleotides in the whole provided sequence If occurences == 0, then select one by random. 
#' For the first dinucleotide: skip point (a).
#' 
#'
#' @param char Original sequence to be sampled
#' @return Sampled sequence
#' @author Lukasz Jan Kielpinski
#' @seealso \code{\link{localDinucleotideSampling}}
#' @export globalDinucleotideSampling


globalDinucleotideSampling <- function(char){
  if (any(!(unlist(strsplit(char, "")) %in% c("A", "C", "G", 
                                              "T")))) 
    stop("Use only A, C, G, T")
  dins <- makeAllKmers(2)
  dinucleotideLength <- nchar(char) - 1
  char_split <- substr(rep(char, dinucleotideLength), 1:dinucleotideLength, 
                       2:nchar(char))
  din_summary <- table(char_split)
  din_prob_mat <- as.integer(din_summary)[match(x = names(din_summary), table = dins)]
  output <- c()
  allowed <- rep(TRUE, 16)
  for (i in 1:dinucleotideLength) {
    output[i] <- sample(x = dins[allowed], prob = din_prob_mat[allowed], size = 1)
    allowed <- substr(dins, 1, 1) == substr(output[i], 2, 
                                            2)
  }
  no_last <- paste(substr(output, 1, 1), collapse = "")
  paste(no_last, substr(output[length(output)], 2, 2), sep = "")
}
