#' Random sampling of the nucleotide sequence with statistical preservation of 
#' the local dinucleotide content
#'
#' Function generates a randomized version of the input DNA sequence, which is 
#' of the same length as the input sequence.
#' For position 'i' the occurences of dinucleotides are counted in the input sequence within a window 
#' from 'i - window_size/2' to 'i + window_size/2' and a new dinucleotide is selected according to:
#' (a) Narrow the set of dinucleotides to those which start from the second nucleotide of a previous dinucleotide.
#' (b) Randomly selected a new dinucleotide from the set defined in (a) with probabilities proportional to number 
#' of occurences of dinucleotides in the window. If occurences == 0, then select one by random. 
#' For the first dinucleotide: skip point (a).
#' 
#'
#' @param char Original sequence to be sampled
#' @param window_size Width of the window from which dinucleotide content should be calculated. 
#' Window is centered around sampled nucleotide.
#' @return Sampled sequence
#' @author Lukasz Jan Kielpinski
#' @seealso \code{\link{globalDinucleotideSampling}}
#' @export localDinucleotideSampling

#Considering dinucleotide frequencies:
localDinucleotideSampling <- function(char, window_size){
  
  if(any(!(unlist(strsplit(char, "")) %in% c("A", "C", "G", "T"))))stop("Use only A, C, G, T")
  
  dins <- makeAllKmers(2) #All dinucleotides
  
  #Calculate probabilities of observing dinucleotides within a window, store them in a matrix
  #(add smalllest possible pseudocount - prevents errors if all dinucleotides starting with a current nt have prob = 0)
  dinucleotideLength <- nchar(char) - 1
  char_split <- substr(rep(char, dinucleotideLength), 1:dinucleotideLength, 2:nchar(char))
  din_prob_mat <- matrix(nrow = 16, ncol = dinucleotideLength)
  for(i in 1:16){
    din <- dins[i]
    din_prob_mat[i,] <- .moving_average(char_split == din, window_size) + .Machine$double.eps
  }
  #end of calculating probabilities
  
  output <- c()
  allowed <- rep(TRUE, 16) #At the beginning all dinucleotides are allowed (with differing probabilities)
  for(i in 1:dinucleotideLength){
    output[i] <- sample(x = dins[allowed], prob = din_prob_mat[allowed,i], size = 1)
    allowed <- substr(dins, 1, 1) == substr(output[i], 2, 2) #in the following cycle the allowed dinucleotides are constrained by the last nt of the preceding dinucleotide
  }
  no_last <- paste(substr(output, 1, 1), collapse = "") #Doesn't have the last nucleotide, which has to be added here:
  paste(no_last, substr(output[length(output)], 2, 2), sep = "")
}

#Aux. functions copied from RNAprobR:
.construct_smoothing_matrix <- function(input_vector, window_size){
  
  vector_length <- length(input_vector)
  #Create empty matrix.
  my_mat <- matrix(nrow=window_size, ncol=(vector_length+window_size-1))
  #Fill the matrix, offsetting by 1 in each cycle.
  for(i in 1:(window_size))
  {
    my_mat[i,i:(vector_length+i-1)] <- input_vector
  }
  
  my_mat
}

.moving_average <- function(input_vector, window_size)
{
  window_side <- window_size/2-0.5
  
  colMeans(.construct_smoothing_matrix(input_vector, window_size),
           na.rm=TRUE)[(window_side+1):(length(input_vector)+window_side)]
}

