#' Generate all kmers of a given length
#'
#' Function generates a character vector with all possible DNA kmers of a given length (sorted)
#' 
#' @param smallestUnit Lenght of kmers to be generated
#' @param nts Nucleotides to use (default = c("A", "C", "G", "T"))
#' @return Sorted vector of kmers (character)
#' @author Lukasz Jan Kielpinski
#' @export makeAllKmers

makeAllKmers <- function(smallestUnit, nts = c("A", "C", "G", "T")){
  seq_list <- list()
  for(i in 1:smallestUnit){seq_list[[i]] <- nts}
  sort(apply(expand.grid(seq_list, stringsAsFactors = F), 1, paste, collapse=""))
}
