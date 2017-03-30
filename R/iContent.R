#' Calculate information content from the changes in position-specific nucleotide abundunce
#'
#' Function calculates information content of different nucleotide positions in a construct.
#' Based on a formula from Gorodkin et al. 1997. 
#' Note that the information content is calculated on the whole set, 
#' hence it is a function of a probing extent.
#' 
#' @param DEresults Output of "GetKmerDE"
#' @param enzyme_source Enzyme of interest
#' @param constr Construct of interest
#' @return i_cont Vector with information content (Bits)
#' @author Lukasz Jan Kielpinski
#' @references Gorodkin, J., L. J. Heyer, S. Brunak, and G. D. Stormo. 1997. 
#' “Displaying the Information Contents of Structural RNA Alignments: 
#' The Structure Logos.” Computer Applications in the Biosciences: 
#' CABIOS 13 (6): 583–86.
#' @export iContent

iContent <- function(DEresults, enzyme_source, constr){
  i_local <- c()
  for(i in 1:length(DEresults[[enzyme_source]][[constr]][[1]])){
    CPM_local <- 2**DEresults[[enzyme_source]][[constr]][[1]][[i]][["logCPM"]]
    FC_local <- 2**DEresults[[enzyme_source]][[constr]][[1]][[i]][["logFC"]]
    before <- CPM_local*(1 + 0.5*(1-FC_local))
    after <- CPM_local*(1 - 0.5*(1-FC_local))
    pT <- after/sum(after)
    pC <- before/sum(before)
    i_local[i] <- sum(pT * log2(pT / pC))
  }
  i_local
}
