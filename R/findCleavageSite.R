#' Assigns preferred cleavage site of HIV-1 RNase H within a heptamer
#'
#' Function assigns where within a given heptamer, the HIV-1 RNase H is expected to cleave.
#' Implementation works only with R7 construct, 
#' based on the observation described in the H-SPA publication.
#'
#' @param fc_results Output of the "enrichmentMatrix" function with smallestUnit set to 7
#' @param kmerLen Length of the cleaved kmer (Currently the function works only for kmerLen = 7)
#' @return cut_site Vector with numbers indicating after which nucleotide of each provided heptamer HIV-1 RNase H preferentially cleaves 
#' @author Lukasz Jan Kielpinski
#' @seealso \code{\link{enrichmentMatrix}}
#' @export findCleavageSite

findCleavageSite <- function(fc_results, kmerLen = 7){
  if(kmerLen != 7)stop("Current implementation works only for heptamers")
  
  bests <- apply(fc_results, MARGIN = 1, FUN = function(scores_line){paste(sort(order(scores_line)[1:2]), collapse = "")})
  cut_site <- match(bests, rev(paste(1:6, 2:7, sep = ""))) + 1
  notNeighbors <- is.na(cut_site)
  cut_site[notNeighbors] <- apply(fc_results[notNeighbors, ], MARGIN = 1, FUN = function(scores_line){8 - which.min(scores_line)})
  cut_site
}
