#' Predict cleavability of bonds within a given sequence by the HIV-1 RNase H
#'
#' Function uses Heptamer Based Predictor ()
#' to assess the cleavability of a longer sequence by HIV-1 RNase H
#'
#' @param rnaSequence Sequence to be scored. Provide the RNA strand, but with T's (not U's)
#' @param hbp Heptamer Based Predictor. See details. 
#' @return Vector with cleavability scores (log2-fold-changes)
#' @author Lukasz Jan Kielpinski
#' @details Heptamer Based Predictor (hbp) is a list consisting of three vectors:
#' @details "kmers_seq" - vector of sequences of heptamers. Each kmer has to be represented.
#' @details "where_cut" - vector of cut sites (after which nt heptamer is cut). Calculated with \code{\link{findCleavageSite}}
#' @details "logFC_cut" - maximum observed extent of the cleavage of a given heptamer. 
#' @details This method will report the best cleavage directed to each bond.
#' NA indicates that no heptamers direct the cleavage towards that bond.
#' @export predictCleavages

predictCleavages <- function(rnaSequence, hbp){
  rnaSequence <- as.character(rnaSequence)
  if(length(rnaSequence) != 1)stop("Input exactly one sequence")
  if(any(!(unlist(strsplit(rnaSequence, "")) %in% c("A", "C", "G", "T"))))stop("Use only A, C, G, T")
  
  seq_len <- nchar(rnaSequence)
  kmer_len <- log(length(hbp$kmers_seq), 4)
  seq_in_kmers <- substr(x = rep(rnaSequence, (seq_len - kmer_len + 1)), 
                         start = 1:(seq_len-kmer_len + 1),
                         stop = kmer_len:seq_len)
  seq_kmers_mapping <- match(x = seq_in_kmers, table = hbp$kmers_seq)
  
  best_cleavage_at_sites <- rep(NA, seq_len)
  
  local_pos <- 1:length(seq_kmers_mapping) + hbp$where_cut[seq_kmers_mapping] - 1
  for(i in 1:length(seq_kmers_mapping)){
    best_cleavage_at_sites[local_pos[i]] <- min(hbp$logFC_cut[seq_kmers_mapping[i]], best_cleavage_at_sites[local_pos[i]], na.rm = T)
  }
  best_cleavage_at_sites
}
