#' Make Heptamer Based Predictor for HIV-1 RNase H
#' 
#' Function makes the hbp (Heptamer Based Predictor) from differential abundances
#' Only for HIV-1 RNase H
#'
#' @param DEresults Output of "GetKmerDE"
#' @param enzyme_source Enzyme of interest. Function intended to be used only with HIV-1 RNase H.
#' @param constr Construct of interest. Function intended to be used only with R7 construct.
#' @return hbp Heptamer Based Predictor (List of kmers_seq, where_cut and logFC_cut)
#' @author Lukasz Jan Kielpinski
#' @export make_hbp

make_hbp <- function(DEresults, enzyme_source = "HIV", constr){
  # Make matrix with fold change scores for each heptamer at each position:
  logFC_pos_mat <- enrichmentMatrix(DEresults, enzyme_source, constr, smallestUnit = 7) 
  # For each heptamer sequence assign where, within the heptamer, HIV-1 RNase H mediated cleavage occurs:
  where_cut <- findCleavageSite(fc_results = logFC_pos_mat, kmerLen = 7)
  # Extract the best cleavage score for each heptamer:
  logFC_max <- apply(logFC_pos_mat, 1, min)
  # Combine into a list:
  list(kmers_seq = makeAllKmers(7), where_cut = where_cut, logFC_cut = logFC_max)
}