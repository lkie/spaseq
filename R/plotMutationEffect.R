#' Plot effects of substitutions
#'
#' Function plots what happens to log2-fold changes when one nucleotide 
#' of a reference sequence is substituted
#' 
#'#' @param stringSet Sequences of all kmers of a given length that have FoldChanges data
#' @param oneString sequence of a reference kmer
#' @param FoldChanges logFC paired with stringSet
#' @return Plot
#' @author Lukasz Jan Kielpinski
#' @export plotMutationEffect

plotMutationEffect <- function(stringSet, oneString, FoldChanges){
  dist2pref <- .hamming_dist(stringSet, oneString)
  distOfOne <- stringSet[dist2pref == 1]
  
  whereD <- .whereDiff(distOfOne, oneString)
  whatD <- .whatDiff(distOfOne, oneString)
  
  plot(0, col =0, xlim = range(whereD), ylim = range(c(FoldChanges[which(dist2pref <= 1)])) + c(-0.2, 0.2), 
       xaxt = "n", ylab = expression(log[2]*(FC)), xlab = "", las = 1, bty = "n")
  for(i in 1:length(whereD)){
    text(x = whereD[i], y = FoldChanges[which(dist2pref == 1)[i]], labels = .T_to_U(whatD[i]), col = 2)
  }
  axis(1, .T_to_U(unlist(strsplit(oneString, ""))), at = 1:nchar(oneString), tick = FALSE, line = -1)
  abline(h = FoldChanges[which(dist2pref == 0)], lty = 2)
}

#Define auxillary functions:
.hamming_dist <- function(stringSet, oneString){
  oneStringSplit <- unlist(strsplit(oneString, ""))
  unlist(
    lapply(strsplit(stringSet, ""), FUN = function(one_string){
      sum(one_string != oneStringSplit)})
  )
}

.whereDiff <- function(stringSet, oneString){
  oneStringSplit <- unlist(strsplit(oneString, ""))
  unlist(
    lapply(strsplit(stringSet, ""), FUN = function(one_string){
      which(one_string != oneStringSplit)})
  )
}

.whatDiff <- function(stringSet, oneString){
  oneStringSplit <- unlist(strsplit(oneString, ""))
  unlist(
    lapply(strsplit(stringSet, ""), FUN = function(one_string){
      one_string[which(one_string != oneStringSplit)]})
  )
}

.T_to_U <- function(char_vector){
  char_vector[char_vector == "T"] <- "U"
  char_vector
}
