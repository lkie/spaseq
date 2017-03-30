#' Plot sequence logo from vector of characters
#'
#' Function to plots sequence logo using seqLogo package, 
#' but using vector of characters (of the same lenth) as input
#' 
#' @param charVector Vector of equally-long characters; only ACGT allowed
#' @param pwm_out Should PWM be returned? (default = FALSE)
#' @param plot Should sequence logo be plotted? (default = TRUE)
#' @param na.rm Remove NA values prior to plotting (default = TRUE)
#' @return pwmMatrix
#' @author Lukasz Jan Kielpinski
#' @import seqLogo
#' @export seqLogoFromChars

seqLogoFromChars <- function(charVector, pwm_out = FALSE, plot = T, na.rm = T, ...){
  #Check the input:
  if(length(unique(nchar(charVector))) != 1)
    stop("Characters have to be of the same length")
  if(!all(unique(unlist(strsplit(charVector, ""))) %in% c("A", "C", "G", "T")))stop("Only A, C, G, T accepted")
  ###
  if(na.rm){charVector <- charVector[!is.na(charVector)]}
  tempMatrix <- matrix(unlist(strsplit(charVector, "")), ncol = nchar(charVector[1]), byrow = T)
  pwmMatrix <- matrix(nrow = 4, ncol = nchar(charVector[1]))
  for(i in 1:nchar(charVector[1])){
    j=1
    for(nt in c("A","C","G","T")){
      pwmMatrix[j,i] <- sum(tempMatrix[,i] == nt)
      j=j+1
    }
  }
  if(pwm_out){
    return(pwmMatrix)
  }
  if(plot)seqLogo(makePWM(pwmMatrix/colSums(pwmMatrix)), ...)
}
