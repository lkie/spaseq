#' Generate variables with all kmers of a given length
#'
#' Function generates variables with character vectors with all possible DNA kmers of given lengths (sorted)
#' Creates variables names dimers, trimers, tetramers, pentamers, hexamers, heptamers, octamers and nonamers
#' for 2 to 9 mers. 
#' 
#' @param whichToMake Variables for which lengths should be created. 
#' @return Creates a variable in a Global Environment (unless variable of the same name exists in a parent environment which is not global environment)
#' @author Lukasz Jan Kielpinski
#' @export makeNamedKmers

makeNamedKmers <- function(whichToMake){
  whichToMake <- as.integer(whichToMake)
  if(!all(whichToMake %in% 2:9)){message("WARNING: Exporting variables only from dimers to nonamers")}
  
  if(any(whichToMake %in% 2)){dimers <<- makeAllKmers(2)}
  if(any(whichToMake %in% 3)){trimers <<- makeAllKmers(3)}
  if(any(whichToMake %in% 4)){tetramers <<- makeAllKmers(4)}
  if(any(whichToMake %in% 5)){pentamers <<- makeAllKmers(5)}
  if(any(whichToMake %in% 6)){hexamers <<- makeAllKmers(6)}
  if(any(whichToMake %in% 7)){heptamers <<- makeAllKmers(7)}
  if(any(whichToMake %in% 8)){octamers <<- makeAllKmers(8)}
  if(any(whichToMake %in% 9)){nonamers <<- makeAllKmers(9)}
}
