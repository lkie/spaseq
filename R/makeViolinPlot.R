#' Make Violin Plot comparing counts in a reference set to counts in background sets
#'
#' Function uses ggplot2 to plot "violins" of background distribution (e.g. sampled content)
#' and prints rombi at the counts of a reference set. 
#' 
#' @param countsInBinsList Result of \code{\link{summarizeBins}}
#' @return Violin plot
#' @author Lukasz Jan Kielpinski
#' @import ggplot2
#' @export makeViolinPlot

makeViolinPlot <- function(countsInBinsList){
  #Reformat the result for easy plotting with ggplot:
  bin_names <- colnames(countsInBinsList[[2]])
  bins_easy_for_ggplot <- .binsForGGPlot(countsInBinsList[[2]], bin_names = bin_names)
  
  ggplot(bins_easy_for_ggplot, aes(factor(bin_names), count)) + 
    geom_violin(fill = "#0000FF88", colour = "#0000FF88") + 
    geom_point(fill = "red", size = 3, shape = 23, 
               data = data.frame(ref=countsInBinsList[[1]]), 
               aes(1:ncol(countsInBinsList[[2]]), ref, colour = "Reference\n(no shuffling)")) + 
    scale_color_manual(values=c("Reference\n(no shuffling)"="black")) +
    labs(x = "Nucleotides between neighboring well-cleaved sites", 
         y = "Occurences in shuffled genomes", colour = "") + 
    theme(panel.background	= element_rect(fill = "transparent"),
          legend.key	= element_rect(fill = "transparent"),
          panel.grid.major.y = element_line(colour = "black"),
          panel.grid.minor.y = element_line(colour = "black"),
          panel.grid.major.x = element_line(colour = "#00000000"),
          axis.ticks.x = element_line(colour = "#00000000"),
          axis.text	= element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15)) 
}

# Function for reformatting matrix into long data frame:
.binsForGGPlot <- function(countsInBinsSampling, bin_names){
  bins_easy_for_ggplot <- data.frame(bin_names = character(), count = integer())
  for(i in 1:ncol(countsInBinsSampling)){
    bins_easy_for_ggplot <- rbind(bins_easy_for_ggplot, data.frame(bin_names[i], countsInBinsSampling[,i]))
  }
  colnames(bins_easy_for_ggplot) <- c("bin_names", "count")
  bins_easy_for_ggplot
}
