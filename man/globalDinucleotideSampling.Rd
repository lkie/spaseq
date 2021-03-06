% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/globalDinucleotideSampling.R
\name{globalDinucleotideSampling}
\alias{globalDinucleotideSampling}
\title{Random sampling of the nucleotide sequence with statistical preservation of 
the global dinucleotide content}
\usage{
globalDinucleotideSampling(char)
}
\arguments{
\item{char}{Original sequence to be sampled}
}
\value{
Sampled sequence
}
\description{
Function generates a randomized version of the input DNA sequence, which is 
of the same length as the input sequence.
A new dinucleotide is selected according to:
(a) Narrow the set of dinucleotides to those which start from the second nucleotide of a previous dinucleotide.
(b) Randomly selected a new dinucleotide from the set defined in (a) with probabilities proportional to number 
of occurences of dinucleotides in the whole provided sequence If occurences == 0, then select one by random. 
For the first dinucleotide: skip point (a).
}
\author{
Lukasz Jan Kielpinski
}
\seealso{
\code{\link{localDinucleotideSampling}}
}

