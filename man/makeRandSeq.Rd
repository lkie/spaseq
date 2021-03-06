% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeRandSeq.R
\name{makeRandSeq}
\alias{makeRandSeq}
\title{Make random sequence}
\usage{
makeRandSeq(length, nts = c("A", "C", "G", "T"), prob = rep(0.25, 4))
}
\arguments{
\item{length}{length of the sequence}

\item{nts}{nucleotides to use}

\item{prob}{at what prob each nucleotide should be sampled}
}
\value{
Character with a randomized sequence
}
\description{
Function plots a heatmap of sequence preferences 
and optionally a barplot with information content
}
\author{
Lukasz Jan Kielpinski
}

