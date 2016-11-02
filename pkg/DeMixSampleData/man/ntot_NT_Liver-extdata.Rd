% $Id$

\name{ntot_NT_Liver-extdata}
\alias{ntot_NT_Liver-extdata}
\docType{data}
\title{DeMix Data}
\description{
  This dataset consists of a matrix of gene expressions after the RMA
  normalization, which includes 211 samples (the first 50 samples are
  normal controls, and the remaining 161 are tumor mixed samples).
  It is available from ref [1].}
}
\references{
  [1]
  Shen-Orr SS, Tibshirani R, Khatri P, Bodian DL, Staedtler F, Perry NM,
  Hastie T, Sarwal MM, Davis MM, & Butte AJ.\cr
  \emph{Cell type-specific gene expression differences in complex tissues.}\cr
  Nature Method (2010) 7(4): 287-289.\cr
  \url{https://doi.org/10.1038/nmeth.1439}
}
\source{
  Where'd it come from?
}
\author{
  P. Roebuck \email{proebuck@mdanderson.org}
}
\examples{
## Path to raw data directory
extdata.dir <- system.file("extdata", package="DeMixSampleData")
extdata.dir

## Directory Listing
files <- list.files(extdata.dir, full.names=TRUE)
files
}
\keyword{datasets}

