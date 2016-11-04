% $Id$

\name{DeMix.Normalization}
\alias{DeMix.Normalization}
\title{
  Normalization for DeMix
}
\usage{
DeMix.Normalization(input,
                    design,
                    method=c("total", "quantile", "median"))
}
\description{
  A normalization function for quantile/total/median normalization.
  The preparation work -- quantile normalization before DeMix (Bayesian
  deconvolution models for both RNA-seq read counts and microarray expressions).} 

\arguments{
  \item{input}{numeric \code{G*S} matrix of gene expressions before
    normalization and selection, where \code{G} is the number of genes and
    \code{S} is the number of samples}
  \item{design}{numeric vector of length \code{S} specifying whether the
    sample type is normal or mixed tumor}
  \item{method}{character string specifying the normalization method}
}
\value{
  Returns same-sized numeric matrix of normalized gene expressions.
}
\author{
  Jaeil Ahn,
  P. Roebuck \email{proebuck@mdanderson.org}
}
\references{
  Ahn J, Yuan Y, Parmigiani G, Suraokar MB, Diao L, Wistuba II, & Wang W.\cr
  \emph{DeMix: deconvolution for mixed cancer transcriptomes using raw measured data}\cr
  Bioinformatics (2013) 29 (15): 1865-1871\cr
  \url{https://doi.org/10.1093/bioinformatics/btt301}\cr
\cr
  Wu H, Wang C, & Wu Z.\cr
  \emph{A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data}\cr
  Biostatistics (2013) 14 (2): 232-243\cr
  \url{https://doi.org/10.1093/biostatistics/kxs033}
}
\examples{
data(ntot_NT_Liver)
input <- ntot_NT_Liver.df
design <- c(rep(0, 60), rep(1, 151))   # 60 normal, then 151 tumor samples
input.norm <- DeMix:::DeMix.Normalization(input, design, method="quantile")
}
\keyword{Normalization}

