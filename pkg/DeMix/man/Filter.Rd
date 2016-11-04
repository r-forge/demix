% $Id$

\name{DeMix.Filter}
\alias{DeMix.Filter}
\title{
  Filters uninformative gene expressions during first stage of deconvolution.
}
\description{
  To avoid unnecessary computational burdens and noises, pre-screening is an
  essential step. This screening process consists of three steps:\cr
  \enumerate{
    \item step 1
    \item step 2
    \item step 3
  }
}
\usage{
DeMix.Filter(input.mat, design, zerofilter=TRUE, conc=0.8, fc=1.2)
}
\arguments{
  \item{input.mat}{numeric matrix of gene expressions after normalization.
    The \code{G*S} matrix where \code{G} is the number of genes and
    \code{S} is the number of samples.}
  \item{design}{integer vector of \code{S} length specifying tissue components,
    where \code{0} represents healthy and \code{1} represents tumor samples}
  \item{zerofilter}{logical scalar. If \code{TRUE}, uses the screening based on    zero counts.}
  \item{conc}{numeric scalar cut-off value used for checking the linearity
    assumption. More than 0.8 or 80\% is recommended.}
  \item{fc}{numeric scalar fold-change cut-off value used for identifying
    informative genes in proportion estimations. 1.2 or 1/1.2 is default.
    This value needs to be set for reducing data-size. Approximately
    2,000-3,000 genes/probes will provide robust estimates.}
}
\value{
  Returns numeric \code{G'*S} matrix of gene expressions after the prescreening,
  where \code{G'} is the number of remaining genes and \code{S} is the number
  of samples.
}
\examples{
data(ntot_NT_Liver)

# First 50 samples normal, remaining 151 tumor for this dataset
input.mat <- as.matrix(ntot_NT_Liver.df)
design <- c(rep(0, 50),
            rep(1, 151))
filtered.mat <- DeMix:::DeMix.Filter(input.mat, design)
}

