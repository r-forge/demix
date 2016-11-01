% $Id$

\name{DeMix.Normalization}
\alias{DeMix.Normalization}
\title{
   Normalization for DeMix
}
\usage{
DeMix.Normalization(input,method, nn,nm)
}
\description{
A normalization function for quantile/total/median normalization.
The preparation work -- quantile normalization before DeMix (Bayesian deconvolution models for both RNA-seq read counts and microarray expressions).} 

\arguments{
\item{input}{
yourDataInput: Matrix of gene expressions before normalization and selection. The G*S matrix where G is the
number of genes and S is the number of samples.}

\item{method}{Method for normalization, "quantile" is a default. "median" and "total" also can be used as alternatives. }


\item{nn}{number of control/normal samples}

\item{nm}{number of mixed tumor samples }
}

\value{
 \item{newtfile}{Matrix of gene expressions. The G*S matrix where G is the number of genes and S is the number of samples.}}

\author{
Jaeil Ahn
}
\references{
Ahn J, Yuan Y, Parmigiani G, Suraokar MB, Diao L, Wistuba II, and Wang W. Bioinformatics 2013 doi: 10.1093/bioinformatics/btt301.

Wu H, Wang C and Wu Z (2013). “A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data.” Biostatistics. doi: 10.1093/biostatistics/kxs033.
}
\keyword{Normalization}
\seealso{
 http://odin.mdacc.tmc.edu/~wwang7/DeMix.html
}
\examples{
data(ntot_NT_Liver)
input <- ntot_NT_Liver
newtfile <- DeMix.Normalization(input, method="quantile" ,60,151)
}

