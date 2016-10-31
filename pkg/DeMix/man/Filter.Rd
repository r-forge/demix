\name{DeMix.Filter}
\alias{DeMix.Filter}
\docType{package}
\title{
    Filtering un-informative gene expressions for the first stage of deconvolution 
}

\description{
To avoid unnecessary computational burdens and noises, 	pre-screening is an essential step. 
This screening process consists of three steps.
}

\usage{
DeMix.Filter(data=newt, ntgroup=c(0,0,  ,1), zero_filter=TRUE, conc=0.8, fc=1.2)
}


\arguments{
\item{data}{
yourDataInput: A Matrix of gene expressions after normalization. The G*S matrix where G is the
number of genes and S is the number of samples.}

\item{ntgroup}{A vector of tissue component 0 for healthy and 1 for tumor samples}

\item{zero_filter}{An indicator of whether the screening based on zero counts will be used or not. (YES=use).
Default is YES. }


\item{conc}{A cut-off value used for checking the linearity assumption. More than 0.8 or 80% is recommended.
 }

\item{fc}{A fold-change Cut-off value used for identifying informative genes in proportion estimations. 1.2 or 1/1.2 is default. This value needs to be set for reducing data-size. Approximatelby 2,000~3,000 genes/probes will provide robust estimates.
}



}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{RNA}{A Matrix of gene expressions after the prescreening. The G'*S matrix where G' is the number of remained genes and S is the number of samples.}

\item{nused}{The number of remained genes for the first Stage of DeMix.}

\item{useornot}{The vector of whether a gene is selected or not.}

}



\examples{
data(ntot_NT_Liver)
input <- ntot_NT_Liver
newtfile <- DeMix::DeMix.Filter(input,c(rep(0, 50), rep(1, 151)), zero_filter=TRUE, conc=0.8, fc=1.2)
}




