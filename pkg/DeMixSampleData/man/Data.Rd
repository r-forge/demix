% $Id$

\name{DeMix.Data}
\alias{DeMix.Data}
\docType{data}
\title{
  DeMix.Data
}
\usage{
data(DataFile)
}
\description{
 There are four data sets in our Data file (1.ntot_NT_Liver, 2.Prop, 3.simulS1, 4.simulS2). 1.ntot_NT_Liver is the input file for 'DeMix.Normalization' function. 3.simulS1 is the input file for 'DeMix.S1' function. 2&4. are the input files for 'DeMix.S2' function. 
}

\arguments{
  \item{ntot_NT_Liver}{
A matrix of gene expressions after the RMA normalization, which includes 211 samples (the first 50 samples are normal controls, and the remaining 161 are tumor mixed samples). It is available from ref [1].}

  \item{Prop}{
A vector of given true tumor proportion for samples in simulS1 and simulS2. It is input to DeMix.S2 for deconvolved expressions of the whole set of genes as the fixed truth.
}

 \item{simulS1}{
A 100 * 20 matrix of simulated RNA-seq counts data. It contains 50 genes subset from simulS2 for estimation of tumor proportions. These first 10 columns correspond to 10 normal control samples and second 10 columns correspond to these 10 mixed tumor samples to be deconvolved.
}

 \item{simulS2}{
A 200 * 20 matrix of simulated RNA-seq counts data. It contains 200 genes for expression deconvolution. These first 10 columns correspond to 10 normal control samples and second 10 columns correspond to these 10 mixed tumor samples to be deconvolved.

}
}


\author{
Jaeil Ahn
}
\references{[1] Shen-Orr, Shai S and Tibshirani, Robert and Khatri, Purvesh and Bodian, Dale L and Staedtler, Frank and Perry, Nicholas M and Hastie, Trevor and Sarwal, Minnie M and Davis, Mark M and Butte, Atul J. Cell type-specific gene expression differences in complex tissues. Nature Method. 7(4):287-289, 2010.
 }
\seealso{
 http://odin.mdacc.tmc.edu/~wwang7/DeMix.html
}
\examples{
 data(ntot_NT_Liver)
}
\concept{Bayesian Model}
\concept{Deconvolution}
\concept{Mixed Cancer Transcriptomes}

