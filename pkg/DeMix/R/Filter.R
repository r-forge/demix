###
### $Id$
###

##-----------------------------------------------------------------------------
DeMix.Filter <- function(input.mat,
                         design,
                         zerofilter=TRUE,
                         conc=0.8,
                         fc=1.2)
{
    ## Check arguments
    stopifnot(is.matrix(input.mat) &&
              is.numeric(input.mat[, 1]) &&
              !anyNA(input.mat))
    stopifnot(is.numeric(design) && !anyNA(design) && length(design) >= 2)
    stopifnot(is.scalar.logical(zerofilter))
    stopifnot(is.scalar.numeric(conc))
    stopifnot(is.scalar.numeric(fc))

    ## 0 denotes normal cell while 1 denote mixed tumor cell
    ## Filtering step 1.
    ## Input data are previously normalized at the data processing stage
    ## Don't use genes with count 0.
    genes.withoutzero <- apply(input.mat, 1, min) > 0
  
    if (zerofilter) {
        ## :PLR: Which means "design" will no longer correspond to "input.mat" cols
        RNA <- input.mat[genes.withoutzero == TRUE, ]
    }

    ## Filtering step 2
    ## Identify genes that satisfy the linearity assumption with more than
    ## conc% probability for both up-regulated or down-regulated genes in tumor.
    Nsam <- RNA[, design == 0]
    Tsam <- RNA[, design == 1]

## :PLR: apply(x, 1, mean) can be rewritten as rowMeans(x)

    Nmean <- apply(Nsam, 1, mean) ##mean gene expression for normals
  

    Cons <- rep(0, length(Nmean)) #logical vector for genes that are highly expressed in tumor
    Nega <- rep(0, length(Nmean)) #logical vector for genes that are highly expressed in normal

    for (i in 1:length(Nmean)) {
        Cons[i] <- sum(Tsam[i, ] > Nmean[i])
        Nega[i] <- sum(Tsam[i, ] < Nmean[i])
    }
    cutoff_con <- round(sum(design == 1) * conc)

    ## Filtering step 3
    ## Default fold-change is set at 1.2 and 1/1.2
    ## fc values need to be adjusted to make final dataset have around 2000~3000.

    efc <- apply(Tsam, 1, mean) / (apply(Nsam, 1, mean)) ##ratio of the mean expression tumor vs. normal

    secfil <- (apply(RNA, 1, max) < 100000 &
               (efc > fc | efc < 1/fc) &
               ((Nega > cutoff_con) | (Cons > cutoff_con)) )

    nused <- sum(secfil) #numbers of to be filtered
    useornot <- secfil
  
    filtered.mat <- RNA[secfil, ]
    filtered.mat
}

