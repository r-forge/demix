###
### $Id$
###


##-----------------------------------------------------------------------------
DeMix.Normalization <- function(input,
                                design,
                                method=c("total", "quantile", "median")) {
    ## Check arguments
    stopifnot(is.matrix(input) || is.data.frame(input))
    input.mat <- as.matrix(input)
    stopifnot(is.matrix(input.mat) &&
              is.numeric(input.mat[, 1]) &&
              !anyNA(input.mat))
    method <- match.arg(method)
    stopifnot(is.numeric(design) && !anyNA(design))
    stopifnot(ncol(input.mat) == length(design))

    ## :PLR: Why throw away dimnames information?
    colnames(input.mat) <- NULL
    rownames(input.mat) <- NULL

    ## 0 denotes normal cell / 1 denotes tumor cell
    seqData <- DSS::newSeqCountSet(input.mat, design)

    ## Determine normalization factor
    seqData <- DSS::estNormFactors(seqData, method)
    k3 <- seqData@normalizationFactor
    k3.median <- median(k3)
    k3 <- k3 / k3.median

    ## Apply normalization factor
    for (i in 1:ncol(input.mat)) {
        input.mat[, i] <- input.mat[, i] / k3[i]
    }

    input.mat
}

