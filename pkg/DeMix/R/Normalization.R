###
### $Id$
###


##-----------------------------------------------------------------------------
DeMix.Normalization <- function(input,
                                method=c("total", "quantile", "median"),
                                groupid,
                                ...) {
    ## Check arguments
    stopifnot(is.matrix(input) || is.data.frame(input))
    method <- match.arg(method)
    stopifnot(is.numeric(groupid) && !anyNA(groupid))


    newt <- as.matrix(input)
    stopifnot(is.matrix(newt) && is.numeric(newt[, 1]) && !anyNA(newt))
    stopifnot(ncol(newt) == length(groupid))

    ## :PLR: Why throw away dimnames information?
    colnames(newt) <- NULL
    rownames(newt) <- NULL

    ## 0 denotes normal cell / 1 denotes tumor cell
    seqData <- DSS::newSeqCountSet(newt, groupid)

    ## Quantile normalization
    seqData <- DSS::estNormFactors(seqData, method)
    k3 <- seqData@normalizationFactor
    k3.median <- median(k3)
    k3 <- k3 / k3.median

    ## Calculating normalization factor
    for (i in 1:ncol(newt)) {
        newt[, i] <- newt[, i] / k3[i]
    }

    newt
}

