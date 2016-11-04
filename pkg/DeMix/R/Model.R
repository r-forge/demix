###
### $Id$
###


# #####################################################################
# How to use :
#
# inputmat1 : matrix after RMA normalization.  G*S matrix
# where G is the number of genes and S is the number of samples.
#
# design : vector of indicators whether the sample is mixed or not.
# For matched 6 samples ex c(0,1,0,1, 0, 1)
#
# cit=64  which tells the bit of the machine. it can be 32 bit for some machines
#
# matched : if sample are matched, we need to specifiy matched=1 otherwise matched =0
# nref : if we have reference genes, just indicate the row numbers of those genes.
# #####################################################################

# inputmat1=as.matrix(inputmat1);design=cnvgroup; nhavepi= 1; givenpi=as.vector(testr2$pi[1,]) ;
#nPoi= 1;ninteg= 28


##-----------------------------------------------------------------------------
DeMix.model <- function(input,
                        design,
                        method=c("total", "quantile", "median"),
                        nhavepi,
                        givenpi,
                        ninteg,
                        ncore=detectCores()) {
    ## Check arguments
    stopifnot(is.matrix(input) && is.numeric(input[, 1]) && !anyNA(input))
    stopifnot(is.numeric(design) && !anyNA(design) && length(design) >= 2)
# :PLRL: "nhavepi" seems logical, but numeric?
    stopifnot(is.scalar.numeric(nhavepi))
# :PLR: Maybe default "givenpi" to null and check for that condition - eliminate "nhavepi" altogether
# :PLR: Precondition checking missing... (givenpi)
    stopifnot(is.scalar.numeric(ninteg))
    stopifnot(is.scalar.numeric(ncore))
    method <- match.arg(method)
# :PLR: Arguments should be ordered by likelihood of needing to specify them


    ## 

    norm.mat <- DeMix.Normalization(input, design, method)
    filtered.mat  <- DeMix.Filter(norm.mat, design)

    nsub   <- ncol(filtered.mat)
    ngenes <- nrow(filtered.mat)

    input.arr <- as.array(matrix(filtered.mat, nrow=1, byrow=FALSE))
    ntumor <- sum(design)
    nnormal <- nsub - ntumor

    rnan <- filtered.mat[, design == 0] # 0 denotes normal cell
    rnat <- filtered.mat[, design == 1] # 1 denotes tumor cell
    ovsn <- ((apply(rnan, 1, sd)^2) - apply(rnan, 1, mean)+1) / (apply(rnan, 1, mean)+1)^2
    ovst <- ((apply(rnat, 1, sd)^2) - apply(rnat, 1, mean)+1) / (apply(rnat, 1, mean)+1)^2

    newovs <-c(abs(ovsn), abs(ovst))

    if (nhavepi == 1) {
        if (!is.vector(givenpi)) {
             stop(sprintf("argument %s must be a numeric vector if pi is known",
                          sQuote("givenpi")))
        }
        givenpi <- as.array(givenpi)
    } else {
        givenpi <- as.array(rep(0, ntumor))
    }

    design <- as.array(design)

    min.ninteg <- 10
    if (ninteg < min.ninteg) {
        ninteg <- min.ninteg
    }

    nPoi <- as.integer(1)    # :TODO: What is 'nPoi' param for?
                             # :TODO: Remove from Bdemix()?
    seeds <- c(629555906, 921927245, 1265635378)


    rres <- .C("Bdemix",
               input.arr,
               as.integer(ncore),
               as.integer(design),
               as.integer(nsub),
               as.integer(ngenes),
               as.integer(nhavepi),
               givenpi,
               as.integer(nPoi),
               as.integer(ninteg),
               newovs,
               rep(0, ntumor*3),
               rep(0, ntumor*3),
               rep(0, nsub*ngenes),
               rep(0, nsub*ngenes),
               rep(0, 500*ngenes),
               rep(0, 2*ngenes),
               seeds)

    outcome2 <- matrix(rres[[14]], ncol=nsub, nrow=ngenes, byrow=TRUE)
    outcome2 <- outcome2[, ((nnormal+1):nsub)]

    outcome3 <- matrix(rres[[15]], ncol=nsub, nrow=ngenes, byrow=TRUE)
    outcome3 <- outcome3[, ((nnormal+1):nsub)]

    outcome1 <- matrix(rres[[12]], ncol=ntumor, nrow=3, byrow=TRUE)
    outcomePoi <- matrix(rres[[13]], ncol=ntumor, nrow=3, byrow=TRUE)
    post <- matrix(rres[[16]], ncol=500, nrow=ngenes, byrow=FALSE)
    mung <- matrix(rres[[17]], ncol=2, nrow=ngenes, byrow=FALSE)

    list(pi=outcome1,
         Poipi=outcomePoi,
         decov=round(outcome2, 0),
         decovn=outcome3,
         munt=mung)
}

