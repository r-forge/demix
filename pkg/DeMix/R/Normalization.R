###
### %Id%
###


##-----------------------------------------------------------------------------
DeMix.Normalization <- function(input,
                                method=c("total", "quantile", "median"),
                                groupid,
                                ...) {
  newt <- as.matrix(input)
  
  method <- match.arg(method)

  colnames(newt) <- NULL
  rownames(newt) <- NULL
    
  # 0 denotes normal cell / 1 denotes tumor cell
  seqData <- newSeqCountSet(as.matrix(newt), groupid)
  
  # Quantile normalization
   seqData <- estNormFactors(seqData, method)
   k3 <- seqData@normalizationFactor
   mk3 <- median(k3)
   k3 <- k3/mk3

  # Calculating normalization factor
   for(i in 1:ncol(newt))
     newt[,i] = newt[,i]/k3[i]

  newt
}

