###
### $Id$
###


# #####################################################################
# How to use :
#
# inputmat1 : matrix after RMA normalization.  G*S matrix
# where G is the number of genes and S is the number of samples.
#
# groupid : vector of indicators whether the sample is mixed or not.
# For matched 6 samples ex c(0,1,0,1, 0, 1)
#
# cit=64  which tells the bit of the machine. it can be 32 bit for some machines
#
# matched : if sample are matched, we need to specifiy matched=1 otherwise matched =0
# nref : if we have reference genes, just indicate the row numbers of those genes.
# #####################################################################

# inputmat1=as.matrix(inputmat1);groupid=cnvgroup; nhavepi= 1; givenpi=as.vector(testr2$pi[1,]) ;
#nPoi= 1;ninteg= 28


##-----------------------------------------------------------------------------
DeMix.model <- function(input,
                        nnormal,
                        ntumor,
                        groupid=c(rep(0, nnormal),
                                  rep(1, ntumor)),
                        nhavepi,
                        givenpi,
                        ninteg,
                        ncore,
                        methods) {


input.norm <- DeMix.Normalization(input, method, groupid)
input.filt <- DeMix.Filter(input.norm, nnormal,ntumor, groupid=c(rep(0,nnormal),rep(1,ntumor))


    if(!is.matrix(inputmat1)) stop(sprintf("argument %s must be a matrix",
                                           sQuote("inputmat1")))
    if(!is.vector(groupid)) stop(sprintf("argument %s must be a vector",sQuote("groupid")))

seeds <- c(629555906, 921927245, 1265635378)



  nsub <- as.integer(dim(inputmat1)[2])
  wgenes <- as.integer(dim(inputmat1)[1])

  dataarray1 <- (as.array(matrix(inputmat1, nrow = 1, byrow = F)))
  intx=sum(groupid);
  intn=nsub-intx;


  rnan<-inputmat1[, groupid==0] # 0 denotes normal cell
  rnat<-inputmat1[, groupid==1] # 1 denotes tumor cell
  ovsn=((apply(rnan,1,  sd)^2)-apply(rnan,1,  mean)+1) / (apply(rnan,1,  mean)+1)^2
  ovst=((apply(rnat,1,  sd)^2)-apply(rnat,1,  mean)+1) / (apply(rnat,1,  mean)+1)^2

  newovs<-c(abs(ovsn), abs(ovst))

    if(nhavepi==1){
        if(!is.vector(givenpi)) stop("argument 5 must be a vector if pi is known");
    givenpi=(as.array(givenpi))
    } else {
    givenpi=(as.array(rep(0, intx)))
  }

    groupid<-(as.array(groupid))

  if(ninteg<10) ninteg=10;

    #dyn.load("main_firstparallel.so")

    rres <- .C("Bdemix", dataarray1, as.integer(ncore), as.integer(groupid), as.integer(nsub),
                         as.integer(wgenes),  as.integer(cbit), as.integer(nhavepi ), givenpi,
                         as.integer(nPoi), as.integer(ninteg), newovs, rep(0, intx*3), rep(0, intx*3),
                         rep(0, nsub*wgenes), rep(0, nsub*wgenes), rep(0, 500*wgenes), rep(0, 2*wgenes),
                         seeds)

    outcome2<-matrix(rres[[14]], ncol=(nsub) , nrow=wgenes, byrow = T)

  outcome2<-outcome2[, ((intn+1):nsub)];

  outcome3<-matrix(rres[[15]], ncol=(nsub) , nrow=wgenes, byrow = T)
  outcome3<-outcome3[, ((intn+1):nsub)];

    outcome1<-matrix(rres[[12]], ncol=intx, nrow=3, byrow=T)
  outcomePoi<-matrix(rres[[13]], ncol=intx, nrow=3, byrow=T)
  post<-matrix(rres[[16]], ncol=500, nrow=wgenes, byrow = F)
  mung<-matrix(rres[[17]], ncol=2, nrow=wgenes, byrow = F)

    list(pi=outcome1,
         Poipi=outcomePoi,
         decov=round(outcome2, 0),
         decovn=outcome3,
         munt=mung)
}

