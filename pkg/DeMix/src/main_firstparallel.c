
/***************************************************************************************************************
* RNAseq modeling
* Gene expression analysis under Three group is first done
* By Jaeil Ahn
* Initial Date 12.04.2011
* Last Update  08.04.2016
* input type : number of sample, phenotype, genes
* ex ./tp 3 1500
****************************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "bayes_para.h"
#include <pthread.h>
//  R CMD SHLIB main_firstparallel.c


// Newlly added for multiple thread
typedef struct threadInfoStruct {
  int offset;
  int opt;      // option
  unsigned long *seed;
  pthread_t threadId;
} threadInfo;
int numThreads;

//# 17 arguments
void Bdemix(double *data, int *ncore,int *nGroup, int *nsamp, int *ngenes, int *ct, int *npi, double *fixpi,
            int *nmodel, int *ninteg, double *newovs, double *output1, double *output2, double *output3,
            double *output32,  double *output4,  double *outputmu, double* seedv)
{

  seed = (unsigned long *)calloc(3,sizeof(unsigned long));
  for (int i=0; i < 3; ++i) {
    seed[i] = seedv[i];
  }

  printf("starting DeMixbayes... \n");
	int i, j, k, l ,s, q;
  int iteration2=iteration;
  int burn2;
  int icount, isave;
  int ii;

  double **tmppi, **tmppiPoi;

	//FILE  *fseed;
  clock_t start_t, end_t;

  numThreads=*ncore;
  nS=*nsamp;            // Number of Samples
	nG=*ngenes;           // Number of Genes
	nPoi=*nmodel;
  integ=*ninteg;
  nHavepi=*npi;

  // printf("Starting random number generator...\n");

  // Random number read
  //seed = (unsigned long *)calloc(3,sizeof(unsigned long));
	//fseed = fopen("/Users/davidgibbs/Code/DeMix/data/seed.dat","r");
	//for (i=0;i<3;i++) {
	//	fscanf(fseed,"%lu ",&(seed[i]));
  //}
	//fclose(fseed);

  // printf("....done with random number generator setup\n");

  // If we have pi, we do not need to go lots of iterations.
  if(nHavepi==1)
  {
    iteration2=500;
    burn2=490;
    // Average of the last ten iterations will be used for deconvolved expressions
  } else
    iteration2=iteration;

	M = exp(-*ct*log(2.0));  /* for *ct bit machines/compilers */

  printf(" ******************************************************************* \n");
  printf(" *   DeMix-Bayes Version 1.01                                      * \n");
  printf(" * Under %d bit computer machine with %d threads                   * \n", *ct, numThreads);
	printf(" * The input data consist of %d subjects with %d genes          * \n", nS, nG);

  p = (PARAM *)calloc(1,sizeof(PARAM));

	// Parameter initialized
	FD     =calloc(nS ,sizeof(double *));       // Where input data are saved
	CD     =calloc(nS ,sizeof(double **));       // Where output data are saved after deconvolution
	tmp_pos= calloc(integ, sizeof(double ));
  like     =calloc(nS ,sizeof(double));       // Liklihood values corresponding to each sample
  Dgroup   =calloc(nS ,sizeof(int));          // Where group ID is save, 1 for tumor 0 for normal
  fake     =calloc(nG ,sizeof(double));

	for(j=0;j<nS;j++)
	{
    FD[j]= calloc(nG, sizeof(double ));
    CD[j]= calloc(nG, sizeof(double  *));

    for(k=0;k<nG;k++)
      CD[j][k]=calloc(2, sizeof(double));

	}

	initialSet(p);
	ReSet(p, newovs, seed);
	load_data(data, nGroup, fixpi,  seed);
  printf(" * A total of %d pure tissue and  %d mixed tissue samples        *  \n", fNorm,intx );
  printf(" ******************************************************************* \n");

  icount=(int)((iteration2-burn)/10);

  tmppi     =calloc(icount ,sizeof(double *));       // Pi save NB
  tmppiPoi  =calloc(icount ,sizeof(double *));       // Pi save Poi

	for(j=0;j<icount;j++){
    tmppi[j]= calloc(intx, sizeof(double ));
    tmppiPoi[j]= calloc(intx, sizeof(double ));
  }

  /*
   Big notes
   IF pis are the same, the mle of the mean expressions for NB or Poisson will be the same
   In this program, we make use of that fact to unburden computational cost.
  */

  isave=0;
  printf("Step 1: MCMC %d iteration will begin \n", iteration2);
  q=0;
  setnormalmean();
  threadInfo *threads = (threadInfo*)calloc(numThreads,sizeof(threadInfo));


	for(i=0;i<iteration2;i++)
	{
    if(i==6) start_t=clock();

    // If we do not have information on the proportions.
    if(nHavepi==0)
    {
      // For a few iteration, we do not update Pi
      if(i>0) //CHNGE BACK TO FIVE
      {
        for(ii=0;ii<numThreads;++ii) {
           threads[ii].offset = ii;
           threads[ii].seed = seed;
           pthread_create(&((threads+ii)->threadId), 0, getpi2_parallel,(void*)(threads+ii));
        }
        for(ii=0;ii<numThreads;++ii) {
          if(pthread_join(threads[ii].threadId, 0)!=0) {
             perror("thread failure");
             exit(-1);
           }
        }

      } else {
        for(l=0;l<nS;l++)
          if(Dgroup[l]==1)
            p->pi[l]= 0.5;
      }
    }

    // Mu_t, sigma_t, Mu_compute
    for(k=0;k<3;k++)
    {

      // K0, k1 for tumor 2 for normal
      if(k<2)
      {
        for(ii=0;ii<numThreads;++ii) {
          threads[ii].offset = ii;
          threads[ii].seed = seed;
          threads[ii].opt = k;
          pthread_create(&((threads+ii)->threadId), 0, gettumor_parallel,(void*)(threads+ii));
        }
      } else {
        for(ii=0;ii<numThreads;++ii) {
          threads[ii].offset = ii;
          threads[ii].seed = seed;
          threads[ii].opt = 1;
          pthread_create(&((threads+ii)->threadId), 0, getnormal_parallel,(void*)(threads+ii));
        }
      }

      for(ii=0;ii<numThreads;++ii) {
        if(pthread_join(threads[ii].threadId, 0)!=0) {
          perror("thread failure");
          exit(-1);
        }
      }
    }

    // Posteriors will be saved at every 10th iterations
  	if(i%10==0 && i >10)
		{
      printf("%d th iteration \n", i);
      if(nHavepi==0)
      {
        // Discard the samples before 'burn-in'
        if(i>=burn)
        {
          s=0;

          /* This process is not necessary. Just give insigiht if Poisson distribution is used*/
          for(l=0;l<nS;l++)
          {
            if(Dgroup[l]==1 )
            {
              // NB
              tmppi[q][s] =p->pi[l];

              // Poisson
              like[l]=golden3(l, getpi3);
              tmppiPoi[q][s] =like[l];

              s=s+1;
            }
          }
          q+=1;
        }

        for(l=0;l<nS;l++)
        {
          if(Dgroup[l]==1 )
          {
            if(nPoi==0 || i<burn)
              printf(" %15.3f \t",  p->pi[l]);
            else
              printf(" %15.3f \t",  like[l]);
          }
        }
        printf("\n");
      }
		}

    // Expressions in the last 10 samples will be averaged
    if(i>iteration2-10)
    {
      if(i==(iteration2-9)) printf("Step 2: Deconvolution starts\n");
      for(l=0;l<nG;l++)
      {
        for(k=0;k<nS;k++)
        {
          if(Dgroup[k]==1)
          {
            CD[k][l][1] +=golden2(k, l, 1, MutliLike);
          }

        }
      }
    }

    if(i==10)
    {
      end_t=clock();
      printf("================================================================\n");
      printf("Total time : %lf hours\n", (double)(end_t-start_t)/CLOCKS_PER_SEC*((double)iteration2/numThreads/5.0/3600.0));
      printf("================================================================\n");
    }
	}

	// Given the mean/variance of normal and tumor, we optimize the deconvoled expression using the Golden
  // search algorithm, which is one-dimensional question.

  printf("Step 2: Deconvolution starts\n");
	for(l=0;l<nG;l++)
	{
		for(i=0;i<nS;i++)
		{
			if(Dgroup[i]==1)
			{
				CD[i][l][1] = CD[i][l][1]/10.0;
				CD[i][l][0] = ((double)FD[i][l]- p->pi[i]*CD[i][l][1])/(1- p->pi[i]);
   		}

		}
	}
  free(threads);

  printf("Step 3: Output is being saved \n");

  // Summarize the outcome
  minimax(tmppi, icount, intx, output1);
  minimax(tmppiPoi, icount, intx, output2);

  saveFiles(output3, output32);

  // pure tumor expression save.
  for(j=0;j<nG;j++)
  {
    outputmu[j*2]=  p->Navg[j];
    outputmu[j*2+1]=  p->Tavg[j];
  }
  printf("Step 4: Saving outputs  : \n");

  /*

  // Release memory
  for(j=0;j<nS;j++)
	{
    free(FD[j]);
    free(CD[j]);
	}
  free(FD); free(CD);

 	for(j=0;j<intx;j++)
  {
    free(tmppi[j]);
    free(tmppiPoi[j]);

  }
  free(tmppi);
  free(tmppiPoi);
   */

}


// Reading the expressions
void load_data(double *mat1, int *group,  double *fixpi, unsigned long *seed)
{
	int  j,k;

  for(k=0;k<nS;k++)
	{
		for(j=0;j<nG;j++)
		{
			FD[k][j]=  mat1[nG*k+j];
	  }
	}

  // 0 for normal 1 for tumor,, group data
	fNorm=0;
	for(j=0;j<nS;j++)
	{
		Dgroup[j]=group[j];
		if(Dgroup[j]==0) fNorm++;

	  //printf("%d : %d \n", j, Dgroup[j]);
	}

	// Determine the number of tumor sample
	intx=nS-fNorm;


  if(nHavepi==1)
  {
    for(k=0;k<intx;k++)
      p->pi[fNorm+k]= fixpi[k];
  }

}

void saveFiles( double *put3, double *put4)
{
	int i, j,k;

	// pure tumor expression save.
	for(j=0;j<nG;j++)
	{
		k=0;
		for(i=0;i<nS;i++)
		{
			if(Dgroup[i]==1)
      {
				put3[j*nS+k]=  CD[i][j][1];
        put4[j*nS+k]=  CD[i][j][0];
      }
			k++;
		}
	}
}


void initialSet(PARAM *qq)
{
	qq->Navg =calloc(nG ,sizeof(double ));
	qq->Tavg =calloc(nG ,sizeof(double ));
  qq->Favg =calloc(nG ,sizeof(double ));

	qq->Neta =calloc(nG ,sizeof(double ));
	qq->Teta = calloc(nG ,sizeof(double ));

	qq->pi =calloc(nS ,sizeof(double ));
}


// Setting the initial values of all parameters
void ReSet(PARAM *qq, double *obs, unsigned long *seed)
{
	int i, j;
	for(j=0;j<nG;j++)
	{
		qq->Navg[j]=1000.0;
		qq->Tavg[j]=1000.0;

    qq->Neta[j]=obs[j];
    qq->Teta[j]=obs[j+nG];

  }

  if(nHavepi==0)
  {
    for(i=0;i<nS;i++)
    {
      if(Dgroup[i]==1)
      {
        p->pi[i]=0.7;
      }
    }
  }
}

/*************************************************************
 Normal  mean , variance calculation in log or original scale
**************************************************************/

// Setting the MLE of the mean
void setnormalmean()
{
	int i, j;

	for(j=0;j<nG;j++)
	{
    p->Navg[j]=0.0;
    p->Tavg[j]=0.0;

    for(i=0;i<fNorm;i++)
    {
      p->Navg[j]+=FD[i][j];

    }
    for(i=fNorm;i<nS;i++)
    {
      p->Tavg[j]+=FD[i][j];

    }
    p->Navg[j]=p->Navg[j]/(double)fNorm;
    p->Tavg[j]=p->Tavg[j]/(double)intx;

    p->Tavg[j]=p->Navg[j] +  1.5*(p->Tavg[j]-p->Navg[j]);


    p->Favg[j]=p->Tavg[j];
    fake[j] = p->Tavg[j];
  }
}

// For normal sample, NB
double PoiGamma(int Sid, int Gid, int opt)
{
	double rtv=0;
	double muig=0;

	muig=(p->Navg[Gid]);

	if(opt==1)
		rtv=gammaln(FD[Sid][Gid]+1.0/p->Neta[Gid])-gammaln(1.0/p->Neta[Gid]);

  rtv=rtv -1.0/p->Neta[Gid]*log(1+muig*p->Neta[Gid]) + FD[Sid][Gid]*(log(muig*p->Neta[Gid]) -log(1+muig*p->Neta[Gid]));

	return rtv;
}


// If we optimize pi w.r.t. the Poisson model
double PoiGammaPi(int Sid, int Gid, double tmppi)
{
	double rtv=0;
	double muig=0;

	muig=(p->Navg[Gid]*(1-tmppi)+ p->Tavg[Gid]*tmppi );

  rtv=FD[Sid][Gid]* log(muig) -muig;

	return rtv;
}

/* For update ratio calculation, we need to fix position at the same number*/
double ft_y(int ntmp, double *tmppos, int samp, int genes)
{
	int i;
	long double tmp=0.0;

	double def_t;

  // Integrations
	for(i=0;i<integ;i++)
	{
		tmppos[i]=(double)FD[samp][genes]/(double)integ*(i+0.5);
		def_t=FD[samp][genes]-tmppos[i];

    // Trapozodal  O(1/N^3)
    if(i==0 || i==(integ-1))
      tmp+=0.4167*NB_nt3(tmppos[i], def_t,samp, genes);
    else if (i==1 || i==(integ-2))
      tmp+=1.0833*NB_nt3(tmppos[i], def_t,samp, genes);
    else
      tmp+=NB_nt3(tmppos[i], def_t,samp, genes);

	}

	if(tmp<=0) tmp=4.94066e-311;

	return (log(tmp*FD[samp][genes] ));
}


long double NB_nt3(double yt, double xt, int samp, int genes)
{
	double tmpn;
  double tmpy;

	if(xt==0 || yt==0) return 0.0;

	yt=yt/(p->pi[samp]);
	xt=xt/(1-p->pi[samp]);

	tmpn=  gammaln(xt +1.0/p->Neta[genes])-gammaln(xt)-gammaln(1.0/p->Neta[genes]);
	tmpn+= gammaln(yt +1.0/p->Teta[genes])-gammaln(yt)-gammaln(1.0/p->Teta[genes]);

  tmpy= p->Navg[genes]*p->Neta[genes];

  tmpn+= 1.0/p->Neta[genes]*log(1.0/(1+tmpy));
	tmpn+=xt *log(tmpy/(1+tmpy));


  tmpy= p->Tavg[genes]*p->Teta[genes];

  tmpn+=1.0/p->Teta[genes]*log(1.0/(1+tmpy));
	tmpn+=yt*log(tmpy/(1+tmpy));

	if(tmpn>700) printf("Warning : %lf  \n",tmpn);

	return exp(tmpn);
}



long double Poisson_nt3(double yt, double xt, int samp, int genes)
{
	double tmpn;

	if(xt==0 || yt==0) return 0.0;

	yt=yt/(p->pi[samp]);
	xt=xt/(1-p->pi[samp]);

	tmpn=  -gammaln(xt);
	tmpn+= -gammaln(yt);

	tmpn+=xt *log(p->Navg[genes])  + p->Navg[genes] - log(xt*(1-p->pi[samp]) );
	tmpn+=yt *log(p->Tavg[genes])  + p->Tavg[genes] - log(yt*(p->pi[samp]) );

	if(tmpn>700) printf("Warning : %lf  \n",tmpn);

	return exp(tmpn);
}


void getnormal(int genes, int opt, unsigned long *seed)		// opt 0 for Mu_n, 1 for sigma_n
{
	int i;
	//double tmp=0.0;
	double old_like, new_like;
	double oldval, newval;
	double priora, priorb;

  // Allow different random walk sd for the faster convergence
	if(opt==0)		// mean
	{
		oldval=p->Navg[genes];
		if(oldval<10)
			newval=rnorm(oldval, 0.5, seed);
		else if(oldval<20)
			newval=rnorm(oldval, 2, seed);
		else if(oldval<50)
			newval=rnorm(oldval, 4.0, seed);
		else if(oldval<100)
			newval=rnorm(oldval, 10.0, seed);
		else if(oldval<500)
			newval=rnorm(oldval, 30, seed);
		else if(oldval<1000)
			newval=rnorm(oldval, 80, seed);
		else if(oldval<2000)
			newval=rnorm(oldval, 200, seed);
		else if(oldval<10000)
			newval=rnorm(oldval, 500, seed);
		else if(oldval<40000)
			newval=rnorm(oldval, 2000, seed);
    else
      newval=rnorm(oldval, 4000, seed);

		if(newval <0.05)
			return;
	} else 	{			//variance
		oldval=p->Neta[genes];
		newval=rnorm(oldval, 0.0005, seed);
		if(newval <0.000005)
			return;
	}
	priora=priorb=0;

	new_like=0.0;old_like=0.0;
	if(opt==0)			// 0 for mean
	{
		for(i=0;i<nS;i++)
		{
			//printf("%d \t %d : %lf \t %lf \t  %lf \t %lf  \t %lf  \n", genes, i, NIG[i][genes], FD[i][genes], nig, tig, p->pi[i]);
			if(Dgroup[i]==0)
			{
        /*
				p->Navg[genes]=newval;
				tmp= PoiGamma(i, genes, 0);
				new_like+= tmp;

				p->Navg[genes]=oldval;
				tmp= PoiGamma(i, genes, 0);
				old_like+= tmp;
         */
			} else {
				/*
				p->Navg[genes]=newval;
				new_like+=  ft_y(integ, tmp_pos, i,  genes);

				p->Navg[genes]=oldval;
				old_like+= ft_y(integ, tmp_pos, i,  genes);
				*/
			}

		}
	} else {   // eta
		for(i=0;i<nS;i++)
		{

			if(Dgroup[i]==0)
			{
        // Used before 02.09.2014, PoiGamma is using the NB
        /*
  			p->Neta[genes]=newval;
				tmp= PoiGamma(i, genes, 1);
				new_like+= tmp;

				p->Neta[genes]=oldval;
				tmp= PoiGamma(i, genes, 1);
				old_like+= tmp;
        */

        p->Neta[genes]=newval;
				new_like+= ft_y(integ, tmp_pos, i,  genes);

				p->Neta[genes]=oldval;
				old_like+= ft_y(integ, tmp_pos, i,  genes);
			} else {
        /*
        p->Neta[genes]=newval;
				new_like+= ft_y(integ, tmp_pos, i,  genes);

				p->Neta[genes]=oldval;
				old_like+= ft_y(integ, tmp_pos, i,  genes);
        */
        /*
				p->Neta[genes]=newval;
				new_like+=  ft_y(integ, tmp_pos, i,  genes);

				p->Neta[genes]=oldval;
				old_like+= ft_y(integ, tmp_pos, i,  genes);
				*/
			}
		}
	}
  if(opt==1)
  {
    //0.05 is the previous mean
    priora=-0.05*(newval-0.05)*(newval-0.05);
    priorb=-0.05*(oldval-0.05)*(oldval-0.05);

  }
	if(log(kiss(seed))<new_like-old_like+priora-priorb)
	{
		if(opt==0)
			p->Navg[genes]=newval;
		else
			p->Neta[genes]=newval;

	} else {
		if(opt==0)
			p->Navg[genes]=oldval;
		else
			p->Neta[genes]=oldval;
	}
}



void* gettumor_parallel(void* args) {
  threadInfo *info = (threadInfo*)args;
  int offset = info->offset;
  int nstart, nend, nsect;

  int l;

  nsect=(int)(nG/numThreads);
  nstart=(offset)*nsect;
  nend=(offset+1)*nsect;

  if(offset==(numThreads-1))
    nend=nG;

  //printf("%d : %d  \t %d \n",info->opt, nstart, nend );
  for(l=nstart;l<nend;l++)
  {
      gettumor(l, info->opt, info->seed)   ;
  }
  return 0;
}

void* getnormal_parallel(void* args) {
  threadInfo *info = (threadInfo*)args;
  int offset = info->offset;
  int nstart, nend, nsect;

  int l;

  nsect=(int)(nG/numThreads);
  nstart=(offset)*nsect;
  nend=(offset+1)*nsect;

  if(offset==(numThreads-1))
    nend=nG;


  for(l=nstart;l<nend;l++)
  {
    getnormal(l, info->opt, info->seed)   ;
  }
  //free(info);
  return 0;
}


void gettumor(int genes, int opt, unsigned long *seed)		// opt 0 for Mu_n, 1 for sigma_n
{
	int i;
	double old_like, new_like;
	double oldval, newval;
	double priora, priorb;
	double tmpprior;

  // Allow different random walk sd for the faster convergence
	if(opt==0)		// mean
	{
		oldval=p->Tavg[genes];
		if(oldval<10)
			newval=rnorm(oldval, 0.5, seed);
		else if(oldval<20)
			newval=rnorm(oldval, 2, seed);
		else if(oldval<50)
			newval=rnorm(oldval, 4.0, seed);
		else if(oldval<100)
			newval=rnorm(oldval, 10.0, seed);
		else if(oldval<500)
			newval=rnorm(oldval, 30, seed);
		else if(oldval<1000)
			newval=rnorm(oldval, 100, seed);
		else if(oldval<2000)
			newval=rnorm(oldval, 200, seed);
		else if(oldval<10000)
			newval=rnorm(oldval, 700, seed);
    else if(oldval<40000)
      newval=rnorm(oldval, 2000, seed);
    else
      newval=rnorm(oldval, 4000, seed);

		if(newval <0.005)
			return;

	} else 	{			//variance
		oldval=p->Teta[genes];
		newval=rnorm(oldval, 0.0005, seed);
		if(newval <0.00005)
			return;
	}
	priora=priorb=0;

	tmpprior=0.05;

	// Having prior on variance is better or not // Not done
	if(opt==1)
	{
		tmpprior=-2.9;
		priora = -log(oldval) - 2*pow(log(oldval)-tmpprior, 2)*0.05;
		priorb = -log(newval) - 2*pow(log(newval)-tmpprior, 2)*0.05;
	}

	new_like=0.0;old_like=0.0;
	if(opt==0)			// 0 for mean
	{
		for(i=0;i<nS;i++)
		{
			if(Dgroup[i]==1)
			{
				p->Tavg[genes]=newval;
				new_like+= ft_y(integ, tmp_pos, i,  genes);

				p->Tavg[genes]=oldval;
				old_like+= ft_y(integ, tmp_pos, i,  genes);
			}
		}
	} else {
		for(i=0;i<nS;i++)
		{
			if(Dgroup[i]==1)
			{
				p->Teta[genes]=newval;
				new_like+= ft_y(integ, tmp_pos, i,  genes);

				p->Teta[genes]=oldval;
				old_like+= ft_y(integ, tmp_pos, i,  genes);
			}
		}
	}

  // Prior ratio is ignored
	if(log(kiss(seed))<new_like-old_like+priora-priorb)
	{
		if(opt==0)
			p->Tavg[genes]=newval;
		else
			p->Teta[genes]=newval;

	} else {
		if(opt==0)
			p->Tavg[genes]=oldval;
		else
			p->Teta[genes]=oldval;
	}
}

void* getpi2_parallel(void* args) {
  threadInfo *info = (threadInfo*)args;
  int offset = info->offset;
  int l;
  for(l=0;l<nS;l++)
    if((Dgroup[l]==1)&&( (l%numThreads)==offset) )
    {
     like[l]=getpi2(l, info->seed)   ;
    }

  //free(info);
  return 0;
}


// ADD filtering step
double getpi2(int samp,  unsigned long *seed)		// opt 0 for Mu_n, 1 for sigma_n
{
	int i;
	double tmp=0.0, tmp2=0.0;
	double old_like, new_like;
	double oldval, newval;
	double priora, priorb;
	int tmpcnt;
  int ncnt, ocnt;

	oldval=p->pi[samp];
	newval=rnorm(oldval, 0.05, seed);
  new_like=0.0;old_like=0.0;



  // If too low or too high, then return the old.
  if(newval> 0.95 || newval <0.05)  return old_like;

	priora=priorb=0;


	tmpcnt=0;
  ncnt=ocnt=0;

	for(i=0;i<nG;i++)
	{
    // Filter out genes with less accuracy
    // Hear filtering is applied to sort out genes that do not have appropriate charateristics
    if(p->Neta[i] < 0.4 && p->Teta[i] < 0.4 &&   (p->Favg[i]/p->Tavg[i]) < 1.5 &&   (p->Favg[i]/p->Tavg[i]) > 0.5)
    {
      // Second filter : IF M>N then T>N,  IF M<N then T < N or the computation is wrong
      if ( (fake[i] > p->Navg[i] && p->Tavg[i]> fake[i]) || (fake[i] < p->Navg[i] && p->Tavg[i]< fake[i]) )
      {
        p->pi[samp]=newval;
        tmp= ft_y(integ, tmp_pos, samp,  i);
        new_like+=tmp;

        p->pi[samp]=oldval;
        tmp2= ft_y(integ, tmp_pos, samp,  i);
        old_like+=tmp2;

        if(tmp>=tmp2)
          ncnt++;
        else
          ocnt++;
      }
    } else {
      // This is a mandatory condtion to be qualified.
      if ( (fake[i] > p->Navg[i] && p->Tavg[i]> fake[i]) || (fake[i] < p->Navg[i] && p->Tavg[i]< fake[i]) )
      {
        p->pi[samp]=newval;
        tmp= PoiGammaPi(samp, i, newval);
        new_like+=tmp;

        p->pi[samp]=oldval;
        tmp2= PoiGammaPi(samp, i, oldval);
        old_like+=tmp2;

        if(tmp>=tmp2)
          ncnt++;
        else
          ocnt++;
      }
    }
  }

  if(isnan(new_like)== 1.0)
    printf("%d error %lf \n", samp, p->pi[samp]);

  if(isnan(old_like)== 1.0)
  printf("%d error %lf \n", samp, p->pi[samp]);

  if((ncnt > 1.25*(ocnt+0.5)) || (log(kiss(seed))<new_like-old_like))
  {
		p->pi[samp]=newval;
    for(i=0;i<nG;i++)
      p->Favg[i]=Fake_avg(i);

    return new_like;
	} else {
      p->pi[samp]=oldval;
    return old_like;
	}
}


// Compute fake mean, which should be very close to the real mean
double Fake_avg(int gene)
{
  int i;
  double tmp=0.0;

  for(i=fNorm;i<nS;i++)
    tmp+=(FD[i][gene]-p->Navg[gene]*(1-p->pi[i])) /p->pi[i];

  return tmp/(double)intx;
}


// Change it to Poisson
double getpi3(int samp, double tmppi)		// opt 0 for Mu_n, 1 for sigma_n
{
	int i;
	double new_like=0.0;

  for(i=0;i<nG;i++)
	{
    new_like+=PoiGammaPi(samp, i, tmppi );
  }
  return new_like*(-1.0);
}

/*************************************************************
 Normal  mean , variance calculation in log or original scale
 **************************************************************/

/* For update ratio calculation, we need to fix position at the same number*/
double ft_logy(int ntmp, double *tmppos, double y, double mung, double mutg, double sng, double stg, double pi1)
{
	int i;
	double rtval;
	double tmp;

	rtval=0;
	for(i=0;i<integ;i++)
		tmp_pos[i] = (y)/(double)integ*(i+0.5);

	for(i=0;i<ntmp;i++)
	{
		tmp= pow(sng, -0.5)* pow(stg, -0.5)*1.0/(y-tmppos[i]) * 1.0/tmppos[i];
		tmp=tmp*exp(-0.5 *	pow(log(y-tmppos[i])-(log(1-pi1) + mung), 2.0)/sng);
		tmp=tmp*exp(-0.5 *	pow(log(tmppos[i])-(log(pi1) + mutg), 2.0)/stg);

		if(tmp<=0)
		{
			tmp=0.000000000000001;
		}

		rtval +=tmp;
	}
	return (log(rtval/(double)ntmp));
}




/*************************************************************************************************
 Minimax Function
 *************************************************************************************************/

// Return value computations
void minimax(double **mat, int nrow, int ncol, double *rtmat)
{
  int i, j;
  double tmpvar, tmpmean, lowc, highc;
  
  for(i=0;i<ncol;i++)
  {
    
    tmpmean=0.0;
    tmpvar=0.0;
    
    for(j=0;j<nrow;j++)
    tmpmean+=mat[j][i];
    
    tmpmean=(double)tmpmean/(double)(nrow);
    
    
    for(j=0;j<nrow;j++)
    tmpvar+=(tmpmean- mat[j][i])*(tmpmean- mat[j][i]);
    
    tmpvar= (double)tmpvar/(double)(nrow);
    
    
    
    lowc=( tmpmean-1.96*pow(tmpvar, 0.5) < 0) ? 0:( tmpmean-1.96*pow(tmpvar, 0.5) );
    highc=( tmpmean+1.96*pow(tmpvar, 0.5) > 1) ? 1:( tmpmean+1.96*pow(tmpvar, 0.5));
    
    
    
    rtmat[i+ncol*1]= lowc;
    rtmat[i+ncol*2]=highc;
    rtmat[i+ncol*0]=tmpmean;
  }
}



double golden3(int v3, double (*fun)(int, double))
{
	int xcount;
	double xmin, fmin;
	const double tol=0.0001;
	const double R=0.61803399;
	const double C = 1.0-R;

	double x1, x2;
	double x0, x3;

	x0=  0.05;
	x3= 0.95;

	x1=x0 + 0.333 * x3;
	x2=x1 + C*(x3-x1);

	double f1=(*fun)(v3, x1);
	double f2=(*fun)(v3, x2);
	xcount=0;
	while(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2)))
	{
		xcount++;
		if(f2<f1)
		{
			shft3(&x0, &x1, &x2, R*x2+C*x3);
			shft2(&f1, &f2, ((*fun)(v3, x2)));
		} else {
			shft3(&x3, &x2, &x1, R*x1+C*x0);
			shft2(&f2, &f1, ((*fun)(v3, x1)));
		}

		if(xcount>100) return x1;
	}

	if(f1<f2){
		xmin=x1;
		fmin=f1;
	} else {
		xmin=x2;
		fmin=f2;
	}

	return xmin;
}


double golden2(int v3, int v4, int which, double (*fun)(int, int, double))
{
	int xcount;
	double xmin, fmin;
	const double tol=0.0001;
	const double R=0.61803399;
	const double C = 1.0-R;

	double x1, x2;
	double x0, x3;

	x0=  0.0;
  //x3= p->Tavg[v4]*10;
	x3=(double)FD[v3][v4]/p->pi[v3];


	x1=x0 + 0.333 * x3;
	x2=x1 + C*(x3-x1);

	double f1=(*fun)(v3, v4, x1);
	double f2=(*fun)(v3, v4, x2);
	xcount=0;
	while(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2)))
	{
		xcount++;
		if(f2<f1)
		{
			shft3(&x0, &x1, &x2, R*x2+C*x3);
			shft2(&f1, &f2, ((*fun)(v3,v4,  x2)));
		} else {
			shft3(&x3, &x2, &x1, R*x1+C*x0);
			shft2(&f2, &f1, ((*fun)(v3, v4, x1)));
		}

		if(xcount>100) return x1;
	}

	if(f1<f2){
		xmin=x1;
		fmin=f1;
	} else {
		xmin=x2;
		fmin=f2;
	}

	return xmin;
}

// NB model density multiplication
double MutliLike(int opt, int opt2,  double trpi)
{
	double tmpn;
	double xt, yt;

	yt= trpi;
	xt = (FD[opt][opt2] - p->pi[opt]*trpi)/(1-p->pi[opt]);

	tmpn=  gammaln(xt +1.0/p->Neta[opt2])-gammaln(xt)-gammaln(1.0/p->Neta[opt2]);
	tmpn+= gammaln(yt +1.0/p->Teta[opt2])-gammaln(yt)-gammaln(1.0/p->Teta[opt2]);

	tmpn+= 1.0/p->Neta[opt2]*log(1.0/(1+p->Navg[opt2]*p->Neta[opt2]));
	tmpn+=xt *log(p->Navg[opt2]*p->Neta[opt2]/(1+p->Navg[opt2]*p->Neta[opt2]));

	tmpn+=1.0/p->Teta[opt2]*log(1.0/(1+p->Tavg[opt2]*p->Teta[opt2]));
	tmpn+=yt*log(p->Tavg[opt2]*p->Teta[opt2]/(1+p->Tavg[opt2]*p->Teta[opt2]));



	// We want to maximize the multiplicaton of two likelihoood
	return (tmpn*(-1.0));
}



void shft2(double *a, double *b, double c)
{
	*a=*b;
	*b=c;
}

void shft3(double *a, double *b, double *c, const double d)
{
	*a=*b;
	*b=*c;
	*c=d;
}



/* Old NB distribution likelihood
// Removed
double NB_nt(double xt, double yt, int samp, int genes)
{
  double tmpn;
  double x;

  x= yt-xt;

  tmpn= gammaln(x/(1-p->pi[samp])+1.0/p->Neta[genes])-gammaln(x/(1-p->pi[samp])+1.0)-gammaln(1.0/p->Neta[genes]);
  tmpn+=gammaln(xt/(p->pi[samp]) +1.0/p->Teta[genes])-gammaln(xt/(p->pi[samp])+1.0)-gammaln(1.0/p->Teta[genes]);

  if(x==0) return 0.0;

  tmpn+=  1.0/p->Neta[genes]*1.0/(1-p->pi[samp]) *log(1.0/(1+p->Navg[genes]*p->Neta[genes]));
  tmpn+=  x/(1-p->pi[samp])* log(p->Navg[genes]*p->Neta[genes]/(1+p->Navg[genes]*p->Neta[genes]));

  tmpn+= 1.0/p->Teta[genes]*1.0/(p->pi[samp])*log(1.0/(1+p->Tavg[genes]*p->Teta[genes]));
  tmpn+=  xt/p->pi[samp]*log(p->Tavg[genes]*p->Teta[genes]/(1+p->Tavg[genes]*p->Teta[genes]));

  return exp(tmpn);
}



// Pi is reflected before calling this function
double NB_nt2(double xt, double yt, int samp, int genes)
{
  double tmpn;

  if(xt==0 || yt==0) return 0.0;

  tmpn=  gammaln(xt +1.0/p->Neta[genes])-gammaln(xt)-gammaln(1.0/p->Neta[genes]);
  tmpn+= gammaln(yt +1.0/p->Teta[genes])-gammaln(yt)-gammaln(1.0/p->Teta[genes]);

  tmpn+=1.0/p->Neta[genes]*log(1.0/(1+p->Navg[genes]*p->Neta[genes]));
  tmpn+=xt*log(p->Navg[genes]*p->Neta[genes]/(1+p->Navg[genes]*p->Neta[genes]));

  tmpn+=1.0/p->Teta[genes]*log(1.0/(1+p->Tavg[genes]*p->Teta[genes]));
  tmpn+=yt*log(p->Tavg[genes]*p->Teta[genes]/(1+p->Tavg[genes]*p->Teta[genes]));

  return exp(tmpn);
}
 */
