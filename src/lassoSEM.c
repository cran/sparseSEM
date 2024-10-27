#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <stdlib.h>
#ifndef FCONE
# define FCONE
#endif

#include "lassoSEM.h" 


/*
void printMat(double *a, int M, int N) //MxN
{
  int i,j;
  Rprintf("Printing the matrix\n\n");
  for(i=0;i<M;i++) 
  {
    for(j=0;j<N;j++)
    {
      Rprintf("%f\t", a[j*M +i]); //a[i,j]
    }
    Rprintf("\n");
  }
}
*/




void centerYX(double *Y,double *X, double *meanY, double *meanX,int M, int N) //M genes; N samples
{

  int i,index;	
  double *Xptr;
  double *Yptr;
  
  int inci = 1;
  int incj = 1;
  int inc0 = 0;
  int lda  = M; //leading dimension
  double *eye;
  eye = (double* ) R_Calloc(N, double);
  double alpha = 1;
  double beta = 0;
  F77_CALL(dcopy)(&N,&alpha,&inc0,eye,&inci);
  char transa = 'N';
  
  F77_CALL(dgemv)(&transa, &M, &N,&alpha, X, &lda, eye, &inci, &beta,meanX, &incj FCONE);
  F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &lda, eye, &inci, &beta,meanY, &incj FCONE);
  double scale;
  scale = 1.0/N;
  F77_CALL(dscal)(&M,&scale,meanY,&inci);
  F77_CALL(dscal)(&M,&scale,meanX,&inci);
  // OUTPUT Y X, set missing values to zero
  scale = -1;
  for(i=0;i<N;i++)
  {
    index = i*M;
    Xptr = &X[index];
    Yptr = &Y[index];
    F77_CALL(daxpy)(&M,&scale,meanY,&inci,Yptr,&incj);
    F77_CALL(daxpy)(&M,&scale,meanX,&inci,Xptr,&incj);
  }
  R_Free(eye);
}	



double constrained_ridge_cff(double *Ycopy, double *Xcopy, double rho_factor, int M, int N,
                             double *B, double *f, double *mue, int verbose)
{
  
  int i,j,k,lda,ldb,ldc,ldk;
  // center Y, X
  double *meanY, *meanX;
  meanY = (double* ) R_Calloc(M, double);
  meanX = (double* ) R_Calloc(M, double);
  
  //copy Y, X; 
  double *Y, *X;
  int MN = M*N;
  Y = (double* ) R_Calloc(MN, double);
  X = (double* ) R_Calloc(MN, double);
  
  int inci = 1;
  int incj = 1;
  F77_CALL(dcopy)(&MN,Ycopy,&inci,Y,&incj);
  F77_CALL(dcopy)(&MN,Xcopy,&inci,X,&incj);
  
  centerYX(Y,X,meanY, meanX,M, N);
  
  if(verbose>7) Rprintf("\t\t\t\t\t\t\t\tEnter Function: Ridge Regression. Shrinkage ratio rho is: %f.\n\n",rho_factor);
  
  int Mi = M -1;
  //for usage in loop
  double *YiPi; //Yi'*Pi
  YiPi =(double* ) R_Calloc(Mi*N, double);
  double xixi,xixiInv; //xi'*xi;
  int jj,index; //jj = 1:(M-1) index of YiPi
  double normYiPi,rho;
  double *bi,*YiPi2Norm; 	//YiPi2Norm: first term of biInv;
  
  double *Hi,*Yi,*xi,*yi,*xii;//xii for Hi calculation Hi= xi*xi'
  Hi = (double* ) R_Calloc(N*N, double);
  Yi =(double* ) R_Calloc(Mi*N, double);
  xi = (double* ) R_Calloc(N, double);
  xii = (double* ) R_Calloc(N, double);
  yi = (double* ) R_Calloc(N, double);
  double alpha, beta;
  char transa = 'N';
  char transb = 'N';
  
  //
  int MiMi = Mi*Mi;
  int NN = N*N;
  YiPi2Norm 	= (double* ) R_Calloc(MiMi, double);	
  bi 			= (double* ) R_Calloc(Mi, double);

  double *xiYi; //xi*Yi
  xiYi = (double* ) R_Calloc(Mi, double);
  double xiYibi, xiyi;

  alpha = 1;
  beta = 0;
  
  //largest Eigenvalue
  double *biInv;
  biInv 		= (double* ) R_Calloc(MiMi, double); //copy of YiPi2Norm
  char jobz = 'N'; // yes for eigenvectors
  char uplo = 'U'; //both ok
  double *w, *work;
  w = (double *) R_Calloc(Mi,double);
  int lwork = 5*Mi + 10;
  work  = (double *) R_Calloc(lwork,double);	
  int liwork = 10;
  int *iwork;
  iwork = (int *) R_Calloc(liwork,int);
  int info = 0;

  int *ipiv;
  ipiv = (int *) R_Calloc(Mi,int);
  double *readPtr,*readPtr2;
  //loop starts here
  for(i=0;i<M;i++)
  {
    readPtr = &X[i];
    F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
    F77_CALL(dcopy)(&N,xi,&inci,xii,&incj);
    readPtr = &Y[i];
    F77_CALL(dcopy)(&N,readPtr,&M,yi,&inci);

    xixi = F77_CALL(ddot)(&N, xi, &inci,xi, &incj);
    xixiInv = -1/(xixi + 1e-10);

    transb = 'N';
    lda = N;
    ldb = N;
    ldc = N;
    F77_CALL(dgemm)(&transa, &transb,&N, &ldb, &inci,&alpha, xi,&lda, xii, &incj, &beta,Hi, &ldc FCONE FCONE);

    F77_CALL(dscal)(&NN,&xixiInv,Hi,&inci); // Hi = -xi*xi'/(xi'*xi);
    for(j=0;j<N;j++) 
    {	index = j*N + j;
      Hi[index] = Hi[index] + 1;
    }

    readPtr2 = &Yi[0];
    jj = 0;
    for(j=0;j<M;j++)
    {	if(j!=i)
    {

      readPtr = &Y[j];
      F77_CALL(dcopy)(&N,readPtr,&M,readPtr2,&Mi);
      jj = jj + 1;
      readPtr2 = &Yi[jj];
    }
    }

    lda = Mi;
    ldb = N;
    ldc = Mi;
    ldk = N; //b copy
    F77_CALL(dgemm)(&transa, &transb,&Mi, &N, &ldk,&alpha, Yi, &lda, Hi, &ldb, &beta, YiPi, &ldc FCONE FCONE);

    transb = 'T';
    ldk = Mi;
    lda = Mi;
    ldb = Mi;
    ldc = Mi;
    F77_CALL(dgemm)(&transa, &transb,&Mi, &ldk, &N,&alpha, YiPi, &lda, Yi, &ldb, &beta, YiPi2Norm, &ldc FCONE FCONE);


    transb = 'N';
    F77_CALL(dcopy)(&MiMi,YiPi2Norm,&inci,biInv,&incj);

    lda = Mi;
    F77_CALL(dsyevd)(&jobz, &uplo,&Mi, biInv, &lda, w, work, &lwork, iwork, &liwork,&info FCONE FCONE);
    normYiPi = w[Mi -1]; //largestEigVal
	
    rho = rho_factor*normYiPi; // 2Norm = sqrt(lambda_Max)
    
    if(verbose>8) Rprintf("\t\t\t\t\t\t\t\t\t Gene number: %d,\t shrinkage rho: %f\n",i,rho);

    for(j=0;j<Mi;j++) 
    {
      index = j*Mi + j;
      YiPi2Norm[index] = YiPi2Norm[index] + rho;
    }

    lda = Mi;
    F77_CALL(dgemv)(&transa, &Mi, &N,&alpha, YiPi, &lda, yi, &inci, &beta,bi, &incj FCONE);

    lda = Mi;
    ldb = Mi;
    F77_CALL(dgesv)(&Mi, &inci, YiPi2Norm, &lda, ipiv, bi, &ldb, &info);
    lda = Mi;
    
    F77_CALL(dgemv)(&transa, &Mi, &N,&alpha, Yi, &lda, xi, &inci, &beta,xiYi, &incj FCONE);

    xiyi = F77_CALL(ddot)(&N, xi, &inci,yi, &incj);

    xiYibi = F77_CALL(ddot)(&Mi, xiYi, &inci,bi, &incj);
    
    f[i] = (xiyi-xiYibi)/xixi;

    jj = 0;
    for(j = 0;j<M;j++)
    {
      if(j!=i)
      {
        B[j*M+i] = bi[jj];
        jj = jj +1;
      }
    }
  }//i = 1:M
  
  double *ImB;
  k = M*M;
  ImB = (double* ) R_Calloc(k, double);
  F77_CALL(dcopy)(&k,B,&inci,ImB,&incj);
  xixiInv = -1;
  F77_CALL(dscal)(&k,&xixiInv,ImB,&inci);
  for(i=0;i<M;i++) 
  {
    index = i*M + i;
    ImB[index] = 1 + ImB[index];
  }
  

  //noise, sigma2learnt,mue;
  double * NOISE; 	//MxN
  NOISE =(double* ) R_Calloc(MN, double);
  transb = 'N';
  ldk = M;
  lda = M;
  ldb = M;
  ldc = M;
  F77_CALL(dgemm)(&transa, &transb,&M, &N, &ldk,&alpha, ImB, &lda, Y, &ldb, &beta, NOISE, &ldc FCONE FCONE);//(I-B)*Y - fX
  for(i=0;i<M;i++)
  {
    // row i of X
    readPtr2 = &X[i];
    readPtr = &NOISE[i];
    alpha = -f[i];
    F77_CALL(daxpy)(&N, &alpha,readPtr2, &ldk,readPtr, &M);
  }//row i = 1:M
  
  double noiseNorm, sigma2learnt;

  noiseNorm = F77_CALL(ddot)(&MN, NOISE, &inci,NOISE, &incj);
  sigma2learnt = noiseNorm/(MN -1);
  
  for(i=0;i<M;i++)
  {
    mue[i] = -f[i]*meanX[i];
  }
  beta = 1;
  ldk = M;
  lda = M;
  alpha = 1;
  F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, ImB, &lda, meanY, &inci, &beta,mue, &incj FCONE);
  
  
  if(verbose>7) Rprintf("\t\t\t\t\t\t\t\tExit function: Ridge Regression. sigma^2 is: %f.\n\n",sigma2learnt);
  
  R_Free(meanY);
  R_Free(meanX);
  R_Free(Y);
  R_Free(X);
  R_Free(YiPi);
  R_Free(YiPi2Norm);
  R_Free(bi);	
  R_Free(xiYi);
  R_Free(NOISE);
  //
  R_Free(Hi);
  R_Free(Yi);
  R_Free(xi);
  R_Free(yi);
  R_Free(xii);
  //
  R_Free(ImB);
  
  //
  R_Free(biInv);
  
  R_Free(w);
  R_Free(iwork);
  R_Free(work);
  
  R_Free(ipiv);
  return sigma2learnt;
  
}


double lambdaMax(double *Y,double *X,double * W,int M, int N)
{	
  double *dxx, *rxy, *DxxRxy,*readPtr1,*readPtr2;
  double lambda_max = 0;		
  dxx				= (double* ) R_Calloc(M, double);
  rxy				= (double* ) R_Calloc(M, double);
  DxxRxy			= (double* ) R_Calloc(M, double);
  int i,k,index,lda;
  int inci = 1;
  int incj = 1; 
  lda = M;
  for(i=0;i<M;i++)
  {
    readPtr1 	= &X[i]; //ith row
    readPtr2 	= &Y[i];

    dxx[i] = F77_CALL(ddot)(&N,readPtr1,&lda,readPtr1,&M);
    rxy[i] 		= F77_CALL(ddot)(&N,readPtr1,&lda,readPtr2,&M);
    DxxRxy[i] 	= rxy[i]/dxx[i];		
  }
  
  double * XDxxRxy;
  int MN = M*N;
  XDxxRxy = (double* ) R_Calloc(MN, double);

  F77_CALL(dcopy)(&MN,X,&inci,XDxxRxy,&incj);
  double alpha;	
  for(i=0;i<M;i++)
  {
    alpha  = -DxxRxy[i];
    readPtr1 = &XDxxRxy[i]; //ith row
    F77_CALL(dscal)(&N,&alpha, readPtr1,&M);//	(n, a, x, incx)
  }
  
  alpha  = 1.0;
  F77_CALL(daxpy)(&MN,&alpha,Y,&inci,XDxxRxy,&inci);
  double *YYXDR; //= Y*XDxxRxy'
  int MM = M*M;
  YYXDR = (double* ) R_Calloc(MM, double);	

  double beta;
  char transa = 'N';
  char transb = 'T';
  alpha = -1;
  beta = 0;
  F77_CALL(dgemm)(&transa, &transb,&M, &M, &N,&alpha, Y,&M, XDxxRxy, &M, &beta,YYXDR, &M FCONE FCONE); //M xK, K xN  --> MxN, N xM --> M <-M, N<-M, k<-N

  for(i=0;i<M;i++)
  {		
    for(k=0;k<M;k++)
    {
      index  = k*M + i;
      if(i==k)
      {
        YYXDR[index] = 0;
      }else
      {
        YYXDR[index] = YYXDR[index]/W[index];
      }
    }
  }

  index = F77_CALL(idamax)(&MM,YYXDR,&inci);
  lambda_max = fabs(YYXDR[index-1]);
  
  R_Free(dxx);
  R_Free(rxy);
  R_Free(DxxRxy);
  //R_Free(XX);
  R_Free(XDxxRxy);
  R_Free(YYXDR);
  
  return lambda_max;	
}

void QlambdaStart(double *Y,double *X, double *Q, double sigma2,int M, int N)
{	
  double *dxx, *rxy, *DxxRxy,*readPtr1,*readPtr2;
  
  dxx				= (double* ) R_Calloc(M, double);
  rxy				= (double* ) R_Calloc(M, double);
  DxxRxy			= (double* ) R_Calloc(M, double);
  int i,index,ldk,lda,ldb,ldc;
  int inci = 1;
  int incj = 1; 
  //double norm;
  lda = M;
  for(i=0;i<M;i++)
  {
    readPtr1 	= &X[i]; //ith row
    readPtr2 	= &Y[i];
    
    dxx[i] = F77_CALL(ddot)(&N,readPtr1,&lda,readPtr1,&M);
    rxy[i] 		= F77_CALL(ddot)(&N,readPtr1,&lda,readPtr2,&M);
    DxxRxy[i] 	= rxy[i]/dxx[i];		
  }
  double Nsigma2  = N*sigma2; 			// int * double --> double

  double * XDxxRxy;
  int MN = M*N;
  XDxxRxy = (double* ) R_Calloc(MN, double);
  F77_CALL(dcopy)(&MN,X,&inci,XDxxRxy,&incj);
  double alpha;	
  for(i=0;i<M;i++)
  {
    alpha  = -DxxRxy[i];
    readPtr1 = &XDxxRxy[i]; //ith row
    F77_CALL(dscal)(&N,&alpha, readPtr1,&M);//	(n, a, x, incx)
  }

  alpha  = 1.0;
  F77_CALL(daxpy)(&MN,&alpha,Y,&inci,XDxxRxy,&incj);
  
  double beta;
  char transa = 'N';
  char transb = 'T';
  alpha = -1;
  beta = 0;
  
  ldb = M;
  ldc = M;
  ldk = M;
  F77_CALL(dgemm)(&transa, &transb,&M, &lda, &N,&alpha, XDxxRxy,&ldb, Y, &ldc, &beta,Q, &ldk FCONE FCONE); //M xK, K xN  --> MxN, N xM --> M <-M, N<-M, k<-N	

  for(i=0;i<M;i++)
  {
    index = i*M + i;
    Q[index]= Q[index] + Nsigma2;
  }	
  
  R_Free(dxx);
  R_Free(rxy);
  R_Free(DxxRxy);
  R_Free(XDxxRxy);
  
  
}


void QlambdaMiddle(double *Y,double *X, double *Q,double *B,double *f, double *mue, double sigma2,int M, int N)
{	

  double *IB, *IBinv,*IBcopy;
  int MM = M*M;
  int MN = M*N;
  IB = (double* ) R_Calloc(MM, double);
  IBinv = (double* ) R_Calloc(MM, double);
  IBcopy = (double* ) R_Calloc(MM, double);
  int inci = 1;
  int incj = 1;
  F77_CALL(dcopy)(&MM,B,&inci,IB,&incj);	
  int i,index;
  double alpha;
  double beta = 0;
  alpha = -1;
  F77_CALL(dscal)(&MM,&alpha,IB,&inci);
  alpha = 0;
  int inc0 = 0;
  F77_CALL(dcopy)(&MM,&alpha,&inc0,IBinv,&inci);
  
  for(i=0;i<M;i++) 
  {
    index = i*M + i;
    IB[index] = 1 + IB[index];
    IBinv[index] = 1;
  }
  F77_CALL(dcopy)(&MM,IB,&inci,IBcopy,&incj);	
  
  int info = 0;
  int *ipiv;
  ipiv = (int *) R_Calloc(M,int);
  int lda = M;
  int ldb = M;
  int ldc = M;
  int ldk = M;
  F77_CALL(dgesv)(&M, &ldk, IBcopy, &lda, ipiv, IBinv, &ldb, &info);
  
  double Nsigma2  = N*sigma2; 			// int * double --> double
  double *Noise;
  Noise = (double* ) R_Calloc(MN, double);	

  char transa = 'N';
  char transb = 'N';
  alpha = 1;
  F77_CALL(dgemm)(&transa, &transb,&M, &N, &ldk,&alpha, IB, &lda, Y, &ldb, &beta, Noise, &ldc FCONE FCONE);
  double *readPtr1, *readPtr2;
  for(i=0;i<M;i++)
  {
    readPtr1 = &X[i];
    readPtr2 = &Noise[i];
    alpha = -f[i]; // y= alpha x + y
    F77_CALL(daxpy)(&N, &alpha,readPtr1, &lda,readPtr2, &M);
  }//row i = 1:M
  
  alpha = -1;
  for(i=0;i<N;i++)
  {
    readPtr1 = &Noise[i*M];
    F77_CALL(daxpy)(&M, &alpha,mue, &inci,readPtr1, &incj);
  }	

  transb = 'T';
  F77_CALL(dgemm)(&transa, &transb,&M, &ldk, &N,&alpha, Noise, &lda, Y, &ldb, &beta, Q, &ldc FCONE FCONE);

  alpha = Nsigma2;
  F77_CALL(daxpy)(&MM, &alpha,IBinv, &inci,Q, &incj);
  
  R_Free(IB);
  R_Free(IBinv);
  R_Free(IBcopy);
  R_Free(Noise);
  R_Free(ipiv);
  
}


void QlambdaMiddleCenter(double *Y,double *X, double *Q,double *B,double *f, double sigma2,int M, int N,
                         double *IBinv)
{	

  double *IB; 	//, *IBinv,*IBcopy
  int MM = M*M;
  int MN = M*N;
  IB = (double* ) R_Calloc(MM, double);

  int inci = 1;
  int incj = 1;

  F77_CALL(dcopy)(&MM,B,&inci,IB,&incj);	
  int i,index;
  double alpha;
  double beta = 0;
  alpha = -1;
  F77_CALL(dscal)(&MM,&alpha,IB,&inci);

  for(i=0;i<M;i++) 
  {
    index = i*M + i;
    IB[index] = 1 + IB[index];

  }

  int lda = M;
  int ldb = M;
  int ldc = M;
  int ldk = M;

  double Nsigma2  = N*sigma2; 			// int * double --> double
  double *Noise;
  Noise = (double* ) R_Calloc(MN, double);	

  char transa = 'N';
  char transb = 'N';
  alpha = 1;
  F77_CALL(dgemm)(&transa, &transb,&M, &N, &ldk,&alpha, IB, &lda, Y, &ldb, &beta, Noise, &ldc FCONE FCONE);
  double *readPtr1, *readPtr2;
  for(i=0;i<M;i++)
  {
    readPtr1 = &X[i];
    readPtr2 = &Noise[i];
    alpha = -f[i]; // y= alpha x + y
    F77_CALL(daxpy)(&N, &alpha,readPtr1, &lda,readPtr2, &M);
  }//row i = 1:M
  
  alpha = -1;
  transb = 'T';
  F77_CALL(dgemm)(&transa, &transb,&M, &ldk, &N,&alpha, Noise, &lda, Y, &ldb, &beta, Q, &ldc FCONE FCONE);

  alpha = Nsigma2;
  F77_CALL(daxpy)(&MM, &alpha,IBinv, &inci,Q, &incj);
  
  R_Free(IB);

  R_Free(Noise);
  
}

void UpdateIBinvPermute(double *QIBinv, double *B, int M)
{
  double *IB,*IBinv;	//, *IBinv,*IBcopy;
  int MM = M*M;
  int lda = M;
  int ldb = M;
  int ldk = M;
  IB = (double* ) R_Calloc(MM, double);
  IBinv = (double* ) R_Calloc(MM, double);
  int inci = 1;
  int incj = 1;
  int inc0 = 0;
  F77_CALL(dcopy)(&MM,B,&inci,IB,&incj);	
  int i,index;
  double alpha;
  //double beta = 0;
  alpha = -1;
  F77_CALL(dscal)(&MM,&alpha,IB,&inci);
  alpha = 0;
  F77_CALL(dcopy)(&MM,&alpha,&inc0,IBinv,&inci);
  for(i=0;i<M;i++) 
  {
    index = i*M + i;
    IB[index] = 1 + IB[index];
    IBinv[index] = 1;
  }
  
  int info = 0;
  int *ipiv;
  ipiv = (int *) R_Calloc(M,int);
  F77_CALL(dgesv)(&M, &ldk, IB, &lda, ipiv, IBinv, &ldb, &info);
  double *ptr1,*ptr2;
  
  for(i=0;i<M;i++)
  {
    index = ipiv[i] -1;
    ptr1 = &QIBinv[index*M];
    ptr2 = &IBinv[i*M];
    F77_CALL(dcopy)(&M,ptr2,&inci,ptr1,&incj);
    
  }
  
  R_Free(IB);
  R_Free(ipiv);
  R_Free(IBinv);
}

void UpdateIBinv(double *QIBinv, double *B, int M)
{
  double *IB;	//, *IBinv,*IBcopy;
  int MM = M*M;
  int lda = M;
  int ldb = M;
  int ldk = M;
  IB = (double* ) R_Calloc(MM, double);
  
  int inci = 1;
  int incj = 1;
  int inc0 = 0;
  F77_CALL(dcopy)(&MM,B,&inci,IB,&incj);	
  int i,index;
  double alpha;
  //double beta = 0;
  alpha = -1;
  F77_CALL(dscal)(&MM,&alpha,IB,&inci);
  alpha = 0;
  F77_CALL(dcopy)(&MM,&alpha,&inc0,QIBinv,&inci);
  for(i=0;i<M;i++) 
  {
    index = i*M + i;
    IB[index] = 1 + IB[index];
    QIBinv[index] = 1;
  }

  int info = 0;
  int *ipiv;
  ipiv = (int *) R_Calloc(M,int);
  F77_CALL(dgesv)(&M, &ldk, IB, &lda, ipiv, QIBinv, &ldb, &info);
  
  R_Free(IB);
  R_Free(ipiv);
}


double Weighted_LassoSf(double * W, double *B, double *f, double *Ycopy,double *Xcopy,
                        double *Q, double lambda_factor, double lambda_factor_prev, double sigma2, int max_iter,
                        int M, int N, int verbose,double *QIBinv,double lambda_max)			//double * mue,
{
  int i,j,index,ldM;
  ldM = M;//fixed
  double *meanY, *meanX;
  meanY = (double* ) R_Calloc(M, double);
  meanX = (double* ) R_Calloc(M, double);
  
  //copy Y, X; 
  double *Y, *X;
  int MN = M*N;
  int MM = M*M;
  Y = (double* ) R_Calloc(MN, double);
  X = (double* ) R_Calloc(MN, double);

  int inci,incj, inc0;
  inci	= 1;
  incj 	= 1;
  inc0 	= 0;
  F77_CALL(dcopy)(&MN,Ycopy,&inci,Y,&incj);
  F77_CALL(dcopy)(&MN,Xcopy,&inci,X,&incj);
  centerYX(Y,X, meanY, meanX,M, N);
  
  double lambda;//lambda_max,

  if(verbose>4) Rprintf("\t\t\t\tEnter Function: weighted_LassoSf. The maximum lambda is: %f\n\n",lambda_max);
  lambda 					= lambda_factor*lambda_max;
  
  //none zeros
  double alpha,beta;
  beta = 0;
  double deltaLambda;
  double *s, *S,*Wcopy;
  S = (double* ) R_Calloc(MM, double);
  s = (double* ) R_Calloc(M, double);
  Wcopy = (double* ) R_Calloc(MM, double);
  F77_CALL(dcopy)(&MM,W,&inci,Wcopy,&incj);
  
  deltaLambda 			= (2*lambda_factor - lambda_factor_prev)*lambda_max;	
  F77_CALL(dscal)(&MM,&deltaLambda,Wcopy,&inci); //wcopy = deltaLambda*W
  
  double *ei,toyZero;
  toyZero= 0;
  ei = (double* ) R_Calloc(M, double);
  F77_CALL(dcopy)(&M,&toyZero,&inc0,ei,&inci);

  double *readPtr,*readPtr2;
  for(i=0;i<M;i++)
  {
    for(j=0;j<M;j++)
    {
      //W[i,j]
      index = j*M  +i;
      if(fabs(Q[index])>= Wcopy[index] && i!= j)
      {
        S[index] 	= 1;
      }else
      {
        S[index] 	= 0;
        B[index] 	= 0;
      }	
    }
    readPtr = &S[i]; //S[i,];
    s[i] = F77_CALL(dasum)(&M, readPtr, &ldM);
  }
  char transa = 'N'; 

  double *f0,*F1;
  f0 	= (double* ) R_Calloc(M, double);
  F1 	= (double* ) R_Calloc(MM, double);
  
  double *y_j;
  y_j 	= (double* ) R_Calloc(N, double);
  double *F1ptr;
  
  
  double XYi, XXi;
  for(i=0;i<M;i++)
  {
    readPtr = &X[i];
    readPtr2 = &Y[i];

    XYi = F77_CALL(ddot)(&N, readPtr, &M,readPtr2, &M);

    XXi = F77_CALL(ddot)(&N, readPtr, &M,readPtr, &M);
    f0[i] 	= XYi/XXi;
    F1ptr	= &F1[M*i];//start from ith column
    alpha = 1/XXi;
    F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, readPtr, &M, &beta,F1ptr, &incj FCONE);
  }
	
  double *IBinv,*zi,*a_iT;		// y_j: one row of Y: Nx1
  IBinv 	= (double* ) R_Calloc(MM, double);
  a_iT 	= (double* ) R_Calloc(N, double);
  
  
  //loop starts here
  int iter = 0;
  double js_i, m_ij,B_old, lambdaW,beta_ij,r_ij, Bij;
  double *eiB;
  eiB = (double* ) R_Calloc(M, double);
  double *BiT;
  BiT = (double* ) R_Calloc(M, double);

  double d_ij, theta_ijp,k_ijp,q_ijp,Bijpp, Bijpm; //case (14)
  double q_ijm, theta_ijm, Bijmm, Bijmp,Lss,candsBij,LssCands;
  
  double dB,ziDb,BF1;

  double delta_BF,FnormOld, FnormChange;
  double *BfOld,*BfNew,*BfChange;
  index = M*(M  +1);
  BfOld = (double* ) R_Calloc(index, double);
  BfNew = (double* ) R_Calloc(index, double);
  BfChange = (double* ) R_Calloc(index, double);
  
  while(iter < max_iter)
  {
    iter = iter + 1;
    F77_CALL(dcopy)(&MM,B,&inci,BfOld,&incj);
    //last column
    F1ptr = &BfOld[MM];
    F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);

    for(i=0;i<M;i++)
    {
      if(s[i] >0)
      { 	//
        if(verbose>6) Rprintf("\t\t\t\t\t updating gene %d \n",i);
        //
        ei[i] = 1;

        zi = &QIBinv[i*M];
 
        for(j=0;j<M;j++) 
        {
          js_i = S[j*M + i]; 		//ith row
          if(js_i >0)
          {
            
            m_ij 	= zi[j];
            B_old 	= B[j*M + i]; //B[i,j]
            if(j!=i)
            {
              readPtr = &Y[j];
              F77_CALL(dcopy)(&N,readPtr,&M,y_j,&inci);
              
              lambdaW 	= lambda*W[j*M + i]; 	//W[i,j];
              readPtr = &B[i];
              
              F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);							
              alpha = -1;
              F77_CALL(dscal)(&M,&alpha,BiT,&inci);
              BiT[j] = 0;
              F77_CALL(dcopy)(&M,ei,&inci,eiB,&incj);
              alpha = 1;
              F77_CALL(daxpy)(&M, &alpha,BiT, &inci,eiB, &incj);
              readPtr = &X[i];
              F77_CALL(dcopy)(&N,readPtr,&M,a_iT,&inci);	

              alpha = -f[i];
              F77_CALL(dscal)(&N,&alpha,a_iT,&inci);							
              
              transa='T'; //y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
              beta = 1;
              alpha = 1;
              F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, eiB, &inci, &beta,a_iT, &incj FCONE);

              r_ij = F77_CALL(ddot)(&N, y_j, &inci,y_j, &incj);

              beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
              
              if (fabs(m_ij)<1e-10) //go to the linear equation 
              {
                if(verbose>7) Rprintf("\t\t\t\t\t\t\t gene %d \t interact with gene %d.\tLinear equation\n",i,j);
                //
                Bij = (beta_ij-lambdaW)/r_ij;
                if(Bij>0) 
                {
                  B[j*M+i] = Bij;//B(i,j)      = Bij;
                }else
                {
                  Bij         = (beta_ij+lambdaW)/r_ij;
                  if(Bij<0)
                  {
                    B[j*M+i] = Bij;//B(i,j)      = Bij;
                  }else
                  {
                    B[j*M+i] = 0;
                  }
                }//B_ij>0 
              }else //m_ij ~=0 go to the quadratic equation
              {
                //
                if(verbose>7) Rprintf("\t\t\t\t\t\t\t gene %d \t interact with gene %d.\tQuadratic equation\n",i,j);
                //
                //assume Bij >0
                d_ij = 1/m_ij + B[j*M+i];
                theta_ijp = r_ij*d_ij + beta_ij - lambdaW;
                k_ijp = d_ij*(beta_ij - lambdaW) - N*sigma2;
                
                q_ijp = theta_ijp*theta_ijp - 4*r_ij * k_ijp;
                Bijpp = (1/(2*r_ij))*(theta_ijp + sqrt(q_ijp));
                Bijpm = (1/(2*r_ij))*(theta_ijp - sqrt(q_ijp));
                
                //assume Bij<0
                q_ijm = q_ijp + 4*lambdaW *(beta_ij - r_ij *d_ij);
                theta_ijm = theta_ijp + 2*lambdaW;
                Bijmm = (1/(2*r_ij))*(theta_ijm - sqrt(q_ijm));
                Bijmp = (1/(2*r_ij))*(theta_ijm + sqrt(q_ijm));
                candsBij = 0;
 
                Lss = sigma2*N*log(fabs(d_ij)+1e-16);
						
                if (Bijpp>0)
                {
                  LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - r_ij*pow(Bijpp,2)/2 + beta_ij*Bijpp -lambdaW*fabs(Bijpp); 
                  
                  if(LssCands>Lss) 
                  {
                    candsBij = Bijpp;
                    Lss 	= LssCands;
                  }	
                }
                if (Bijpm>0)
                {
                  LssCands = sigma2*N*log(fabs(d_ij - Bijpm)+1e-16) - r_ij*pow(Bijpm,2)/2 + beta_ij*Bijpm -lambdaW*fabs(Bijpm); 
                  if(LssCands>Lss) 
                  {
                    candsBij = Bijpm;
                    Lss 	= LssCands;
                  }	
                }								
                //
                if (Bijmm<0)
                {
                  LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - r_ij*pow(Bijmm,2)/2 + beta_ij*Bijmm -lambdaW*fabs(Bijmm);  
                  if(LssCands>Lss) 
                  {
                    candsBij = Bijmm;
                    Lss 	= LssCands;
                  }	
                }
                if (Bijmp<0)
                {
                  LssCands = sigma2*N*log(fabs(d_ij - Bijmp)+1e-16) - r_ij*pow(Bijmp,2)/2 + beta_ij*Bijmp -lambdaW*fabs(Bijmp); 
                  if(LssCands>Lss) 
                  {
                    candsBij = Bijmp;
                    Lss 	= LssCands;
                  }	
                }
                B[j*M+i] = candsBij;
              }//m_ij
            }
            dB = B_old - B[j*M +i];
            ziDb = 1/(1 + dB*m_ij);
            F77_CALL(dscal)(&M,&ziDb,zi,&inci);
            
          }//js_i >0
        }//j = 1:M	

        readPtr = &B[i];
        F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);
        
        F1ptr = &F1[M*i];
        BF1 = F77_CALL(ddot)(&M, BiT, &inci,F1ptr, &incj);
        
        f[i] = f0[i] - BF1;
        ei[i] = 0; // re-set ei for next i
      }else
      {
        readPtr = &B[i];
        F77_CALL(dcopy)(&M,&toyZero,&inc0,readPtr,&ldM);
        f[i] = f0[i];
      } // s[i]
    }//i= 1:M
    
    F77_CALL(dcopy)(&MM,B,&inci,BfNew,&incj);
    F1ptr = &BfNew[MM];
    F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
    //convergence 		
    index = (M+1)*M;			//daxpy(n, a, x, incx, y, incy) 	y := a*x + y
    alpha = -1;
    F77_CALL(dcopy)(&index,BfOld,&inci,BfChange,&incj);
    F77_CALL(daxpy)(&index, &alpha,BfNew, &inci,BfChange, &incj);
    
    FnormOld = F77_CALL(dnrm2)(&index,BfOld,&inci);	
    FnormChange = F77_CALL(dnrm2)(&index,BfChange,&inci);	
    //
    delta_BF = FnormChange/(FnormOld + 1e-10);
    UpdateIBinv(QIBinv, B,M);	
    
    if(verbose>5) Rprintf("\t\t\t\t\t\tdelta_BF: %f\n",delta_BF);
    if(delta_BF<1e-3)		//break out
    {
      break;
    }
    
  }
  
  if(verbose>4) Rprintf("\t\t\t\tCurrent lambda: %f;\t number of iteration is: %d.\tExiting Weighted_LassoSf\n\n",lambda, iter);
  
  R_Free(meanY);
  R_Free(meanX);
  R_Free(Y);
  R_Free(X);
  
  R_Free(S);
  R_Free(s);
  R_Free(f0);
  R_Free(F1);
  R_Free(Wcopy);
  
  //R_Free(xi);
  R_Free(y_j);
  
  R_Free(ei);
  R_Free(IBinv);
  //R_Free(zi);
  R_Free(a_iT);
  
  R_Free(eiB);
  R_Free(BiT);
  R_Free(BfOld);
  R_Free(BfNew);
  R_Free(BfChange);
  
  return lambda;
}//weighted_LassoSf

double Weighted_LassoSf_MLf(double * W, double *BL, double *fL, double *Ycopy,double *Xcopy,
                            double *Q, double lambda_factor, double lambda_factor_prev, double sigma2, int max_iter,
                            int M, int N, int verbose,
                            double *BC, double *fC, double *mue,double *QIBinv,double *IBinvZero,double lambda_max)
{
  //SET TO PART1: LASSO
  double *B, *f;
  B = &BL[0];
  f = &fL[0];
  
  int i,j,index,ldk,ldM;
  ldM = M;//fixed
  double *meanY, *meanX;
  meanY = (double* ) R_Calloc(M, double);
  meanX = (double* ) R_Calloc(M, double);
  
  //copy Y, X; 
  double *Y, *X;
  int MN = M*N;
  int MM = M*M;
  Y = (double* ) R_Calloc(MN, double);
  X = (double* ) R_Calloc(MN, double);

  int inci,incj,inc0;
  inci	= 1;
  incj 	= 1;
  inc0 	= 0;
  F77_CALL(dcopy)(&MN,Ycopy,&inci,Y,&incj);
  F77_CALL(dcopy)(&MN,Xcopy,&inci,X,&incj);
  centerYX(Y,X, meanY, meanX,M, N);

  double lambda;//lambda_max,

  if(verbose>4) Rprintf("\t\t\t\tEnter Function: weighted_LassoSf. The maximum lambda is: %f\n\n",lambda_max);
  lambda 					= lambda_factor*lambda_max;
  
  //none zeros
  double alpha,beta;
  beta = 0;
  double deltaLambda;
  double *s, *S,*Wcopy;
  S = (double* ) R_Calloc(MM, double);
  s = (double* ) R_Calloc(M, double);
  Wcopy = (double* ) R_Calloc(MM, double);
  F77_CALL(dcopy)(&MM,W,&inci,Wcopy,&incj);
  
  deltaLambda 			= (2*lambda_factor - lambda_factor_prev)*lambda_max;	
  F77_CALL(dscal)(&MM,&deltaLambda,Wcopy,&inci); //wcopy = deltaLambda*W
  
  //ei = 0
  double *ei,toyZero;
  toyZero= 0;
  ei = (double* ) R_Calloc(M, double);
  F77_CALL(dcopy)(&M,&toyZero,&inc0,ei,&inci);

  double *readPtr,*readPtr2;
  for(i=0;i<M;i++)
  {
    for(j=0;j<M;j++)
    {
      //W[i,j]
      index = j*M  +i;
      if(fabs(Q[index])>= Wcopy[index] && i!= j)
      {
        S[index] 	= 1;
      }else
      {
        S[index] 	= 0;
        B[index] 	= 0;
      }	
    }
    readPtr = &S[i]; //S[i,];
    s[i] = F77_CALL(dasum)(&M, readPtr, &ldM);
  }
  char transa = 'N'; 

  double *f0,*F1;

  f0 	= (double* ) R_Calloc(M, double);
  F1 	= (double* ) R_Calloc(MM, double);
  
  //double *xi, *y_j;
  double *y_j;

  y_j 	= (double* ) R_Calloc(N, double);
  double *F1ptr;
  
  
  double XYi, XXi;
  for(i=0;i<M;i++)
  {
    readPtr = &X[i];
    readPtr2 = &Y[i];

    XYi = F77_CALL(ddot)(&N, readPtr, &M,readPtr2, &M);
    XXi = F77_CALL(ddot)(&N, readPtr, &M,readPtr, &M);
    f0[i] 	= XYi/XXi;
    F1ptr	= &F1[M*i];//start from ith column
    alpha = 1/XXi;
    F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, readPtr, &M, &beta,F1ptr, &incj FCONE);
  }

  double *IBinv,*zi,*a_iT;// y_j: one row of Y: Nx1
  IBinv 	= (double* ) R_Calloc(MM, double);
  a_iT 	= (double* ) R_Calloc(N, double);

  //loop starts here
  int iter = 0;
  double js_i, m_ij,B_old, lambdaW,beta_ij,r_ij, Bij;
  //dynamic variable keep intermidiate values 
  double *eiB;
  eiB = (double* ) R_Calloc(M, double);
  double *BiT;
  BiT = (double* ) R_Calloc(M, double);
  //quadratic function
  double d_ij, theta_ijp,k_ijp,q_ijp,Bijpp, Bijpm; //case (14)
  double q_ijm, theta_ijm, Bijmm, Bijmp,Lss,candsBij,LssCands;
  
  //converge of gene i
  double dB,ziDb,BF1;
  
  //converge of while
  double delta_BF,FnormOld, FnormChange;
  double *BfOld,*BfNew,*BfChange;
  index = M*(M  +1);
  BfOld = (double* ) R_Calloc(index, double);
  BfNew = (double* ) R_Calloc(index, double);
  BfChange = (double* ) R_Calloc(index, double);
  
  while(iter < max_iter)
  {
    iter = iter + 1;
    F77_CALL(dcopy)(&MM,B,&inci,BfOld,&incj);
    //last column
    F1ptr = &BfOld[MM];
    F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);

    for(i=0;i<M;i++)
    {
      if(s[i] >0)
      { 	//
        if(verbose>6) Rprintf("\t\t\t\t\t updating gene %d \n",i);
        ei[i] = 1;
        zi = &QIBinv[i*M];
        for(j=0;j<M;j++) 
        {
          js_i = S[j*M + i]; 		//ith row
          if(js_i >0)
          {
            
            m_ij 	= zi[j];
            B_old 	= B[j*M + i]; //B[i,j]
            if(j!=i)
            {
              readPtr = &Y[j];
              F77_CALL(dcopy)(&N,readPtr,&M,y_j,&inci);
              //Y[j,:]
              
              lambdaW 	= lambda*W[j*M + i]; 	//W[i,j];
              //BiT = -B[i:]
              readPtr = &B[i];
              
              F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);							
              alpha = -1;
              F77_CALL(dscal)(&M,&alpha,BiT,&inci);
              BiT[j] = 0;
              //eiB
              F77_CALL(dcopy)(&M,ei,&inci,eiB,&incj);
              alpha = 1;
              F77_CALL(daxpy)(&M, &alpha,BiT, &inci,eiB, &incj);
              readPtr = &X[i];
              F77_CALL(dcopy)(&N,readPtr,&M,a_iT,&inci);	
              alpha = -f[i];
              F77_CALL(dscal)(&N,&alpha,a_iT,&inci);							
              
              transa='T'; //y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
              beta = 1;
              alpha = 1;
              F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, eiB, &inci, &beta,a_iT, &incj FCONE);

              r_ij = F77_CALL(ddot)(&N, y_j, &inci,y_j, &incj);

              beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
              
              if (fabs(m_ij)<1e-10) //go to the linear equation 
              {
                //
                if(verbose>7) Rprintf("\t\t\t\t\t\t\t gene %d \t interact with gene %d.\tLinear equation\n",i,j);
                //
                Bij = (beta_ij-lambdaW)/r_ij;
                if(Bij>0) 
                {
                  B[j*M+i] = Bij;//B(i,j)      = Bij;
                }else
                {
                  Bij         = (beta_ij+lambdaW)/r_ij;
                  if(Bij<0)
                  {
                    B[j*M+i] = Bij;//B(i,j)      = Bij;
                  }else
                  {
                    B[j*M+i] = 0;
                  }
                }//B_ij>0 
              }else //m_ij ~=0 go to the quadratic equation
              {
                //
                if(verbose>7) Rprintf("\t\t\t\t\t\t\t gene %d \t interact with gene %d.\tQuadratic equation\n",i,j);

                d_ij = 1/m_ij + B[j*M+i];
                theta_ijp = r_ij*d_ij + beta_ij - lambdaW;
                k_ijp = d_ij*(beta_ij - lambdaW) - N*sigma2;
                
                q_ijp = theta_ijp*theta_ijp - 4*r_ij * k_ijp;
                Bijpp = (1/(2*r_ij))*(theta_ijp + sqrt(q_ijp));
                Bijpm = (1/(2*r_ij))*(theta_ijp - sqrt(q_ijp));
                
                //assume Bij<0
                q_ijm = q_ijp + 4*lambdaW *(beta_ij - r_ij *d_ij);
                theta_ijm = theta_ijp + 2*lambdaW;
                Bijmm = (1/(2*r_ij))*(theta_ijm - sqrt(q_ijm));
                Bijmp = (1/(2*r_ij))*(theta_ijm + sqrt(q_ijm));
                candsBij = 0;
 
                Lss = sigma2*N*log(fabs(d_ij)+1e-16);
							
                if (Bijpp>0)
                {
                  LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - r_ij*pow(Bijpp,2)/2 + beta_ij*Bijpp -lambdaW*fabs(Bijpp); 
                  
                  if(LssCands>Lss) 
                  {
                    candsBij = Bijpp;
                    Lss 	= LssCands;
                  }	
                }
                if (Bijpm>0)
                {
                  LssCands = sigma2*N*log(fabs(d_ij - Bijpm)+1e-16) - r_ij*pow(Bijpm,2)/2 + beta_ij*Bijpm -lambdaW*fabs(Bijpm); 
                  if(LssCands>Lss) 
                  {
                    candsBij = Bijpm;
                    Lss 	= LssCands;
                  }	
                }								
                //
                if (Bijmm<0)
                {
                  LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - r_ij*pow(Bijmm,2)/2 + beta_ij*Bijmm -lambdaW*fabs(Bijmm);  
                  if(LssCands>Lss) 
                  {
                    candsBij = Bijmm;
                    Lss 	= LssCands;
                  }	
                }
                if (Bijmp<0)
                {
                  LssCands = sigma2*N*log(fabs(d_ij - Bijmp)+1e-16) - r_ij*pow(Bijmp,2)/2 + beta_ij*Bijmp -lambdaW*fabs(Bijmp); 
                  if(LssCands>Lss) 
                  {
                    candsBij = Bijmp;
                    Lss 	= LssCands;
                  }	
                }
                B[j*M+i] = candsBij;
              }//m_ij
            }//if(j!= i)
            dB = B_old - B[j*M +i];
            //update c_ij
            ziDb = 1/(1 + dB*m_ij);
            F77_CALL(dscal)(&M,&ziDb,zi,&inci);
          }//js_i >0
        }//j = 1:M	
        //f
        readPtr = &B[i];
        F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);
        
        F1ptr = &F1[M*i];
        BF1 = F77_CALL(ddot)(&M, BiT, &inci,F1ptr, &incj);
        
        f[i] = f0[i] - BF1;
        ei[i] = 0; // re-set ei for next i
      }else
      {
        readPtr = &B[i];
        F77_CALL(dcopy)(&M,&toyZero,&inc0,readPtr,&ldM);
        f[i] = f0[i];
      } // s[i]
    }//i= 1:M
    
    F77_CALL(dcopy)(&MM,B,&inci,BfNew,&incj);
    F1ptr = &BfNew[MM];
    F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
    //convergence 		
    index = (M+1)*M;			//daxpy(n, a, x, incx, y, incy) 	y := a*x + y
    alpha = -1;
    F77_CALL(dcopy)(&index,BfOld,&inci,BfChange,&incj);
    F77_CALL(daxpy)(&index, &alpha,BfNew, &inci,BfChange, &incj);
    
    FnormOld = F77_CALL(dnrm2)(&index,BfOld,&inci);	
    FnormChange = F77_CALL(dnrm2)(&index,BfChange,&inci);	
    //
    delta_BF = FnormChange/(FnormOld + 1e-10);
    UpdateIBinv(QIBinv, B,M);
    
    if(verbose>5) Rprintf("\t\t\t\t\t\tdelta_BF: %f\n",delta_BF);
    if(delta_BF<1e-3)		//break out
    {
      break;
    }
    
  }
  
  if(verbose>3) Rprintf("\t\t\t\tCurrent lambda: %f;\t number of iteration is: %d.\tExiting Weighted_LassoSf\n\n",lambda, iter);
  
  F77_CALL(dcopy)(&MM,BL,&inci,BC,&incj);
  F77_CALL(dcopy)(&M,fL,&inci,fC,&incj);
  B = &BC[0];
  f = &fC[0];
  
  if(verbose>4) Rprintf("Enter Function: constrained-MLf. Shrinkage lambda is: 0. \n");
  // SET SL
  for(i=0;i<MM;i++)
  {
    if(BL[i]==0)
    {
      S[i] = 0;
    }else
    {
      S[i] = 1;
    }
  }	

  for(i=0;i<M;i++)
  {
    readPtr = &S[i]; //S[i,];
    s[i] = F77_CALL(dasum)(&M, readPtr, &ldM);
  }
  
  beta = 1;
  iter = 0;
  double theta_ij,k_ij,q_ij,Bijp, Bijm; //case (14)
  max_iter = max_iter/5;

  while(iter < max_iter)
  {
    iter = iter + 1;
    //converge Bfold = [B f];
    F77_CALL(dcopy)(&MM,B,&inci,BfOld,&incj);
    //last column
    F1ptr = &BfOld[MM];
    F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);

    for(i=0;i<M;i++)
    {
      if(s[i] >0)
      {
        //
        if(verbose>6) Rprintf("\t\t updating gene %d \n",i);
        //
        ei[i] = 1;
        zi = &IBinvZero[i*M];

        for(j=0;j<M;j++)
        {
          js_i = S[j*M + i]; 		//ith row
          if(js_i >0)
          {
            m_ij 	= zi[j];
            B_old 	= B[j*M + i]; //B[i,j]
            
            //y_j
            readPtr = &Y[j];
            F77_CALL(dcopy)(&N,readPtr,&M,y_j,&inci);
            //BiT = -B[i:]
            readPtr = &B[i];
            F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);							
            alpha = -1;
            F77_CALL(dscal)(&M,&alpha,BiT,&inci);
            BiT[j] = 0;
            F77_CALL(dcopy)(&M,ei,&inci,eiB,&incj);						
            alpha = 1;
            F77_CALL(daxpy)(&M, &alpha,BiT, &inci,eiB, &incj);
            readPtr = &X[i];
            F77_CALL(dcopy)(&N,readPtr,&M,a_iT,&inci);
            alpha = -f[i];
            F77_CALL(dscal)(&N,&alpha,a_iT,&inci);							
            
            transa='T'; //y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            alpha = 1;
            F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, eiB, &inci, &beta,a_iT, &incj FCONE);

            r_ij = F77_CALL(ddot)(&N, y_j, &inci,y_j, &inci);

            beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
            
            if (fabs(m_ij)<1e-10) //go to the linear equation 
            {
              //
              if(verbose>7) Rprintf("\t\t\t gene %d \t interact with gene %d.\tLinear equation\n",i,j);
              B[j*M+i] = beta_ij/r_ij;
            }else //m_ij ~=0 go to the quadratic equation
            {
              //
              if(verbose>7) Rprintf("\t\t\t gene %d \t interact with gene %d.\tQuadratic equation\n",i,j);

              d_ij = 1/m_ij + B[j*M+i];
              theta_ij = r_ij*d_ij + beta_ij;
              k_ij = d_ij*beta_ij - N*sigma2;
              
              q_ij = theta_ij*theta_ij - 4*r_ij* k_ij;
              Bijp = (1/(2*r_ij))*(theta_ij + sqrt(q_ij));
              Bijm = (1/(2*r_ij))*(theta_ij - sqrt(q_ij));
              
              candsBij = 0;
              Lss = sigma2*N*log(fabs(d_ij)+1e-16);
              LssCands = sigma2*N*log(fabs(d_ij - Bijp)+1e-16) - r_ij*pow(Bijp,2)/2 + beta_ij*Bijp;
              if(LssCands>Lss) 
              {
                candsBij = Bijp;
                Lss 	= LssCands;
              }	
							
              LssCands = sigma2*N*log(fabs(d_ij - Bijm)+1e-16) - r_ij*pow(Bijm,2)/2 + beta_ij*Bijm;
              if(LssCands>Lss) 
              {
                candsBij = Bijm;
                Lss 	= LssCands;
              }	
              
              B[j*M+i] = candsBij;
            }//m_ij
            
            dB = B_old - B[j*M +i];
            //update c_ij
            ziDb = 1/(1 + dB*m_ij);
            F77_CALL(dscal)(&M,&ziDb,zi,&inci);	
            

            
          }//js_i >0
        }//j = 1:M	

        readPtr = &B[i];
        F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);
        F1ptr = &F1[M*i];
        BF1 = F77_CALL(ddot)(&M, BiT, &inci,F1ptr, &incj);
        
        f[i] = f0[i] - BF1;
        ei[i] = 0; // re-set ei for next i
        
      }//[si]
      
      
    }//i= 1:M
    
    //convergence 
    F77_CALL(dcopy)(&MM,B,&inci,BfNew,&incj);
    F1ptr = &BfNew[MM];
    F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
    index = (M+1)*M;			//daxpy(n, a, x, incx, y, incy) 	y := a*x + y
    alpha = -1;
    F77_CALL(dcopy)(&index,BfOld,&inci,BfChange,&incj);
    F77_CALL(daxpy)(&index, &alpha,BfNew, &inci,BfChange, &incj);
    
    FnormOld = F77_CALL(dnrm2)(&index,BfOld,&inci);	
    FnormChange = F77_CALL(dnrm2)(&index,BfChange,&inci);	
    
    delta_BF = FnormChange/(FnormOld + 1e-10);
    if(verbose>5) Rprintf("\t\tdelta_BF: %f\n",delta_BF);

    UpdateIBinv(IBinvZero, B,M);
    
    if(delta_BF<1e-2)		//break out
    {
      break;
    }
    
  }//while
  
  
  if(verbose>3) Rprintf("\t number of iteration is: %d.\nExiting constrained_MLf\n",iter);

  F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
  alpha = -1; 
  F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
  //diagonal + 1
  for(j=0;j<M;j++) 
  {
    index = j*M + j;
    IBinv[index] = 1 + IBinv[index];
    mue[j] = -f[j]*meanX[j];
  }

  transa = 'N';
  alpha = 1;
  //beta = 1;
  ldk = M;
  F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, IBinv, &ldM, meanY, &inci, &beta,mue, &incj FCONE);

  
  R_Free(meanY);
  R_Free(meanX);
  R_Free(Y);
  R_Free(X);
  
  R_Free(S);
  R_Free(s);
  R_Free(f0);
  R_Free(F1);
  R_Free(Wcopy);

  R_Free(y_j);
  
  R_Free(ei);
  R_Free(IBinv);
  R_Free(a_iT);
  
  R_Free(eiB);
  R_Free(BiT);
  R_Free(BfOld);
  R_Free(BfNew);
  R_Free(BfChange);
  

  return lambda;

}//weighted_LassoSf









