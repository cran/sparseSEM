#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <stdlib.h>

//version Note: block coordinate ascent: after each block updated: update IBinv before next iteration
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

void centerYX(double *Y,double *X, double *meanY, double *meanX,int M, int N) //M genes; N samples
{
	//matrix is vectorized by column: 1-M is the first column of Y
	//missing values are from X; set corresponding Y to zero.  in main_SMLX.m

	int i,index;	
	double *Xptr;
	double *Yptr;

	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	int lda  = M; //leading dimension
	double *eye;
	eye = (double* ) Calloc(N, double);
	double alpha = 1;
	double beta = 0;
	F77_CALL(dcopy)(&N,&alpha,&inc0,eye,&inci);
	char transa = 'N';

	F77_CALL(dgemv)(&transa, &M, &N,&alpha, X, &lda, eye, &inci, &beta,meanX, &incj);
	F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &lda, eye, &inci, &beta,meanY, &incj);
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
	Free(eye);
}	

//--------------------- LINEAR SYSTEM SOLVER END ---------

//ridge regression; return sigma2learnt
double constrained_ridge_cff(double *Ycopy, double *Xcopy, double rho_factor, int M, int N,
		double *B, double *f, double *mue, int verbose)
{
	
	int i,j,k,lda,ldb,ldc,ldk;
	// center Y, X
	double *meanY, *meanX;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	
	//copy Y, X; 
	double *Y, *X;
	int MN = M*N;
	Y = (double* ) Calloc(MN, double);
	X = (double* ) Calloc(MN, double);
	
	//F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
		//double *dy, const int *incy);
	int inci = 1;
	int incj = 1;
	F77_CALL(dcopy)(&MN,Ycopy,&inci,Y,&incj);
	F77_CALL(dcopy)(&MN,Xcopy,&inci,X,&incj);
	
	centerYX(Y,X,meanY, meanX,M, N);
	
//	double NPresent = 0;
//	for(i=0;i<M;i++)
//	{
//		for(j=0;j<N;j++)
//		{
//			if(Y[j*M + i]!=0) NPresent = NPresent + 1; //Y[i,j]
//		}
//	}	
	if(verbose>7) Rprintf("\t\t\t\t\t\t\t\tEnter Function: Ridge Regression. Shrinkage ratio rho is: %f.\n\n",rho_factor);
	
	int Mi = M -1;
	//for usage in loop
	double *YiPi; //Yi'*Pi
	YiPi =(double* ) Calloc(Mi*N, double);
	double xixi,xixiInv; //xi'*xi;
	int jj,index; //jj = 1:(M-1) index of YiPi
	double normYiPi,rho;
	double *bi,*YiPi2Norm; 	//YiPi2Norm: first term of biInv;
	
	double *Hi,*Yi,*xi,*yi,*xii;//xii for Hi calculation Hi= xi*xi'
	Hi = (double* ) Calloc(N*N, double);
	Yi =(double* ) Calloc(Mi*N, double);
	xi = (double* ) Calloc(N, double);
	xii = (double* ) Calloc(N, double);
	yi = (double* ) Calloc(N, double);
	double alpha, beta;
	char transa = 'N';
	char transb = 'N';
	
	//
	int MiMi = Mi*Mi;
	int NN = N*N;
	YiPi2Norm 	= (double* ) Calloc(MiMi, double);	
	bi 			= (double* ) Calloc(Mi, double);
	//YiPiyi 		= (double* ) Calloc(Mi, double);
	
	//bi,fi
	double *xiYi; //xi*Yi
	xiYi = (double* ) Calloc(Mi, double);
	double xiYibi, xiyi;
	//main loop:
//Rprintf("check point 1: before loop\n");
	alpha = 1;
	beta = 0;
		
//largest Eigenvalue
	double *biInv;
	biInv 		= (double* ) Calloc(MiMi, double); //copy of YiPi2Norm
	//dsyevd
	char jobz = 'N'; // yes for eigenvectors
	char uplo = 'U'; //both ok
	double *w, *work;
	w = (double *) Calloc(Mi,double);
	int lwork = 5*Mi + 10;
	work  = (double *) Calloc(lwork,double);	
	int liwork = 10;
	int *iwork;
	iwork = (int *) Calloc(liwork,int);
	int info = 0;
	//dsyevd function 

	//linear system
	int *ipiv;
	//ipiv = (int *) R_alloc(N,sizeof(int));
	ipiv = (int *) Calloc(Mi,int);
	double *readPtr,*readPtr2;
	//loop starts here
	for(i=0;i<M;i++)
	{
		//xi = X[i,:]
		readPtr = &X[i];
		F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
		F77_CALL(dcopy)(&N,xi,&inci,xii,&incj);
		readPtr = &Y[i];
		F77_CALL(dcopy)(&N,readPtr,&M,yi,&inci);

		//xixi = F77_CALL(dnrm2)(&N,xi,&inci);
		//xixi = pow(xixi,2);
		xixi = F77_CALL(ddot)(&N, xi, &inci,xi, &incj);
		xixiInv = -1/xixi;
//printMat(xi,1,N);		
		//xi'*xi
		
		//YiPi
		//Hi          = xi*xi'/(xi'*xi);
        //Pi          = eye(N)-Hi;
//Rprintf("check point 2: xixi: %f.\n",xixi);

		//MatrixMult(xi,xi, Hi,alpha, beta, N, k, N);
		transb = 'N';
		lda = N;
		ldb = N;
		ldc = N;
		F77_CALL(dgemm)(&transa, &transb,&N, &ldb, &inci,&alpha, xi,&lda, xii, &incj, &beta,Hi, &ldc);
		
		
		//F77_NAME(dscal)(const int *n, const double *alpha, double *dx, const int *incx);
		//k= N*N;
		F77_CALL(dscal)(&NN,&xixiInv,Hi,&inci); // Hi = -xi*xi'/(xi'*xi);
		for(j=0;j<N;j++) 
		{	index = j*N + j;
			Hi[index] = Hi[index] + 1;
		}//Pi
//printMat(Hi,N,N);	
		
	
		//Yi
		readPtr2 = &Yi[0];
		jj = 0;
		for(j=0;j<M;j++)
		{	if(j!=i)
			{
				//copy one j row
				readPtr = &Y[j];
				F77_CALL(dcopy)(&N,readPtr,&M,readPtr2,&Mi);
				jj = jj + 1;
				readPtr2 = &Yi[jj];
			}
		}//jj 1:(M-1), YiPi[jj,:]
		//YiPi=Yi*Pi
//printMat(Yi,Mi,N);		
		//transb = 'N';
		lda = Mi;
		ldb = N;
		ldc = Mi;
		ldk = N; //b copy
		F77_CALL(dgemm)(&transa, &transb,&Mi, &N, &ldk,&alpha, Yi, &lda, Hi, &ldb, &beta, YiPi, &ldc);
		//F77_CALL(dgemm)(&transa, &transb,&Mi, &N, &N,&alpha, Yi, &Mi, Hi, &N, &beta, YiPi, &Mi);
		
		//Rprintf("check point 3: YiPi\n");
		
		//YiPi*Yi' --> MixMi
		transb = 'T';
		ldk = Mi;
		lda = Mi;
		ldb = Mi;
		ldc = Mi;
		F77_CALL(dgemm)(&transa, &transb,&Mi, &ldk, &N,&alpha, YiPi, &lda, Yi, &ldb, &beta, YiPi2Norm, &ldc);
		//F77_CALL(dgemm)(&transa, &transb,&Mi, &Mi, &N,&alpha, YiPi, &Mi, Yi, &Mi, &beta, YiPi2Norm, &Mi); //M-> N -> K
		//2-norm YiPi this is the largest eigenvalue of (Yi'*Pi)*YiPi
		//Matrix 2-norm;
		//Repeat compution; use biInv;
		//normYiPi = Mat2NormSq(YiPi,Mi,N); //MixN
		//YiPi2Norm
		
		//Rprintf("check point 4: YiPi2Norm\n");
		
//normYiPi = largestEigVal(YiPi2Norm, Mi,verbose);  //YiPi2Norm make a copy inside largestEigVal, YiPi2Norm will be intermediate value of biInv;
		//transa = 'N';
		transb = 'N';
		//j = Mi*Mi;
		F77_CALL(dcopy)(&MiMi,YiPi2Norm,&inci,biInv,&incj);
		
		//F77_CALL(dgeev)(&transa, &transb,&Mi, biInv, &Mi, wr, wi, vl, &ldvl,vr, &ldvr, work, &lwork, &info);
		lda = Mi;
		F77_CALL(dsyevd)(&jobz, &uplo,&Mi, biInv, &lda, w, work, &lwork, iwork, &liwork,&info);
		normYiPi = w[Mi -1];
//printMat(w,1,Mi);		
//Rprintf("Eigenvalue: %f.\n",normYiPi);		
		rho = rho_factor*normYiPi; // 2Norm = sqrt(lambda_Max)
		
		if(verbose>8) Rprintf("\t\t\t\t\t\t\t\t\t Gene number: %d,\t shrinkage rho: %f\n",i,rho);
		//biInv = (YiPi*Yi+rho*I) ; (M-1) x (M-1)

		//
		for(j=0;j<Mi;j++) 
		{
			index = j*Mi + j;
			YiPi2Norm[index] = YiPi2Norm[index] + rho;
		}
		//biInv;

		//Inverse
		//MatrixInverse(biInv,Mi);
		//Rprintf("check point 5: biInv inversed\n");
		//YiPiyi = Yi'*Pi*yi  = YiPi *yi;
		//F77_NAME(dgemv)(const char *trans, const int *m, const int *n,
		//const double *alpha, const double *a, const int *lda,
		//const double *x, const int *incx, const double *beta,
		//double *y, const int *incy);
		lda = Mi;
		F77_CALL(dgemv)(&transa, &Mi, &N,&alpha, YiPi, &lda, yi, &inci, &beta,bi, &incj);
		

//linearSystem(YiPi2Norm,Mi,bi);//A NxN matrix
		lda = Mi;
		ldb = Mi;
		F77_CALL(dgesv)(&Mi, &inci, YiPi2Norm, &lda, ipiv, bi, &ldb, &info);
		//Rprintf("check point 6: bi updated\n");
		//------------------------------------------------Ridge coefficient beta obtained for row i
		// f(i)        = (xi'*yi-xi'*Yi*bi)/(xi'*xi);
		//xiYi (M-1) x1
		lda = Mi;
		
		F77_CALL(dgemv)(&transa, &Mi, &N,&alpha, Yi, &lda, xi, &inci, &beta,xiYi, &incj);
		
		//xiyi = xi*yi 	= X[i,j]*Y[i,j]
		//dot product 
		xiyi = F77_CALL(ddot)(&N, xi, &inci,yi, &incj);
		
		//xiYibi = xiYi*bi
		xiYibi = F77_CALL(ddot)(&Mi, xiYi, &inci,bi, &incj);

		f[i] = (xiyi-xiYibi)/xixi;
		//Rprintf("check point 7: fi calculated\n");
		//update B
		jj = 0;
		for(j = 0;j<M;j++)
		{
			if(j!=i)
			{
				//B[i,j] = bi[jj];
				B[j*M+i] = bi[jj];
				jj = jj +1;
			}
		}
		
		//end of 1:M		
	}//i = 1:M

	//I -B
	double *ImB;
	k = M*M;
	ImB = (double* ) Calloc(k, double);
	F77_CALL(dcopy)(&k,B,&inci,ImB,&incj);
	//F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
	//	double *dy, const int *incy);
	xixiInv = -1;
	F77_CALL(dscal)(&k,&xixiInv,ImB,&inci);
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		ImB[index] = 1 + ImB[index];
	}
	
	
	
	
	//noise, sigma2learnt,mue;
	double * NOISE; 	//MxN
	NOISE =(double* ) Calloc(MN, double);
	transb = 'N';
	ldk = M;
	lda = M;
	ldb = M;
	ldc = M;
	F77_CALL(dgemm)(&transa, &transb,&M, &N, &ldk,&alpha, ImB, &lda, Y, &ldb, &beta, NOISE, &ldc);//(I-B)*Y - fX
	for(i=0;i<M;i++)
	{
		// row i of X
		readPtr2 = &X[i];
		//F77_CALL(dcopy)(&N,readPtr2,&M,xi,&inci);
		//row i of noise
		readPtr = &NOISE[i];
		alpha = -f[i];
		//F77_CALL(daxpy)(&N, &alpha,xi, &inci,readPtr, &M);
		F77_CALL(daxpy)(&N, &alpha,readPtr2, &ldk,readPtr, &M);
		//NOISE[i,1:N]	
	}//row i = 1:M
	
	
	//NoiseMatF77(NOISE,ImB,Y, f, X, M, N);
		
	double noiseNorm, sigma2learnt;
	//noiseNorm = FrobeniusNorm(NOISE, M, N); //NOISE  = sparse(speye(M)-B)*Y-bsxfun(@times,f,X);
	//noiseNorm = F77_CALL(dnrm2)(&MN,NOISE,&inci);
	//sigma2learnt = noiseNorm*noiseNorm/(MN -1); //sigma2learnt    = sum(sum(NOISE.^2))/(sum(NPresent)-1);
	noiseNorm = F77_CALL(ddot)(&MN, NOISE, &inci,NOISE, &incj);
	sigma2learnt = noiseNorm/(MN -1);
	
	//mue          = (IM-B)*meanY-bsxfun(@times,f,meanX);
	//dgemv mue = Ax + beta*mue
	

	
	
	
	//mue[i] = -f[i]*meanX[i];
	for(i=0;i<M;i++)
	{
		mue[i] = -f[i]*meanX[i];
	}
	beta = 1;
	ldk = M;
	lda = M;
	alpha = 1;
	F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, ImB, &lda, meanY, &inci, &beta,mue, &incj);
	
	
	if(verbose>7) Rprintf("\t\t\t\t\t\t\t\tExit function: Ridge Regression. sigma^2 is: %f.\n\n",sigma2learnt);
	
	Free(meanY);
	Free(meanX);
	Free(Y);
	Free(X);
	Free(YiPi);
	//Free(biInv);
	Free(YiPi2Norm);
	Free(bi);	
	//Free(YiPiyi);
	Free(xiYi);
	Free(NOISE);
	//
	Free(Hi);
	Free(Yi);
	Free(xi);
	Free(yi);
	//
	Free(ImB);
	
	//
	Free(biInv);

	Free(w);
	Free(iwork);
	Free(work);
	
	Free(ipiv);
	return sigma2learnt;

}


//by Weighted_LassoSf xi
//lambda_max          = max(max(abs(N*sigma2*IM - (Y*(Y'-X'*DxxRxy)))./W  ));
double lambdaMax(double *Y,double *X,double * W,int M, int N)
{	
	// Oct 08, 2012: assume one eQTL for each gene; This fucntion needs significant revision if this assumption doesnot hold
	double *dxx, *rxy, *DxxRxy,*readPtr1,*readPtr2;
	double lambda_max = 0;		
	dxx				= (double* ) Calloc(M, double);
	rxy				= (double* ) Calloc(M, double);
	DxxRxy			= (double* ) Calloc(M, double);
	int i,k,index,lda;
	int inci = 1;
	int incj = 1; 
	lda = M;
	//for(i=0;i<M;i++)
	//{
	//	dxx[i] 		= 0;
	//	rxy[i] 		= 0;
	//	for(j=0;j<N;j++)
	//	{
	//		index = j*M+i;
	//		dxx[i] 	= dxx[i] + pow(X[index],2);	// sum of each row X[i,j]
	//		rxy[i] 	= rxy[i] + X[index]*Y[index];			
	//	}
	//	DxxRxy[i] 	= rxy[i]/dxx[i];		
	//}
	for(i=0;i<M;i++)
	{
		readPtr1 	= &X[i]; //ith row
		readPtr2 	= &Y[i];

		//norm  		= F77_CALL(dnrm2)(&N,readPtr1,&M);	
		//dxx[i] 		= pow(norm,2);
		dxx[i] = F77_CALL(ddot)(&N,readPtr1,&lda,readPtr1,&M);
		//res = ddot(n, x, incx, y, incy)
		rxy[i] 		= F77_CALL(ddot)(&N,readPtr1,&lda,readPtr2,&M);
		DxxRxy[i] 	= rxy[i]/dxx[i];		
	}
	
	
	//abs(N*sigma2*IM - (Y*(Y'-X'*DxxRxy)))./W ; W[i,i] = inf.
	//printMat(DxxRxy,M,1);

	//cache X[k,:]*DxxRxy[k]
	double * XDxxRxy;
	int MN = M*N;
	XDxxRxy = (double* ) Calloc(MN, double);
	//for(i=0;i<M;i++)
	//{
	//	for(j=0;j<N;j++)
	//	{
			//X[i,j] * DxxRxy[i];
	//		index = j*M + i;
	//		XDxxRxy[index] = -X[index]*DxxRxy[i];
	//	}
	//}
	F77_CALL(dcopy)(&MN,X,&inci,XDxxRxy,&incj);
	double alpha;	
	for(i=0;i<M;i++)
	{
		alpha  = -DxxRxy[i];
		readPtr1 = &XDxxRxy[i]; //ith row
		F77_CALL(dscal)(&N,&alpha, readPtr1,&M);//	(n, a, x, incx)
	}
	
	
	//printMat(XDxxRxy,M,1);
	// Y- XDxxRxy  			daxpy(n, a, x, incx, y, incy) y= ax + y
	alpha  = 1.0;
	// XDxxRxy <- alpha*Y + XDxxRxy
	F77_CALL(daxpy)(&MN,&alpha,Y,&inci,XDxxRxy,&inci);
	//printMat(XDxxRxy,M,1);
	double *YYXDR; //= Y*XDxxRxy'
	int MM = M*M;
	YYXDR = (double* ) Calloc(MM, double);	
	
//C := alpha*op(A)*op(B) + beta*C,
	double beta;
	char transa = 'N';
	char transb = 'T';
	alpha = -1;
	beta = 0;
	F77_CALL(dgemm)(&transa, &transb,&M, &M, &N,&alpha, Y,&M, XDxxRxy, &M, &beta,YYXDR, &M); //M xK, K xN  --> MxN, N xM --> M <-M, N<-M, k<-N
	//diagonal -->0; other element /wij
	
	//printMat(YYXDR,M,3);
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
	//printMat(YYXDR,M,3);
	//BLAS_extern int    /* IDAMAX - return the index of the element with max abs value */
	//F77_NAME(idamax)(const int *n, const double *dx, const int *incx);
	index = F77_CALL(idamax)(&MM,YYXDR,&inci);
	//Rprintf("index: %d\n",index);
	lambda_max = fabs(YYXDR[index-1]);

	Free(dxx);
	Free(rxy);
	Free(DxxRxy);
	//Free(XX);
	Free(XDxxRxy);
	Free(YYXDR);
	
	return lambda_max;	
}

//Q[i,k] =	N*sigma2*IM - (Y*(Y'-X'*DxxRxy)))
void QlambdaStart(double *Y,double *X, double *Q, double sigma2,int M, int N)
{	
	// Oct 08, 2012: assume one eQTL for each gene; This fucntion needs significant revision if this assumption doesnot hold
	double *dxx, *rxy, *DxxRxy,*readPtr1,*readPtr2;
	
	dxx				= (double* ) Calloc(M, double);
	rxy				= (double* ) Calloc(M, double);
	DxxRxy			= (double* ) Calloc(M, double);
	int i,index,ldk,lda,ldb,ldc;
	int inci = 1;
	int incj = 1; 
	//double norm;
	lda = M;
	for(i=0;i<M;i++)
	{
		readPtr1 	= &X[i]; //ith row
		readPtr2 	= &Y[i];

		//norm  		= F77_CALL(dnrm2)(&N,readPtr1,&M);	
		//dxx[i] 		= pow(norm,2);
		dxx[i] = F77_CALL(ddot)(&N,readPtr1,&lda,readPtr1,&M);
		//res = ddot(n, x, incx, y, incy)
		rxy[i] 		= F77_CALL(ddot)(&N,readPtr1,&lda,readPtr2,&M);
		DxxRxy[i] 	= rxy[i]/dxx[i];		
	}
	//abs(N*sigma2*IM - (Y*(Y'-X'*DxxRxy)))./W ; W[i,i] = inf.
	double Nsigma2  = N*sigma2; 			// int * double --> double

	//cache X[k,:]*DxxRxy[k]
	double * XDxxRxy;
	int MN = M*N;
	XDxxRxy = (double* ) Calloc(MN, double);
	F77_CALL(dcopy)(&MN,X,&inci,XDxxRxy,&incj);
	double alpha;	
	for(i=0;i<M;i++)
	{
		alpha  = -DxxRxy[i];
		readPtr1 = &XDxxRxy[i]; //ith row
		F77_CALL(dscal)(&N,&alpha, readPtr1,&M);//	(n, a, x, incx)
	}
	
	// Y- XDxxRxy  			daxpy(n, a, x, incx, y, incy) y= ax + y
	alpha  = 1.0;
	// XDxxRxy <- alpha*Y + XDxxRxy
	F77_CALL(daxpy)(&MN,&alpha,Y,&inci,XDxxRxy,&incj);
	//double *YYXDR; //= Y*XDxxRxy' 		--> Q

	double beta;
	char transa = 'N';
	char transb = 'T';
	alpha = -1;
	beta = 0;
	//F77_CALL(dgemm)(&transa, &transb,&M, &M, &N,&alpha, Y,&M, XDxxRxy, &M, &beta,Q, &M); //M xK, K xN  --> MxN, N xM --> M <-M, N<-M, k<-N
	//transpose 	(Y-X*DxxRxy)*Y'

	ldb = M;
	ldc = M;
	ldk = M;
	F77_CALL(dgemm)(&transa, &transb,&M, &lda, &N,&alpha, XDxxRxy,&ldb, Y, &ldc, &beta,Q, &ldk); //M xK, K xN  --> MxN, N xM --> M <-M, N<-M, k<-N	
	
	//diagonal -->0; other element /wij
	for(i=0;i<M;i++)
	{
		index = i*M + i;
		Q[index]= Q[index] + Nsigma2;
	}	
	
	Free(dxx);
	Free(rxy);
	Free(DxxRxy);
	//Free(XX);
	Free(XDxxRxy);

	
}

// 8888888888888888888888888888888888888888888888888888888888888888888888
//Q = N*sigma2*inv(I-B)-(Y-B*Y-fL*X)-mueL*ones(1,N))*Y';
void QlambdaMiddle(double *Y,double *X, double *Q,double *B,double *f, double *mue, double sigma2,int M, int N)
{	
	// Oct 08, 2012: assume one eQTL for each gene; This fucntion needs significant revision if this assumption doesnot hold
	//I - B; copy of IB for inverse
	double *IB, *IBinv,*IBcopy;
	int MM = M*M;
	int MN = M*N;
	IB = (double* ) Calloc(MM, double);
	IBinv = (double* ) Calloc(MM, double);
	IBcopy = (double* ) Calloc(MM, double);
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
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci);//initialized
	F77_CALL(dcopy)(&MM,&alpha,&inc0,IBinv,&inci);
	
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		IB[index] = 1 + IB[index];
		IBinv[index] = 1;
	}
	F77_CALL(dcopy)(&MM,IB,&inci,IBcopy,&incj);	

	
	//MatrixInverse(IBinv,M);
	//By linear solver: not inverse (IB*x = IM) result stored in IM; 
	//multiLinearSystem(IB, M,IBinv,M);
	int info = 0;
	int *ipiv;
	ipiv = (int *) Calloc(M,int);
	int lda = M;
	int ldb = M;
	int ldc = M;
	int ldk = M;
	F77_CALL(dgesv)(&M, &ldk, IBcopy, &lda, ipiv, IBinv, &ldb, &info);

	
	
	//abs(N*sigma2*inv(I-B) - NOISE*Y'.
	double Nsigma2  = N*sigma2; 			// int * double --> double
	double *Noise;
	Noise = (double* ) Calloc(MN, double);	
	//(I-B)*Y-bsxfun(@times,f,X);
	char transa = 'N';
	char transb = 'N';
	alpha = 1;
	F77_CALL(dgemm)(&transa, &transb,&M, &N, &ldk,&alpha, IB, &lda, Y, &ldb, &beta, Noise, &ldc);
	double *readPtr1, *readPtr2;
	for(i=0;i<M;i++)
	{
		readPtr1 = &X[i];
		readPtr2 = &Noise[i];
		alpha = -f[i]; // y= alpha x + y
		F77_CALL(daxpy)(&N, &alpha,readPtr1, &lda,readPtr2, &M);
	}//row i = 1:M

	//NoiseMat(Noise,B,Y, f, X, M, N);
	//Noise - Mue
	//Errs(irho,cv)   = norm((1-Missing_test).*(A*Ytest-bsxfun(@times,fR,Xtest)-mueR*ones(1,Ntest)),'fro')^2;
	//Noise - mue: Mx N
	alpha = -1;
	for(i=0;i<N;i++)
	{
		readPtr1 = &Noise[i*M];
		F77_CALL(daxpy)(&M, &alpha,mue, &inci,readPtr1, &incj);
	}	
	//Nsigma2*IBinv -  Noise *Y'
	//-Noise*Y' -->Q  	C := alpha*op(A)*op(B) + beta*C,
	//alpha = -1;
	transb = 'T';
	F77_CALL(dgemm)(&transa, &transb,&M, &ldk, &N,&alpha, Noise, &lda, Y, &ldb, &beta, Q, &ldc);
	//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
	alpha = Nsigma2;
	F77_CALL(daxpy)(&MM, &alpha,IBinv, &inci,Q, &incj);
	
	Free(IB);
	Free(IBinv);
	Free(IBcopy);
	Free(Noise);
	Free(ipiv);
	
}


void QlambdaMiddleCenter(double *Y,double *X, double *Q,double *B,double *f, double sigma2,int M, int N,
					double *IBinv)
{	
	// Oct 08, 2012: assume one eQTL for each gene; This fucntion needs significant revision if this assumption doesnot hold
	//I - B; copy of IB for inverse
	double *IB; 	//, *IBinv,*IBcopy
	int MM = M*M;
	int MN = M*N;
	IB = (double* ) Calloc(MM, double);
	//IBinv = (double* ) Calloc(MM, double);
	//IBcopy = (double* ) Calloc(MM, double);
	int inci = 1;
	int incj = 1;
	//int inc0 = 0;
	F77_CALL(dcopy)(&MM,B,&inci,IB,&incj);	
	int i,index;
	double alpha;
	double beta = 0;
	alpha = -1;
	F77_CALL(dscal)(&MM,&alpha,IB,&inci);
	//alpha = 0;
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci);//initialized
	//F77_CALL(dcopy)(&MM,&alpha,&inc0,IBinv,&inci);
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		IB[index] = 1 + IB[index];
		//IBinv[index] = 1;
	}
	//F77_CALL(dcopy)(&MM,IB,&inci,IBcopy,&incj);	

	
	//MatrixInverse(IBinv,M);
	//By linear solver: not inverse (IB*x = IM) result stored in IM; 
	//multiLinearSystem(IB, M,IBinv,M);
	//int info = 0;
	//int *ipiv;
	//ipiv = (int *) Calloc(M,int);
	int lda = M;
	int ldb = M;
	int ldc = M;
	int ldk = M;
	//F77_CALL(dgesv)(&M, &ldk, IBcopy, &lda, ipiv, IBinv, &ldb, &info);

	
	
	//abs(N*sigma2*inv(I-B) - NOISE*Y'.
	double Nsigma2  = N*sigma2; 			// int * double --> double
	double *Noise;
	Noise = (double* ) Calloc(MN, double);	
	//(I-B)*Y-bsxfun(@times,f,X);
	char transa = 'N';
	char transb = 'N';
	alpha = 1;
	F77_CALL(dgemm)(&transa, &transb,&M, &N, &ldk,&alpha, IB, &lda, Y, &ldb, &beta, Noise, &ldc);
	double *readPtr1, *readPtr2;
	for(i=0;i<M;i++)
	{
		readPtr1 = &X[i];
		readPtr2 = &Noise[i];
		alpha = -f[i]; // y= alpha x + y
		F77_CALL(daxpy)(&N, &alpha,readPtr1, &lda,readPtr2, &M);
	}//row i = 1:M

	//NoiseMat(Noise,B,Y, f, X, M, N);
	//Noise - Mue
	//Errs(irho,cv)   = norm((1-Missing_test).*(A*Ytest-bsxfun(@times,fR,Xtest)-mueR*ones(1,Ntest)),'fro')^2;
	//Noise - mue: Mx N

	//Nsigma2*IBinv -  Noise *Y'
	//-Noise*Y' -->Q  	C := alpha*op(A)*op(B) + beta*C,
	alpha = -1;
	transb = 'T';
	F77_CALL(dgemm)(&transa, &transb,&M, &ldk, &N,&alpha, Noise, &lda, Y, &ldb, &beta, Q, &ldc);
	//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
	alpha = Nsigma2;
	F77_CALL(daxpy)(&MM, &alpha,IBinv, &inci,Q, &incj);
	
	Free(IB);
	//Free(IBinv);
	//Free(IBcopy);
	Free(Noise);
	//Free(ipiv);
	
}


// 8888888888888888888888888888888888888888888888888888888888888888888888
//BLOCK COORDINATE ASCENSION: QIBinv = inv(I-B): by multi linear system
void UpdateIBinvPermute(double *QIBinv, double *B, int M)
{
	//I - B; copy of IB for inverse
	double *IB,*IBinv;	//, *IBinv,*IBcopy;
	int MM = M*M;
	int lda = M;
	int ldb = M;
	int ldk = M;
	IB = (double* ) Calloc(MM, double);
	IBinv = (double* ) Calloc(MM, double);
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
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci);//initialized
	F77_CALL(dcopy)(&MM,&alpha,&inc0,IBinv,&inci);
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		IB[index] = 1 + IB[index];
		IBinv[index] = 1;
	}
	
	//MatrixInverse(IBinv,M);
	//By linear solver: not inverse (IB*x = IM) result stored in IM; 
	//multiLinearSystem(IB, M,IBinv,M);
	int info = 0;
	int *ipiv;
	ipiv = (int *) Calloc(M,int);
	F77_CALL(dgesv)(&M, &ldk, IB, &lda, ipiv, IBinv, &ldb, &info);
	double *ptr1,*ptr2;
	//for(i=0;i<M;i++) Rprintf("IPIV: \n: %d \t",ipiv[i]);
	//Rprintf("\n");
	
	
	for(i=0;i<M;i++)
	{
		index = ipiv[i] -1;
		ptr1 = &QIBinv[index*M];
		ptr2 = &IBinv[i*M];
		F77_CALL(dcopy)(&M,ptr2,&inci,ptr1,&incj);
		
	}
	
	Free(IB);
	Free(ipiv);
	Free(IBinv);
}


// 8888888888888888888888888888888888888888888888888888888888888888888888
//BLOCK COORDINATE ASCENSION: QIBinv = inv(I-B): by multi linear system
void UpdateIBinv(double *QIBinv, double *B, int M)
{
	//I - B; copy of IB for inverse
	double *IB;	//, *IBinv,*IBcopy;
	int MM = M*M;
	int lda = M;
	int ldb = M;
	int ldk = M;
	IB = (double* ) Calloc(MM, double);

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
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci);//initialized
	F77_CALL(dcopy)(&MM,&alpha,&inc0,QIBinv,&inci);
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		IB[index] = 1 + IB[index];
		QIBinv[index] = 1;
	}
	
	//MatrixInverse(IBinv,M);
	//By linear solver: not inverse (IB*x = IM) result stored in IM; 
	//multiLinearSystem(IB, M,IBinv,M);
	int info = 0;
	int *ipiv;
	ipiv = (int *) Calloc(M,int);
	F77_CALL(dgesv)(&M, &ldk, IB, &lda, ipiv, QIBinv, &ldb, &info);

	Free(IB);
	Free(ipiv);
}

//no Missing
//--------------------------------------------------------------------------------  WEIGHTED_LASSOSF
double Weighted_LassoSf(double * W, double *B, double *f, double *Ycopy,double *Xcopy,
		double *Q, double lambda_factor, double lambda_factor_prev, double sigma2, int max_iter,
		int M, int N, int verbose,double *QIBinv,double lambda_max)			//double * mue,
{
	int i,j,index,ldM;
	//lda = M;
	//ldb = M;ldb,
	ldM = M;//fixed
	// return lambda;
	double *meanY, *meanX;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	
	//copy Y, X; 
	double *Y, *X;
	int MN = M*N;
	int MM = M*M;
	Y = (double* ) Calloc(MN, double);
	X = (double* ) Calloc(MN, double);
	
	//F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
		//double *dy, const int *incy);
	int inci,incj, inc0;
	inci	= 1;
	incj 	= 1;
	inc0 	= 0;
	F77_CALL(dcopy)(&MN,Ycopy,&inci,Y,&incj);
	F77_CALL(dcopy)(&MN,Xcopy,&inci,X,&incj);
	centerYX(Y,X, meanY, meanX,M, N);
	
	//return value
	//double sigma2 			= SIGMA2[0];
	double lambda;//lambda_max,
	//lambdaMax
	//lambda_max 				= lambdaMax(Y,X,W,M, N);
	if(verbose>4) Rprintf("\t\t\t\tEnter Function: weighted_LassoSf. The maximum lambda is: %f\n\n",lambda_max);
	lambda 					= lambda_factor*lambda_max;
	
	//none zeros
	double alpha,beta;
	beta = 0;
	double deltaLambda;
	double *s, *S,*Wcopy;
	S = (double* ) Calloc(MM, double);
	s = (double* ) Calloc(M, double);
	Wcopy = (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,W,&inci,Wcopy,&incj);

	deltaLambda 			= (2*lambda_factor - lambda_factor_prev)*lambda_max;	
	F77_CALL(dscal)(&MM,&deltaLambda,Wcopy,&inci); //wcopy = deltaLambda*W
	
	//ei = 0
	double *ei,toyZero;
	toyZero= 0;
	ei = (double* ) Calloc(M, double);
	//F77_CALL(dscal)(&M,&toyZero,ei,&inci);
	F77_CALL(dcopy)(&M,&toyZero,&inc0,ei,&inci);
/*	double *eye;
	eye = (double* ) Calloc(M, double);
	alpha = 1;
	F77_CALL(dcopy)(&M,&alpha,&inc0,eye,&inci);
*/
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
/*	ldk = M;
	//lda = M;
	F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, S, &ldM, eye, &inci, &beta,s, &incj);
*/	
//printMat(W,M,M);
	//f0, F1
	double *f0,*F1;
	//int qdif = M*M;
	f0 	= (double* ) Calloc(M, double);
	F1 	= (double* ) Calloc(MM, double);
	
	double *y_j;
	//xi 	= (double* ) Calloc(N, double);
	y_j 	= (double* ) Calloc(N, double);
	double *F1ptr;


	double XYi, XXi;
	for(i=0;i<M;i++)
	{
		readPtr = &X[i];
		//F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
		readPtr2 = &Y[i];
		//F77_CALL(dcopy)(&N,readPtr2,&M,y_j,&inci);

		//dot product
		//XYi = F77_CALL(ddot)(&N, xi, &inci,y_j, &incj);
		XYi = F77_CALL(ddot)(&N, readPtr, &M,readPtr2, &M);
		//XXi = F77_CALL(ddot)(&N, xi, &inci,xi, &incj);
		//norm2 = F77_CALL(dnrm2)(&N,xi,&inci);
		//XXi 	= pow(norm2,2);
		XXi = F77_CALL(ddot)(&N, readPtr, &M,readPtr, &M);
		f0[i] 	= XYi/XXi;
		F1ptr	= &F1[M*i];//start from ith column
		//Y*X(i,:)' 		y := alpha*A*x + beta*y,
		alpha = 1/XXi;
		F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, readPtr, &M, &beta,F1ptr, &incj);
	}
	
	//printMat(f0,M,1);

	// entering loop		
	double *IBinv,*zi,*a_iT;		// y_j: one row of Y: Nx1
	IBinv 	= (double* ) Calloc(MM, double);
	//zi 		= (double* ) Calloc(M, double);
	a_iT 	= (double* ) Calloc(N, double);


	
	//loop starts here
	int iter = 0;
	double js_i, m_ij,B_old, lambdaW,beta_ij,r_ij, Bij;
	//dynamic variable keep intermidiate values 
	double *eiB;
	eiB = (double* ) Calloc(M, double);
	double *BiT;
	BiT = (double* ) Calloc(M, double);
	//quadratic function
	double d_ij, theta_ijp,k_ijp,q_ijp,Bijpp, Bijpm; //case (14)
	double q_ijm, theta_ijm, Bijmm, Bijmp,Lss,candsBij,LssCands;
	
	//converge of gene i
	double dB,ziDb,BF1;
	
	//converge of while
	double delta_BF,FnormOld, FnormChange;
	double *BfOld,*BfNew,*BfChange;
	index = M*(M  +1);
	BfOld = (double* ) Calloc(index, double);
	BfNew = (double* ) Calloc(index, double);
	BfChange = (double* ) Calloc(index, double);
	
	//	//linear system
	//int *ipiv;
	//ipiv = (int *) R_alloc(N,sizeof(int));
	//ipiv = (int *) Calloc(M,int);
	//int nrhs = 1;
	//int info = 0;
//for(i=0;i<M;i++) ipiv[i] = i+1;
	//
	//double *aiTQuad; 
	//aiTQuad 	= (double* ) Calloc(N, double);

	//alpha =  F77_CALL(dnrm2)(&index,BfOld,&inci);	
//Rprintf("BfOldnorm:%f\n",alpha);	
	//int anyUpdateRow = 0;
	
	while(iter < max_iter)
	{
		iter = iter + 1;
		//converge Bfold = [B f];
		F77_CALL(dcopy)(&MM,B,&inci,BfOld,&incj);
		//last column
		F1ptr = &BfOld[MM];
		F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
	//alpha =  F77_CALL(dnrm2)(&index,BfOld,&inci);	
//Rprintf("BfOldnorm:%f\n",alpha);
		//printMat(BfOld,M,M+1);
		//
		for(i=0;i<M;i++)
		{
			if(s[i] >0)
			{ 	//
				if(verbose>6) Rprintf("\t\t\t\t\t updating gene %d \n",i);
				//
				ei[i] = 1;
				//zi   IBinv = I -B
				//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
				//alpha = -1; 
				//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
				//diagonal + 1
				//for(j=0;j<M;j++) 
				//{
				//	index = j*M + j;
				//	IBinv[index] = 1 + IBinv[index];
				//}

				//call matrix inverse
				//MatrixInverse(IBinv,M);
				//zi is ith column of IBinv
				//for(j =0;j<M;j++) zi[j] 	= IBinv[i*M+j];
				
				// by linear solver: zi = (I - B)\ei;
				//copy ei to zi;
				zi = &QIBinv[i*M];
				//F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
				//
				//if (i==0) printMat(zi,M,1);
				//
				//call linearSystem(double *a, int N,double *B) a is NxN
//linearSystem(IBinv, M,zi);  //zi is updated.

				//F77_CALL(dgesv)(&M, &nrhs, IBinv, &ldM, ipiv, zi, &ldb, &info);
				//anyUpdateRow = 1;
//j reserved				//for j = js_i
				for(j=0;j<M;j++) 
				{
					js_i = S[j*M + i]; 		//ith row
					if(js_i >0)
					{

						m_ij 	= zi[j];
						B_old 	= B[j*M + i]; //B[i,j]
						if(j!=i)
						{
						
							//y_j; jth row Nx1
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
							//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
							alpha = 1;
							F77_CALL(daxpy)(&M, &alpha,BiT, &inci,eiB, &incj);
							//a_iT      = (ei'-BiT)*Y-f(i)*X(i,:);
							readPtr = &X[i];
							F77_CALL(dcopy)(&N,readPtr,&M,a_iT,&inci);	

							//a_iT = -f[i]*xi 		dscal(n, a, x, incx) 		x = a*x
							alpha = -f[i];
							F77_CALL(dscal)(&N,&alpha,a_iT,&inci);							

							transa='T'; //y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
							beta = 1;
							alpha = 1;
							F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, eiB, &inci, &beta,a_iT, &incj);
							
							//r_ij                = y_j'*y_j;	
r_ij = F77_CALL(ddot)(&N, y_j, &inci,y_j, &incj);
							//norm2 = F77_CALL(dnrm2)(&N,y_j,&inci);
							//r_ij 	= pow(norm2,2);
							
							//beta_ij             = a_iT*y_j;
							beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
							
							if (fabs(m_ij)<1e-10) //go to the linear equation 
							{
								//
								if(verbose>7) Rprintf("\t\t\t\t\t\t\t gene %d \t interact with gene %d.\tLinear equation\n",i,j);
								//
								Bij = (beta_ij-lambdaW)/r_ij;
								//Rprintf("\t\t gene (%d,\t%d): linear Bij %f\n",i,j,Bij);
				//
				//if (i==0) 
				//{
				//	Rprintf("\t\t\t beta_ij: %f;\t r_ij:%f; lambdaW: %f\n ",beta_ij,r_ij,lambdaW);
				//	printMat(a_iT,N,1);
					//printMat(y_j,N,1);
				//	printMat(eiB,M,1);
				//}
				//
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
								//Lss = quadraticLss(sigma2, N, d_ij, candsBij, r_ij, lambdaW,beta_ij);								
								//Lss = quadraticLssJuan(sigma2,N, d_ij, candsBij, a_iT,y_j,lambdaW);
//beta = F77_CALL(ddot)(&N, a_iT, &inci,a_iT, &incj);
								//norm2 = F77_CALL(dnrm2)(&N,a_iT,&inci);
								//beta 	= pow(norm2,2);
							
								//Lss = sigma2*N*log(fabs(d_ij - candsBij)+1e-16) - beta/2 -lambdaW*fabs(candsBij);
								//Lss = sigma2*N*log(fabs(d_ij)+1e-16) - beta/2;
								Lss = sigma2*N*log(fabs(d_ij)+1e-16);
								
								// a_iT = a_iT - candsBij*y_j 		daxpy(n, a, x, incx, y, incy) 	y := a*x + y							
								if (Bijpp>0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijpp, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijpp, a_iT,y_j,lambdaW);	
									//aiTQuad dcopy(n, x, incx, y, incy)  y= x
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijpp;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - beta/2 -lambdaW*fabs(Bijpp);
									LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - r_ij*pow(Bijpp,2)/2 + beta_ij*Bijpp -lambdaW*fabs(Bijpp); 
									
									if(LssCands>Lss) 
									{
										candsBij = Bijpp;
										Lss 	= LssCands;
									}	
								}
								if (Bijpm>0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijpm, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijpm, a_iT,y_j,lambdaW);								
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijpm;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijpm)+1e-16) - beta/2 -lambdaW*fabs(Bijpm);
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
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijmm, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijmm, a_iT,y_j,lambdaW);								
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijmm;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - beta/2 -lambdaW*fabs(Bijmm);
									LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - r_ij*pow(Bijmm,2)/2 + beta_ij*Bijmm -lambdaW*fabs(Bijmm);  
									if(LssCands>Lss) 
									{
										candsBij = Bijmm;
										Lss 	= LssCands;
									}	
								}
								if (Bijmp<0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijmp, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijmp, a_iT,y_j,lambdaW);									
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijmp;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijmp)+1e-16) - beta/2 -lambdaW*fabs(Bijmp);			
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
						
						//update QIBinv ith column updated; jth row: 
						//readPtr = &QIBinv[j];
						//F77_CALL(dscal)(&M,&ziDb,readPtr,&M);
						//QIBinv[j*M + i] = QIBinv[j*M + i]*ziDb;
						// adjugate of QIBinv = ziDb*(QIBinvadj - (-1)^kl) + (-1)^kl
						//for(k=0;k<M;k++) // kth row
						//{
						//	for(l=0;l<M;l++)//lth column
						//	{
						//		if(k!=j && l!=i) 
						//		{
						//			alpha = pow(-1.0,k+l);
						//			QIBinv[l*M+k] = ziDb*(QIBinv[l*M+k] -alpha) + alpha;
						//		}	
						//	}
						//}
						
					}//js_i >0
				}//j = 1:M	
				//f
				//BF1 = B(i,:)*F1(:,i)
				readPtr = &B[i];
				F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);

				F1ptr = &F1[M*i];
				BF1 = F77_CALL(ddot)(&M, BiT, &inci,F1ptr, &incj);

				f[i] = f0[i] - BF1;
				ei[i] = 0; // re-set ei for next i
			}else//s[i]  no un-zero weight in this gene
			{
				readPtr = &B[i];
				//F77_CALL(dscal)(&M,&toyZero,readPtr,&ldM);
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
		
//Rprintf("FnormOld: %f; \t FnormChange: %f\n",FnormOld,FnormChange);
		//MASTER IS HERE: Update before break: IBinv will be used in Qij
		// BLOCK COORDINATE ASCEND: Update IBinv
		UpdateIBinv(QIBinv, B,M);	
		
		if(verbose>5) Rprintf("\t\t\t\t\t\tdelta_BF: %f\n",delta_BF);
		if(delta_BF<1e-3)		//break out
		{
			break;
		}
		

		
	}//while

	//QIBinv: compute not updated ones
	//F77_NAME(dgetrs)(const char* trans, const int* n, const int* nrhs,
	//	 const double* a, const int* lda, const int* ipiv,
	//	 double* b, const int* ldb, int* info);
	//if(anyUpdateRow>0) //QIBinv --> = Cij/det(I-B) != eye(M); ipiv not zeros
	//{
	//	transa = 'N';
	//	for(i=0;i<M;i++)
	//	{
	//		if(s[i] ==0)
	//		{ 	
	//			info = 0;
	//			ei[i] = 1;
	//			zi = &QIBinv[i*M];
	//			F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
				// dgetrs
	//			F77_CALL(dgetrs)(&transa, &M, &nrhs,IBinv,&M,ipiv,zi,&ldb,&info); //column i updated
	//			ei[i] = 0;
	//		}
	//	}
	//}else
	//{
		//initialize IBinv
	//	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	//	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
	//}
	
	
	
	if(verbose>4) Rprintf("\t\t\t\tCurrent lambda: %f;\t number of iteration is: %d.\tExiting Weighted_LassoSf\n\n",lambda, iter);
	
	//mue          = (IM-B)*meanY-bsxfun(@times,f,meanX);

	//IBinv = I -B
	//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
	//alpha = -1; 
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
	//diagonal + 1
	//for(j=0;j<M;j++) 
	//{
	//	index = j*M + j;
	//	IBinv[index] = 1 + IBinv[index];
	//	mue[j] = -f[j]*meanX[j];
	//}
	//mue[i] = -f[i]*meanX[i];
	
	//	//dgemv mue = Ax + beta*mue
	//transa = 'N';
	//alpha = 1;
	//beta = 1;
	//ldk = M;
	//F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, IBinv, &ldM, meanY, &inci, &beta,mue, &incj);

	//mue                     	= (IM-B)*meanY-bsxfun(@times,f,meanX); diag(f)*meanX

	Free(meanY);
	Free(meanX);
	Free(Y);
	Free(X);
	
	Free(S);
	Free(s);
	Free(f0);
	Free(F1);
	Free(Wcopy);
	
	//Free(xi);
	Free(y_j);
	
	Free(ei);
	Free(IBinv);
	//Free(zi);
	Free(a_iT);
	
	Free(eiB);
	Free(BiT);
	Free(BfOld);
	Free(BfNew);
	Free(BfChange);
	
	//Free(ipiv);
	
	//Free(aiTQuad);
	//Free(eye);
	return lambda;
	//sigma2 remains the same
}//weighted_LassoSf

//combine lasso and zero_lasso for CV_selecting lambda
double Weighted_LassoSf_MLf(double * W, double *BL, double *fL, double *Ycopy,double *Xcopy,
		double *Q, double lambda_factor, double lambda_factor_prev, double sigma2, int max_iter,
		int M, int N, int verbose,
		double *BC, double *fC, double *mue,double *QIBinv,double *IBinvZero,double lambda_max)
{
	//SET TO PART1: LASSO
	double *B, *f;
	B = &BL[0];
	f = &fL[0];
	//mue is used for calculation when Y, X are not centered: for testing set;

	int i,j,index,ldk,ldM;
	//lda = M;
	//ldb = M;ldb,
	ldM = M;//fixed
	// return lambda;
	double *meanY, *meanX;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	
	//copy Y, X; 
	double *Y, *X;
	int MN = M*N;
	int MM = M*M;
	Y = (double* ) Calloc(MN, double);
	X = (double* ) Calloc(MN, double);
	
	//F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
		//double *dy, const int *incy);
	int inci,incj,inc0;
	inci	= 1;
	incj 	= 1;
	inc0 	= 0;
	F77_CALL(dcopy)(&MN,Ycopy,&inci,Y,&incj);
	F77_CALL(dcopy)(&MN,Xcopy,&inci,X,&incj);
	centerYX(Y,X, meanY, meanX,M, N);
	
	//return value
	//double sigma2 			= SIGMA2[0];
	double lambda;//lambda_max,
	//lambdaMax
	//lambda_max 				= lambdaMax(Y,X,W,M, N);
	if(verbose>4) Rprintf("\t\t\t\tEnter Function: weighted_LassoSf. The maximum lambda is: %f\n\n",lambda_max);
	lambda 					= lambda_factor*lambda_max;
	
	//none zeros
	double alpha,beta;
	beta = 0;
	double deltaLambda;
	double *s, *S,*Wcopy;
	S = (double* ) Calloc(MM, double);
	s = (double* ) Calloc(M, double);
	Wcopy = (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,W,&inci,Wcopy,&incj);

	deltaLambda 			= (2*lambda_factor - lambda_factor_prev)*lambda_max;	
	F77_CALL(dscal)(&MM,&deltaLambda,Wcopy,&inci); //wcopy = deltaLambda*W
	
	//ei = 0
	double *ei,toyZero;
	toyZero= 0;
	ei = (double* ) Calloc(M, double);
	//F77_CALL(dscal)(&M,&toyZero,ei,&inci);
	F77_CALL(dcopy)(&M,&toyZero,&inc0,ei,&inci);
/*	double *eye;
	eye = (double* ) Calloc(M, double);
	alpha = 1;
	F77_CALL(dcopy)(&M,&alpha,&inc0,eye,&inci);
*/
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
/*	ldk = M;
	//lda = M;
	F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, S, &ldM, eye, &inci, &beta,s, &incj);
*/	
//printMat(W,M,M);
	//f0, F1
	double *f0,*F1;
	//int qdif = M*M;
	f0 	= (double* ) Calloc(M, double);
	F1 	= (double* ) Calloc(MM, double);
	
	//double *xi, *y_j;
	double *y_j;
	//xi 	= (double* ) Calloc(N, double);
	y_j 	= (double* ) Calloc(N, double);
	double *F1ptr;


	double XYi, XXi;
	for(i=0;i<M;i++)
	{
		readPtr = &X[i];
		//F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
		readPtr2 = &Y[i];
		//F77_CALL(dcopy)(&N,readPtr,&M,y_j,&inci);

		//dot product
		//XYi = F77_CALL(ddot)(&N, xi, &inci,y_j, &incj);
		XYi = F77_CALL(ddot)(&N, readPtr, &M,readPtr2, &M);
		//XXi = F77_CALL(ddot)(&N, xi, &inci,xi, &incj);
		//norm2 = F77_CALL(dnrm2)(&N,xi,&inci);
		//XXi 	= pow(norm2,2);
		XXi = F77_CALL(ddot)(&N, readPtr, &M,readPtr, &M);
		f0[i] 	= XYi/XXi;
		F1ptr	= &F1[M*i];//start from ith column
		//Y*X(i,:)' 		y := alpha*A*x + beta*y, alpha*Y *xi + beta*F1
		alpha = 1/XXi;
		F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, readPtr, &M, &beta,F1ptr, &incj);
	}
	
	//printMat(f0,M,1);

	// entering loop
	double *IBinv,*zi,*a_iT;// y_j: one row of Y: Nx1
	IBinv 	= (double* ) Calloc(MM, double);
	//zi 		= (double* ) Calloc(M, double);
	a_iT 	= (double* ) Calloc(N, double);

	
	
	//loop starts here
	int iter = 0;
	double js_i, m_ij,B_old, lambdaW,beta_ij,r_ij, Bij;
	//dynamic variable keep intermidiate values 
	double *eiB;
	eiB = (double* ) Calloc(M, double);
	double *BiT;
	BiT = (double* ) Calloc(M, double);
	//quadratic function
	double d_ij, theta_ijp,k_ijp,q_ijp,Bijpp, Bijpm; //case (14)
	double q_ijm, theta_ijm, Bijmm, Bijmp,Lss,candsBij,LssCands;
	
	//converge of gene i
	double dB,ziDb,BF1;
	
	//converge of while
	double delta_BF,FnormOld, FnormChange;
	double *BfOld,*BfNew,*BfChange;
	index = M*(M  +1);
	BfOld = (double* ) Calloc(index, double);
	BfNew = (double* ) Calloc(index, double);
	BfChange = (double* ) Calloc(index, double);
	
	//	//linear system
	//int *ipiv;
	//ipiv = (int *) R_alloc(N,sizeof(int));
	//ipiv = (int *) Calloc(M,int);
	//int nrhs = 1;
	//int info = 0;
	//for(i=0;i<M;i++) ipiv[i] = i+1;
	//
	//double *aiTQuad; 
	//aiTQuad 	= (double* ) Calloc(N, double);

	//alpha =  F77_CALL(dnrm2)(&index,BfOld,&inci);	
//Rprintf("BfOldnorm:%f\n",alpha);	
	//int anyUpdateRow = 0;
	while(iter < max_iter)
	{
		iter = iter + 1;
		//converge Bfold = [B f];
		F77_CALL(dcopy)(&MM,B,&inci,BfOld,&incj);
		//last column
		F1ptr = &BfOld[MM];
		F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
	//alpha =  F77_CALL(dnrm2)(&index,BfOld,&inci);	
//Rprintf("BfOldnorm:%f\n",alpha);
		//printMat(BfOld,M,M+1);
		//
		for(i=0;i<M;i++)
		{
			if(s[i] >0)
			{ 	//
				if(verbose>6) Rprintf("\t\t\t\t\t updating gene %d \n",i);
				//
				ei[i] = 1;
				//zi   IBinv = I -B
				//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
				//alpha = -1; 
				//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
				//diagonal + 1
				//for(j=0;j<M;j++) 
				//{
				//	index = j*M + j;
				//	IBinv[index] = 1 + IBinv[index];
				//}

				//call matrix inverse
				//MatrixInverse(IBinv,M);
				//zi is ith column of IBinv
				//for(j =0;j<M;j++) zi[j] 	= IBinv[i*M+j];
				
				// by linear solver: zi = (I - B)\ei;
				//copy ei to zi;
				zi = &QIBinv[i*M];
				//F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
				//
				//if (i==0) printMat(zi,M,1);
				//
				//call linearSystem(double *a, int N,double *B) a is NxN
//linearSystem(IBinv, M,zi);  //zi is updated.

				//F77_CALL(dgesv)(&M, &nrhs, IBinv, &ldM, ipiv, zi, &ldb, &info);
				//anyUpdateRow = 1;
//j reserved				//for j = js_i
				for(j=0;j<M;j++) 
				{
					js_i = S[j*M + i]; 		//ith row
					if(js_i >0)
					{

						m_ij 	= zi[j];
						B_old 	= B[j*M + i]; //B[i,j]
						if(j!=i)
						{
						
							//y_j; jth row Nx1
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
							//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
							alpha = 1;
							F77_CALL(daxpy)(&M, &alpha,BiT, &inci,eiB, &incj);
							//a_iT      = (ei'-BiT)*Y-f(i)*X(i,:);
							readPtr = &X[i];
							F77_CALL(dcopy)(&N,readPtr,&M,a_iT,&inci);	

							//a_iT = -f[i]*xi 		dscal(n, a, x, incx) 		x = a*x
							alpha = -f[i];
							F77_CALL(dscal)(&N,&alpha,a_iT,&inci);							

							transa='T'; //y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
							beta = 1;
							alpha = 1;
							F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, eiB, &inci, &beta,a_iT, &incj);
							
							//r_ij                = y_j'*y_j;	
r_ij = F77_CALL(ddot)(&N, y_j, &inci,y_j, &incj);
							//norm2 = F77_CALL(dnrm2)(&N,y_j,&inci);
							//r_ij 	= pow(norm2,2);
							
							//beta_ij             = a_iT*y_j;
							beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
							
							if (fabs(m_ij)<1e-10) //go to the linear equation 
							{
								//
								if(verbose>7) Rprintf("\t\t\t\t\t\t\t gene %d \t interact with gene %d.\tLinear equation\n",i,j);
								//
								Bij = (beta_ij-lambdaW)/r_ij;
								//Rprintf("\t\t gene (%d,\t%d): linear Bij %f\n",i,j,Bij);
				//
				//if (i==0) 
				//{
				//	Rprintf("\t\t\t beta_ij: %f;\t r_ij:%f; lambdaW: %f\n ",beta_ij,r_ij,lambdaW);
				//	printMat(a_iT,N,1);
					//printMat(y_j,N,1);
				//	printMat(eiB,M,1);
				//}
				//
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
								//Lss = quadraticLss(sigma2, N, d_ij, candsBij, r_ij, lambdaW,beta_ij);								
								//Lss = quadraticLssJuan(sigma2,N, d_ij, candsBij, a_iT,y_j,lambdaW);
//beta = F77_CALL(ddot)(&N, a_iT, &inci,a_iT, &incj);
								//norm2 = F77_CALL(dnrm2)(&N,a_iT,&inci);
								//beta 	= pow(norm2,2);
							
								//Lss = sigma2*N*log(fabs(d_ij - candsBij)+1e-16) - beta/2 -lambdaW*fabs(candsBij);
								//Lss = sigma2*N*log(fabs(d_ij)+1e-16) - beta/2;
								Lss = sigma2*N*log(fabs(d_ij)+1e-16);
								
								// a_iT = a_iT - candsBij*y_j 		daxpy(n, a, x, incx, y, incy) 	y := a*x + y							
								if (Bijpp>0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijpp, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijpp, a_iT,y_j,lambdaW);	
									//aiTQuad dcopy(n, x, incx, y, incy)  y= x
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijpp;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - beta/2 -lambdaW*fabs(Bijpp);
									LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - r_ij*pow(Bijpp,2)/2 + beta_ij*Bijpp -lambdaW*fabs(Bijpp); 
									
									if(LssCands>Lss) 
									{
										candsBij = Bijpp;
										Lss 	= LssCands;
									}	
								}
								if (Bijpm>0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijpm, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijpm, a_iT,y_j,lambdaW);								
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijpm;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijpm)+1e-16) - beta/2 -lambdaW*fabs(Bijpm);
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
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijmm, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijmm, a_iT,y_j,lambdaW);								
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijmm;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - beta/2 -lambdaW*fabs(Bijmm);
									LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - r_ij*pow(Bijmm,2)/2 + beta_ij*Bijmm -lambdaW*fabs(Bijmm);  
									if(LssCands>Lss) 
									{
										candsBij = Bijmm;
										Lss 	= LssCands;
									}	
								}
								if (Bijmp<0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijmp, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijmp, a_iT,y_j,lambdaW);									
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijmp;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijmp)+1e-16) - beta/2 -lambdaW*fabs(Bijmp);			
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
						
						//update QIBinv ith column updated; jth row: 
						//readPtr = &QIBinv[j];
						//F77_CALL(dscal)(&M,&ziDb,readPtr,&M);
						//QIBinv[j*M + i] = QIBinv[j*M + i]*ziDb;
						// adjugate of QIBinv = ziDb*(QIBinvadj - (-1)^kl) + (-1)^kl
						//for(k=0;k<M;k++) // kth row
						//{
						//	for(l=0;l<M;l++)//lth column
						//	{
						//		if(k!=j && l!=i) 
						///		{
						//			alpha = pow(-1.0,k+l);
						//			QIBinv[l*M+k] = ziDb*(QIBinv[l*M+k] -alpha) + alpha;
						//		}	
						//	}
						//}
						
					}//js_i >0
				}//j = 1:M	
				//f
				//BF1 = B(i,:)*F1(:,i)
				readPtr = &B[i];
				F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);

				F1ptr = &F1[M*i];
				BF1 = F77_CALL(ddot)(&M, BiT, &inci,F1ptr, &incj);

				f[i] = f0[i] - BF1;
				ei[i] = 0; // re-set ei for next i
			}else//s[i]  no un-zero weight in this gene
			{
				readPtr = &B[i];
				//F77_CALL(dscal)(&M,&toyZero,readPtr,&ldM);
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
		
//Rprintf("FnormOld: %f; \t FnormChange: %f\n",FnormOld,FnormChange);
//MASTER IS HERE	IBinv will be used in Qij	
		// BLOCK COORDINATE ASCEND: Update IBinv
		UpdateIBinv(QIBinv, B,M);
		
		if(verbose>5) Rprintf("\t\t\t\t\t\tdelta_BF: %f\n",delta_BF);
		if(delta_BF<1e-3)		//break out
		{
			break;
		}

		
		
	}//while

	
	//QIBinv: compute not updated ones
	//F77_NAME(dgetrs)(const char* trans, const int* n, const int* nrhs,
	//	 const double* a, const int* lda, const int* ipiv,
	//	 double* b, const int* ldb, int* info);
	//printMat(IBinv,M,M);
	//if(anyUpdateRow>0) //QIBinv --> = Cij/det(I-B) != eye(M); ipiv not zeros
	//{	
	//	transa = 'N';
	//	for(i=0;i<M;i++)
	//	{
	//	//	if(s[i] ==0)
	//		{ 	
	//			info = 0;
	//			ei[i] = 1;
	//			zi = &QIBinv[i*M];
	//			F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
	//			// dgetrs
	//			F77_CALL(dgetrs)(&transa, &M, &nrhs,IBinv,&M,ipiv,zi,&ldb,&info); //column i updated
	//			ei[i] = 0;
	//		}
	//	}
	//}else
	//{
		//initialize IBinv
	//	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	//	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
	//}
	
	//printMat(QIBinv,M,M);
	
	
	if(verbose>3) Rprintf("\t\t\t\tCurrent lambda: %f;\t number of iteration is: %d.\tExiting Weighted_LassoSf\n\n",lambda, iter);
	
	//mue          = (IM-B)*meanY-bsxfun(@times,f,meanX);

	//IBinv = I -B
	//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
	//alpha = -1; 
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
	//diagonal + 1
	//for(j=0;j<M;j++) 
	//{
	//	index = j*M + j;
	//	IBinv[index] = 1 + IBinv[index];
	//	mueL[j] = -f[j]*meanX[j];
	//}
	//mue[i] = -f[i]*meanX[i];
	
	//	//dgemv mue = Ax + beta*mue
	//transa = 'N';
	//alpha = 1;
	//beta = 1;
	//ldk = M;
	//F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, IBinv, &ldM, meanY, &inci, &beta,mueL, &incj);

	//----------------------------------------------------------------------------------END OF LASSO
		//SET TO PART2: LASSO with lambda zero
	//double *B, double *f;
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
	//none zeros
	//double *s;
	//s = (double* ) Calloc(M, double);
	for(i=0;i<M;i++)
	{
		readPtr = &S[i]; //S[i,];
		s[i] = F77_CALL(dasum)(&M, readPtr, &ldM);
	}

	beta = 1;
	//f0, F1
	//ei = 0	
	// entering loop
	//loop starts here
	iter = 0;
	//double js_i, m_ij,B_old, beta_ij,r_ij;
	//int *ipiv;
	//ipiv = (int *) R_alloc(N,sizeof(int));
	//info = 0;
	
	//dynamic variable keep intermidiate values eiB = ei-BiT

	//quadratic function
	double theta_ij,k_ij,q_ij,Bijp, Bijm; //case (14)
	
	//converge of gene i
	
	//converge of while
	max_iter = max_iter/5;
	//linear system

	//
	//double *aiTQuad; 
	//aiTQuad 	= (double* ) Calloc(N, double);
	//double *zzi;
	//zzi 		= (double* ) Calloc(M, double);
	//zi = &zzi[0];
	while(iter < max_iter)
	{
		iter = iter + 1;
		//converge Bfold = [B f];
		F77_CALL(dcopy)(&MM,B,&inci,BfOld,&incj);
		//last column
		F1ptr = &BfOld[MM];
		F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
		//
//	printMat(ei,M,1);
//	printMat(F1,M,M);
		//
		
		// inner loop
		for(i=0;i<M;i++)
		{
			if(s[i] >0)
			{
				//
				if(verbose>6) Rprintf("\t\t updating gene %d \n",i);
				//
			ei[i] = 1;
				//zi
				//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
				//alpha = -1; 
				//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
				//diagonal + 1
				//for(j=0;j<M;j++) 
				//{
				//	index = j*M + j;
				//	IBinv[index] = 1 + IBinv[index];
				//}
				//call matrix inverse
				//MatrixInverse(IBinv,M);
				//zi is ith column of IBinv
				//for(j =0;j<M;j++) zi[j] 	= IBinv[i*M+j];
				
				//by linear solver: 
				//copy ei to zi;
				zi = &IBinvZero[i*M];
				//F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
				//call linearSystem(double *a, int N,double *B) a is NxN
				//linearSystem(IBinv, M,zi);  //zi is updated.

//Rprintf("M: %d; irhs: %d; ldM: %d; ldb: %d; info: %d\n",M,nrhs,ldM,ldb,info);
				//F77_CALL(dgesv)(&M, &nrhs, IBinv, &ldM, ipiv, zi, &ldb, &info);
								
				//printMat(zi,M,1);
				//for j = js_i
				for(j=0;j<M;j++)
				{
					js_i = S[j*M + i]; 		//ith row
					if(js_i >0)
					{
						//if(verbose>4) Rprintf("\t\t\t gene %d \t interact with gene %d\n",i,j);
					
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
						//a_iT      = (ei'-BiT)*Y-f(i)*X(i,:);
						//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
						alpha = 1;
						F77_CALL(daxpy)(&M, &alpha,BiT, &inci,eiB, &incj);
						//a_iT      = (ei'-BiT)*Y-f(i)*X(i,:);
						readPtr = &X[i];
						F77_CALL(dcopy)(&N,readPtr,&M,a_iT,&inci);
						//a_iT = -f[i]*xi 		dscal(n, a, x, incx) 		x = a*x
						alpha = -f[i];
						F77_CALL(dscal)(&N,&alpha,a_iT,&inci);							

						transa='T'; //y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
						//beta = 1;
						alpha = 1;
						F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, eiB, &inci, &beta,a_iT, &incj);

						//r_ij                = y_j'*y_j;	
						r_ij = F77_CALL(ddot)(&N, y_j, &inci,y_j, &inci);
						//norm2 = F77_CALL(dnrm2)(&N,y_j,&inci);
						//r_ij 	= pow(norm2,2);	
						//beta_ij             = a_iT*y_j;
						beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
						
						if (fabs(m_ij)<1e-10) //go to the linear equation 
						{
							//
							if(verbose>7) Rprintf("\t\t\t gene %d \t interact with gene %d.\tLinear equation\n",i,j);
							//
							//Bij = (beta_ij+lambdaW)/r_ij; 
							B[j*M+i] = beta_ij/r_ij;
							//Rprintf("\t\t\t beta_ij: %f;\t r_ij:%f\n ",beta_ij,r_ij);
							//printMat(a_iT,M,1);
						}else //m_ij ~=0 go to the quadratic equation
						{
							//
							if(verbose>7) Rprintf("\t\t\t gene %d \t interact with gene %d.\tQuadratic equation\n",i,j);
							//					
						
							//assume Bij >0
							d_ij = 1/m_ij + B[j*M+i];
							theta_ij = r_ij*d_ij + beta_ij;
							k_ij = d_ij*beta_ij - N*sigma2;
								
							q_ij = theta_ij*theta_ij - 4*r_ij* k_ij;
							Bijp = (1/(2*r_ij))*(theta_ij + sqrt(q_ij));
							Bijm = (1/(2*r_ij))*(theta_ij - sqrt(q_ij));
								
							candsBij = 0;
							//Lss = quadraticLss(sigma2, N, d_ij, candsBij, r_ij, lambdaW,beta_ij);
							//Lss = quadraticLssJuan(sigma2,N, d_ij, candsBij, a_iT,y_j,lambdaW);
							//beta = F77_CALL(ddot)(&N, a_iT, &inci,a_iT, &inci);
							//norm2 = F77_CALL(dnrm2)(&N,a_iT,&inci);
							//beta 	= pow(norm2,2);
							//Lss = sigma2*N*log(fabs(d_ij - candsBij)+1e-16) - beta/2 -lambdaW*fabs(candsBij);
							//Lss = sigma2*N*log(fabs(d_ij)+1e-16) - beta/2;															
							Lss = sigma2*N*log(fabs(d_ij)+1e-16);
							//Bijp
							//LssCands = quadraticLss(sigma2, N, d_ij, Bijp, r_ij, lambdaW,beta_ij);
							//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijp, a_iT,y_j,lambdaW);							
							//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
							//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
							//alpha = -Bijp;
							//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
							//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
							//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
							//beta 	= pow(norm2,2);							
							//LssCands = sigma2*N*log(fabs(d_ij - Bijp)+1e-16) - beta/2;
							LssCands = sigma2*N*log(fabs(d_ij - Bijp)+1e-16) - r_ij*pow(Bijp,2)/2 + beta_ij*Bijp;
							if(LssCands>Lss) 
							{
								candsBij = Bijp;
								Lss 	= LssCands;
							}	
							//Bijm>
							//LssCands = quadraticLss(sigma2, N, d_ij, Bijm, r_ij, lambdaW,beta_ij);
							//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijm, a_iT,y_j,lambdaW);							
							//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
							//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
							//alpha = -Bijm;
							//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
							//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &inci);
							//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
							//beta 	= pow(norm2,2);
							//LssCands = sigma2*N*log(fabs(d_ij - Bijm)+1e-16) - beta/2;							
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

						//update IBinvZero ith column updated; jth row: 
						//readPtr = &IBinvZero[j];
						//F77_CALL(dscal)(&M,&ziDb,readPtr,&M);
						//IBinvZero[j*M + i] = IBinvZero[j*M + i]*ziDb;
						// adjugate of QIBinv = ziDb*(QIBinvadj - (-1)^kl) + (-1)^kl
						//for(k=0;k<M;k++) // kth row
						//{
						//	for(l=0;l<M;l++)//lth column
						//	{
						//		if(k!=j && l!=i) 
						//		{
						//			alpha = pow(-1.0,k+l);
						//			IBinvZero[l*M+k] = ziDb*(IBinvZero[l*M+k] -alpha) + alpha;
						//		}	
						//	}
						//}


						
					}//js_i >0
				}//j = 1:M	
			
				//f
				//BF1 = B(i,:)*F1(:,i)
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
		
		
		// BLOCK COORDINATE ASCEND: Update IBinv
		UpdateIBinv(IBinvZero, B,M);
		
		if(delta_BF<1e-2)		//break out
		{
			break;
		}

	}//while
		
	
	if(verbose>3) Rprintf("\t number of iteration is: %d.\nExiting constrained_MLf\n",iter);

	//IBinv = I -B
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
	//mue[i] = -f[i]*meanX[i];
	//	//dgemv mue = Ax + beta*mue
	transa = 'N';
	alpha = 1;
	//beta = 1;
	ldk = M;
	F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, IBinv, &ldM, meanY, &inci, &beta,mue, &incj);
	
	
	//----------------------------------------------------------------------------------END OF LASSO_ZERO
	
	
	
	
	//mue                     	= (IM-B)*meanY-bsxfun(@times,f,meanX); diag(f)*meanX

	Free(meanY);
	Free(meanX);
	Free(Y);
	Free(X);
	
	Free(S);
	Free(s);
	Free(f0);
	Free(F1);
	Free(Wcopy);
	
	//Free(xi);
	Free(y_j);
	
	Free(ei);
	Free(IBinv);
	//Free(zzi);
	Free(a_iT);
	
	Free(eiB);
	Free(BiT);
	Free(BfOld);
	Free(BfNew);
	Free(BfChange);
	
	//Free(ipiv);
	
	//Free(aiTQuad);
	//Free(eye);
	return lambda;
	//sigma2 remains the same
}//weighted_LassoSf

//--------------------------------------------------------------------------------END WEIGHTED_LASSOSF

//----------------------------------------------- TEST FUNCTION
int cv_gene_nets_support_Rdg(double *Y, double *X, int Kcv,double *lambda_factors, double *rho_factors, 
			int maxiter, int M, int N,int Nlambdas, int Nrho,int verbose,double *W, 			//double sigma2learnt,
			double *sigmaLasso)		//,double *IBinv
{	//sigma2_ms ilambda_ms
	// 
	// return sigmaLASSO = sigma_RidgeCVf for lasso in main
	//return ilambda_cv_ms: index of lambda_factors
	//0. lambda(i) path:
	//1. 	cross validation by ridge_cvf: rho_factor 		--> ONLY  for ilambda = 1; save Kcv set;
	//2. 	constrained_ridge_eff: weights W
	//3. 	weighted_lassoSf: none-zeros with lambda(i)
	//4. 	constrained_MLf: likelihood
	//5. return lambda
	int MM = M*M;
	int MN = M*N;
	// to save for each lambda
	double *Q, *BL, *fL,*mueL; 
	int i,j,index;	
	Q = (double* ) Calloc(MM, double); //ptr Q
	BL = (double* ) Calloc(MM, double); //ptr BL
	fL = (double* ) Calloc(M, double); //ptr fL
	mueL = (double* ) Calloc(M, double);
	//
	double *BC, *fC;
	BC = (double* ) Calloc(MM, double);
	fC = (double* ) Calloc(M, double);
	
	int Ntest = N/Kcv; 		//C takes floor
	int Nlearn = N - Ntest; 
	double * Errs, *Sigmas2, *ErrorMean, *ErrorCV;
	
	if(Nrho>Nlambdas)
	{
		i = Nrho*Kcv;
	}else
	{
		i= Nlambdas*Kcv;
	}
	
	Errs = (double* ) Calloc(i, double);
	Sigmas2 = (double* ) Calloc(Nlambdas, double);
	ErrorMean= (double* ) Calloc(Nlambdas, double);
	i = Nlambdas*Kcv;
	ErrorCV = (double* ) Calloc(i, double);
	
	//parameters inside the loop
	double *Ylearn, *Xlearn, *Ytest,*Xtest;
	int NlearnM = Nlearn*M;
	int NtestM = Ntest*M;
	Ylearn = (double* ) Calloc(NlearnM, double);
	Xlearn = (double* ) Calloc(NlearnM, double);
	Ytest = (double* ) Calloc(NtestM, double);
	Xtest = (double* ) Calloc(NtestM, double);

	//centerYX function
	double *Ylearn_centered, *Xlearn_centered, *meanY, *meanX;
	Ylearn_centered = (double* ) Calloc(NlearnM, double);
	Xlearn_centered = (double* ) Calloc(NlearnM, double);
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	

	//main loop
	double err_mean;
	//initialize
	err_mean = 1e5;

	int ilambda, cv,testStart, testEnd;
	//min_attained,min_attained = 0;
	ilambda = 0;
	
	//Weighted_LassoSf
	double lambda,lambda_factor,lambda_factor_prev;
	lambda_factor_prev = 1;

	//SL
//	double sigma2Ridge; 	//check which one is better: sigma2Ridge of sigma_hat(lambda)
	double *SL;				//zero or not-zero
	SL = (double* ) Calloc(MM, double);

	//convergence 
	
	if(verbose>1) Rprintf("\t\tEnter Function: cv_support. Nlambdas: %d; \t %d-fold cross validation.\n", Nlambdas,Kcv);
	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	int lda,ldb,ldc,ldk;
	double *Xptr, *XsubPtr; //submatrix
	double alpha, beta;
//Rprintf("Initialized.\n");		


//-------------------------------------------------------------ridge cv
	double sigma2learnt; // sigma of ridge regression for lasso_cv_support
	if(verbose>4) Rprintf("\n\t\t\t\t\tEnter Function: ridge_cvf. %d-fold cross validation.\n", Kcv);
	int irho = 0;
	int min_attained = 0;
	double err_mean_prev = 0;
	double sigma2R,i_rho_factor;
	double *BR, *fR, *mueR;
	BR =&BL[0];
	fR = &fL[0];
	mueR= &mueL[0];//BORROWED; no return: reset in lambda CV
	
	//testNoise, ImB: shared with Ridge & Lasso CV
	double *testNoise,*ImB,*NOISE;
	NOISE =(double* ) Calloc(MN, double);
	//testNoise =(double* ) Calloc(NtestM, double);
	testNoise = &NOISE[0];
	ImB = (double* ) Calloc(MM, double);
	char transa = 'N'; 
	char transb = 'N';
	lda = M;
	ldb = M;
	ldc = M;
	ldk = M;
	double testNorm;
	//result of CV 
	double rho_factor; // no pointers needed
	
	//Ridge CV main loop
	while(irho<Nrho && min_attained ==0)
	{

		i_rho_factor = rho_factors[irho];
		err_mean_prev = err_mean;
		err_mean = 0;
		//cross validation
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>6) Rprintf("\t\t\t\t\t\t\t crossValidation %d/Kcv\n\n",cv);
			// start and end point
			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);
//Rprintf("\t\t\t testStart: %d \t End %d\n",testStart,testEnd);
			//assign submatrices
			//SubMatrix(X,Xtest,Xlearn,M,N,testStart,testEnd);
			Xptr = &X[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Xtest,&incj);
			if(testStart ==1)
			{
				Xptr = &X[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Xlearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &X[testEnd*M]; //index
				XsubPtr = &Xlearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
			}
			
			//SubMatrix(Y,Ytest,Ylearn,M,N,testStart,testEnd);
			Xptr = &Y[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Ytest,&incj);
			if(testStart ==1)
			{
				Xptr = &Y[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Ylearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &Y[testEnd*M]; //index
				XsubPtr = &Ylearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
			}
			
			
			//SubMatrixInt(Missing,Missing_test,Missing_learn,M,N,testStart,testEnd);
		
			
			// ridge SEM			
			sigma2R = constrained_ridge_cff(Ylearn, Xlearn, i_rho_factor, M, Nlearn,BR,fR,mueR,verbose);

			//noise error
			F77_CALL(dcopy)(&MM,BR,&inci,ImB,&incj);
			alpha = -1;
			F77_CALL(dscal)(&MM,&alpha,ImB,&inci);
			for(i=0;i<M;i++) 
			{
				index = i*M + i;
				ImB[index] = 1 + ImB[index];
			} //I-BR
			alpha = 1; 
			beta = 0;			
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				XsubPtr = &testNoise[i];
				//row i of noise
				alpha = -fR[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &ldk,XsubPtr, &M);
				//NOISE[i,1:N]	
			}//row i = 1:M
			
			// y = ax + y  F77_CALL(daxpy)(&MNtest, &alpha,Xptr, &inci,XsubPtr, &incj);
			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				Xptr = &testNoise[i*M];
				F77_CALL(daxpy)(&M, &alpha,mueR, &inci,Xptr, &incj);
			}

			//testNorm = FrobeniusNorm(testNoise, M, N);
			//i = M*Ntest;
			//testNorm = F77_CALL(dnrm2)(&i,testNoise,&inci);
			//testNorm = testNorm*testNorm;
			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);
			//testNorm = F77_CALL(dnrm2)(&NtestM,testNoise,&inci);
			//testNorm = testNorm*testNorm;
			Errs[cv*Nrho+irho] = testNorm;
	
			
			
			err_mean = err_mean + testNorm;
//Rprintf("Errmean: %f; \t sigma2R: %f\n",err_mean, sigma2R);
		}//cv = 0: Kcv
		//err_mean = err_mean/Kcv;   	calculate sum instead
		//check convergence
		
		irho = irho + 1;
		if(verbose>5) Rprintf("\t\t\t\t\t\t Test rho ratio %f; (%d/Nrho)\t error: %f.\n",i_rho_factor,irho, err_mean);	
		if(irho>1)
		{
			if(err_mean_prev<=err_mean)
			{
				min_attained = 1;
				irho = irho -1;
			}
		}
			
	}//while
	
	//sum(sum(1-Missing))
	//int Nmiss = 0;
	//for(i=0;i<MN;i++) Nmiss = Nmiss + Missing[i];
	if(verbose>4) Rprintf("\t\t\t\t\tExit RidgeCV. sigma2R: %f\t",sigma2R);	
	
	//int Npresent = MN- Nmiss;
	if(irho == Nrho || irho == 1) //nonstop in while loop
	{	
		//sigma2ridge             = sum(Errs(irho,:))/(sum(sum(1-Missing))-1);
		sigma2R = err_mean/(MN -1);
	}else
	{
		sigma2R = err_mean_prev/(MN -1);
	}
	//rho_factor_m            = rho_factors(irho)*N/(N-Ntest);
	rho_factor = rho_factors[irho-1]*N/(N-Ntest);
	if(verbose>4) Rprintf("sigma2learnt: %f\n",sigma2R);	
	//sigma2R is for LASSO CV;
	//rho_factor_m: for ridge regression: use the whole dataset
	if(verbose==0) Rprintf("Step 2: ridge CV; find rho : %f\n", rho_factor);
	sigma2learnt = constrained_ridge_cff(Y, X, rho_factor, M, N,BR,fR,mueR,verbose);	
	if(verbose==0) Rprintf("Step 3: ridge; calculate weights.\n");
	// weight W: for CV(this function) and selection path (Main)
	for(i=0;i<MM;i++) W[i] = 1/fabs(BL[i]+ 1e-10);
	
	//sigma2learnt = sigma2cr; //for lambda CV
	//sigma2R is used for path in Main
	sigmaLasso[0] = sigma2R;
	
//------------------------------------------------------------ridge cv end return sigma2R, weight	
	//IBinv: update in path
	double *IBinv,*IBinvZero,*lambda_Max;
	double *IBinvPath,*IBinvPathZero;
	irho  = MM*Kcv;
	IBinvPath = (double* ) Calloc(irho, double);
	IBinvPathZero = (double* ) Calloc(irho, double);
	lambda_Max = (double* ) Calloc(Kcv, double);
	double lambda_max_cv;
	beta = 0;
	//initialized
	
	for(cv=0;cv<Kcv;cv++)
	{
		IBinv = &IBinvPath[cv*MM];
		F77_CALL(dcopy)(&MM,&beta,&inc0,IBinv,&incj);
		for(i=0;i<M;i++) IBinv[i*M+i] =1;
	}
	F77_CALL(dcopy)(&irho,IBinvPath,&inci,IBinvPathZero,&incj);
	
	while(ilambda < Nlambdas)
	{	
		//ilambda = ilambda  +1;
		//err_mean_prev = err_mean;
		//err_sigma_prev = err_sigma;
		err_mean = 0;
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>3) Rprintf("\t\t %d/Kcv cross validation.\n", cv);
			// test, learn
			// start and end point
			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);
			//assign submatrices
			//SubMatrix(X,Xtest,Xlearn,M,N,testStart,testEnd);
			//SubMatrix(Y,Ytest,Ylearn,M,N,testStart,testEnd);
			//SubMatrixInt(Missing,Missing_test,Missing_learn,M,N,testStart,testEnd);
			//SubMatrix(X,Xtest,Xlearn,M,N,testStart,testEnd);
			Xptr = &X[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Xtest,&incj);
			if(testStart ==1)
			{
				Xptr = &X[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Xlearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &X[testEnd*M]; //index
				XsubPtr = &Xlearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
			}
			
			//SubMatrix(Y,Ytest,Ylearn,M,N,testStart,testEnd);
			Xptr = &Y[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Ytest,&incj);
			if(testStart ==1)
			{
				Xptr = &Y[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Ylearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &Y[testEnd*M]; //index
				XsubPtr = &Ylearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
			}
			
			
			//SubMatrixInt(Missing,Missing_test,Missing_learn,M,N,testStart,testEnd);
			//Missing_test = &Missing[(testStart-1)*M];
			//Learn matrix
			
			//Ylearn_centered
			F77_CALL(dcopy)(&NlearnM,Xlearn,&inci,Xlearn_centered,&incj);
			F77_CALL(dcopy)(&NlearnM,Ylearn,&inci,Ylearn_centered,&incj);

			centerYX(Ylearn_centered,Xlearn_centered,meanY, meanX,M, Nlearn);
			//first ilambda

//Rprintf("Initialized.\n");			
			if(ilambda == 0)
			{
				//if(verbose>3) Rprintf("\t\t\t step 0 for first lambda: Ridge cross validation; \tRidge weights.\n");

				//for(i=0;i<M;i++) fL[i] = 1;
				alpha = 1; 
				F77_CALL(dcopy)(&M,&alpha,&inc0,fL,&incj); // call dcopy(n, x, inci, y, incy)
				alpha = 0;
				F77_CALL(dcopy)(&MM,&alpha,&inc0,BL,&incj); 
				//F77_CALL(dscal)(&MM,&alpha,BL,&inci);
				// Q(lambda)S4
				QlambdaStart(Ylearn_centered,Xlearn_centered, Q, sigma2learnt,M, Nlearn);
				// set BL, fL to  zeros: not necessary; will test this after Debug
			
				lambda_Max[cv]	= lambdaMax(Ylearn_centered,Xlearn_centered,W,M, Nlearn);
			}//ilambda ==0
			lambda_max_cv = lambda_Max[cv];
			//if(ilambda > 0) sigma2learnt = Sigmas2[cv*Nlambdas];//in Juan: 5 sigma2learnt saved: one for each fold
			//Weighted_LassoSf
			lambda_factor = lambda_factors[ilambda];
			//SIGMA2[0] = sigma2learnt;
			
			//lambda = Weighted_LassoSf(W, BL, fL, Ylearn, Xlearn, Q, lambda_factor,
			//				lambda_factor_prev, sigma2learnt, maxiter, M, Nlearn, verbose,IBinv);
							
			IBinv = &IBinvPath[cv*MM];
			IBinvZero= &IBinvPathZero[cv*MM];
			lambda = Weighted_LassoSf_MLf(W, BL, fL, Ylearn,Xlearn, Q, lambda_factor, 
							lambda_factor_prev, sigma2learnt, maxiter,M, Nlearn, verbose,
							BC, fC, mueL,IBinv,IBinvZero,lambda_max_cv);
//if(ilambda==0 && cv==0)
//{
//	Rprintf("In cv support: \n");
//	printMat(IBinv,M,M);
//	printMat(fL,1,M);
//printMat(BL,M,M);
	//printMat(fC,1,M);
	//printMat(mueL,1,M);
//}

							
			if(verbose>3) Rprintf("\t\t\t step 1 SML lasso regression, lambda: %f.\n",lambda);
			//sign of B SL		//SL int:	
			//F77_CALL(dscal)(&MM,&alpha,SL,&inci);

			//Q(lambda)
			//QlambdaMiddle(Ylearn,Xlearn, Q,BL,fL, mueL, sigma2learnt,M, Nlearn);
			QlambdaMiddleCenter(Ylearn_centered,Xlearn_centered, Q,BL,fL,sigma2learnt,M, Nlearn,IBinv);
//if(ilambda==0 && cv==0)
//{
	//Rprintf("In cv support: \n");
	//printMat(Q,M,M);
//	printMat(fL,1,M);
	//printMat(BC,M,M);
	//printMat(fC,1,M);
	//printMat(mueL,1,M);
//}			
			
			// constrained_MLf
			if(verbose>3) Rprintf("\t\t\t step 2 SML ZeroRegression.\n");
			

			//constrained_MLf(BC, fC, SL, Ylearn,Xlearn, sigma2learnt, maxiter, mueL,M, Nlearn, verbose);

			//noise error Errs[ilambda, cv]
			
			//error function: norm((1-Missing).*(A*Y-fC*X)-mueC),'fro')^2;
//Problem in JUAN'S CODE: FOR CV mean: need to normalize, ow. it is not a standard error.
//Errs[index] = errorFunc(BC,Ytest, fC, Xtest, mueL, M, Ntest);
			
			F77_CALL(dcopy)(&MM,BC,&inci,ImB,&incj);
			alpha = -1;
			F77_CALL(dscal)(&MM,&alpha,ImB,&inci);
			for(i=0;i<M;i++) 
			{
				index = i*M + i;
				ImB[index] = 1 + ImB[index];
			} //I-BR
			alpha = 1; 
			beta = 0;
			testNoise = &NOISE[NtestM*cv];
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				//XsubPtr = &testNoise[i];
				XsubPtr = &NOISE[NtestM*cv+i];
				//F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
				//row i of noise ldM
				alpha = -fC[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &M,XsubPtr, &ldk);
				//NOISE[i,1:N]	
			}//row i = 1:M
			
			// y = ax + y  F77_CALL(daxpy)(&MNtest, &alpha,Xptr, &inci,XsubPtr, &incj);
			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				//Xptr = &testNoise[i*M];
				Xptr = &NOISE[NtestM*cv+i*M];
				F77_CALL(daxpy)(&M, &alpha,mueL, &inci,Xptr, &incj);
			}

			//testNorm = FrobeniusNorm(testNoise, M, N);
			// TEST_NOISE computed;
			
			
			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);
			//testNorm = F77_CALL(dnrm2)(&NtestM,testNoise,&inci); //just dotproduct will work
			//testNorm = testNorm*testNorm;
			index = cv*Nlambdas+ilambda;
			
			Errs[index] = testNorm;
	
	
//			Sigmas2[index] = sigma2learnt;
			ErrorCV[index] = Errs[index]/NtestM;
			err_mean = err_mean + Errs[index];
			if(verbose>3) Rprintf("\t\t\t cv: %d \tend; err_mean = %f.\n", cv, Errs[index]);
				
		}//cv
		//err
		err_mean = err_mean/MN;
		ErrorMean[ilambda] = err_mean;
		
		//mean, sigma
		for(i=0;i<MN;i++)
		{
			NOISE[i] = pow(NOISE[i],2);
		}
		alpha = -1;
		F77_CALL(daxpy)(&MN,&alpha,&err_mean,&inc0,NOISE,&inci);
		testNorm = F77_CALL(dnrm2)(&MN,NOISE,&inci);
		Sigmas2[ilambda] = testNorm/sqrt(Kcv*(MN -1));
		
		
/*		if(minimumErr>err_mean) 
		{
			minimumErr = err_mean;
			ilambda_min = ilambda;
		}
*/		//check for convergence
		//ilambda = ilambda + 1; //don't move this line
		//if(ilambda>3 && min_attained==0)
		//{
		//	if ((err_mean_prev + err_sigma_prev) < err_mean)
		//	{
		//		min_attained = 1;
		//		ilambda_min = ilambda; //number; not for indexing
		//	}
		//}
		lambda_factor_prev = lambda_factor;
		
		if(verbose>2) Rprintf("\t\t\t %d/Nlambdas. %d fold cv; \t Err_Mean: %f; std:%f; \t sigma2learnt:%f.\n", ilambda,Kcv,err_mean,Sigmas2[ilambda],sigma2learnt);
		ilambda = ilambda + 1; 
	}//while ilambda

	// select ilambda_ms
	
	// step1. calculate sigma2
	//ErrorCV = ErrorCV - ErrorMean
//printMat(Sigmas2,Nlambdas,1);	
	// step2. index of minimal ErrorMean
	double *errorMeanCpy;
	errorMeanCpy = (double* ) Calloc(Nlambdas, double);
	//int inc0  = 0;
	double minimumErr = -1e5 - MN;
	F77_CALL(dcopy)(&Nlambdas,&minimumErr,&inc0,errorMeanCpy,&inci);
	alpha = 1;
	F77_CALL(daxpy)(&Nlambdas,&alpha,ErrorMean,&inci,errorMeanCpy,&incj); // y = ax + y
	int ilambda_ms;	
	ilambda_ms = F77_CALL(idamax)(&Nlambdas, errorMeanCpy, &inci);//index of the max(errs_mean)<--min(mean)
	index = ilambda_ms - 1;
//	minimumErr = ErrorMean[index] + Sigmas2[index]; //actually max
//printMat(ErrorMean,1,Nlambdas);
//printMat(Sigmas2,1, Nlambdas);


	//double lowBound;
//	for(i=index-1;i>0;i--) 
//	{
//		if(ErrorMean[i] < minimumErr)
//		{
//			ilambda_ms = i + 1;
//		}
//	}
	//return index of lambda_factors
	if(verbose>1) Rprintf("\t\tExit Function: cv_support. optimal lambda index: %d.\n\n", ilambda_ms);

//	Free(Ws);
	Free(NOISE);
	Free(ImB);
	
	
	Free(Q);
	Free(BL);
	Free(fL);
	Free(mueL);
	Free(BC);
	Free(fC);
	
	Free(Errs);
	Free(Sigmas2);
	Free(ErrorMean);
	Free(ErrorCV);
	Free(errorMeanCpy);
	
	Free(Ylearn);
	Free(Xlearn);
	Free(Ytest);
	Free(Xtest);	

	//Free(Missing_learn);
	//Free(Missing_test);
	Free(Ylearn_centered);
	Free(Xlearn_centered);	
	Free(meanY);
	Free(meanX);
	Free(SL);
//	Free(errs_mean);	
//	Free(errs_sigma);	

//	Free(rho_factor_ptr);


	Free(IBinvPath);
	Free(IBinvPathZero);
	
	Free(lambda_Max);
	
	return ilambda_ms;
	//index + 1 -->position; return position
	
	
	
} //end of cv_gene_nets_support		

			
			
			
//cv_gene_nets_support: return ilambda_cv_ms Missing


// ---------------------------- Main function that calls by R (interface)
//-----------------------------main function calls:
// centerYX
// cv_gene_nets_support
// ridge_cvf
// constraint - ridge_cff
// weighted_lassoSf


//--------------------------------not elastic at this moment
void mainSML(double *Y, double *X, int *m, int *n, int *Missing,double*B, double *f,double*stat,int*VB)
{
	//stat contains: correct postive, true positve; false postive, positve detected; power; fdr 6x1
	// assume B is the true imput
	int M, N, i, j,index,verbose;
	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	M 			= m[0];
	N 			= n[0];	
	verbose 	= VB[0];
	int MN 		= M*N;
	int MM 		= M*M;
	double *Strue;
	Strue 		= (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,B,&inci,Strue,&incj);
	stat[1] 	= 0;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index = j*M  +i;
			//Strue[index] = B[index];
			if(i!=j && B[index]!=0)
			{	stat[1] = stat[1] + 1;} //stat[1] total positive
		}
	}
	double alpha = 1;	
	F77_CALL(dcopy)(&M,&alpha,&inc0,f,&inci);
	//for(i=0;i<M;i++) f[i] = 1;

	alpha = 0;
	//F77_CALL(dscal)(&MM,&alpha,B,&inci);
	F77_CALL(dcopy)(&MM,&alpha,&inc0,B,&inci);
	//assume missing values are from X		//set missing counterpart in Y to zero

	for(i=0;i<MN;i++)
	{
		if(Missing[i] == 1) Y[i] = 0;
	}
	
	//call cv_gene_nets_support ------------------------SYSTEM PARAMETERS
	int maxiter 	= 500;
	int Kcv 		= 5;
	int L_lambda 	= 20; // number of lambdas in lambda_factors	stop at 0.001
	double *lambda_factors;
	lambda_factors 	= (double* ) Calloc(L_lambda, double);
	double step 	= -0.2;
	for(i=0;i<L_lambda;i++) 
	{
		lambda_factors[i] 	= pow(10.0,step);
		step 				= step - 0.2;
	}
	
	//rho factor
	step 					= -6;
	double * rho_factors;
	int L 					= 31; // number of rho_factors
	rho_factors 			= (double* ) Calloc(L, double);
	for(i=0;i<L;i++) 
	{
		rho_factors[i] 		= pow(10.0,step);
		step 				= step + 0.2;
	}
	
	//---------------------------------------------END  SYSTEM PARAMETERS
	int Nlambdas,Nrho;
	Nlambdas 				= L_lambda;
	Nrho 					= L;
//call ridge_cvf
	double sigma2; //return value;

	//double *mueL;
	//mueL = (double* ) Calloc(M, double); 

	// weight W
	double *W;  //weight on diagonal?
	W = (double* ) Calloc(MM, double);
	
	// IBinv: SAVE COMPUTATION: get IBinv from lasso: 1) CV_support; 2) selection path
	double *QIBinv;
	//IBinv = (double* ) Calloc(MM, double);
	QIBinv = (double* ) Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1
	
	//ilambda_cv_ms = cv_gene_nets_support(Y, X, Kcv,lambda_factors, 
	//				rho_factors, maxiter, M, N,Nlambdas, Nrho,verbose,W,sigma2cr);
	
	ilambda_cv_ms = cv_gene_nets_support_Rdg(Y, X, Kcv,lambda_factors, rho_factors, 
			maxiter, M, N,Nlambdas, Nrho,verbose,W, &sigma2);
	
	//ilambda_cv_ms = 17;
	//sigma2 = 0.0156827270598132;
	//step = 1;
	//F77_CALL(dcopy)(&MM,&step,&inc0,W,&inci);
	
	
	
	if(verbose==0) Rprintf("Step 1: CV support; return number of lambda needed: %d\n", ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	Ycopy = (double* ) Calloc(MN, double);
	Xcopy = (double* ) Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda; 
	double lambda_max;
	lambda_max = lambdaMax(Ycopy,Xcopy,W,M, N);
	//
	//for(i=0;i<M;i++) f[i] = 1;
	//F77_CALL(dscal)(&MM,&alpha,B,&inci);
	if(verbose==0) Rprintf("Step 4: lasso selection path.\n");
	//printMat(B,M,M);
	//ilambda_cv_ms = 0;
//printMat(IBinv,M,M);	
	for(ilambda = 0;ilambda<ilambda_cv_ms;ilambda++)
	//for(ilambda = 0;ilambda<1;ilambda++)
	{
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, ilambda_cv_ms,lambda_factor);
		// call Weighted_LassoSf		Missing,
		lambda = Weighted_LassoSf(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max); 	// mueL not calculated
//printMat(IBinv,M,M);	
if(verbose>0) Rprintf("\tlambda: %f\n", lambda);						
		// call QlambdaMiddle
		//QlambdaMiddle(Y,X, Q,B,f, mueL, sigma2,M,N);
		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];
//printMat(IBinv,M,M);
	}//ilambda; selection path
	//return B,f; 
	//correct positive
	//printMat(B,M,M);
	stat[0] = 0;// correct positive
	stat[2] = 0;//false positive
	stat[3] = 0;//positive detected
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index = j*M + i;
			if(Strue[index]==0 && B[index]!=0) stat[2] = stat[2] + 1;
			if(i!=j)
			{
				//stat[3]
				if(B[index]!=0) 
				{
					stat[3] = stat[3] + 1;
					//stat[0]
					if(Strue[index]!=0) stat[0] = stat[0] + 1;
				}
			}
		}
	}
	//power
	stat[4] = stat[0]/stat[1];
	stat[5] = stat[2]/stat[3];
	if(verbose==0) Rprintf("Step 5: Finish calculation; detection power in stat vector.\n");
	Free(Strue);
	Free(meanY);
	Free(meanX);
	Free(lambda_factors);
	Free(rho_factors);
	Free(Ycopy);
	Free(Xcopy);
	//Free(rho_factor_m);
	//Free(mueL);
	Free(W);
	//Free(IBinv);
	Free(QIBinv);
	Free(Q);
	//-------------------------------- some function changes Y, X permantly.
}


//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------
//version v1: elastic net penalty: converted from lassoSMLv10beta.c
/*#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <stdio.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <stdlib.h>

*/
//version Note: block coordinate ascent: after each block updated: update IBinv before next iteration
//package SEM_SML: comment out same functions in two C files
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

void centerYX(double *Y,double *X, double *meanY, double *meanX,int M, int N) //M genes; N samples
{
	//matrix is vectorized by column: 1-M is the first column of Y
	//missing values are from X; set corresponding Y to zero.  in main_SMLX.m

	int i,index;	
	double *Xptr;
	double *Yptr;

	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	int lda  = M; //leading dimension
	double *eye;
	eye = (double* ) Calloc(N, double);
	double alpha = 1;
	double beta = 0;
	F77_CALL(dcopy)(&N,&alpha,&inc0,eye,&inci);
	char transa = 'N';

	F77_CALL(dgemv)(&transa, &M, &N,&alpha, X, &lda, eye, &inci, &beta,meanX, &incj);
	F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &lda, eye, &inci, &beta,meanY, &incj);
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
	Free(eye);
}	

//--------------------- LINEAR SYSTEM SOLVER END ---------

//ridge regression; return sigma2learnt
double constrained_ridge_cff(double *Ycopy, double *Xcopy, double rho_factor, int M, int N,
		double *B, double *f, double *mue, int verbose)
{
	
	int i,j,k,lda,ldb,ldc,ldk;
	// center Y, X
	double *meanY, *meanX;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	
	//copy Y, X; 
	double *Y, *X;
	int MN = M*N;
	Y = (double* ) Calloc(MN, double);
	X = (double* ) Calloc(MN, double);
	
	//F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
		//double *dy, const int *incy);
	int inci = 1;
	int incj = 1;
	F77_CALL(dcopy)(&MN,Ycopy,&inci,Y,&incj);
	F77_CALL(dcopy)(&MN,Xcopy,&inci,X,&incj);
	
	centerYX(Y,X,meanY, meanX,M, N);
	
//	double NPresent = 0;
//	for(i=0;i<M;i++)
//	{
//		for(j=0;j<N;j++)
//		{
//			if(Y[j*M + i]!=0) NPresent = NPresent + 1; //Y[i,j]
//		}
//	}	
	if(verbose>7) Rprintf("\t\t\t\t\t\t\t\tEnter Function: Ridge Regression. Shrinkage ratio rho is: %f.\n\n",rho_factor);
	
	int Mi = M -1;
	//for usage in loop
	double *YiPi; //Yi'*Pi
	YiPi =(double* ) Calloc(Mi*N, double);
	double xixi,xixiInv; //xi'*xi;
	int jj,index; //jj = 1:(M-1) index of YiPi
	double normYiPi,rho;
	double *bi,*YiPi2Norm; 	//YiPi2Norm: first term of biInv;
	
	double *Hi,*Yi,*xi,*yi,*xii;//xii for Hi calculation Hi= xi*xi'
	Hi = (double* ) Calloc(N*N, double);
	Yi =(double* ) Calloc(Mi*N, double);
	xi = (double* ) Calloc(N, double);
	xii = (double* ) Calloc(N, double);
	yi = (double* ) Calloc(N, double);
	double alpha, beta;
	char transa = 'N';
	char transb = 'N';
	
	//
	int MiMi = Mi*Mi;
	int NN = N*N;
	YiPi2Norm 	= (double* ) Calloc(MiMi, double);	
	bi 			= (double* ) Calloc(Mi, double);
	//YiPiyi 		= (double* ) Calloc(Mi, double);
	
	//bi,fi
	double *xiYi; //xi*Yi
	xiYi = (double* ) Calloc(Mi, double);
	double xiYibi, xiyi;
	//main loop:
//Rprintf("check point 1: before loop\n");
	alpha = 1;
	beta = 0;
		
//largest Eigenvalue
	double *biInv;
	biInv 		= (double* ) Calloc(MiMi, double); //copy of YiPi2Norm
	//dsyevd
	char jobz = 'N'; // yes for eigenvectors
	char uplo = 'U'; //both ok
	double *w, *work;
	w = (double *) Calloc(Mi,double);
	int lwork = 5*Mi + 10;
	work  = (double *) Calloc(lwork,double);	
	int liwork = 10;
	int *iwork;
	iwork = (int *) Calloc(liwork,int);
	int info = 0;
	//dsyevd function 

	//linear system
	int *ipiv;
	//ipiv = (int *) R_alloc(N,sizeof(int));
	ipiv = (int *) Calloc(Mi,int);
	double *readPtr,*readPtr2;
	//loop starts here
	for(i=0;i<M;i++)
	{
		//xi = X[i,:]
		readPtr = &X[i];
		F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
		F77_CALL(dcopy)(&N,xi,&inci,xii,&incj);
		readPtr = &Y[i];
		F77_CALL(dcopy)(&N,readPtr,&M,yi,&inci);

		//xixi = F77_CALL(dnrm2)(&N,xi,&inci);
		//xixi = pow(xixi,2);
		xixi = F77_CALL(ddot)(&N, xi, &inci,xi, &incj);
		xixiInv = -1/xixi;
//printMat(xi,1,N);		
		//xi'*xi
		
		//YiPi
		//Hi          = xi*xi'/(xi'*xi);
        //Pi          = eye(N)-Hi;
//Rprintf("check point 2: xixi: %f.\n",xixi);

		//MatrixMult(xi,xi, Hi,alpha, beta, N, k, N);
		transb = 'N';
		lda = N;
		ldb = N;
		ldc = N;
		F77_CALL(dgemm)(&transa, &transb,&N, &ldb, &inci,&alpha, xi,&lda, xii, &incj, &beta,Hi, &ldc);
		
		
		//F77_NAME(dscal)(const int *n, const double *alpha, double *dx, const int *incx);
		//k= N*N;
		F77_CALL(dscal)(&NN,&xixiInv,Hi,&inci); // Hi = -xi*xi'/(xi'*xi);
		for(j=0;j<N;j++) 
		{	index = j*N + j;
			Hi[index] = Hi[index] + 1;
		}//Pi
//printMat(Hi,N,N);	
		
	
		//Yi
		readPtr2 = &Yi[0];
		jj = 0;
		for(j=0;j<M;j++)
		{	if(j!=i)
			{
				//copy one j row
				readPtr = &Y[j];
				F77_CALL(dcopy)(&N,readPtr,&M,readPtr2,&Mi);
				jj = jj + 1;
				readPtr2 = &Yi[jj];
			}
		}//jj 1:(M-1), YiPi[jj,:]
		//YiPi=Yi*Pi
//printMat(Yi,Mi,N);		
		//transb = 'N';
		lda = Mi;
		ldb = N;
		ldc = Mi;
		ldk = N; //b copy
		F77_CALL(dgemm)(&transa, &transb,&Mi, &N, &ldk,&alpha, Yi, &lda, Hi, &ldb, &beta, YiPi, &ldc);
		//F77_CALL(dgemm)(&transa, &transb,&Mi, &N, &N,&alpha, Yi, &Mi, Hi, &N, &beta, YiPi, &Mi);
		
		//Rprintf("check point 3: YiPi\n");
		
		//YiPi*Yi' --> MixMi
		transb = 'T';
		ldk = Mi;
		lda = Mi;
		ldb = Mi;
		ldc = Mi;
		F77_CALL(dgemm)(&transa, &transb,&Mi, &ldk, &N,&alpha, YiPi, &lda, Yi, &ldb, &beta, YiPi2Norm, &ldc);
		//F77_CALL(dgemm)(&transa, &transb,&Mi, &Mi, &N,&alpha, YiPi, &Mi, Yi, &Mi, &beta, YiPi2Norm, &Mi); //M-> N -> K
		//2-norm YiPi this is the largest eigenvalue of (Yi'*Pi)*YiPi
		//Matrix 2-norm;
		//Repeat compution; use biInv;
		//normYiPi = Mat2NormSq(YiPi,Mi,N); //MixN
		//YiPi2Norm
		
		//Rprintf("check point 4: YiPi2Norm\n");
		
//normYiPi = largestEigVal(YiPi2Norm, Mi,verbose);  //YiPi2Norm make a copy inside largestEigVal, YiPi2Norm will be intermediate value of biInv;
		//transa = 'N';
		transb = 'N';
		//j = Mi*Mi;
		F77_CALL(dcopy)(&MiMi,YiPi2Norm,&inci,biInv,&incj);
		
		//F77_CALL(dgeev)(&transa, &transb,&Mi, biInv, &Mi, wr, wi, vl, &ldvl,vr, &ldvr, work, &lwork, &info);
		lda = Mi;
		F77_CALL(dsyevd)(&jobz, &uplo,&Mi, biInv, &lda, w, work, &lwork, iwork, &liwork,&info);
		normYiPi = w[Mi -1];
//printMat(w,1,Mi);		
//Rprintf("Eigenvalue: %f.\n",normYiPi);		
		rho = rho_factor*normYiPi; // 2Norm = sqrt(lambda_Max)
		
		if(verbose>8) Rprintf("\t\t\t\t\t\t\t\t\t Gene number: %d,\t shrinkage rho: %f\n",i,rho);
		//biInv = (YiPi*Yi+rho*I) ; (M-1) x (M-1)

		//
		for(j=0;j<Mi;j++) 
		{
			index = j*Mi + j;
			YiPi2Norm[index] = YiPi2Norm[index] + rho;
		}
		//biInv;

		//Inverse
		//MatrixInverse(biInv,Mi);
		//Rprintf("check point 5: biInv inversed\n");
		//YiPiyi = Yi'*Pi*yi  = YiPi *yi;
		//F77_NAME(dgemv)(const char *trans, const int *m, const int *n,
		//const double *alpha, const double *a, const int *lda,
		//const double *x, const int *incx, const double *beta,
		//double *y, const int *incy);
		lda = Mi;
		F77_CALL(dgemv)(&transa, &Mi, &N,&alpha, YiPi, &lda, yi, &inci, &beta,bi, &incj);
		

//linearSystem(YiPi2Norm,Mi,bi);//A NxN matrix
		lda = Mi;
		ldb = Mi;
		F77_CALL(dgesv)(&Mi, &inci, YiPi2Norm, &lda, ipiv, bi, &ldb, &info);
		//Rprintf("check point 6: bi updated\n");
		//------------------------------------------------Ridge coefficient beta obtained for row i
		// f(i)        = (xi'*yi-xi'*Yi*bi)/(xi'*xi);
		//xiYi (M-1) x1
		lda = Mi;
		
		F77_CALL(dgemv)(&transa, &Mi, &N,&alpha, Yi, &lda, xi, &inci, &beta,xiYi, &incj);
		
		//xiyi = xi*yi 	= X[i,j]*Y[i,j]
		//dot product 
		xiyi = F77_CALL(ddot)(&N, xi, &inci,yi, &incj);
		
		//xiYibi = xiYi*bi
		xiYibi = F77_CALL(ddot)(&Mi, xiYi, &inci,bi, &incj);

		f[i] = (xiyi-xiYibi)/xixi;
		//Rprintf("check point 7: fi calculated\n");
		//update B
		jj = 0;
		for(j = 0;j<M;j++)
		{
			if(j!=i)
			{
				//B[i,j] = bi[jj];
				B[j*M+i] = bi[jj];
				jj = jj +1;
			}
		}
		
		//end of 1:M		
	}//i = 1:M

	//I -B
	double *ImB;
	k = M*M;
	ImB = (double* ) Calloc(k, double);
	F77_CALL(dcopy)(&k,B,&inci,ImB,&incj);
	//F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
	//	double *dy, const int *incy);
	xixiInv = -1;
	F77_CALL(dscal)(&k,&xixiInv,ImB,&inci);
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		ImB[index] = 1 + ImB[index];
	}
	
	
	
	
	//noise, sigma2learnt,mue;
	double * NOISE; 	//MxN
	NOISE =(double* ) Calloc(MN, double);
	transb = 'N';
	ldk = M;
	lda = M;
	ldb = M;
	ldc = M;
	F77_CALL(dgemm)(&transa, &transb,&M, &N, &ldk,&alpha, ImB, &lda, Y, &ldb, &beta, NOISE, &ldc);//(I-B)*Y - fX
	for(i=0;i<M;i++)
	{
		// row i of X
		readPtr2 = &X[i];
		//F77_CALL(dcopy)(&N,readPtr2,&M,xi,&inci);
		//row i of noise
		readPtr = &NOISE[i];
		alpha = -f[i];
		//F77_CALL(daxpy)(&N, &alpha,xi, &inci,readPtr, &M);
		F77_CALL(daxpy)(&N, &alpha,readPtr2, &ldk,readPtr, &M);
		//NOISE[i,1:N]	
	}//row i = 1:M
	
	
	//NoiseMatF77(NOISE,ImB,Y, f, X, M, N);
		
	double noiseNorm, sigma2learnt;
	//noiseNorm = FrobeniusNorm(NOISE, M, N); //NOISE  = sparse(speye(M)-B)*Y-bsxfun(@times,f,X);
	//noiseNorm = F77_CALL(dnrm2)(&MN,NOISE,&inci);
	//sigma2learnt = noiseNorm*noiseNorm/(MN -1); //sigma2learnt    = sum(sum(NOISE.^2))/(sum(NPresent)-1);
	noiseNorm = F77_CALL(ddot)(&MN, NOISE, &inci,NOISE, &incj);
	sigma2learnt = noiseNorm/(MN -1);
	
	//mue          = (IM-B)*meanY-bsxfun(@times,f,meanX);
	//dgemv mue = Ax + beta*mue
	

	
	
	
	//mue[i] = -f[i]*meanX[i];
	for(i=0;i<M;i++)
	{
		mue[i] = -f[i]*meanX[i];
	}
	beta = 1;
	ldk = M;
	lda = M;
	alpha = 1;
	F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, ImB, &lda, meanY, &inci, &beta,mue, &incj);
	
	
	if(verbose>7) Rprintf("\t\t\t\t\t\t\t\tExit function: Ridge Regression. sigma^2 is: %f.\n\n",sigma2learnt);
	
	Free(meanY);
	Free(meanX);
	Free(Y);
	Free(X);
	Free(YiPi);
	//Free(biInv);
	Free(YiPi2Norm);
	Free(bi);	
	//Free(YiPiyi);
	Free(xiYi);
	Free(NOISE);
	//
	Free(Hi);
	Free(Yi);
	Free(xi);
	Free(yi);
	//
	Free(ImB);
	
	//
	Free(biInv);

	Free(w);
	Free(iwork);
	Free(work);
	
	Free(ipiv);
	return sigma2learnt;

}

*/
//by Weighted_LassoSf xi
//lambda_max          = max(max(abs(N*sigma2*IM - (Y*(Y'-X'*DxxRxy)))./W  ));
double lambdaMax_adaEN(double *Y,double *X,double * Wori,int M, int N,double alpha_factor) 	//------adaEN Apr. 2013	
{	
	// Oct 08, 2012: assume one eQTL for each gene; This fucntion needs significant revision if this assumption doesnot hold
	double *dxx, *rxy, *DxxRxy,*readPtr1,*readPtr2;
	double lambda_max = 0;		
	dxx				= (double* ) Calloc(M, double);
	rxy				= (double* ) Calloc(M, double);
	DxxRxy			= (double* ) Calloc(M, double);
	int i,k,index,lda;
	int inci = 1;
	int incj = 1; 
	lda = M;
	int MN = M*N;
	int MM = M*M;
	//------adaEN Apr. 2013	
	double *W;
	W 				= (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,Wori,&inci,W,&incj);
	F77_CALL(dscal)(&MM,&alpha_factor,W,&inci); 
	//------adaEN Apr. 2013	
	//for(i=0;i<M;i++)
	//{
	//	dxx[i] 		= 0;
	//	rxy[i] 		= 0;
	//	for(j=0;j<N;j++)
	//	{
	//		index = j*M+i;
	//		dxx[i] 	= dxx[i] + pow(X[index],2);	// sum of each row X[i,j]
	//		rxy[i] 	= rxy[i] + X[index]*Y[index];			
	//	}
	//	DxxRxy[i] 	= rxy[i]/dxx[i];		
	//}
	for(i=0;i<M;i++)
	{
		readPtr1 	= &X[i]; //ith row
		readPtr2 	= &Y[i];

		//norm  		= F77_CALL(dnrm2)(&N,readPtr1,&M);	
		//dxx[i] 		= pow(norm,2);
		dxx[i] = F77_CALL(ddot)(&N,readPtr1,&lda,readPtr1,&M);
		//res = ddot(n, x, incx, y, incy)
		rxy[i] 		= F77_CALL(ddot)(&N,readPtr1,&lda,readPtr2,&M);
		DxxRxy[i] 	= rxy[i]/dxx[i];		
	}
	
	
	//abs(N*sigma2*IM - (Y*(Y'-X'*DxxRxy)))./W ; W[i,i] = inf.
	//printMat(DxxRxy,M,1);

	//cache X[k,:]*DxxRxy[k]
	double * XDxxRxy;

	XDxxRxy = (double* ) Calloc(MN, double);
	//for(i=0;i<M;i++)
	//{
	//	for(j=0;j<N;j++)
	//	{
			//X[i,j] * DxxRxy[i];
	//		index = j*M + i;
	//		XDxxRxy[index] = -X[index]*DxxRxy[i];
	//	}
	//}
	F77_CALL(dcopy)(&MN,X,&inci,XDxxRxy,&incj);
	double alpha;	
	for(i=0;i<M;i++)
	{
		alpha  = -DxxRxy[i];
		readPtr1 = &XDxxRxy[i]; //ith row
		F77_CALL(dscal)(&N,&alpha, readPtr1,&M);//	(n, a, x, incx)
	}
	
	
	//printMat(XDxxRxy,M,1);
	// Y- XDxxRxy  			daxpy(n, a, x, incx, y, incy) y= ax + y
	alpha  = 1.0;
	// XDxxRxy <- alpha*Y + XDxxRxy
	F77_CALL(daxpy)(&MN,&alpha,Y,&inci,XDxxRxy,&inci);
	//printMat(XDxxRxy,M,1);
	double *YYXDR; //= Y*XDxxRxy'

	YYXDR = (double* ) Calloc(MM, double);	
	
//C := alpha*op(A)*op(B) + beta*C,
	double beta;
	char transa = 'N';
	char transb = 'T';
	alpha = -1;
	beta = 0;
	F77_CALL(dgemm)(&transa, &transb,&M, &M, &N,&alpha, Y,&M, XDxxRxy, &M, &beta,YYXDR, &M); //M xK, K xN  --> MxN, N xM --> M <-M, N<-M, k<-N
	//diagonal -->0; other element /wij
	
	//printMat(YYXDR,M,3);
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
	//printMat(YYXDR,M,3);
	//BLAS_extern int    /* IDAMAX - return the index of the element with max abs value */
	//F77_NAME(idamax)(const int *n, const double *dx, const int *incx);
	index = F77_CALL(idamax)(&MM,YYXDR,&inci);
	//Rprintf("index: %d\n",index);
	lambda_max = fabs(YYXDR[index-1]);

	Free(dxx);
	Free(rxy);
	Free(DxxRxy);
	//Free(XX);
	Free(XDxxRxy);
	Free(YYXDR);
	//------adaEN Apr. 2013	
	Free(W);
	//------adaEN Apr. 2013	
	return lambda_max;	
}

/*
//Q[i,k] =	N*sigma2*IM - (Y*(Y'-X'*DxxRxy)))
void QlambdaStart(double *Y,double *X, double *Q, double sigma2,int M, int N)
{	
	// Oct 08, 2012: assume one eQTL for each gene; This fucntion needs significant revision if this assumption doesnot hold
	double *dxx, *rxy, *DxxRxy,*readPtr1,*readPtr2;
	
	dxx				= (double* ) Calloc(M, double);
	rxy				= (double* ) Calloc(M, double);
	DxxRxy			= (double* ) Calloc(M, double);
	int i,index,ldk,lda,ldb,ldc;
	int inci = 1;
	int incj = 1; 
	//double norm;
	lda = M;
	for(i=0;i<M;i++)
	{
		readPtr1 	= &X[i]; //ith row
		readPtr2 	= &Y[i];

		//norm  		= F77_CALL(dnrm2)(&N,readPtr1,&M);	
		//dxx[i] 		= pow(norm,2);
		dxx[i] = F77_CALL(ddot)(&N,readPtr1,&lda,readPtr1,&M);
		//res = ddot(n, x, incx, y, incy)
		rxy[i] 		= F77_CALL(ddot)(&N,readPtr1,&lda,readPtr2,&M);
		DxxRxy[i] 	= rxy[i]/dxx[i];		
	}
	//abs(N*sigma2*IM - (Y*(Y'-X'*DxxRxy)))./W ; W[i,i] = inf.
	double Nsigma2  = N*sigma2; 			// int * double --> double

	//cache X[k,:]*DxxRxy[k]
	double * XDxxRxy;
	int MN = M*N;
	XDxxRxy = (double* ) Calloc(MN, double);
	F77_CALL(dcopy)(&MN,X,&inci,XDxxRxy,&incj);
	double alpha;	
	for(i=0;i<M;i++)
	{
		alpha  = -DxxRxy[i];
		readPtr1 = &XDxxRxy[i]; //ith row
		F77_CALL(dscal)(&N,&alpha, readPtr1,&M);//	(n, a, x, incx)
	}
	
	// Y- XDxxRxy  			daxpy(n, a, x, incx, y, incy) y= ax + y
	alpha  = 1.0;
	// XDxxRxy <- alpha*Y + XDxxRxy
	F77_CALL(daxpy)(&MN,&alpha,Y,&inci,XDxxRxy,&incj);
	//double *YYXDR; //= Y*XDxxRxy' 		--> Q

	double beta;
	char transa = 'N';
	char transb = 'T';
	alpha = -1;
	beta = 0;
	//F77_CALL(dgemm)(&transa, &transb,&M, &M, &N,&alpha, Y,&M, XDxxRxy, &M, &beta,Q, &M); //M xK, K xN  --> MxN, N xM --> M <-M, N<-M, k<-N
	//transpose 	(Y-X*DxxRxy)*Y'

	ldb = M;
	ldc = M;
	ldk = M;
	F77_CALL(dgemm)(&transa, &transb,&M, &lda, &N,&alpha, XDxxRxy,&ldb, Y, &ldc, &beta,Q, &ldk); //M xK, K xN  --> MxN, N xM --> M <-M, N<-M, k<-N	
	
	//diagonal -->0; other element /wij
	for(i=0;i<M;i++)
	{
		index = i*M + i;
		Q[index]= Q[index] + Nsigma2;
	}	
	
	Free(dxx);
	Free(rxy);
	Free(DxxRxy);
	//Free(XX);
	Free(XDxxRxy);

	
}

// 8888888888888888888888888888888888888888888888888888888888888888888888
//Q = N*sigma2*inv(I-B)-(Y-B*Y-fL*X)-mueL*ones(1,N))*Y';
void QlambdaMiddle(double *Y,double *X, double *Q,double *B,double *f, double *mue, double sigma2,int M, int N)
{	
	// Oct 08, 2012: assume one eQTL for each gene; This fucntion needs significant revision if this assumption doesnot hold
	//I - B; copy of IB for inverse
	double *IB, *IBinv,*IBcopy;
	int MM = M*M;
	int MN = M*N;
	IB = (double* ) Calloc(MM, double);
	IBinv = (double* ) Calloc(MM, double);
	IBcopy = (double* ) Calloc(MM, double);
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
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci);//initialized
	F77_CALL(dcopy)(&MM,&alpha,&inc0,IBinv,&inci);
	
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		IB[index] = 1 + IB[index];
		IBinv[index] = 1;
	}
	F77_CALL(dcopy)(&MM,IB,&inci,IBcopy,&incj);	

	
	//MatrixInverse(IBinv,M);
	//By linear solver: not inverse (IB*x = IM) result stored in IM; 
	//multiLinearSystem(IB, M,IBinv,M);
	int info = 0;
	int *ipiv;
	ipiv = (int *) Calloc(M,int);
	int lda = M;
	int ldb = M;
	int ldc = M;
	int ldk = M;
	F77_CALL(dgesv)(&M, &ldk, IBcopy, &lda, ipiv, IBinv, &ldb, &info);

	
	
	//abs(N*sigma2*inv(I-B) - NOISE*Y'.
	double Nsigma2  = N*sigma2; 			// int * double --> double
	double *Noise;
	Noise = (double* ) Calloc(MN, double);	
	//(I-B)*Y-bsxfun(@times,f,X);
	char transa = 'N';
	char transb = 'N';
	alpha = 1;
	F77_CALL(dgemm)(&transa, &transb,&M, &N, &ldk,&alpha, IB, &lda, Y, &ldb, &beta, Noise, &ldc);
	double *readPtr1, *readPtr2;
	for(i=0;i<M;i++)
	{
		readPtr1 = &X[i];
		readPtr2 = &Noise[i];
		alpha = -f[i]; // y= alpha x + y
		F77_CALL(daxpy)(&N, &alpha,readPtr1, &lda,readPtr2, &M);
	}//row i = 1:M

	//NoiseMat(Noise,B,Y, f, X, M, N);
	//Noise - Mue
	//Errs(irho,cv)   = norm((1-Missing_test).*(A*Ytest-bsxfun(@times,fR,Xtest)-mueR*ones(1,Ntest)),'fro')^2;
	//Noise - mue: Mx N
	alpha = -1;
	for(i=0;i<N;i++)
	{
		readPtr1 = &Noise[i*M];
		F77_CALL(daxpy)(&M, &alpha,mue, &inci,readPtr1, &incj);
	}	
	//Nsigma2*IBinv -  Noise *Y'
	//-Noise*Y' -->Q  	C := alpha*op(A)*op(B) + beta*C,
	//alpha = -1;
	transb = 'T';
	F77_CALL(dgemm)(&transa, &transb,&M, &ldk, &N,&alpha, Noise, &lda, Y, &ldb, &beta, Q, &ldc);
	//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
	alpha = Nsigma2;
	F77_CALL(daxpy)(&MM, &alpha,IBinv, &inci,Q, &incj);
	
	Free(IB);
	Free(IBinv);
	Free(IBcopy);
	Free(Noise);
	Free(ipiv);
	
}


void QlambdaMiddleCenter(double *Y,double *X, double *Q,double *B,double *f, double sigma2,int M, int N,
					double *IBinv)
{	
	// Oct 08, 2012: assume one eQTL for each gene; This fucntion needs significant revision if this assumption doesnot hold
	//I - B; copy of IB for inverse
	double *IB; 	//, *IBinv,*IBcopy
	int MM = M*M;
	int MN = M*N;
	IB = (double* ) Calloc(MM, double);
	//IBinv = (double* ) Calloc(MM, double);
	//IBcopy = (double* ) Calloc(MM, double);
	int inci = 1;
	int incj = 1;
	//int inc0 = 0;
	F77_CALL(dcopy)(&MM,B,&inci,IB,&incj);	
	int i,index;
	double alpha;
	double beta = 0;
	alpha = -1;
	F77_CALL(dscal)(&MM,&alpha,IB,&inci);
	//alpha = 0;
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci);//initialized
	//F77_CALL(dcopy)(&MM,&alpha,&inc0,IBinv,&inci);
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		IB[index] = 1 + IB[index];
		//IBinv[index] = 1;
	}
	//F77_CALL(dcopy)(&MM,IB,&inci,IBcopy,&incj);	

	
	//MatrixInverse(IBinv,M);
	//By linear solver: not inverse (IB*x = IM) result stored in IM; 
	//multiLinearSystem(IB, M,IBinv,M);
	//int info = 0;
	//int *ipiv;
	//ipiv = (int *) Calloc(M,int);
	int lda = M;
	int ldb = M;
	int ldc = M;
	int ldk = M;
	//F77_CALL(dgesv)(&M, &ldk, IBcopy, &lda, ipiv, IBinv, &ldb, &info);

	
	
	//abs(N*sigma2*inv(I-B) - NOISE*Y'.
	double Nsigma2  = N*sigma2; 			// int * double --> double
	double *Noise;
	Noise = (double* ) Calloc(MN, double);	
	//(I-B)*Y-bsxfun(@times,f,X);
	char transa = 'N';
	char transb = 'N';
	alpha = 1;
	F77_CALL(dgemm)(&transa, &transb,&M, &N, &ldk,&alpha, IB, &lda, Y, &ldb, &beta, Noise, &ldc);
	double *readPtr1, *readPtr2;
	for(i=0;i<M;i++)
	{
		readPtr1 = &X[i];
		readPtr2 = &Noise[i];
		alpha = -f[i]; // y= alpha x + y
		F77_CALL(daxpy)(&N, &alpha,readPtr1, &lda,readPtr2, &M);
	}//row i = 1:M

	//NoiseMat(Noise,B,Y, f, X, M, N);
	//Noise - Mue
	//Errs(irho,cv)   = norm((1-Missing_test).*(A*Ytest-bsxfun(@times,fR,Xtest)-mueR*ones(1,Ntest)),'fro')^2;
	//Noise - mue: Mx N

	//Nsigma2*IBinv -  Noise *Y'
	//-Noise*Y' -->Q  	C := alpha*op(A)*op(B) + beta*C,
	alpha = -1;
	transb = 'T';
	F77_CALL(dgemm)(&transa, &transb,&M, &ldk, &N,&alpha, Noise, &lda, Y, &ldb, &beta, Q, &ldc);
	//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
	alpha = Nsigma2;
	F77_CALL(daxpy)(&MM, &alpha,IBinv, &inci,Q, &incj);
	
	Free(IB);
	//Free(IBinv);
	//Free(IBcopy);
	Free(Noise);
	//Free(ipiv);
	
}


// 8888888888888888888888888888888888888888888888888888888888888888888888
//BLOCK COORDINATE ASCENSION: QIBinv = inv(I-B): by multi linear system
void UpdateIBinvPermute(double *QIBinv, double *B, int M)
{
	//I - B; copy of IB for inverse
	double *IB,*IBinv;	//, *IBinv,*IBcopy;
	int MM = M*M;
	int lda = M;
	int ldb = M;
	int ldk = M;
	IB = (double* ) Calloc(MM, double);
	IBinv = (double* ) Calloc(MM, double);
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
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci);//initialized
	F77_CALL(dcopy)(&MM,&alpha,&inc0,IBinv,&inci);
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		IB[index] = 1 + IB[index];
		IBinv[index] = 1;
	}
	
	//MatrixInverse(IBinv,M);
	//By linear solver: not inverse (IB*x = IM) result stored in IM; 
	//multiLinearSystem(IB, M,IBinv,M);
	int info = 0;
	int *ipiv;
	ipiv = (int *) Calloc(M,int);
	F77_CALL(dgesv)(&M, &ldk, IB, &lda, ipiv, IBinv, &ldb, &info);
	double *ptr1,*ptr2;
	//for(i=0;i<M;i++) Rprintf("IPIV: \n: %d \t",ipiv[i]);
	//Rprintf("\n");
	
	
	for(i=0;i<M;i++)
	{
		index = ipiv[i] -1;
		ptr1 = &QIBinv[index*M];
		ptr2 = &IBinv[i*M];
		F77_CALL(dcopy)(&M,ptr2,&inci,ptr1,&incj);
		
	}
	
	Free(IB);
	Free(ipiv);
	Free(IBinv);
}


// 8888888888888888888888888888888888888888888888888888888888888888888888
//BLOCK COORDINATE ASCENSION: QIBinv = inv(I-B): by multi linear system
void UpdateIBinv(double *QIBinv, double *B, int M)
{
	//I - B; copy of IB for inverse
	double *IB;	//, *IBinv,*IBcopy;
	int MM = M*M;
	int lda = M;
	int ldb = M;
	int ldk = M;
	IB = (double* ) Calloc(MM, double);

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
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci);//initialized
	F77_CALL(dcopy)(&MM,&alpha,&inc0,QIBinv,&inci);
	for(i=0;i<M;i++) 
	{
		index = i*M + i;
		IB[index] = 1 + IB[index];
		QIBinv[index] = 1;
	}
	
	//MatrixInverse(IBinv,M);
	//By linear solver: not inverse (IB*x = IM) result stored in IM; 
	//multiLinearSystem(IB, M,IBinv,M);
	int info = 0;
	int *ipiv;
	ipiv = (int *) Calloc(M,int);
	F77_CALL(dgesv)(&M, &ldk, IB, &lda, ipiv, QIBinv, &ldb, &info);

	Free(IB);
	Free(ipiv);
}
*/
//no Missing
//--------------------------------------------------------------------------------  WEIGHTED_LASSOSF
double Weighted_LassoSf_adaEN(double * Wori, double *B, double *f, double *Ycopy,double *Xcopy, 				//------adaEN Apr. 2013	
		double *Q, double lambda_factor, double lambda_factor_prev, double sigma2, int max_iter,
		int M, int N, int verbose,double *QIBinv,double lambda_max,			//double * mue,
		double alpha_factor)
{
	int i,j,index,ldM;
	//lda = M;
	//ldb = M;ldb,
	ldM = M;//fixed
	// return lambda;
	double *meanY, *meanX;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	
	//copy Y, X; 
	double *Y, *X;
	int MN = M*N;
	int MM = M*M;
	Y = (double* ) Calloc(MN, double);
	X = (double* ) Calloc(MN, double);
	
	//F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
		//double *dy, const int *incy);
	int inci,incj, inc0;
	inci	= 1;
	incj 	= 1;
	inc0 	= 0;
	F77_CALL(dcopy)(&MN,Ycopy,&inci,Y,&incj);
	F77_CALL(dcopy)(&MN,Xcopy,&inci,X,&incj);
	centerYX(Y,X, meanY, meanX,M, N);
	
	//return value
	//double sigma2 			= SIGMA2[0];
	double lambda;//lambda_max,
	//lambdaMax
	//lambda_max 				= lambdaMax(Y,X,W,M, N);
	if(verbose>4) Rprintf("\t\t\t\tEnter Function: weighted_LassoSf. The maximum lambda is: %f\n\n",lambda_max);
	lambda 					= lambda_factor*lambda_max;
	
	//none zeros
	double alpha,beta;
	beta = 0;
	double deltaLambda;
//------adaEN Apr. 2013		
	double *s, *S,*Wcopy, *W;
	S = (double* ) Calloc(MM, double);
	s = (double* ) Calloc(M, double);
	W 						= (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,Wori,&inci,W,&incj);
	F77_CALL(dscal)(&MM,&alpha_factor,W,&inci); 
//------adaEN Apr. 2013	
	Wcopy = (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,W,&inci,Wcopy,&incj);

	deltaLambda 			= (2*lambda_factor - lambda_factor_prev)*lambda_max;	
	F77_CALL(dscal)(&MM,&deltaLambda,Wcopy,&inci); //wcopy = deltaLambda*W
	
	//ei = 0
	double *ei,toyZero;
	toyZero= 0;
	ei = (double* ) Calloc(M, double);
	//F77_CALL(dscal)(&M,&toyZero,ei,&inci);
	F77_CALL(dcopy)(&M,&toyZero,&inc0,ei,&inci);
/*	double *eye;
	eye = (double* ) Calloc(M, double);
	alpha = 1;
	F77_CALL(dcopy)(&M,&alpha,&inc0,eye,&inci);
*/
	double *readPtr,*readPtr2;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			//W[i,j]
			index = j*M  +i;
//------adaEN Apr. 2013
// version 1_3: Qij			
			//if(fabs(Q[index])>= Wcopy[index] && i!= j)
			if(fabs(Q[index] -(1-alpha_factor)*lambda*B[index])>= Wcopy[index] && i!= j)	
// version 1_3: Qij	
//------adaEN Apr. 2013
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
/*	ldk = M;
	//lda = M;
	F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, S, &ldM, eye, &inci, &beta,s, &incj);
*/	
//printMat(W,M,M);
	//f0, F1
	double *f0,*F1;
	//int qdif = M*M;
	f0 	= (double* ) Calloc(M, double);
	F1 	= (double* ) Calloc(MM, double);
	
	double *y_j;
	//xi 	= (double* ) Calloc(N, double);
	y_j 	= (double* ) Calloc(N, double);
	double *F1ptr;


	double XYi, XXi;
	for(i=0;i<M;i++)
	{
		readPtr = &X[i];
		//F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
		readPtr2 = &Y[i];
		//F77_CALL(dcopy)(&N,readPtr2,&M,y_j,&inci);

		//dot product
		//XYi = F77_CALL(ddot)(&N, xi, &inci,y_j, &incj);
		XYi = F77_CALL(ddot)(&N, readPtr, &M,readPtr2, &M);
		//XXi = F77_CALL(ddot)(&N, xi, &inci,xi, &incj);
		//norm2 = F77_CALL(dnrm2)(&N,xi,&inci);
		//XXi 	= pow(norm2,2);
		XXi = F77_CALL(ddot)(&N, readPtr, &M,readPtr, &M);
		f0[i] 	= XYi/XXi;
		F1ptr	= &F1[M*i];//start from ith column
		//Y*X(i,:)' 		y := alpha*A*x + beta*y,
		alpha = 1/XXi;
		F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, readPtr, &M, &beta,F1ptr, &incj);
	}
	
	//printMat(f0,M,1);

	// entering loop		
	double *IBinv,*zi,*a_iT;		// y_j: one row of Y: Nx1
	IBinv 	= (double* ) Calloc(MM, double);
	//zi 		= (double* ) Calloc(M, double);
	a_iT 	= (double* ) Calloc(N, double);


	
	//loop starts here
	int iter = 0;
	double js_i, m_ij,B_old, lambdaW,beta_ij,r_ij, Bij;
	//dynamic variable keep intermidiate values 
	double *eiB;
	eiB = (double* ) Calloc(M, double);
	double *BiT;
	BiT = (double* ) Calloc(M, double);
	//quadratic function
	double d_ij, theta_ijp,k_ijp,q_ijp,Bijpp, Bijpm; //case (14)
	double q_ijm, theta_ijm, Bijmm, Bijmp,Lss,candsBij,LssCands;
	
	//converge of gene i
	double dB,ziDb,BF1;
	
	//converge of while
	double delta_BF,FnormOld, FnormChange;
	double *BfOld,*BfNew,*BfChange;
	index = M*(M  +1);
	BfOld = (double* ) Calloc(index, double);
	BfNew = (double* ) Calloc(index, double);
	BfChange = (double* ) Calloc(index, double);
	
	//	//linear system
	//int *ipiv;
	//ipiv = (int *) R_alloc(N,sizeof(int));
	//ipiv = (int *) Calloc(M,int);
	//int nrhs = 1;
	//int info = 0;
//for(i=0;i<M;i++) ipiv[i] = i+1;
	//
	//double *aiTQuad; 
	//aiTQuad 	= (double* ) Calloc(N, double);

	//alpha =  F77_CALL(dnrm2)(&index,BfOld,&inci);	
//Rprintf("BfOldnorm:%f\n",alpha);	
	//int anyUpdateRow = 0;
	
	while(iter < max_iter)
	{
		iter = iter + 1;
		//converge Bfold = [B f];
		F77_CALL(dcopy)(&MM,B,&inci,BfOld,&incj);
		//last column
		F1ptr = &BfOld[MM];
		F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
	//alpha =  F77_CALL(dnrm2)(&index,BfOld,&inci);	
//Rprintf("BfOldnorm:%f\n",alpha);
		//printMat(BfOld,M,M+1);
		//
		for(i=0;i<M;i++)
		{
			if(s[i] >0)
			{ 	//
				if(verbose>6) Rprintf("\t\t\t\t\t updating gene %d \n",i);
				//
				ei[i] = 1;
				//zi   IBinv = I -B
				//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
				//alpha = -1; 
				//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
				//diagonal + 1
				//for(j=0;j<M;j++) 
				//{
				//	index = j*M + j;
				//	IBinv[index] = 1 + IBinv[index];
				//}

				//call matrix inverse
				//MatrixInverse(IBinv,M);
				//zi is ith column of IBinv
				//for(j =0;j<M;j++) zi[j] 	= IBinv[i*M+j];
				
				// by linear solver: zi = (I - B)\ei;
				//copy ei to zi;
				zi = &QIBinv[i*M];
				//F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
				//
				//if (i==0) printMat(zi,M,1);
				//
				//call linearSystem(double *a, int N,double *B) a is NxN
//linearSystem(IBinv, M,zi);  //zi is updated.

				//F77_CALL(dgesv)(&M, &nrhs, IBinv, &ldM, ipiv, zi, &ldb, &info);
				//anyUpdateRow = 1;
//j reserved				//for j = js_i
				for(j=0;j<M;j++) 
				{
					js_i = S[j*M + i]; 		//ith row
					if(js_i >0)
					{

						m_ij 	= zi[j];
						B_old 	= B[j*M + i]; //B[i,j]
						if(j!=i)
						{
						
							//y_j; jth row Nx1
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
							//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
							alpha = 1;
							F77_CALL(daxpy)(&M, &alpha,BiT, &inci,eiB, &incj);
							//a_iT      = (ei'-BiT)*Y-f(i)*X(i,:);
							readPtr = &X[i];
							F77_CALL(dcopy)(&N,readPtr,&M,a_iT,&inci);	

							//a_iT = -f[i]*xi 		dscal(n, a, x, incx) 		x = a*x
							alpha = -f[i];
							F77_CALL(dscal)(&N,&alpha,a_iT,&inci);							

							transa='T'; //y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
							beta = 1;
							alpha = 1;
							F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, eiB, &inci, &beta,a_iT, &incj);
							
							//r_ij                = y_j'*y_j;	
r_ij = F77_CALL(ddot)(&N, y_j, &inci,y_j, &incj);
//------adaEN Apr. 2013	
r_ij = r_ij + (1 -alpha_factor)*lambda;
//------adaEN Apr. 2013		
							//norm2 = F77_CALL(dnrm2)(&N,y_j,&inci);
							//r_ij 	= pow(norm2,2);
							
							//beta_ij             = a_iT*y_j;
							beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
							
							if (fabs(m_ij)<1e-10) //go to the linear equation 
							{
								//
								if(verbose>7) Rprintf("\t\t\t\t\t\t\t gene %d \t interact with gene %d.\tLinear equation\n",i,j);
								//
								Bij = (beta_ij-lambdaW)/r_ij;
								//Rprintf("\t\t gene (%d,\t%d): linear Bij %f\n",i,j,Bij);
				//
				//if (i==0) 
				//{
				//	Rprintf("\t\t\t beta_ij: %f;\t r_ij:%f; lambdaW: %f\n ",beta_ij,r_ij,lambdaW);
				//	printMat(a_iT,N,1);
					//printMat(y_j,N,1);
				//	printMat(eiB,M,1);
				//}
				//
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
								//Lss = quadraticLss(sigma2, N, d_ij, candsBij, r_ij, lambdaW,beta_ij);								
								//Lss = quadraticLssJuan(sigma2,N, d_ij, candsBij, a_iT,y_j,lambdaW);
//beta = F77_CALL(ddot)(&N, a_iT, &inci,a_iT, &incj);
								//norm2 = F77_CALL(dnrm2)(&N,a_iT,&inci);
								//beta 	= pow(norm2,2);
							
								//Lss = sigma2*N*log(fabs(d_ij - candsBij)+1e-16) - beta/2 -lambdaW*fabs(candsBij);
								//Lss = sigma2*N*log(fabs(d_ij)+1e-16) - beta/2;
								Lss = sigma2*N*log(fabs(d_ij)+1e-16);
								
								// a_iT = a_iT - candsBij*y_j 		daxpy(n, a, x, incx, y, incy) 	y := a*x + y							
								if (Bijpp>0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijpp, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijpp, a_iT,y_j,lambdaW);	
									//aiTQuad dcopy(n, x, incx, y, incy)  y= x
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijpp;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - beta/2 -lambdaW*fabs(Bijpp);
									LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - r_ij*pow(Bijpp,2)/2 + beta_ij*Bijpp -lambdaW*fabs(Bijpp); 
									
									if(LssCands>Lss) 
									{
										candsBij = Bijpp;
										Lss 	= LssCands;
									}	
								}
								if (Bijpm>0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijpm, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijpm, a_iT,y_j,lambdaW);								
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijpm;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijpm)+1e-16) - beta/2 -lambdaW*fabs(Bijpm);
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
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijmm, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijmm, a_iT,y_j,lambdaW);								
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijmm;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - beta/2 -lambdaW*fabs(Bijmm);
									LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - r_ij*pow(Bijmm,2)/2 + beta_ij*Bijmm -lambdaW*fabs(Bijmm);  
									if(LssCands>Lss) 
									{
										candsBij = Bijmm;
										Lss 	= LssCands;
									}	
								}
								if (Bijmp<0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijmp, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijmp, a_iT,y_j,lambdaW);									
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijmp;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijmp)+1e-16) - beta/2 -lambdaW*fabs(Bijmp);			
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
						
						//update QIBinv ith column updated; jth row: 
						//readPtr = &QIBinv[j];
						//F77_CALL(dscal)(&M,&ziDb,readPtr,&M);
						//QIBinv[j*M + i] = QIBinv[j*M + i]*ziDb;
						// adjugate of QIBinv = ziDb*(QIBinvadj - (-1)^kl) + (-1)^kl
						//for(k=0;k<M;k++) // kth row
						//{
						//	for(l=0;l<M;l++)//lth column
						//	{
						//		if(k!=j && l!=i) 
						//		{
						//			alpha = pow(-1.0,k+l);
						//			QIBinv[l*M+k] = ziDb*(QIBinv[l*M+k] -alpha) + alpha;
						//		}	
						//	}
						//}
						
					}//js_i >0
				}//j = 1:M	
				//f
				//BF1 = B(i,:)*F1(:,i)
				readPtr = &B[i];
				F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);

				F1ptr = &F1[M*i];
				BF1 = F77_CALL(ddot)(&M, BiT, &inci,F1ptr, &incj);

				f[i] = f0[i] - BF1;
				ei[i] = 0; // re-set ei for next i
			}else//s[i]  no un-zero weight in this gene
			{
				readPtr = &B[i];
				//F77_CALL(dscal)(&M,&toyZero,readPtr,&ldM);
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
		
//Rprintf("FnormOld: %f; \t FnormChange: %f\n",FnormOld,FnormChange);
		//MASTER IS HERE: Update before break: IBinv will be used in Qij
		// BLOCK COORDINATE ASCEND: Update IBinv
		UpdateIBinv(QIBinv, B,M);	
		
		if(verbose>5) Rprintf("\t\t\t\t\t\tdelta_BF: %f\n",delta_BF);
		if(delta_BF<1e-3)		//break out
		{
			break;
		}
		

		
	}//while

	//QIBinv: compute not updated ones
	//F77_NAME(dgetrs)(const char* trans, const int* n, const int* nrhs,
	//	 const double* a, const int* lda, const int* ipiv,
	//	 double* b, const int* ldb, int* info);
	//if(anyUpdateRow>0) //QIBinv --> = Cij/det(I-B) != eye(M); ipiv not zeros
	//{
	//	transa = 'N';
	//	for(i=0;i<M;i++)
	//	{
	//		if(s[i] ==0)
	//		{ 	
	//			info = 0;
	//			ei[i] = 1;
	//			zi = &QIBinv[i*M];
	//			F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
				// dgetrs
	//			F77_CALL(dgetrs)(&transa, &M, &nrhs,IBinv,&M,ipiv,zi,&ldb,&info); //column i updated
	//			ei[i] = 0;
	//		}
	//	}
	//}else
	//{
		//initialize IBinv
	//	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	//	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
	//}
	
	
	
	if(verbose>4) Rprintf("\t\t\t\tCurrent lambda: %f;\t number of iteration is: %d.\tExiting Weighted_LassoSf\n\n",lambda, iter);
	
	//mue          = (IM-B)*meanY-bsxfun(@times,f,meanX);

	//IBinv = I -B
	//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
	//alpha = -1; 
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
	//diagonal + 1
	//for(j=0;j<M;j++) 
	//{
	//	index = j*M + j;
	//	IBinv[index] = 1 + IBinv[index];
	//	mue[j] = -f[j]*meanX[j];
	//}
	//mue[i] = -f[i]*meanX[i];
	
	//	//dgemv mue = Ax + beta*mue
	//transa = 'N';
	//alpha = 1;
	//beta = 1;
	//ldk = M;
	//F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, IBinv, &ldM, meanY, &inci, &beta,mue, &incj);

	//mue                     	= (IM-B)*meanY-bsxfun(@times,f,meanX); diag(f)*meanX

	Free(meanY);
	Free(meanX);
	Free(Y);
	Free(X);
	
	Free(S);
	Free(s);
	Free(f0);
	Free(F1);
	Free(Wcopy);
	
	//Free(xi);
	Free(y_j);
	
	Free(ei);
	Free(IBinv);
	//Free(zi);
	Free(a_iT);
	
	Free(eiB);
	Free(BiT);
	Free(BfOld);
	Free(BfNew);
	Free(BfChange);
	
	//Free(ipiv);
	//------adaEN Apr. 2013	
	Free(W);
	//------adaEN Apr. 2013	
	//Free(aiTQuad);
	//Free(eye);
	return lambda;
	//sigma2 remains the same
}//weighted_LassoSf

//combine lasso and zero_lasso for CV_selecting lambda
double Weighted_LassoSf_MLf_adaEN(double * Wori, double *BL, double *fL, double *Ycopy,double *Xcopy,
		double *Q, double lambda_factor, double lambda_factor_prev, double sigma2, int max_iter,
		int M, int N, int verbose,
		double *BC, double *fC, double *mue,double *QIBinv,double *IBinvZero,double lambda_max,
		double alpha_factor)
{
	//SET TO PART1: LASSO
	double *B, *f;
	B = &BL[0];
	f = &fL[0];
	//mue is used for calculation when Y, X are not centered: for testing set;

	int i,j,index,ldk,ldM;
	//lda = M;
	//ldb = M;ldb,
	ldM = M;//fixed
	// return lambda;
	double *meanY, *meanX;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	
	//copy Y, X; 
	double *Y, *X;
	int MN = M*N;
	int MM = M*M;
	Y = (double* ) Calloc(MN, double);
	X = (double* ) Calloc(MN, double);
	
	//F77_NAME(dcopy)(const int *n, const double *dx, const int *incx,
		//double *dy, const int *incy);
	int inci,incj,inc0;
	inci	= 1;
	incj 	= 1;
	inc0 	= 0;
	F77_CALL(dcopy)(&MN,Ycopy,&inci,Y,&incj);
	F77_CALL(dcopy)(&MN,Xcopy,&inci,X,&incj);
	centerYX(Y,X, meanY, meanX,M, N);
	
	//return value
	//double sigma2 			= SIGMA2[0];
	double lambda;//lambda_max,
	//lambdaMax
	//lambda_max 				= lambdaMax(Y,X,W,M, N);
	if(verbose>4) Rprintf("\t\t\t\tEnter Function: weighted_LassoSf. The maximum lambda is: %f\n\n",lambda_max);
	lambda 					= lambda_factor*lambda_max;
	
	//none zeros
	double alpha,beta;
	beta = 0;
	double deltaLambda;
	
//------adaEN Apr. 2013		
	double *s, *S,*Wcopy,*W; //global copy of W
	S = (double* ) Calloc(MM, double);
	s = (double* ) Calloc(M, double);
	W = (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,Wori,&inci,W,&incj);
	F77_CALL(dscal)(&MM,&alpha_factor,W,&inci); //MAKE sure alpha_factor is bcasted
//------adaEN Apr. 2013	

	Wcopy = (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,W,&inci,Wcopy,&incj);

	deltaLambda 			= (2*lambda_factor - lambda_factor_prev)*lambda_max;	
	F77_CALL(dscal)(&MM,&deltaLambda,Wcopy,&inci); //wcopy = deltaLambda*W
	
	//ei = 0
	double *ei,toyZero;
	toyZero= 0;
	ei = (double* ) Calloc(M, double);
	//F77_CALL(dscal)(&M,&toyZero,ei,&inci);
	F77_CALL(dcopy)(&M,&toyZero,&inc0,ei,&inci);
/*	double *eye;
	eye = (double* ) Calloc(M, double);
	alpha = 1;
	F77_CALL(dcopy)(&M,&alpha,&inc0,eye,&inci);
*/
	double *readPtr,*readPtr2;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			//W[i,j]
			index = j*M  +i;
//------adaEN Apr. 2013				
// version 1_3: Qij			
			//if(fabs(Q[index])>= Wcopy[index] && i!= j)
			if(fabs(Q[index] -(1-alpha_factor)*lambda*B[index])>= Wcopy[index] && i!= j)	
// version 1_3: Qij	
//------adaEN Apr. 2013	
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
/*	ldk = M;
	//lda = M;
	F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, S, &ldM, eye, &inci, &beta,s, &incj);
*/	
//printMat(W,M,M);
	//f0, F1
	double *f0,*F1;
	//int qdif = M*M;
	f0 	= (double* ) Calloc(M, double);
	F1 	= (double* ) Calloc(MM, double);
	
	//double *xi, *y_j;
	double *y_j;
	//xi 	= (double* ) Calloc(N, double);
	y_j 	= (double* ) Calloc(N, double);
	double *F1ptr;


	double XYi, XXi;
	for(i=0;i<M;i++)
	{
		readPtr = &X[i];
		//F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
		readPtr2 = &Y[i];
		//F77_CALL(dcopy)(&N,readPtr,&M,y_j,&inci);

		//dot product
		//XYi = F77_CALL(ddot)(&N, xi, &inci,y_j, &incj);
		XYi = F77_CALL(ddot)(&N, readPtr, &M,readPtr2, &M);
		//XXi = F77_CALL(ddot)(&N, xi, &inci,xi, &incj);
		//norm2 = F77_CALL(dnrm2)(&N,xi,&inci);
		//XXi 	= pow(norm2,2);
		XXi = F77_CALL(ddot)(&N, readPtr, &M,readPtr, &M);
		f0[i] 	= XYi/XXi;
		F1ptr	= &F1[M*i];//start from ith column
		//Y*X(i,:)' 		y := alpha*A*x + beta*y, alpha*Y *xi + beta*F1
		alpha = 1/XXi;
		F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, readPtr, &M, &beta,F1ptr, &incj);
	}
	
	//printMat(f0,M,1);

	// entering loop
	double *IBinv,*zi,*a_iT;// y_j: one row of Y: Nx1
	IBinv 	= (double* ) Calloc(MM, double);
	//zi 		= (double* ) Calloc(M, double);
	a_iT 	= (double* ) Calloc(N, double);

	
	
	//loop starts here
	int iter = 0;
	double js_i, m_ij,B_old, lambdaW,beta_ij,r_ij, Bij;
	//dynamic variable keep intermidiate values 
	double *eiB;
	eiB = (double* ) Calloc(M, double);
	double *BiT;
	BiT = (double* ) Calloc(M, double);
	//quadratic function
	double d_ij, theta_ijp,k_ijp,q_ijp,Bijpp, Bijpm; //case (14)
	double q_ijm, theta_ijm, Bijmm, Bijmp,Lss,candsBij,LssCands;
	
	//converge of gene i
	double dB,ziDb,BF1;
	
	//converge of while
	double delta_BF,FnormOld, FnormChange;
	double *BfOld,*BfNew,*BfChange;
	index = M*(M  +1);
	BfOld = (double* ) Calloc(index, double);
	BfNew = (double* ) Calloc(index, double);
	BfChange = (double* ) Calloc(index, double);
	
	//	//linear system
	//int *ipiv;
	//ipiv = (int *) R_alloc(N,sizeof(int));
	//ipiv = (int *) Calloc(M,int);
	//int nrhs = 1;
	//int info = 0;
	//for(i=0;i<M;i++) ipiv[i] = i+1;
	//
	//double *aiTQuad; 
	//aiTQuad 	= (double* ) Calloc(N, double);

	//alpha =  F77_CALL(dnrm2)(&index,BfOld,&inci);	
//Rprintf("BfOldnorm:%f\n",alpha);	
	//int anyUpdateRow = 0;
	while(iter < max_iter)
	{
		iter = iter + 1;
		//converge Bfold = [B f];
		F77_CALL(dcopy)(&MM,B,&inci,BfOld,&incj);
		//last column
		F1ptr = &BfOld[MM];
		F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
	//alpha =  F77_CALL(dnrm2)(&index,BfOld,&inci);	
//Rprintf("BfOldnorm:%f\n",alpha);
		//printMat(BfOld,M,M+1);
		//
		for(i=0;i<M;i++)
		{
			if(s[i] >0)
			{ 	//
				if(verbose>6) Rprintf("\t\t\t\t\t updating gene %d \n",i);
				//
				ei[i] = 1;
				//zi   IBinv = I -B
				//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
				//alpha = -1; 
				//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
				//diagonal + 1
				//for(j=0;j<M;j++) 
				//{
				//	index = j*M + j;
				//	IBinv[index] = 1 + IBinv[index];
				//}

				//call matrix inverse
				//MatrixInverse(IBinv,M);
				//zi is ith column of IBinv
				//for(j =0;j<M;j++) zi[j] 	= IBinv[i*M+j];
				
				// by linear solver: zi = (I - B)\ei;
				//copy ei to zi;
				zi = &QIBinv[i*M];
				//F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
				//
				//if (i==0) printMat(zi,M,1);
				//
				//call linearSystem(double *a, int N,double *B) a is NxN
//linearSystem(IBinv, M,zi);  //zi is updated.

				//F77_CALL(dgesv)(&M, &nrhs, IBinv, &ldM, ipiv, zi, &ldb, &info);
				//anyUpdateRow = 1;
//j reserved				//for j = js_i
				for(j=0;j<M;j++) 
				{
					js_i = S[j*M + i]; 		//ith row
					if(js_i >0)
					{

						m_ij 	= zi[j];
						B_old 	= B[j*M + i]; //B[i,j]
						if(j!=i)
						{
						
							//y_j; jth row Nx1
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
							//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
							alpha = 1;
							F77_CALL(daxpy)(&M, &alpha,BiT, &inci,eiB, &incj);
							//a_iT      = (ei'-BiT)*Y-f(i)*X(i,:);
							readPtr = &X[i];
							F77_CALL(dcopy)(&N,readPtr,&M,a_iT,&inci);	

							//a_iT = -f[i]*xi 		dscal(n, a, x, incx) 		x = a*x
							alpha = -f[i];
							F77_CALL(dscal)(&N,&alpha,a_iT,&inci);							

							transa='T'; //y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
							beta = 1;
							alpha = 1;
							F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, eiB, &inci, &beta,a_iT, &incj);
							
							//r_ij                = y_j'*y_j;	
r_ij = F77_CALL(ddot)(&N, y_j, &inci,y_j, &incj);
//------adaEN Apr. 2013	
r_ij = r_ij + (1 -alpha_factor)*lambda;
//------adaEN Apr. 2013	
							//norm2 = F77_CALL(dnrm2)(&N,y_j,&inci);
							//r_ij 	= pow(norm2,2);
							
							//beta_ij             = a_iT*y_j;
							beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
							
							if (fabs(m_ij)<1e-10) //go to the linear equation 
							{
								//
								if(verbose>7) Rprintf("\t\t\t\t\t\t\t gene %d \t interact with gene %d.\tLinear equation\n",i,j);
								//
								Bij = (beta_ij-lambdaW)/r_ij;
								//Rprintf("\t\t gene (%d,\t%d): linear Bij %f\n",i,j,Bij);
				//
				//if (i==0) 
				//{
				//	Rprintf("\t\t\t beta_ij: %f;\t r_ij:%f; lambdaW: %f\n ",beta_ij,r_ij,lambdaW);
				//	printMat(a_iT,N,1);
					//printMat(y_j,N,1);
				//	printMat(eiB,M,1);
				//}
				//
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
								//Lss = quadraticLss(sigma2, N, d_ij, candsBij, r_ij, lambdaW,beta_ij);								
								//Lss = quadraticLssJuan(sigma2,N, d_ij, candsBij, a_iT,y_j,lambdaW);
//beta = F77_CALL(ddot)(&N, a_iT, &inci,a_iT, &incj);
								//norm2 = F77_CALL(dnrm2)(&N,a_iT,&inci);
								//beta 	= pow(norm2,2);
							
								//Lss = sigma2*N*log(fabs(d_ij - candsBij)+1e-16) - beta/2 -lambdaW*fabs(candsBij);
								//Lss = sigma2*N*log(fabs(d_ij)+1e-16) - beta/2;
								Lss = sigma2*N*log(fabs(d_ij)+1e-16);
								
								// a_iT = a_iT - candsBij*y_j 		daxpy(n, a, x, incx, y, incy) 	y := a*x + y							
								if (Bijpp>0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijpp, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijpp, a_iT,y_j,lambdaW);	
									//aiTQuad dcopy(n, x, incx, y, incy)  y= x
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijpp;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - beta/2 -lambdaW*fabs(Bijpp);
									LssCands = sigma2*N*log(fabs(d_ij - Bijpp)+1e-16) - r_ij*pow(Bijpp,2)/2 + beta_ij*Bijpp -lambdaW*fabs(Bijpp); 
									
									if(LssCands>Lss) 
									{
										candsBij = Bijpp;
										Lss 	= LssCands;
									}	
								}
								if (Bijpm>0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijpm, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijpm, a_iT,y_j,lambdaW);								
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijpm;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijpm)+1e-16) - beta/2 -lambdaW*fabs(Bijpm);
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
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijmm, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijmm, a_iT,y_j,lambdaW);								
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijmm;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - beta/2 -lambdaW*fabs(Bijmm);
									LssCands = sigma2*N*log(fabs(d_ij - Bijmm)+1e-16) - r_ij*pow(Bijmm,2)/2 + beta_ij*Bijmm -lambdaW*fabs(Bijmm);  
									if(LssCands>Lss) 
									{
										candsBij = Bijmm;
										Lss 	= LssCands;
									}	
								}
								if (Bijmp<0)
								{
									//LssCands = quadraticLss(sigma2, N, d_ij, Bijmp, r_ij, lambdaW,beta_ij);
									//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijmp, a_iT,y_j,lambdaW);									
									//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
									//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
									//alpha = -Bijmp;
									//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
									//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
									//beta 	= pow(norm2,2);
									//LssCands = sigma2*N*log(fabs(d_ij - Bijmp)+1e-16) - beta/2 -lambdaW*fabs(Bijmp);			
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
						
						//update QIBinv ith column updated; jth row: 
						//readPtr = &QIBinv[j];
						//F77_CALL(dscal)(&M,&ziDb,readPtr,&M);
						//QIBinv[j*M + i] = QIBinv[j*M + i]*ziDb;
						// adjugate of QIBinv = ziDb*(QIBinvadj - (-1)^kl) + (-1)^kl
						//for(k=0;k<M;k++) // kth row
						//{
						//	for(l=0;l<M;l++)//lth column
						//	{
						//		if(k!=j && l!=i) 
						///		{
						//			alpha = pow(-1.0,k+l);
						//			QIBinv[l*M+k] = ziDb*(QIBinv[l*M+k] -alpha) + alpha;
						//		}	
						//	}
						//}
						
					}//js_i >0
				}//j = 1:M	
				//f
				//BF1 = B(i,:)*F1(:,i)
				readPtr = &B[i];
				F77_CALL(dcopy)(&M,readPtr,&ldM,BiT,&inci);

				F1ptr = &F1[M*i];
				BF1 = F77_CALL(ddot)(&M, BiT, &inci,F1ptr, &incj);

				f[i] = f0[i] - BF1;
				ei[i] = 0; // re-set ei for next i
			}else//s[i]  no un-zero weight in this gene
			{
				readPtr = &B[i];
				//F77_CALL(dscal)(&M,&toyZero,readPtr,&ldM);
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
		
//Rprintf("FnormOld: %f; \t FnormChange: %f\n",FnormOld,FnormChange);
//MASTER IS HERE	IBinv will be used in Qij	
		// BLOCK COORDINATE ASCEND: Update IBinv
		UpdateIBinv(QIBinv, B,M);
		
		if(verbose>5) Rprintf("\t\t\t\t\t\tdelta_BF: %f\n",delta_BF);
		if(delta_BF<1e-3)		//break out
		{
			break;
		}

		
		
	}//while

	
	//QIBinv: compute not updated ones
	//F77_NAME(dgetrs)(const char* trans, const int* n, const int* nrhs,
	//	 const double* a, const int* lda, const int* ipiv,
	//	 double* b, const int* ldb, int* info);
	//printMat(IBinv,M,M);
	//if(anyUpdateRow>0) //QIBinv --> = Cij/det(I-B) != eye(M); ipiv not zeros
	//{	
	//	transa = 'N';
	//	for(i=0;i<M;i++)
	//	{
	//	//	if(s[i] ==0)
	//		{ 	
	//			info = 0;
	//			ei[i] = 1;
	//			zi = &QIBinv[i*M];
	//			F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
	//			// dgetrs
	//			F77_CALL(dgetrs)(&transa, &M, &nrhs,IBinv,&M,ipiv,zi,&ldb,&info); //column i updated
	//			ei[i] = 0;
	//		}
	//	}
	//}else
	//{
		//initialize IBinv
	//	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	//	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
	//}
	
	//printMat(QIBinv,M,M);
	
	
	if(verbose>3) Rprintf("\t\t\t\tCurrent lambda: %f;\t number of iteration is: %d.\tExiting Weighted_LassoSf\n\n",lambda, iter);
	
	//mue          = (IM-B)*meanY-bsxfun(@times,f,meanX);

	//IBinv = I -B
	//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
	//alpha = -1; 
	//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
	//diagonal + 1
	//for(j=0;j<M;j++) 
	//{
	//	index = j*M + j;
	//	IBinv[index] = 1 + IBinv[index];
	//	mueL[j] = -f[j]*meanX[j];
	//}
	//mue[i] = -f[i]*meanX[i];
	
	//	//dgemv mue = Ax + beta*mue
	//transa = 'N';
	//alpha = 1;
	//beta = 1;
	//ldk = M;
	//F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, IBinv, &ldM, meanY, &inci, &beta,mueL, &incj);

	//----------------------------------------------------------------------------------END OF LASSO
		//SET TO PART2: LASSO with lambda zero
	//double *B, double *f;
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
	//none zeros
	//double *s;
	//s = (double* ) Calloc(M, double);
	for(i=0;i<M;i++)
	{
		readPtr = &S[i]; //S[i,];
		s[i] = F77_CALL(dasum)(&M, readPtr, &ldM);
	}

	beta = 1;
	//f0, F1
	//ei = 0	
	// entering loop
	//loop starts here
	iter = 0;
	//double js_i, m_ij,B_old, beta_ij,r_ij;
	//int *ipiv;
	//ipiv = (int *) R_alloc(N,sizeof(int));
	//info = 0;
	
	//dynamic variable keep intermidiate values eiB = ei-BiT

	//quadratic function
	double theta_ij,k_ij,q_ij,Bijp, Bijm; //case (14)
	
	//converge of gene i
	
	//converge of while
	max_iter = max_iter/5;
	//linear system

	//
	//double *aiTQuad; 
	//aiTQuad 	= (double* ) Calloc(N, double);
	//double *zzi;
	//zzi 		= (double* ) Calloc(M, double);
	//zi = &zzi[0];
	while(iter < max_iter)
	{
		iter = iter + 1;
		//converge Bfold = [B f];
		F77_CALL(dcopy)(&MM,B,&inci,BfOld,&incj);
		//last column
		F1ptr = &BfOld[MM];
		F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
		//
//	printMat(ei,M,1);
//	printMat(F1,M,M);
		//
		
		// inner loop
		for(i=0;i<M;i++)
		{
			if(s[i] >0)
			{
				//
				if(verbose>6) Rprintf("\t\t updating gene %d \n",i);
				//
			ei[i] = 1;
				//zi
				//F77_CALL(dcopy)(&MM,B,&inci,IBinv,&incj);
				//alpha = -1; 
				//F77_CALL(dscal)(&MM,&alpha,IBinv,&inci); // dscal(n, a, x, incx) x = a*x
				//diagonal + 1
				//for(j=0;j<M;j++) 
				//{
				//	index = j*M + j;
				//	IBinv[index] = 1 + IBinv[index];
				//}
				//call matrix inverse
				//MatrixInverse(IBinv,M);
				//zi is ith column of IBinv
				//for(j =0;j<M;j++) zi[j] 	= IBinv[i*M+j];
				
				//by linear solver: 
				//copy ei to zi;
				zi = &IBinvZero[i*M];
				//F77_CALL(dcopy)(&M,ei,&inci,zi,&incj);
				//call linearSystem(double *a, int N,double *B) a is NxN
				//linearSystem(IBinv, M,zi);  //zi is updated.

//Rprintf("M: %d; irhs: %d; ldM: %d; ldb: %d; info: %d\n",M,nrhs,ldM,ldb,info);
				//F77_CALL(dgesv)(&M, &nrhs, IBinv, &ldM, ipiv, zi, &ldb, &info);
								
				//printMat(zi,M,1);
				//for j = js_i
				for(j=0;j<M;j++)
				{
					js_i = S[j*M + i]; 		//ith row
					if(js_i >0)
					{
						//if(verbose>4) Rprintf("\t\t\t gene %d \t interact with gene %d\n",i,j);
					
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
						//a_iT      = (ei'-BiT)*Y-f(i)*X(i,:);
						//eiB = ei-BiT			//daxpy(n, a, x, incx, y, incy) 		y := a*x + y
						alpha = 1;
						F77_CALL(daxpy)(&M, &alpha,BiT, &inci,eiB, &incj);
						//a_iT      = (ei'-BiT)*Y-f(i)*X(i,:);
						readPtr = &X[i];
						F77_CALL(dcopy)(&N,readPtr,&M,a_iT,&inci);
						//a_iT = -f[i]*xi 		dscal(n, a, x, incx) 		x = a*x
						alpha = -f[i];
						F77_CALL(dscal)(&N,&alpha,a_iT,&inci);							

						transa='T'; //y := alpha*A**T*x + beta*y, 		 dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
						//beta = 1;
						alpha = 1;
						F77_CALL(dgemv)(&transa, &M, &N,&alpha, Y, &ldM, eiB, &inci, &beta,a_iT, &incj);

						//r_ij                = y_j'*y_j;	
						r_ij = F77_CALL(ddot)(&N, y_j, &inci,y_j, &inci);
						//norm2 = F77_CALL(dnrm2)(&N,y_j,&inci);
						//r_ij 	= pow(norm2,2);	
						//beta_ij             = a_iT*y_j;
						beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
						
						if (fabs(m_ij)<1e-10) //go to the linear equation 
						{
							//
							if(verbose>7) Rprintf("\t\t\t gene %d \t interact with gene %d.\tLinear equation\n",i,j);
							//
							//Bij = (beta_ij+lambdaW)/r_ij; 
							B[j*M+i] = beta_ij/r_ij;
							//Rprintf("\t\t\t beta_ij: %f;\t r_ij:%f\n ",beta_ij,r_ij);
							//printMat(a_iT,M,1);
						}else //m_ij ~=0 go to the quadratic equation
						{
							//
							if(verbose>7) Rprintf("\t\t\t gene %d \t interact with gene %d.\tQuadratic equation\n",i,j);
							//					
						
							//assume Bij >0
							d_ij = 1/m_ij + B[j*M+i];
							theta_ij = r_ij*d_ij + beta_ij;
							k_ij = d_ij*beta_ij - N*sigma2;
								
							q_ij = theta_ij*theta_ij - 4*r_ij* k_ij;
							Bijp = (1/(2*r_ij))*(theta_ij + sqrt(q_ij));
							Bijm = (1/(2*r_ij))*(theta_ij - sqrt(q_ij));
								
							candsBij = 0;
							//Lss = quadraticLss(sigma2, N, d_ij, candsBij, r_ij, lambdaW,beta_ij);
							//Lss = quadraticLssJuan(sigma2,N, d_ij, candsBij, a_iT,y_j,lambdaW);
							//beta = F77_CALL(ddot)(&N, a_iT, &inci,a_iT, &inci);
							//norm2 = F77_CALL(dnrm2)(&N,a_iT,&inci);
							//beta 	= pow(norm2,2);
							//Lss = sigma2*N*log(fabs(d_ij - candsBij)+1e-16) - beta/2 -lambdaW*fabs(candsBij);
							//Lss = sigma2*N*log(fabs(d_ij)+1e-16) - beta/2;															
							Lss = sigma2*N*log(fabs(d_ij)+1e-16);
							//Bijp
							//LssCands = quadraticLss(sigma2, N, d_ij, Bijp, r_ij, lambdaW,beta_ij);
							//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijp, a_iT,y_j,lambdaW);							
							//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
							//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
							//alpha = -Bijp;
							//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
							//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &incj);
							//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
							//beta 	= pow(norm2,2);							
							//LssCands = sigma2*N*log(fabs(d_ij - Bijp)+1e-16) - beta/2;
							LssCands = sigma2*N*log(fabs(d_ij - Bijp)+1e-16) - r_ij*pow(Bijp,2)/2 + beta_ij*Bijp;
							if(LssCands>Lss) 
							{
								candsBij = Bijp;
								Lss 	= LssCands;
							}	
							//Bijm>
							//LssCands = quadraticLss(sigma2, N, d_ij, Bijm, r_ij, lambdaW,beta_ij);
							//LssCands = quadraticLssJuan(sigma2,N, d_ij, Bijm, a_iT,y_j,lambdaW);							
							//F77_CALL(dcopy)(&N,a_iT,&inci,aiTQuad,&incj);
							//aiTQuad = aiTQuad - candsBij*y_j  daxpy(n, a, x, incx, y, incy) 	y := a*x + y
							//alpha = -Bijm;
							//F77_CALL(daxpy)(&N, &alpha,y_j, &inci,aiTQuad, &incj);
							//beta = F77_CALL(ddot)(&N, aiTQuad, &inci,aiTQuad, &inci);
							//norm2 = F77_CALL(dnrm2)(&N,aiTQuad,&inci);
							//beta 	= pow(norm2,2);
							//LssCands = sigma2*N*log(fabs(d_ij - Bijm)+1e-16) - beta/2;							
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

						//update IBinvZero ith column updated; jth row: 
						//readPtr = &IBinvZero[j];
						//F77_CALL(dscal)(&M,&ziDb,readPtr,&M);
						//IBinvZero[j*M + i] = IBinvZero[j*M + i]*ziDb;
						// adjugate of QIBinv = ziDb*(QIBinvadj - (-1)^kl) + (-1)^kl
						//for(k=0;k<M;k++) // kth row
						//{
						//	for(l=0;l<M;l++)//lth column
						//	{
						//		if(k!=j && l!=i) 
						//		{
						//			alpha = pow(-1.0,k+l);
						//			IBinvZero[l*M+k] = ziDb*(IBinvZero[l*M+k] -alpha) + alpha;
						//		}	
						//	}
						//}


						
					}//js_i >0
				}//j = 1:M	
			
				//f
				//BF1 = B(i,:)*F1(:,i)
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
		
		
		// BLOCK COORDINATE ASCEND: Update IBinv
		UpdateIBinv(IBinvZero, B,M);
		
		if(delta_BF<1e-2)		//break out
		{
			break;
		}

	}//while
		
	
	if(verbose>3) Rprintf("\t number of iteration is: %d.\nExiting constrained_MLf\n",iter);

	//IBinv = I -B
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
	//mue[i] = -f[i]*meanX[i];
	//	//dgemv mue = Ax + beta*mue
	transa = 'N';
	alpha = 1;
	//beta = 1;
	ldk = M;
	F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, IBinv, &ldM, meanY, &inci, &beta,mue, &incj);
	
	
	//----------------------------------------------------------------------------------END OF LASSO_ZERO
	
	
	
	
	//mue                     	= (IM-B)*meanY-bsxfun(@times,f,meanX); diag(f)*meanX

	Free(meanY);
	Free(meanX);
	Free(Y);
	Free(X);
	
	Free(S);
	Free(s);
	Free(f0);
	Free(F1);
	Free(Wcopy);
	
	//Free(xi);
	Free(y_j);
	
	Free(ei);
	Free(IBinv);
	//Free(zzi);
	Free(a_iT);
	
	Free(eiB);
	Free(BiT);
	Free(BfOld);
	Free(BfNew);
	Free(BfChange);
	
	//Free(ipiv);
	//------adaEN Apr. 2013	
	Free(W);
	//------adaEN Apr. 2013		
	//Free(aiTQuad);
	//Free(eye);
	return lambda;
	//sigma2 remains the same
}//weighted_LassoSf

//--------------------------------------------------------------------------------END WEIGHTED_LASSOSF

//----------------------------------------------- TEST FUNCTION
int cv_gene_nets_support_adaEN(double *Y, double *X, int Kcv,double *lambda_factors, double *rho_factors, 
			int maxiter, int M, int N,int Nlambdas, int Nrho,int verbose,double *W, 			//double sigma2learnt,
			double *sigmaLasso,
			int i_alpha, double alpha_factor,double * ErrorEN,double *sigma2learnt_EN,
						double *ErrorEN_min, double* steEN_min)	 //version V1_1delta.c			//,double *IBinv	//------adaEN Apr. 2013	
{	//sigma2_ms ilambda_ms
	// 
	// return sigmaLASSO = sigma_RidgeCVf for lasso in main
	//return ilambda_cv_ms: index of lambda_factors
	//0. lambda(i) path:
	//1. 	cross validation by ridge_cvf: rho_factor 		--> ONLY  for ilambda = 1; save Kcv set;
	//2. 	constrained_ridge_eff: weights W
	//3. 	weighted_lassoSf: none-zeros with lambda(i)
	//4. 	constrained_MLf: likelihood
	//5. return lambda
	int MM = M*M;
	int MN = M*N;
	// to save for each lambda
	double *Q, *BL, *fL,*mueL; 
	int i,j,index;	
	Q = (double* ) Calloc(MM, double); //ptr Q
	BL = (double* ) Calloc(MM, double); //ptr BL
	fL = (double* ) Calloc(M, double); //ptr fL
	mueL = (double* ) Calloc(M, double);
	//
	double *BC, *fC;
	BC = (double* ) Calloc(MM, double);
	fC = (double* ) Calloc(M, double);
	
	int Ntest = N/Kcv; 		//C takes floor
	int Nlearn = N - Ntest; 
	double * Errs, *Sigmas2, *ErrorMean, *ErrorCV;
	
	if(Nrho>Nlambdas)
	{
		i = Nrho*Kcv;
	}else
	{
		i= Nlambdas*Kcv;
	}
	
	Errs = (double* ) Calloc(i, double);
	Sigmas2 = (double* ) Calloc(Nlambdas, double);
	ErrorMean= (double* ) Calloc(Nlambdas, double);
	i = Nlambdas*Kcv;
	ErrorCV = (double* ) Calloc(i, double);
	
	//parameters inside the loop
	double *Ylearn, *Xlearn, *Ytest,*Xtest;
	int NlearnM = Nlearn*M;
	int NtestM = Ntest*M;
	Ylearn = (double* ) Calloc(NlearnM, double);
	Xlearn = (double* ) Calloc(NlearnM, double);
	Ytest = (double* ) Calloc(NtestM, double);
	Xtest = (double* ) Calloc(NtestM, double);

	//centerYX function
	double *Ylearn_centered, *Xlearn_centered, *meanY, *meanX;
	Ylearn_centered = (double* ) Calloc(NlearnM, double);
	Xlearn_centered = (double* ) Calloc(NlearnM, double);
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	

	//main loop
	double err_mean;
	//initialize
	err_mean = 1e5;

	int ilambda, cv,testStart, testEnd;
	//min_attained,min_attained = 0;
	ilambda = 0;
	
	//Weighted_LassoSf
	double lambda,lambda_factor,lambda_factor_prev;
	lambda_factor_prev = 1;

	//SL
//	double sigma2Ridge; 	//check which one is better: sigma2Ridge of sigma_hat(lambda)
	double *SL;				//zero or not-zero
	SL = (double* ) Calloc(MM, double);

	//convergence
	//------adaEN Apr. 2013	
	if(verbose>1) Rprintf("\t\ti_alpha %d; \t Enter Function: cv_support. Nlambdas: %d; \t %d-fold cross validation.\n", i_alpha, Nlambdas,Kcv);
	if(verbose>4) Rprintf("\n\t\t\t\t\tEnter Function: ridge_cvf. %d-fold cross validation.\n", Kcv);
	//if(verbose>1) Rprintf("\t\tEnter Function: cv_support. Nlambdas: %d; \t %d-fold cross validation.\n", Nlambdas,Kcv);
	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	int lda,ldb,ldc,ldk;
	double *Xptr, *XsubPtr; //submatrix
	double alpha, beta;
//Rprintf("Initialized.\n");		


//-------------------------------------------------------------ridge cv
	double sigma2learnt; // sigma of ridge regression for lasso_cv_support
	if(verbose>4) Rprintf("\n\t\t\t\t\tEnter Function: ridge_cvf. %d-fold cross validation.\n", Kcv);
	int irho = 0;
	int min_attained = 0;
	double err_mean_prev = 0;
	double sigma2R,i_rho_factor;
	double *BR, *fR, *mueR;
	BR =&BL[0];
	fR = &fL[0];
	mueR= &mueL[0];//BORROWED; no return: reset in lambda CV
	
	//testNoise, ImB: shared with Ridge & Lasso CV
	double *testNoise,*ImB,*NOISE;
	NOISE =(double* ) Calloc(MN, double);
	//testNoise =(double* ) Calloc(NtestM, double);
	testNoise = &NOISE[0];
	ImB = (double* ) Calloc(MM, double);
	char transa = 'N'; 
	char transb = 'N';
	lda = M;
	ldb = M;
	ldc = M;
	ldk = M;
	double testNorm;
	//result of CV 
	double rho_factor; // no pointers needed

	
//---------------------adaEN Apr. 2013	
if(i_alpha ==0)	
{	
	//Ridge CV main loop
	while(irho<Nrho && min_attained ==0)
	{

		i_rho_factor = rho_factors[irho];
		err_mean_prev = err_mean;
		err_mean = 0;
		//cross validation
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>6) Rprintf("\t\t\t\t\t\t\t crossValidation %d/Kcv\n\n",cv);
			// start and end point
			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);
//Rprintf("\t\t\t testStart: %d \t End %d\n",testStart,testEnd);
			//assign submatrices
			//SubMatrix(X,Xtest,Xlearn,M,N,testStart,testEnd);
			Xptr = &X[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Xtest,&incj);
			if(testStart ==1)
			{
				Xptr = &X[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Xlearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &X[testEnd*M]; //index
				XsubPtr = &Xlearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
			}
			
			//SubMatrix(Y,Ytest,Ylearn,M,N,testStart,testEnd);
			Xptr = &Y[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Ytest,&incj);
			if(testStart ==1)
			{
				Xptr = &Y[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Ylearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &Y[testEnd*M]; //index
				XsubPtr = &Ylearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
			}
			
			
			//SubMatrixInt(Missing,Missing_test,Missing_learn,M,N,testStart,testEnd);
		
			
			// ridge SEM			
			sigma2R = constrained_ridge_cff(Ylearn, Xlearn, i_rho_factor, M, Nlearn,BR,fR,mueR,verbose);

			//noise error
			F77_CALL(dcopy)(&MM,BR,&inci,ImB,&incj);
			alpha = -1;
			F77_CALL(dscal)(&MM,&alpha,ImB,&inci);
			for(i=0;i<M;i++) 
			{
				index = i*M + i;
				ImB[index] = 1 + ImB[index];
			} //I-BR
			alpha = 1; 
			beta = 0;			
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				XsubPtr = &testNoise[i];
				//row i of noise
				alpha = -fR[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &ldk,XsubPtr, &M);
				//NOISE[i,1:N]	
			}//row i = 1:M
			
			// y = ax + y  F77_CALL(daxpy)(&MNtest, &alpha,Xptr, &inci,XsubPtr, &incj);
			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				Xptr = &testNoise[i*M];
				F77_CALL(daxpy)(&M, &alpha,mueR, &inci,Xptr, &incj);
			}

			//testNorm = FrobeniusNorm(testNoise, M, N);
			//i = M*Ntest;
			//testNorm = F77_CALL(dnrm2)(&i,testNoise,&inci);
			//testNorm = testNorm*testNorm;
			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);
			//testNorm = F77_CALL(dnrm2)(&NtestM,testNoise,&inci);
			//testNorm = testNorm*testNorm;
			Errs[cv*Nrho+irho] = testNorm;
	
			
			
			err_mean = err_mean + testNorm;
//Rprintf("Errmean: %f; \t sigma2R: %f\n",err_mean, sigma2R);
		}//cv = 0: Kcv
		//err_mean = err_mean/Kcv;   	calculate sum instead
		//check convergence
		
		irho = irho + 1;
		if(verbose>5) Rprintf("\t\t\t\t\t\t Test rho ratio %f; (%d/Nrho)\t error: %f.\n",i_rho_factor,irho, err_mean);	
		if(irho>1)
		{
			if(err_mean_prev<=err_mean)
			{
				min_attained = 1;
				irho = irho -1;
			}
		}
			
	}//while
	
	//sum(sum(1-Missing))
	//int Nmiss = 0;
	//for(i=0;i<MN;i++) Nmiss = Nmiss + Missing[i];
	if(verbose>4) Rprintf("\t\t\t\t\tExit RidgeCV. sigma2R: %f\t",sigma2R);	
	
	//int Npresent = MN- Nmiss;
	if(irho == Nrho || irho == 1) //nonstop in while loop
	{	
		//sigma2ridge             = sum(Errs(irho,:))/(sum(sum(1-Missing))-1);
		sigma2R = err_mean/(MN -1);
	}else
	{
		sigma2R = err_mean_prev/(MN -1);
	}
	//rho_factor_m            = rho_factors(irho)*N/(N-Ntest);
	rho_factor = rho_factors[irho-1]*N/(N-Ntest);
	if(verbose>4) Rprintf("sigma2learnt: %f\n",sigma2R);	
	//sigma2R is for LASSO CV;
	//rho_factor_m: for ridge regression: use the whole dataset
	if(verbose==0) Rprintf("Step 1: ridge CV; find rho : %f\n", rho_factor);
	sigma2learnt = constrained_ridge_cff(Y, X, rho_factor, M, N,BR,fR,mueR,verbose);	
	if(verbose==0) Rprintf("Step 2: ridge; calculate weights.\n");
	// weight W: for CV(this function) and selection path (Main)
	for(i=0;i<MM;i++) W[i] = 1/fabs(BL[i]+ 1e-10);
	
	//sigma2learnt = sigma2cr; //for lambda CV
	//sigma2R is used for path in Main
	sigmaLasso[0] = sigma2R;
	sigma2learnt_EN[0] = sigma2learnt;
	}else	//i_alpha ==0
	{
		sigma2R 		= sigmaLasso[0];
		sigma2learnt 	= sigma2learnt_EN[0];
	}
//-----------------------------------------------adaEN Apr. 2013	
//------------------------------------------------------------ridge cv end return sigma2R, weight	
	//IBinv: update in path
	double *IBinv,*IBinvZero,*lambda_Max;
	double *IBinvPath,*IBinvPathZero;
	irho  = MM*Kcv;
	IBinvPath = (double* ) Calloc(irho, double);
	IBinvPathZero = (double* ) Calloc(irho, double);
	lambda_Max = (double* ) Calloc(Kcv, double);
	double lambda_max_cv;
	beta = 0;
	//initialized
	
	for(cv=0;cv<Kcv;cv++)
	{
		IBinv = &IBinvPath[cv*MM];
		F77_CALL(dcopy)(&MM,&beta,&inc0,IBinv,&incj);
		for(i=0;i<M;i++) IBinv[i*M+i] =1;
	}
	F77_CALL(dcopy)(&irho,IBinvPath,&inci,IBinvPathZero,&incj);
	
	while(ilambda < Nlambdas)
	{	
		//ilambda = ilambda  +1;
		//err_mean_prev = err_mean;
		//err_sigma_prev = err_sigma;
		err_mean = 0;
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>3) Rprintf("\t\t %d/Kcv cross validation.\n", cv);
			// test, learn
			// start and end point
			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);
			//assign submatrices
			//SubMatrix(X,Xtest,Xlearn,M,N,testStart,testEnd);
			//SubMatrix(Y,Ytest,Ylearn,M,N,testStart,testEnd);
			//SubMatrixInt(Missing,Missing_test,Missing_learn,M,N,testStart,testEnd);
			//SubMatrix(X,Xtest,Xlearn,M,N,testStart,testEnd);
			Xptr = &X[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Xtest,&incj);
			if(testStart ==1)
			{
				Xptr = &X[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Xlearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &X[testEnd*M]; //index
				XsubPtr = &Xlearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
			}
			
			//SubMatrix(Y,Ytest,Ylearn,M,N,testStart,testEnd);
			Xptr = &Y[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Ytest,&incj);
			if(testStart ==1)
			{
				Xptr = &Y[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Ylearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &Y[testEnd*M]; //index
				XsubPtr = &Ylearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
			}
			
			
			//SubMatrixInt(Missing,Missing_test,Missing_learn,M,N,testStart,testEnd);
			//Missing_test = &Missing[(testStart-1)*M];
			//Learn matrix
			
			//Ylearn_centered
			F77_CALL(dcopy)(&NlearnM,Xlearn,&inci,Xlearn_centered,&incj);
			F77_CALL(dcopy)(&NlearnM,Ylearn,&inci,Ylearn_centered,&incj);

			centerYX(Ylearn_centered,Xlearn_centered,meanY, meanX,M, Nlearn);
			//first ilambda

//Rprintf("Initialized.\n");			
			if(ilambda == 0)
			{
				//if(verbose>3) Rprintf("\t\t\t step 0 for first lambda: Ridge cross validation; \tRidge weights.\n");

				//for(i=0;i<M;i++) fL[i] = 1;
				alpha = 1; 
				F77_CALL(dcopy)(&M,&alpha,&inc0,fL,&incj); // call dcopy(n, x, inci, y, incy)
				alpha = 0;
				F77_CALL(dcopy)(&MM,&alpha,&inc0,BL,&incj); 
				//F77_CALL(dscal)(&MM,&alpha,BL,&inci);
				// Q(lambda)S4
				QlambdaStart(Ylearn_centered,Xlearn_centered, Q, sigma2learnt,M, Nlearn);
				// set BL, fL to  zeros: not necessary; will test this after Debug
			
				//lambda_Max[cv]	= lambdaMax(Ylearn_centered,Xlearn_centered,W,M, Nlearn);
				lambda_Max[cv]	= lambdaMax_adaEN(Ylearn_centered,Xlearn_centered,W,M, Nlearn,alpha_factor); 	//------adaEN Apr. 2013	
			}//ilambda ==0
			lambda_max_cv = lambda_Max[cv];
			//if(ilambda > 0) sigma2learnt = Sigmas2[cv*Nlambdas];//in Juan: 5 sigma2learnt saved: one for each fold
			//Weighted_LassoSf
			lambda_factor = lambda_factors[ilambda];
			//SIGMA2[0] = sigma2learnt;
			
			//lambda = Weighted_LassoSf(W, BL, fL, Ylearn, Xlearn, Q, lambda_factor,
			//				lambda_factor_prev, sigma2learnt, maxiter, M, Nlearn, verbose,IBinv);
							
			IBinv = &IBinvPath[cv*MM];
			IBinvZero= &IBinvPathZero[cv*MM];
			lambda = Weighted_LassoSf_MLf_adaEN(W, BL, fL, Ylearn,Xlearn, Q, lambda_factor, 
							lambda_factor_prev, sigma2learnt, maxiter,M, Nlearn, verbose,
							BC, fC, mueL,IBinv,IBinvZero,lambda_max_cv,
							alpha_factor); 								//------adaEN Apr. 2013	
//if(ilambda==0 && cv==0)
//{
//	Rprintf("In cv support: \n");
//	printMat(IBinv,M,M);
//	printMat(fL,1,M);
//printMat(BL,M,M);
	//printMat(fC,1,M);
	//printMat(mueL,1,M);
//}

							
			if(verbose>3) Rprintf("\t\t\t step 1 SML lasso regression, lambda: %f.\n",lambda);
			//sign of B SL		//SL int:	
			//F77_CALL(dscal)(&MM,&alpha,SL,&inci);

			//Q(lambda)
			//QlambdaMiddle(Ylearn,Xlearn, Q,BL,fL, mueL, sigma2learnt,M, Nlearn);
			QlambdaMiddleCenter(Ylearn_centered,Xlearn_centered, Q,BL,fL,sigma2learnt,M, Nlearn,IBinv);
//if(ilambda==0 && cv==0)
//{
	//Rprintf("In cv support: \n");
	//printMat(Q,M,M);
//	printMat(fL,1,M);
	//printMat(BC,M,M);
	//printMat(fC,1,M);
	//printMat(mueL,1,M);
//}			
			
			// constrained_MLf
			if(verbose>3) Rprintf("\t\t\t step 2 SML ZeroRegression.\n");
			

			//constrained_MLf(BC, fC, SL, Ylearn,Xlearn, sigma2learnt, maxiter, mueL,M, Nlearn, verbose);

			//noise error Errs[ilambda, cv]
			
			//error function: norm((1-Missing).*(A*Y-fC*X)-mueC),'fro')^2;
//Problem in JUAN'S CODE: FOR CV mean: need to normalize, ow. it is not a standard error.
//Errs[index] = errorFunc(BC,Ytest, fC, Xtest, mueL, M, Ntest);
			
			F77_CALL(dcopy)(&MM,BC,&inci,ImB,&incj);
			alpha = -1;
			F77_CALL(dscal)(&MM,&alpha,ImB,&inci);
			for(i=0;i<M;i++) 
			{
				index = i*M + i;
				ImB[index] = 1 + ImB[index];
			} //I-BR
			alpha = 1; 
			beta = 0;
			testNoise = &NOISE[NtestM*cv];
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				//XsubPtr = &testNoise[i];
				XsubPtr = &NOISE[NtestM*cv+i];
				//F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
				//row i of noise ldM
				alpha = -fC[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &M,XsubPtr, &ldk);
				//NOISE[i,1:N]	
			}//row i = 1:M
			
			// y = ax + y  F77_CALL(daxpy)(&MNtest, &alpha,Xptr, &inci,XsubPtr, &incj);
			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				//Xptr = &testNoise[i*M];
				Xptr = &NOISE[NtestM*cv+i*M];
				F77_CALL(daxpy)(&M, &alpha,mueL, &inci,Xptr, &incj);
			}

			//testNorm = FrobeniusNorm(testNoise, M, N);
			// TEST_NOISE computed;
			
			
			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);
			//testNorm = F77_CALL(dnrm2)(&NtestM,testNoise,&inci); //just dotproduct will work
			//testNorm = testNorm*testNorm;
			index = cv*Nlambdas+ilambda;
			
			Errs[index] = testNorm;
	
	
//			Sigmas2[index] = sigma2learnt;
			ErrorCV[index] = Errs[index]/NtestM;
			err_mean = err_mean + Errs[index];
			if(verbose>3) Rprintf("\t\t\t cv: %d \tend; err_mean = %f.\n", cv, Errs[index]);
				
		}//cv
		//err
		err_mean = err_mean/MN;
		ErrorMean[ilambda] = err_mean;
		
		//mean, sigma
		for(i=0;i<MN;i++)
		{
			NOISE[i] = pow(NOISE[i],2);
		}
		alpha = -1;
		F77_CALL(daxpy)(&MN,&alpha,&err_mean,&inc0,NOISE,&inci);
		testNorm = F77_CALL(dnrm2)(&MN,NOISE,&inci);
		Sigmas2[ilambda] = testNorm/sqrt(Kcv*(MN -1));
		
		
/*		if(minimumErr>err_mean) 
		{
			minimumErr = err_mean;
			ilambda_min = ilambda;
		}
*/		//check for convergence
		//ilambda = ilambda + 1; //don't move this line
		//if(ilambda>3 && min_attained==0)
		//{
		//	if ((err_mean_prev + err_sigma_prev) < err_mean)
		//	{
		//		min_attained = 1;
		//		ilambda_min = ilambda; //number; not for indexing
		//	}
		//}
		lambda_factor_prev = lambda_factor;
		
		if(verbose>2) Rprintf("\t\t\t %d/Nlambdas. %d fold cv; \t Err_Mean: %f; std:%f; \t sigma2learnt:%f.\n", ilambda,Kcv,err_mean,Sigmas2[ilambda],sigma2learnt);
		ilambda = ilambda + 1; 
	}//while ilambda

	// select ilambda_ms
	
	// step1. calculate sigma2
	//ErrorCV = ErrorCV - ErrorMean
//printMat(Sigmas2,Nlambdas,1);	
	// step2. index of minimal ErrorMean
	double *errorMeanCpy;
	errorMeanCpy = (double* ) Calloc(Nlambdas, double);
	//int inc0  = 0;
	double minimumErr = -1e5 - MN;
	F77_CALL(dcopy)(&Nlambdas,&minimumErr,&inc0,errorMeanCpy,&inci);
	alpha = 1;
	F77_CALL(daxpy)(&Nlambdas,&alpha,ErrorMean,&inci,errorMeanCpy,&incj); // y = ax + y
	int ilambda_ms;	
	ilambda_ms = F77_CALL(idamax)(&Nlambdas, errorMeanCpy, &inci);//index of the max(errs_mean)<--min(mean)
	index = ilambda_ms - 1;
	minimumErr = ErrorMean[index] + Sigmas2[index]; //actually max
//printMat(ErrorMean,1,Nlambdas);
//printMat(Sigmas2,1, Nlambdas);
//--------------------------------------------------------------R_package changes
//	double *readPtr1;

//	readPtr1 = &mseStd[0];
//	F77_CALL(dcopy)(&Nlambdas,ErrorMean,&inci,readPtr1,&incj);
//	readPtr1 = &mseStd[Nlambdas];
//	F77_CALL(dcopy)(&Nlambdas,Sigmas2,&inci,readPtr1,&incj);
// version 2 of R_package:ｆollow EN_MPI_epsilon
//version 3.1
	double minDist = fabs(ErrorMean[index -1] - minimumErr);
	double distance;
	int tempIlambda_ms = index;
	for(i=index-1;i>0;i--) 
	{
		distance = fabs(ErrorMean[i] - minimumErr);
		if(distance <minDist)
		{
			minDist = distance;
			tempIlambda_ms = i + 1;
		}
	}
	ilambda_ms = tempIlambda_ms;
	
	//version 3.1: find the ilambda with ErrorMean closest to ErrorMean[ilambda_ms - 1]
	index = ilambda_ms - 1;
	ErrorEN[i_alpha] = ErrorMean[index];

//------adaEN Apr. 2013		
		//version V1_1epsilon.c
	ErrorEN_min[i_alpha] = ErrorMean[index];
	steEN_min[i_alpha] = Sigmas2[index];
	//version V1_1epsilon.c

//--------------------------------------------------------------R_package changes
//-------------------------------------------------------------adaEN Apr. 2013	
	//double lowBound;
//	for(i=index-1;i>0;i--) 
//	{
//		if(ErrorMean[i] < minimumErr)
//		{
//			ilambda_ms = i + 1;
//		}
//	}
	//return index of lambda_factors
	if(verbose>1) Rprintf("\t\tExit Function: cv_support. optimal lambda index: %d.\n\n", ilambda_ms);

//	Free(Ws);
	Free(NOISE);
	Free(ImB);
	
	
	Free(Q);
	Free(BL);
	Free(fL);
	Free(mueL);
	Free(BC);
	Free(fC);
	
	Free(Errs);
	Free(Sigmas2);
	Free(ErrorMean);
	Free(ErrorCV);
	Free(errorMeanCpy);
	
	Free(Ylearn);
	Free(Xlearn);
	Free(Ytest);
	Free(Xtest);	

	//Free(Missing_learn);
	//Free(Missing_test);
	Free(Ylearn_centered);
	Free(Xlearn_centered);	
	Free(meanY);
	Free(meanX);
	Free(SL);
//	Free(errs_mean);	
//	Free(errs_sigma);	

//	Free(rho_factor_ptr);


	Free(IBinvPath);
	Free(IBinvPathZero);
	
	Free(lambda_Max);
	
	return ilambda_ms;
	//index + 1 -->position; return position
	
	
	
} //end of cv_gene_nets_support		

			
			
			
//cv_gene_nets_support: return ilambda_cv_ms Missing


// ---------------------------- Main function that calls by R (interface)
//-----------------------------main function calls:
// centerYX
// cv_gene_nets_support
// ridge_cvf
// constraint - ridge_cff
// weighted_lassoSf



void mainSML_adaEN(double *Y, double *X, int *m, int *n, int *Missing,double*B, double *f,double*stat,int*VB)
{
	//stat contains: correct postive, true positve; false postive, positve detected; power; fdr 6x1
	// assume B is the true imput
	int M, N, i, j,index,verbose;
	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	M 			= m[0];
	N 			= n[0];	
	verbose 	= VB[0];
	int MN 		= M*N;
	int MM 		= M*M;
	double *Strue;
	Strue 		= (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,B,&inci,Strue,&incj);
	stat[1] 	= 0;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index = j*M  +i;
			//Strue[index] = B[index];
			if(i!=j && B[index]!=0)
			{	stat[1] = stat[1] + 1;} //stat[1] total positive
		}
	}
	double alpha = 1;	
	F77_CALL(dcopy)(&M,&alpha,&inc0,f,&inci);
	//for(i=0;i<M;i++) f[i] = 1;

	alpha = 0;
	//F77_CALL(dscal)(&MM,&alpha,B,&inci);
	F77_CALL(dcopy)(&MM,&alpha,&inc0,B,&inci);
	//assume missing values are from X		//set missing counterpart in Y to zero

	for(i=0;i<MN;i++)
	{
		if(Missing[i] == 1) Y[i] = 0;
	}
	
	//call cv_gene_nets_support ------------------------SYSTEM PARAMETERS
	int maxiter 	= 500;
	int Kcv 		= 5;
	int L_lambda 	= 20; // number of lambdas in lambda_factors	stop at 0.001
	double *lambda_factors;
	lambda_factors 	= (double* ) Calloc(L_lambda, double);
	double step 	= -0.2;
	for(i=0;i<L_lambda;i++) 
	{
		lambda_factors[i] 	= pow(10.0,step);
		step 				= step - 0.2;
	}
	
	//rho factor
	step 					= -6;
	double * rho_factors;
	int L 					= 31; // number of rho_factors
	rho_factors 			= (double* ) Calloc(L, double);
	for(i=0;i<L;i++) 
	{
		rho_factors[i] 		= pow(10.0,step);
		step 				= step + 0.2;
	}
	//------adaEN Apr. 2013	
	//adaEN parameter
	double *alpha_factors, *ErrorEN,*lambdaEN, sigma2learnt;
	int L_alpha  		= 19; 
	alpha_factors 		= (double* ) Calloc(L_alpha, double); //variable in all ranks
	ErrorEN 			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
	lambdaEN 			= (double* ) Calloc(L_alpha, double);
		//R_package Version 2: EN_MPI_epsilon
//version V1_1delta.c: record the minimum lambda point of Error +/- ste.
	double *ErrorEN_min, *steEN_min;
	ErrorEN_min			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
	steEN_min 			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
//version V1_1delta.c
	//R_package Version 2: EN_MPI_epsilon
	
	step 					= 0.05;
	for(i=0;i<L_alpha;i++) 
	{
		alpha_factors[i]	= 0.95 - step*i;
	}

//------adaEN Apr. 2013	
	
	
	
	
	
	//---------------------------------------------END  SYSTEM PARAMETERS
	int Nlambdas,Nrho;
	Nlambdas 				= L_lambda;
	Nrho 					= L;
//call ridge_cvf
	double sigma2; //return value;

	//double *mueL;
	//mueL = (double* ) Calloc(M, double); 

	// weight W
	double *W;  //weight on diagonal?
	W = (double* ) Calloc(MM, double);
	
	// IBinv: SAVE COMPUTATION: get IBinv from lasso: 1) CV_support; 2) selection path
	double *QIBinv;
	//IBinv = (double* ) Calloc(MM, double);
	QIBinv = (double* ) Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1
	
	//ilambda_cv_ms = cv_gene_nets_support(Y, X, Kcv,lambda_factors, 
	//				rho_factors, maxiter, M, N,Nlambdas, Nrho,verbose,W,sigma2cr);
	
	
	//------------------------------adaEN Apr. 2013		

	int i_alpha;
	double alpha_factor;
	for(i_alpha=0;i_alpha<L_alpha;i_alpha++)
	{	
		alpha_factor 		= alpha_factors[i_alpha];		
		
		ilambda_cv_ms = cv_gene_nets_support_adaEN(Y, X, Kcv,lambda_factors, rho_factors, 
			maxiter, M, N,Nlambdas, Nrho,verbose,W, &sigma2,
			i_alpha,alpha_factor,ErrorEN, &sigma2learnt,
			ErrorEN_min,steEN_min);	 //version V1_1delta.c	);	
			
		lambdaEN[i_alpha] 	= ilambda_cv_ms;
	}
	
	//find the min of ErrorEN;
	double minEN 			= ErrorEN[0];
	int ind_minEN 			= 0;
	for(i_alpha = 1;i_alpha<L_alpha;i_alpha++)
	{
		if(minEN > ErrorEN[i_alpha])
		{
			minEN 			= ErrorEN[i_alpha];
			ind_minEN 		= i_alpha;
		}	
	}
	
//R_package Version 2: EN_MPI_epsilon	
		//version V1_1delta.c
	double minErr_ = ErrorEN_min[ind_minEN] + steEN_min[ind_minEN];
	double distance;
	int temIndex = ind_minEN; //index
	for(i_alpha = ind_minEN-1;i_alpha>=0;i_alpha--)
	{
		distance = ErrorEN[i_alpha] - minErr_;
		if(distance <=0)
		{
			temIndex = i_alpha;
		}
	}
	ind_minEN = temIndex;
	//version V1_1delta.c
//R_package Version 2: EN_MPI_epsilon
	
	
	
	//sigma2 					= sigmaEN[ind_minEN];
	ilambda_cv_ms 			= lambdaEN[ind_minEN];
	alpha_factor  			= alpha_factors[ind_minEN];
	Rprintf("\tAdaptive_EN %d-fold CV, alpha: %f.\n", Kcv,alpha_factor);
	//printMat(ErrorEN,L_alpha,1);
	//--------------------------------adaEN Apr. 2013	
	
	//ilambda_cv_ms = 17;
	//sigma2 = 0.0156827270598132;
	//step = 1;
	//F77_CALL(dcopy)(&MM,&step,&inc0,W,&inci);
	
	
	
	if(verbose==0) Rprintf("Step 3: CV support; alpha: %f, number of lambda needed: %d\n", alpha_factor,ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	Ycopy = (double* ) Calloc(MN, double);
	Xcopy = (double* ) Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda; 
	double lambda_max;
	//lambda_max = lambdaMax(Ycopy,Xcopy,W,M, N);
	
	lambda_max = lambdaMax_adaEN(Ycopy,Xcopy,W,M, N,alpha_factor);
	
	//
	//for(i=0;i<M;i++) f[i] = 1;
	//F77_CALL(dscal)(&MM,&alpha,B,&inci);
	if(verbose==0) Rprintf("Step 4: lasso selection path.\n");
	//printMat(B,M,M);
	//ilambda_cv_ms = 0;
//printMat(IBinv,M,M);	
	for(ilambda = 0;ilambda<ilambda_cv_ms;ilambda++)
	//for(ilambda = 0;ilambda<1;ilambda++)
	{
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, ilambda_cv_ms,lambda_factor);
		// call Weighted_LassoSf		Missing,
		lambda = Weighted_LassoSf_adaEN(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max,
					alpha_factor); 	// mueL not calculated
//printMat(IBinv,M,M);			
if(verbose>0) Rprintf("\tlambda: %f\n", lambda);				
		// call QlambdaMiddle
		//QlambdaMiddle(Y,X, Q,B,f, mueL, sigma2,M,N);
		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];
//printMat(IBinv,M,M);
	}//ilambda; selection path
	//return B,f; 
	//correct positive
	//printMat(B,M,M);
	stat[0] = 0;// correct positive
	stat[2] = 0;//false positive
	stat[3] = 0;//positive detected
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index = j*M + i;
			if(Strue[index]==0 && B[index]!=0) stat[2] = stat[2] + 1;
			if(i!=j)
			{
				//stat[3]
				if(B[index]!=0) 
				{
					stat[3] = stat[3] + 1;
					//stat[0]
					if(Strue[index]!=0) stat[0] = stat[0] + 1;
				}
			}
		}
	}
	//power
	stat[4] = stat[0]/stat[1];
	stat[5] = stat[2]/stat[3];
	if(verbose==0) Rprintf("Step 5: Finish calculation; detection power in stat vector.\n");
	Free(Strue);
	Free(meanY);
	Free(meanX);
	Free(lambda_factors);
	Free(rho_factors);
	Free(Ycopy);
	Free(Xcopy);
	//Free(rho_factor_m);
	//Free(mueL);
	Free(W);
	//Free(IBinv);
	Free(QIBinv);
	Free(Q);
	
	//------adaEN Apr. 2013		
	Free(alpha_factors);
	Free(ErrorEN);
	Free(lambdaEN);
	
	Free(ErrorEN_min);
	Free(steEN_min);
//------adaEN Apr. 2013
	//-------------------------------- some function changes Y, X permantly.
}






//for R_package only
//ONLY difference is: errorMean, sigma2 is copied to mseStd to return to main.
//Donot combine this utility into code above: the above code is kept to be the same as plain C code
int cv_gene_nets_support_adaENcv(double *Y, double *X, int Kcv,double *lambda_factors, double *rho_factors, 
			int maxiter, int M, int N,int Nlambdas, int Nrho,int verbose,double *W, 			//double sigma2learnt,
			double *sigmaLasso,
			int i_alpha, double alpha_factor,double * ErrorEN,double *sigma2learnt_EN,//------adaEN Apr. 2013	
			double *mseStd, 		//
			double *ErrorEN_min, double* steEN_min)	 //version V1_1delta.c			
//------adaEN Apr. 2013	
// if i_alpha ==0: calculate W
// keep record of ErrorEN = errs_mean_min (the error in selecting lambda)
// i_alpha goes global, like the function call of cv_gene_nets_support_adaEN
//------adaEN Apr. 2013	
{	//sigma2_ms ilambda_ms
	// 
	// return sigmaLASSO = sigma_RidgeCVf for lasso in main
	//return ilambda_cv_ms: index of lambda_factors
	//0. lambda(i) path:
	//1. 	cross validation by ridge_cvf: rho_factor 		--> ONLY  for ilambda = 1; save Kcv set;
	//2. 	constrained_ridge_eff: weights W
	//3. 	weighted_lassoSf: none-zeros with lambda(i)
	//4. 	constrained_MLf: likelihood
	//5. return lambda
	int MM = M*M;
	int MN = M*N;
	// to save for each lambda
	double *Q, *BL, *fL,*mueL; 
	int i,j,index;	
	Q = (double* ) Calloc(MM, double); //ptr Q
	BL = (double* ) Calloc(MM, double); //ptr BL
	fL = (double* ) Calloc(M, double); //ptr fL
	mueL = (double* ) Calloc(M, double);
	//
	double *BC, *fC;
	BC = (double* ) Calloc(MM, double);
	fC = (double* ) Calloc(M, double);
	
	int Ntest = N/Kcv; 		//C takes floor
	int Nlearn = N - Ntest; 
	double * Errs, *Sigmas2, *ErrorMean, *ErrorCV;
	
	if(Nrho>Nlambdas)
	{
		i = Nrho*Kcv;
	}else
	{
		i= Nlambdas*Kcv;
	}
	
	Errs = (double* ) Calloc(i, double);
	Sigmas2 = (double* ) Calloc(Nlambdas, double);
	ErrorMean= (double* ) Calloc(Nlambdas, double);
	i = Nlambdas*Kcv;
	ErrorCV = (double* ) Calloc(i, double);
	
	//parameters inside the loop
	double *Ylearn, *Xlearn, *Ytest,*Xtest;
	int NlearnM = Nlearn*M;
	int NtestM = Ntest*M;
	Ylearn = (double* ) Calloc(NlearnM, double);
	Xlearn = (double* ) Calloc(NlearnM, double);
	Ytest = (double* ) Calloc(NtestM, double);
	Xtest = (double* ) Calloc(NtestM, double);

	//centerYX function
	double *Ylearn_centered, *Xlearn_centered, *meanY, *meanX;
	Ylearn_centered = (double* ) Calloc(NlearnM, double);
	Xlearn_centered = (double* ) Calloc(NlearnM, double);
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	

	//main loop
	double err_mean;
	//initialize
	err_mean = 1e5;

	int ilambda, cv,testStart, testEnd;
	//min_attained,min_attained = 0;
	ilambda = 0;
	
	//Weighted_LassoSf
	double lambda,lambda_factor,lambda_factor_prev;
	lambda_factor_prev = 1;

	//SL
//	double sigma2Ridge; 	//check which one is better: sigma2Ridge of sigma_hat(lambda)
	double *SL;				//zero or not-zero
	SL = (double* ) Calloc(MM, double);

	//convergence
	//------adaEN Apr. 2013	
	if(verbose>1) Rprintf("\t\ti_alpha %d; \t Enter Function: cv_support. Nlambdas: %d; \t %d-fold cross validation.\n", i_alpha, Nlambdas,Kcv);
	if(verbose>4) Rprintf("\n\t\t\t\t\tEnter Function: ridge_cvf. %d-fold cross validation.\n", Kcv);
	//if(verbose>1) Rprintf("\t\tEnter Function: cv_support. Nlambdas: %d; \t %d-fold cross validation.\n", Nlambdas,Kcv);
	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	int lda,ldb,ldc,ldk;
	double *Xptr, *XsubPtr; //submatrix
	double alpha, beta;
//Rprintf("Initialized.\n");		


//-------------------------------------------------------------ridge cv
	double sigma2learnt; // sigma of ridge regression for lasso_cv_support
	if(verbose>4) Rprintf("\n\t\t\t\t\tEnter Function: ridge_cvf. %d-fold cross validation.\n", Kcv);
	int irho = 0;
	int min_attained = 0;
	double err_mean_prev = 0;
	double sigma2R,i_rho_factor;
	double *BR, *fR, *mueR;
	BR =&BL[0];
	fR = &fL[0];
	mueR= &mueL[0];//BORROWED; no return: reset in lambda CV
	
	//testNoise, ImB: shared with Ridge & Lasso CV
	double *testNoise,*ImB,*NOISE;
	NOISE =(double* ) Calloc(MN, double);
	//testNoise =(double* ) Calloc(NtestM, double);
	testNoise = &NOISE[0];
	ImB = (double* ) Calloc(MM, double);
	char transa = 'N'; 
	char transb = 'N';
	lda = M;
	ldb = M;
	ldc = M;
	ldk = M;
	double testNorm;
	//result of CV 
	double rho_factor; // no pointers needed

	
//---------------------adaEN Apr. 2013	
if(i_alpha ==0)	
{	
	//Ridge CV main loop
	while(irho<Nrho && min_attained ==0)
	{

		i_rho_factor = rho_factors[irho];
		err_mean_prev = err_mean;
		err_mean = 0;
		//cross validation
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>6) Rprintf("\t\t\t\t\t\t\t crossValidation %d/Kcv\n\n",cv);
			// start and end point
			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);
//Rprintf("\t\t\t testStart: %d \t End %d\n",testStart,testEnd);
			//assign submatrices
			//SubMatrix(X,Xtest,Xlearn,M,N,testStart,testEnd);
			Xptr = &X[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Xtest,&incj);
			if(testStart ==1)
			{
				Xptr = &X[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Xlearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &X[testEnd*M]; //index
				XsubPtr = &Xlearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
			}
			
			//SubMatrix(Y,Ytest,Ylearn,M,N,testStart,testEnd);
			Xptr = &Y[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Ytest,&incj);
			if(testStart ==1)
			{
				Xptr = &Y[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Ylearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &Y[testEnd*M]; //index
				XsubPtr = &Ylearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
			}
			
			
			//SubMatrixInt(Missing,Missing_test,Missing_learn,M,N,testStart,testEnd);
		
			
			// ridge SEM			
			sigma2R = constrained_ridge_cff(Ylearn, Xlearn, i_rho_factor, M, Nlearn,BR,fR,mueR,verbose);

			//noise error
			F77_CALL(dcopy)(&MM,BR,&inci,ImB,&incj);
			alpha = -1;
			F77_CALL(dscal)(&MM,&alpha,ImB,&inci);
			for(i=0;i<M;i++) 
			{
				index = i*M + i;
				ImB[index] = 1 + ImB[index];
			} //I-BR
			alpha = 1; 
			beta = 0;			
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				XsubPtr = &testNoise[i];
				//row i of noise
				alpha = -fR[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &ldk,XsubPtr, &M);
				//NOISE[i,1:N]	
			}//row i = 1:M
			
			// y = ax + y  F77_CALL(daxpy)(&MNtest, &alpha,Xptr, &inci,XsubPtr, &incj);
			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				Xptr = &testNoise[i*M];
				F77_CALL(daxpy)(&M, &alpha,mueR, &inci,Xptr, &incj);
			}

			//testNorm = FrobeniusNorm(testNoise, M, N);
			//i = M*Ntest;
			//testNorm = F77_CALL(dnrm2)(&i,testNoise,&inci);
			//testNorm = testNorm*testNorm;
			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);
			//testNorm = F77_CALL(dnrm2)(&NtestM,testNoise,&inci);
			//testNorm = testNorm*testNorm;
			Errs[cv*Nrho+irho] = testNorm;
	
			
			
			err_mean = err_mean + testNorm;
//Rprintf("Errmean: %f; \t sigma2R: %f\n",err_mean, sigma2R);
		}//cv = 0: Kcv
		//err_mean = err_mean/Kcv;   	calculate sum instead
		//check convergence
		
		irho = irho + 1;
		if(verbose>5) Rprintf("\t\t\t\t\t\t Test rho ratio %f; (%d/Nrho)\t error: %f.\n",i_rho_factor,irho, err_mean);	
		if(irho>1)
		{
			if(err_mean_prev<=err_mean)
			{
				min_attained = 1;
				irho = irho -1;
			}
		}
			
	}//while
	
	//sum(sum(1-Missing))
	//int Nmiss = 0;
	//for(i=0;i<MN;i++) Nmiss = Nmiss + Missing[i];
	if(verbose>4) Rprintf("\t\t\t\t\tExit RidgeCV. sigma2R: %f\t",sigma2R);	
	
	//int Npresent = MN- Nmiss;
	if(irho == Nrho || irho == 1) //nonstop in while loop
	{	
		//sigma2ridge             = sum(Errs(irho,:))/(sum(sum(1-Missing))-1);
		sigma2R = err_mean/(MN -1);
	}else
	{
		sigma2R = err_mean_prev/(MN -1);
	}
	//rho_factor_m            = rho_factors(irho)*N/(N-Ntest);
	rho_factor = rho_factors[irho-1]*N/(N-Ntest);
	if(verbose>4) Rprintf("sigma2learnt: %f\n",sigma2R);	
	//sigma2R is for LASSO CV;
	//rho_factor_m: for ridge regression: use the whole dataset
	if(verbose==0) Rprintf("Step 1: ridge CV; find rho : %f\n", rho_factor);
	sigma2learnt = constrained_ridge_cff(Y, X, rho_factor, M, N,BR,fR,mueR,verbose);	
	if(verbose==0) Rprintf("Step 2: ridge; calculate weights.\n");
	// weight W: for CV(this function) and selection path (Main)
	for(i=0;i<MM;i++) W[i] = 1/fabs(BL[i]+ 1e-10);
	
	//sigma2learnt = sigma2cr; //for lambda CV
	//sigma2R is used for path in Main
	sigmaLasso[0] = sigma2R;
	sigma2learnt_EN[0] = sigma2learnt;
	}else	//i_alpha ==0
	{
		sigma2R 		= sigmaLasso[0];
		sigma2learnt 	= sigma2learnt_EN[0];
	}
//-----------------------------------------------adaEN Apr. 2013	
//------------------------------------------------------------ridge cv end return sigma2R, weight	
	//IBinv: update in path
	double *IBinv,*IBinvZero,*lambda_Max;
	double *IBinvPath,*IBinvPathZero;
	irho  = MM*Kcv;
	IBinvPath = (double* ) Calloc(irho, double);
	IBinvPathZero = (double* ) Calloc(irho, double);
	lambda_Max = (double* ) Calloc(Kcv, double);
	double lambda_max_cv;
	beta = 0;
	//initialized
	
	for(cv=0;cv<Kcv;cv++)
	{
		IBinv = &IBinvPath[cv*MM];
		F77_CALL(dcopy)(&MM,&beta,&inc0,IBinv,&incj);
		for(i=0;i<M;i++) IBinv[i*M+i] =1;
	}
	F77_CALL(dcopy)(&irho,IBinvPath,&inci,IBinvPathZero,&incj);
	
	while(ilambda < Nlambdas)
	{	
		//ilambda = ilambda  +1;
		//err_mean_prev = err_mean;
		//err_sigma_prev = err_sigma;
		err_mean = 0;
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>3) Rprintf("\t\t %d/Kcv cross validation.\n", cv);
			// test, learn
			// start and end point
			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);
			//assign submatrices
			//SubMatrix(X,Xtest,Xlearn,M,N,testStart,testEnd);
			//SubMatrix(Y,Ytest,Ylearn,M,N,testStart,testEnd);
			//SubMatrixInt(Missing,Missing_test,Missing_learn,M,N,testStart,testEnd);
			//SubMatrix(X,Xtest,Xlearn,M,N,testStart,testEnd);
			Xptr = &X[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Xtest,&incj);
			if(testStart ==1)
			{
				Xptr = &X[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Xlearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &X[testEnd*M]; //index
				XsubPtr = &Xlearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,X,&inci,Xlearn,&incj);
			}
			
			//SubMatrix(Y,Ytest,Ylearn,M,N,testStart,testEnd);
			Xptr = &Y[(testStart-1)*M];//index
			i = (testEnd - testStart + 1)*M;//length
			F77_CALL(dcopy)(&i,Xptr,&inci,Ytest,&incj);
			if(testStart ==1)
			{
				Xptr = &Y[testEnd*M]; //index
				i = (N - testEnd)*M;//length
				F77_CALL(dcopy)(&i,Xptr,&inci,Ylearn,&incj);
			}else if(testEnd!=N) //two segments
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
				j = (N - testEnd)*M;//length
				Xptr = &Y[testEnd*M]; //index
				XsubPtr = &Ylearn[i];//index
				F77_CALL(dcopy)(&j,Xptr,&inci,XsubPtr,&incj);
			}else // 1- start is learn
			{
				i = (testStart-1)*M;//length
				F77_CALL(dcopy)(&i,Y,&inci,Ylearn,&incj);
			}
			
			
			//SubMatrixInt(Missing,Missing_test,Missing_learn,M,N,testStart,testEnd);
			//Missing_test = &Missing[(testStart-1)*M];
			//Learn matrix
			
			//Ylearn_centered
			F77_CALL(dcopy)(&NlearnM,Xlearn,&inci,Xlearn_centered,&incj);
			F77_CALL(dcopy)(&NlearnM,Ylearn,&inci,Ylearn_centered,&incj);

			centerYX(Ylearn_centered,Xlearn_centered,meanY, meanX,M, Nlearn);
			//first ilambda

//Rprintf("Initialized.\n");			
			if(ilambda == 0)
			{
				//if(verbose>3) Rprintf("\t\t\t step 0 for first lambda: Ridge cross validation; \tRidge weights.\n");

				//for(i=0;i<M;i++) fL[i] = 1;
				alpha = 1; 
				F77_CALL(dcopy)(&M,&alpha,&inc0,fL,&incj); // call dcopy(n, x, inci, y, incy)
				alpha = 0;
				F77_CALL(dcopy)(&MM,&alpha,&inc0,BL,&incj); 
				//F77_CALL(dscal)(&MM,&alpha,BL,&inci);
				// Q(lambda)S4
				QlambdaStart(Ylearn_centered,Xlearn_centered, Q, sigma2learnt,M, Nlearn);
				// set BL, fL to  zeros: not necessary; will test this after Debug
			
				//lambda_Max[cv]	= lambdaMax(Ylearn_centered,Xlearn_centered,W,M, Nlearn);
				lambda_Max[cv]	= lambdaMax_adaEN(Ylearn_centered,Xlearn_centered,W,M, Nlearn,alpha_factor); 	//------adaEN Apr. 2013	
			}//ilambda ==0
			lambda_max_cv = lambda_Max[cv];
			//if(ilambda > 0) sigma2learnt = Sigmas2[cv*Nlambdas];//in Juan: 5 sigma2learnt saved: one for each fold
			//Weighted_LassoSf
			lambda_factor = lambda_factors[ilambda];
			//SIGMA2[0] = sigma2learnt;
			
			//lambda = Weighted_LassoSf(W, BL, fL, Ylearn, Xlearn, Q, lambda_factor,
			//				lambda_factor_prev, sigma2learnt, maxiter, M, Nlearn, verbose,IBinv);
							
			IBinv = &IBinvPath[cv*MM];
			IBinvZero= &IBinvPathZero[cv*MM];
			lambda = Weighted_LassoSf_MLf_adaEN(W, BL, fL, Ylearn,Xlearn, Q, lambda_factor, 
							lambda_factor_prev, sigma2learnt, maxiter,M, Nlearn, verbose,
							BC, fC, mueL,IBinv,IBinvZero,lambda_max_cv,
							alpha_factor); 								//------adaEN Apr. 2013	
//if(ilambda==0 && cv==0)
//{
//	Rprintf("In cv support: \n");
//	printMat(IBinv,M,M);
//	printMat(fL,1,M);
//printMat(BL,M,M);
	//printMat(fC,1,M);
	//printMat(mueL,1,M);
//}

							
			if(verbose>3) Rprintf("\t\t\t step 1 SML lasso regression, lambda: %f.\n",lambda);
			//sign of B SL		//SL int:	
			//F77_CALL(dscal)(&MM,&alpha,SL,&inci);

			//Q(lambda)
			//QlambdaMiddle(Ylearn,Xlearn, Q,BL,fL, mueL, sigma2learnt,M, Nlearn);
			QlambdaMiddleCenter(Ylearn_centered,Xlearn_centered, Q,BL,fL,sigma2learnt,M, Nlearn,IBinv);
//if(ilambda==0 && cv==0)
//{
	//Rprintf("In cv support: \n");
	//printMat(Q,M,M);
//	printMat(fL,1,M);
	//printMat(BC,M,M);
	//printMat(fC,1,M);
	//printMat(mueL,1,M);
//}			
			
			// constrained_MLf
			if(verbose>3) Rprintf("\t\t\t step 2 SML ZeroRegression.\n");
			

			//constrained_MLf(BC, fC, SL, Ylearn,Xlearn, sigma2learnt, maxiter, mueL,M, Nlearn, verbose);

			//noise error Errs[ilambda, cv]
			
			//error function: norm((1-Missing).*(A*Y-fC*X)-mueC),'fro')^2;
//Problem in JUAN'S CODE: FOR CV mean: need to normalize, ow. it is not a standard error.
//Errs[index] = errorFunc(BC,Ytest, fC, Xtest, mueL, M, Ntest);
			
			F77_CALL(dcopy)(&MM,BC,&inci,ImB,&incj);
			alpha = -1;
			F77_CALL(dscal)(&MM,&alpha,ImB,&inci);
			for(i=0;i<M;i++) 
			{
				index = i*M + i;
				ImB[index] = 1 + ImB[index];
			} //I-BR
			alpha = 1; 
			beta = 0;
			testNoise = &NOISE[NtestM*cv];
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				//XsubPtr = &testNoise[i];
				XsubPtr = &NOISE[NtestM*cv+i];
				//F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
				//row i of noise ldM
				alpha = -fC[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &M,XsubPtr, &ldk);
				//NOISE[i,1:N]	
			}//row i = 1:M
			
			// y = ax + y  F77_CALL(daxpy)(&MNtest, &alpha,Xptr, &inci,XsubPtr, &incj);
			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				//Xptr = &testNoise[i*M];
				Xptr = &NOISE[NtestM*cv+i*M];
				F77_CALL(daxpy)(&M, &alpha,mueL, &inci,Xptr, &incj);
			}

			//testNorm = FrobeniusNorm(testNoise, M, N);
			// TEST_NOISE computed;
			
			
			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);
			//testNorm = F77_CALL(dnrm2)(&NtestM,testNoise,&inci); //just dotproduct will work
			//testNorm = testNorm*testNorm;
			index = cv*Nlambdas+ilambda;
			
			Errs[index] = testNorm;
	
	
//			Sigmas2[index] = sigma2learnt;
			ErrorCV[index] = Errs[index]/NtestM;
			err_mean = err_mean + Errs[index];
			if(verbose>3) Rprintf("\t\t\t cv: %d \tend; err_mean = %f.\n", cv, Errs[index]);
				
		}//cv
		//err
		err_mean = err_mean/MN;
		ErrorMean[ilambda] = err_mean;
		
		//mean, sigma
		for(i=0;i<MN;i++)
		{
			NOISE[i] = pow(NOISE[i],2);
		}
		alpha = -1;
		F77_CALL(daxpy)(&MN,&alpha,&err_mean,&inc0,NOISE,&inci);
		testNorm = F77_CALL(dnrm2)(&MN,NOISE,&inci);
		Sigmas2[ilambda] = testNorm/sqrt(Kcv*(MN -1));
		
		
/*		if(minimumErr>err_mean) 
		{
			minimumErr = err_mean;
			ilambda_min = ilambda;
		}
*/		//check for convergence
		//ilambda = ilambda + 1; //don't move this line
		//if(ilambda>3 && min_attained==0)
		//{
		//	if ((err_mean_prev + err_sigma_prev) < err_mean)
		//	{
		//		min_attained = 1;
		//		ilambda_min = ilambda; //number; not for indexing
		//	}
		//}
		lambda_factor_prev = lambda_factor;
		
		if(verbose>2) Rprintf("\t\t\t %d/Nlambdas. %d fold cv; \t Err_Mean: %f; std:%f; \t sigma2learnt:%f.\n", ilambda,Kcv,err_mean,Sigmas2[ilambda],sigma2learnt);
		ilambda = ilambda + 1; 
	}//while ilambda

	// select ilambda_ms
	
	// step1. calculate sigma2
	//ErrorCV = ErrorCV - ErrorMean
//printMat(Sigmas2,Nlambdas,1);	
	// step2. index of minimal ErrorMean
	double *errorMeanCpy;
	errorMeanCpy = (double* ) Calloc(Nlambdas, double);
	//int inc0  = 0;
	double minimumErr = -1e5 - MN;
	F77_CALL(dcopy)(&Nlambdas,&minimumErr,&inc0,errorMeanCpy,&inci);
	alpha = 1;
	F77_CALL(daxpy)(&Nlambdas,&alpha,ErrorMean,&inci,errorMeanCpy,&incj); // y = ax + y
	int ilambda_ms;	
	ilambda_ms = F77_CALL(idamax)(&Nlambdas, errorMeanCpy, &inci);//index of the max(errs_mean)<--min(mean)
	index = ilambda_ms - 1;
	minimumErr = ErrorMean[index] + Sigmas2[index]; //actually max
//printMat(ErrorMean,1,Nlambdas);
//printMat(Sigmas2,1, Nlambdas);

//-------------------------------------------------------------adaEN Apr. 2013	
//	
//--------------------------------------------------------------R_package changes
	double *readPtr1;

	readPtr1 = &mseStd[0];
	F77_CALL(dcopy)(&Nlambdas,ErrorMean,&inci,readPtr1,&incj);
	readPtr1 = &mseStd[Nlambdas];
	F77_CALL(dcopy)(&Nlambdas,Sigmas2,&inci,readPtr1,&incj);
// version 2 of R_package:ｆollow EN_MPI_epsilon
//version 3.1
	double minDist = fabs(ErrorMean[index -1] - minimumErr);
	double distance;
	int tempIlambda_ms = index;
	for(i=index-1;i>0;i--) 
	{
		distance = fabs(ErrorMean[i] - minimumErr);
		if(distance <minDist)
		{
			minDist = distance;
			tempIlambda_ms = i + 1;
		}
	}
	ilambda_ms = tempIlambda_ms;
	
	//version 3.1: find the ilambda with ErrorMean closest to ErrorMean[ilambda_ms - 1]
	index = ilambda_ms - 1;
	ErrorEN[i_alpha] = ErrorMean[index];

//------adaEN Apr. 2013		
		//version V1_1epsilon.c
	ErrorEN_min[i_alpha] = ErrorMean[index];
	steEN_min[i_alpha] = Sigmas2[index];
	//version V1_1epsilon.c

//--------------------------------------------------------------R_package changes	
	//double lowBound;
//	for(i=index-1;i>0;i--) 
//	{
//		if(ErrorMean[i] < minimumErr)
//		{
//			ilambda_ms = i + 1;
//		}
//	}
	//return index of lambda_factors
	if(verbose>1) Rprintf("\t\tExit Function: cv_support. optimal lambda index: %d.\n\n", ilambda_ms);

//	Free(Ws);
	Free(NOISE);
	Free(ImB);
	
	
	Free(Q);
	Free(BL);
	Free(fL);
	Free(mueL);
	Free(BC);
	Free(fC);
	
	Free(Errs);
	Free(Sigmas2);
	Free(ErrorMean);
	Free(ErrorCV);
	Free(errorMeanCpy);
	
	Free(Ylearn);
	Free(Xlearn);
	Free(Ytest);
	Free(Xtest);	

	//Free(Missing_learn);
	//Free(Missing_test);
	Free(Ylearn_centered);
	Free(Xlearn_centered);	
	Free(meanY);
	Free(meanX);
	Free(SL);
//	Free(errs_mean);	
//	Free(errs_sigma);	

//	Free(rho_factor_ptr);


	Free(IBinvPath);
	Free(IBinvPathZero);
	
	Free(lambda_Max);
	
	return ilambda_ms;
	//index + 1 -->position; return position
	
	
	
} //end of cv_gene_nets_support		

	












void mainSML_adaENcv(double *Y, double *X, int *m, int *n, int *Missing, double*B, double *f,double*stat,
			double*alpha_factors,int *nAlpha, 	// must be scalar
			double *lambda_factors, int *nLambda, double *mseStd, int*VB) 						// mseStd: nLmabda x 2 matrix, keep mse + std 
{
	//stat contains: correct postive, true positve; false postive, positve detected; power; fdr 6x1
	// assume B is the true imput
	int M, N, i, j,index,verbose;
	int inci 				= 1;
	int incj 				= 1;
	int inc0 				= 0;
	M 						= m[0];
	N 						= n[0];	
	verbose 				= VB[0];
	int MN 					= M*N;
	int MM 					= M*M;
	double *Strue;
	Strue 					= (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,B,&inci,Strue,&incj);
	stat[1] 				= 0;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index 			= j*M  +i;
			if(i!=j && B[index]!=0)
			{	stat[1] 	= stat[1] + 1;} //stat[1] total positive
		}
	}
	double alpha 			= 1;	
	F77_CALL(dcopy)(&M,&alpha,&inc0,f,&inci);

	alpha 					= 0;
	F77_CALL(dcopy)(&MM,&alpha,&inc0,B,&inci);

	for(i=0;i<MN;i++)
	{
		if(Missing[i] == 1) Y[i] = 0;
	}
	
	//call cv_gene_nets_support ------------------------SYSTEM PARAMETERS
	int maxiter 			= 500;
	int Kcv 				= 5;

	double step 	= -0.2;
	//rho factor
	step 					= -6;
	double * rho_factors;
	int L 					= 31; // number of rho_factors
	rho_factors 			= (double* ) Calloc(L, double);
	for(i=0;i<L;i++) 
	{
		rho_factors[i] 		= pow(10.0,step);
		step 				= step + 0.2;
	}
	//------adaEN Apr. 2013	
	//adaEN parameter
	double  *ErrorEN,*lambdaEN, sigma2learnt;	//*alpha_factors,
	int L_alpha  		= nAlpha[0]; 
	//alpha_factors 		= (double* ) Calloc(L_alpha, double); //variable in all ranks
	ErrorEN 			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
	lambdaEN 			= (double* ) Calloc(L_alpha, double);
	//R_package Version 2: EN_MPI_epsilon
//version V1_1delta.c: record the minimum lambda point of Error +/- ste.
	double *ErrorEN_min, *steEN_min;
	ErrorEN_min			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
	steEN_min 			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
//version V1_1delta.c
	//R_package Version 2: EN_MPI_epsilon
	
	step 					= 0.05;
//------adaEN Apr. 2013	
	
	
	
	
	
	//---------------------------------------------END  SYSTEM PARAMETERS
	int Nlambdas,Nrho;
	//----------------------------------------------------------R-package
	Nlambdas 				= nLambda[0]; 			

	//----------------------------------------------------------R-package	
	Nrho 					= L;
//call ridge_cvf
	double sigma2; //return value;

	//double *mueL;
	//mueL = (double* ) Calloc(M, double); 

	// weight W
	double *W;  //weight on diagonal?
	W = (double* ) Calloc(MM, double);
	
	// IBinv: SAVE COMPUTATION: get IBinv from lasso: 1) CV_support; 2) selection path
	double *QIBinv;
	//IBinv = (double* ) Calloc(MM, double);
	QIBinv = (double* ) Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1
	
	//ilambda_cv_ms = cv_gene_nets_support(Y, X, Kcv,lambda_factors, 
	//				rho_factors, maxiter, M, N,Nlambdas, Nrho,verbose,W,sigma2cr);
	
	
	//------------------------------adaEN Apr. 2013		

	int i_alpha;
	double alpha_factor;
	for(i_alpha=0;i_alpha<L_alpha;i_alpha++)
	{	
		alpha_factor 		= alpha_factors[i_alpha];		
		
		ilambda_cv_ms = cv_gene_nets_support_adaENcv(Y, X, Kcv,lambda_factors, rho_factors, 
			maxiter, M, N,Nlambdas, Nrho,verbose,W, &sigma2,
			i_alpha,alpha_factor,ErrorEN, &sigma2learnt,mseStd,
			ErrorEN_min,steEN_min);	 //version V1_1delta.c			
			
		lambdaEN[i_alpha] 	= ilambda_cv_ms;
	}
	
	//find the min of ErrorEN;
	double minEN 			= ErrorEN[0];
	int ind_minEN 			= 0;
	for(i_alpha = 1;i_alpha<L_alpha;i_alpha++)
	{
		if(minEN > ErrorEN[i_alpha])
		{
			minEN 			= ErrorEN[i_alpha];
			ind_minEN 		= i_alpha;
		}	
	}
	
	//R_package Version 2: EN_MPI_epsilon	
		//version V1_1delta.c
	double minErr_ = ErrorEN_min[ind_minEN] + steEN_min[ind_minEN];
	double distance;
	int temIndex = ind_minEN; //index
	for(i_alpha = ind_minEN-1;i_alpha>=0;i_alpha--)
	{
		distance = ErrorEN[i_alpha] - minErr_;
		if(distance <=0)
		{
			temIndex = i_alpha;
		}
	}
	ind_minEN = temIndex;
	//version V1_1delta.c
	//R_package Version 2: EN_MPI_epsilon	
	
	
	//sigma2 					= sigmaEN[ind_minEN];
	ilambda_cv_ms 			= lambdaEN[ind_minEN];
	alpha_factor  			= alpha_factors[ind_minEN];
	Rprintf("\tAdaptive_EN %d-fold CV, alpha: %f.\n", Kcv,alpha_factor);
	//printMat(ErrorEN,L_alpha,1);
	//--------------------------------adaEN Apr. 2013	

	
	if(verbose==0) Rprintf("Step 3: CV support; alpha: %f, number of lambda needed: %d\n", alpha_factor,ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	Ycopy = (double* ) Calloc(MN, double);
	Xcopy = (double* ) Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda;
	double lambda_max;
	//lambda_max = lambdaMax(Ycopy,Xcopy,W,M, N);
	
	lambda_max = lambdaMax_adaEN(Ycopy,Xcopy,W,M, N,alpha_factor);
	
	//
	//for(i=0;i<M;i++) f[i] = 1;
	//F77_CALL(dscal)(&MM,&alpha,B,&inci);
	if(verbose==0) Rprintf("Step 4: lasso selection path.\n");
	//printMat(B,M,M);
	//ilambda_cv_ms = 0;
//printMat(IBinv,M,M);	
	for(ilambda = 0;ilambda<ilambda_cv_ms;ilambda++)
	//for(ilambda = 0;ilambda<1;ilambda++)
	{
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, ilambda_cv_ms,lambda_factor);
		// call Weighted_LassoSf		Missing,
		lambda = Weighted_LassoSf_adaEN(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max,
					alpha_factor); 	// mueL not calculated
					if(verbose>0) Rprintf("\tlambda: %f\n", lambda);
//printMat(IBinv,M,M);							
		// call QlambdaMiddle
		//QlambdaMiddle(Y,X, Q,B,f, mueL, sigma2,M,N);
		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];
//printMat(IBinv,M,M);
	}//ilambda; selection path
	//return B,f; 
	//correct positive
	//printMat(B,M,M);
	stat[0] = 0;// correct positive
	stat[2] = 0;//false positive
	stat[3] = 0;//positive detected
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index = j*M + i;
			if(Strue[index]==0 && B[index]!=0) stat[2] = stat[2] + 1;
			if(i!=j)
			{
				//stat[3]
				if(B[index]!=0) 
				{
					stat[3] = stat[3] + 1;
					//stat[0]
					if(Strue[index]!=0) stat[0] = stat[0] + 1;
				}
			}
		}
	}
	//power
	stat[4] = stat[0]/stat[1];
	stat[5] = stat[2]/stat[3];
	if(verbose==0) Rprintf("Step 5: Finish calculation; detection power in stat vector.\n");
	Free(Strue);
	Free(meanY);
	Free(meanX);
	Free(rho_factors);
	Free(Ycopy);
	Free(Xcopy);
	Free(W);
	Free(QIBinv);
	Free(Q);
	
	//------adaEN Apr. 2013		
	Free(ErrorEN);
	Free(lambdaEN);
	
	Free(ErrorEN_min);
	Free(steEN_min);
//------adaEN Apr. 2013
	//-------------------------------- some function changes Y, X permantly.
}


//point estimate given 1 alpha, 1lambda, estimate B,F
//facilitate further computation such as:
//1.  	User decided lambda
//2.	Stability selection

void mainSML_adaENpointLmabda(double *Y, double *X, int *m, int *n, int *Missing, double*B, double *f,double*stat,
			double*alpha_factors,//int *nAlpha, 	 must be scalar
			double *lambda_factors, int *nLambda, double *mseStd, int*VB) 						// mseStd: nLmabda x 2 matrix, keep mse + std 
{
	//stat contains: correct postive, true positve; false postive, positve detected; power; fdr 6x1
	// assume B is the true imput
	int M, N, i, j,index,verbose;
	int inci 				= 1;
	int incj 				= 1;
	int inc0 				= 0;
	M 						= m[0];
	N 						= n[0];	
	verbose 				= VB[0];
	int MN 					= M*N;
	int MM 					= M*M;
	double *Strue;
	Strue 					= (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,B,&inci,Strue,&incj);
	stat[1] 				= 0;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index 			= j*M  +i;
			if(i!=j && B[index]!=0)
			{	stat[1] 	= stat[1] + 1;} //stat[1] total positive
		}
	}
	double alpha 			= 1;	
	F77_CALL(dcopy)(&M,&alpha,&inc0,f,&inci);

	alpha 					= 0;
	F77_CALL(dcopy)(&MM,&alpha,&inc0,B,&inci);

	for(i=0;i<MN;i++)
	{
		if(Missing[i] == 1) Y[i] = 0;
	}
	
	//call cv_gene_nets_support ------------------------SYSTEM PARAMETERS
	int maxiter 			= 500;
	int Kcv 				= 5;

	double step 	= -0.2;
	//rho factor
	step 					= -6;
	double * rho_factors;
	int L 					= 31; // number of rho_factors
	rho_factors 			= (double* ) Calloc(L, double);
	for(i=0;i<L;i++) 
	{
		rho_factors[i] 		= pow(10.0,step);
		step 				= step + 0.2;
	}
	//------adaEN Apr. 2013	
	//adaEN parameter
	double  *ErrorEN,*lambdaEN, sigma2learnt;	//*alpha_factors,
	int L_alpha  		= 1; 
	//alpha_factors 		= (double* ) Calloc(L_alpha, double); //variable in all ranks
	ErrorEN 			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
	lambdaEN 			= (double* ) Calloc(L_alpha, double);
	//R_package Version 2: EN_MPI_epsilon
//version V1_1delta.c: record the minimum lambda point of Error +/- ste.
	double *ErrorEN_min, *steEN_min;
	ErrorEN_min			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
	steEN_min 			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
//version V1_1delta.c
	//R_package Version 2: EN_MPI_epsilon
	step 					= 0.05;
//------adaEN Apr. 2013	
	
	
	
	
	
	//---------------------------------------------END  SYSTEM PARAMETERS
	int Nlambdas,Nrho;
	//----------------------------------------------------------R-package
	Nlambdas 				= 1; 			

	//----------------------------------------------------------R-package	
	Nrho 					= L;
//call ridge_cvf
	double sigma2; //return value;

	//double *mueL;
	//mueL = (double* ) Calloc(M, double); 

	// weight W
	double *W;  //weight on diagonal?
	W = (double* ) Calloc(MM, double);
	
	// IBinv: SAVE COMPUTATION: get IBinv from lasso: 1) CV_support; 2) selection path
	double *QIBinv;
	//IBinv = (double* ) Calloc(MM, double);
	QIBinv = (double* ) Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1
	
	//ilambda_cv_ms = cv_gene_nets_support(Y, X, Kcv,lambda_factors, 
	//				rho_factors, maxiter, M, N,Nlambdas, Nrho,verbose,W,sigma2cr);
	
	
	//------------------------------adaEN Apr. 2013		

	int i_alpha;
	double alpha_factor;
	for(i_alpha=0;i_alpha<L_alpha;i_alpha++)
	{	
		alpha_factor 		= alpha_factors[i_alpha];		
		
		ilambda_cv_ms = cv_gene_nets_support_adaENcv(Y, X, Kcv,lambda_factors, rho_factors, 
			maxiter, M, N,Nlambdas, Nrho,verbose,W, &sigma2,
			i_alpha,alpha_factor,ErrorEN, &sigma2learnt,mseStd, ErrorEN_min,steEN_min);	
			
		lambdaEN[i_alpha] 	= ilambda_cv_ms;
	}
	
	//find the min of ErrorEN;
	double minEN 			= ErrorEN[0];
	int ind_minEN 			= 0;
	for(i_alpha = 1;i_alpha<L_alpha;i_alpha++)
	{
		if(minEN > ErrorEN[i_alpha])
		{
			minEN 			= ErrorEN[i_alpha];
			ind_minEN 		= i_alpha;
		}	
	}
	
		//R_package Version 2: EN_MPI_epsilon	
		//version V1_1delta.c
	double minErr_ = ErrorEN_min[ind_minEN] + steEN_min[ind_minEN];
	double distance;
	int temIndex = ind_minEN; //index
	for(i_alpha = ind_minEN-1;i_alpha>=0;i_alpha--)
	{
		distance = ErrorEN[i_alpha] - minErr_;
		if(distance <=0)
		{
			temIndex = i_alpha;
		}
	}
	ind_minEN = temIndex;
	//version V1_1delta.c
	//R_package Version 2: EN_MPI_epsilon	
	
	
	//sigma2 					= sigmaEN[ind_minEN];
	ilambda_cv_ms 			= lambdaEN[ind_minEN];
	alpha_factor  			= alpha_factors[ind_minEN];
	Rprintf("\tAdaptive_EN %d-fold CV, alpha: %f.\n", Kcv,alpha_factor);
	//printMat(ErrorEN,L_alpha,1);
	//--------------------------------adaEN Apr. 2013	

	//----------------------------------------------------------R-package
	Nlambdas 				= nLambda[0]; 			

	//----------------------------------------------------------R-package	
	if(verbose==0) Rprintf("Step 3: CV support; alpha: %f, number of lambda needed: %d\n", alpha_factor,ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	Ycopy = (double* ) Calloc(MN, double);
	Xcopy = (double* ) Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda; 
	double lambda_max;
	//lambda_max = lambdaMax(Ycopy,Xcopy,W,M, N);
	
	lambda_max = lambdaMax_adaEN(Ycopy,Xcopy,W,M, N,alpha_factor);
	
	//
	//for(i=0;i<M;i++) f[i] = 1;
	//F77_CALL(dscal)(&MM,&alpha,B,&inci);
	if(verbose==0) Rprintf("Step 4: lasso selection path.\n");
	//printMat(B,M,M);
	//ilambda_cv_ms = 0;
//printMat(IBinv,M,M);	
	for(ilambda = 0;ilambda<Nlambdas;ilambda++)
	//for(ilambda = 0;ilambda<1;ilambda++)
	{
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, Nlambdas,lambda_factor);
		// call Weighted_LassoSf		Missing,
		lambda = Weighted_LassoSf_adaEN(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max,
					alpha_factor); 	// mueL not calculated
					if(verbose>0) Rprintf("\tlambda: %f\n", lambda);
//printMat(IBinv,M,M);							
		// call QlambdaMiddle
		//QlambdaMiddle(Y,X, Q,B,f, mueL, sigma2,M,N);
		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];
//printMat(IBinv,M,M);
	}//ilambda; selection path
	//return B,f; 
	//correct positive
	//printMat(B,M,M);
	stat[0] = 0;// correct positive
	stat[2] = 0;//false positive
	stat[3] = 0;//positive detected
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index = j*M + i;
			if(Strue[index]==0 && B[index]!=0) stat[2] = stat[2] + 1;
			if(i!=j)
			{
				//stat[3]
				if(B[index]!=0) 
				{
					stat[3] = stat[3] + 1;
					//stat[0]
					if(Strue[index]!=0) stat[0] = stat[0] + 1;
				}
			}
		}
	}
	//power
	stat[4] = stat[0]/stat[1];
	stat[5] = stat[2]/stat[3];
	if(verbose==0) Rprintf("Step 5: Finish calculation; detection power in stat vector.\n");
	Free(Strue);
	Free(meanY);
	Free(meanX);
	Free(rho_factors);
	Free(Ycopy);
	Free(Xcopy);
	Free(W);
	Free(QIBinv);
	Free(Q);
	
	//------adaEN Apr. 2013		
	Free(ErrorEN);
	Free(lambdaEN);
		Free(ErrorEN_min);
	Free(steEN_min);
	
//------adaEN Apr. 2013
	//-------------------------------- some function changes Y, X permantly.
}


//-------------------------------------------------------------------
//----------------------- Stability Selection -----------------------
//-----------------------     June 2013       -----------------------
//-------------------------------------------------------------------

// compute Bs for all {alpha, lambda} return to R
// R will repeat calling this functions for 100times, compute FDR Ev for all setups.
//from input: {alpha, lambda} are in descending order with among of shrinkage 
void mainSML_adaENstabilitySelection(double *Y, double *X, int *m, int *n, int *Missing, 
			double*B, double *f,double*stat,double*alpha_factors,int *nAlpha, 
			double *lambda_factors, int *nLambda, double *mseStd, int*VB,
			double *Bout) 						// mseStd: nLmabda x 2 matrix, keep mse + std 
{
	//stat contains: correct postive, true positve; false postive, positve detected; power; fdr 6x1
	// assume B is the true imput
	int M, N, i, j,index,verbose;
	int inci 				= 1;
	int incj 				= 1;
	int inc0 				= 0;
	M 						= m[0];
	N 						= n[0];	
	verbose 				= VB[0];
	int MN 					= M*N;
	int MM 					= M*M;
	double *Strue;
	Strue 					= (double* ) Calloc(MM, double);
	F77_CALL(dcopy)(&MM,B,&inci,Strue,&incj);
	stat[1] 				= 0;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index 			= j*M  +i;
			if(i!=j && B[index]!=0)
			{	stat[1] 	= stat[1] + 1;} //stat[1] total positive
		}
	}
	double alpha 			= 1;	
	F77_CALL(dcopy)(&M,&alpha,&inc0,f,&inci);

	alpha 					= 0;
	F77_CALL(dcopy)(&MM,&alpha,&inc0,B,&inci);

	for(i=0;i<MN;i++)
	{
		if(Missing[i] == 1) Y[i] = 0;
	}
	
	//call cv_gene_nets_support ------------------------SYSTEM PARAMETERS
	int maxiter 			= 500;
	int Kcv 				= 5;

	double step 	= -0.2;
	//rho factor
	step 					= -6;
	double * rho_factors;
	int L 					= 31; // number of rho_factors
	rho_factors 			= (double* ) Calloc(L, double);
	for(i=0;i<L;i++) 
	{
		rho_factors[i] 		= pow(10.0,step);
		step 				= step + 0.2;
	}
	//------adaEN Apr. 2013	
	//adaEN parameter
	double  *ErrorEN,*lambdaEN, sigma2learnt;	//*alpha_factors,
	int L_alpha  		= 1; 
	//alpha_factors 		= (double* ) Calloc(L_alpha, double); //variable in all ranks
	ErrorEN 			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
	lambdaEN 			= (double* ) Calloc(L_alpha, double);
	//R_package Version 2: EN_MPI_epsilon
//version V1_1delta.c: record the minimum lambda point of Error +/- ste.
	double *ErrorEN_min, *steEN_min;
	ErrorEN_min			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
	steEN_min 			= (double* ) Calloc(L_alpha, double);//rank0 variable actually.
//version V1_1delta.c
	//R_package Version 2: EN_MPI_epsilon
	
	
	step 					= 0.05;
//------adaEN Apr. 2013	
	
	
	
	
	
	//---------------------------------------------END  SYSTEM PARAMETERS
	int Nlambdas,Nrho;
	//----------------------------------------------------------R-package_1alpha point lambda
	Nlambdas 				= 1; 			

	//----------------------------------------------------------R-package_1alpha point lambda	
	Nrho 					= L;
//call ridge_cvf
	double sigma2; //return value;

	//double *mueL;
	//mueL = (double* ) Calloc(M, double); 

	// weight W
	double *W;  //weight on diagonal?
	W = (double* ) Calloc(MM, double);
	
	// IBinv: SAVE COMPUTATION: get IBinv from lasso: 1) CV_support; 2) selection path
	double *QIBinv;
	//IBinv = (double* ) Calloc(MM, double);
	QIBinv = (double* ) Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1
	
	//ilambda_cv_ms = cv_gene_nets_support(Y, X, Kcv,lambda_factors, 
	//				rho_factors, maxiter, M, N,Nlambdas, Nrho,verbose,W,sigma2cr);
	
	
	//------------------------------adaEN Apr. 2013		

	int i_alpha;
	double alpha_factor;
	for(i_alpha=0;i_alpha<L_alpha;i_alpha++)
	{	
		alpha_factor 		= alpha_factors[i_alpha];		
		
		ilambda_cv_ms = cv_gene_nets_support_adaENcv(Y, X, Kcv,lambda_factors, rho_factors, 
			maxiter, M, N,Nlambdas, Nrho,verbose,W, &sigma2,
			i_alpha,alpha_factor,ErrorEN, &sigma2learnt,mseStd, ErrorEN_min,steEN_min);	
			
		lambdaEN[i_alpha] 	= ilambda_cv_ms;
	}
	
	//find the min of ErrorEN;
	double minEN 			= ErrorEN[0];
	int ind_minEN 			= 0;
	for(i_alpha = 1;i_alpha<L_alpha;i_alpha++)
	{
		if(minEN > ErrorEN[i_alpha])
		{
			minEN 			= ErrorEN[i_alpha];
			ind_minEN 		= i_alpha;
		}	
	}
	
		//R_package Version 2: EN_MPI_epsilon	
		//version V1_1delta.c
	double minErr_ = ErrorEN_min[ind_minEN] + steEN_min[ind_minEN];
	double distance;
	int temIndex = ind_minEN; //index
	for(i_alpha = ind_minEN-1;i_alpha>=0;i_alpha--)
	{
		distance = ErrorEN[i_alpha] - minErr_;
		if(distance <=0)
		{
			temIndex = i_alpha;
		}
	}
	ind_minEN = temIndex;
	//version V1_1delta.c
	//R_package Version 2: EN_MPI_epsilon	
	
	
	//sigma2 					= sigmaEN[ind_minEN];
	ilambda_cv_ms 			= lambdaEN[ind_minEN];
	alpha_factor  			= alpha_factors[ind_minEN];
	Rprintf("\tAdaptive_EN %d-fold CV for ridge-weight, alpha: %f.\n", Kcv,alpha_factor);
	//printMat(ErrorEN,L_alpha,1);
	//--------------------------------adaEN Apr. 2013	

	//----------------------------------------------------------R-package
	Nlambdas 				= nLambda[0]; 			

	//----------------------------------------------------------R-package	
	if(verbose==0) Rprintf("Step 3: CV support; alpha: %f, number of lambda needed: %d\n", alpha_factor,ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) Calloc(M, double);
	meanX = (double* ) Calloc(M, double);
	Ycopy = (double* ) Calloc(MN, double);
	Xcopy = (double* ) Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda; lambda = 1.0;
	double lambda_max;
	//lambda_max = lambdaMax(Ycopy,Xcopy,W,M, N);
	//----------------------------------------------------------R-package_stability selection
	alpha_factor = 1;
	double *readPtr;
	//----------------------------------------------------------R-package_stability selection		
	lambda_max = lambdaMax_adaEN(Ycopy,Xcopy,W,M, N,alpha_factor);

	//F77_CALL(dscal)(&MM,&alpha,B,&inci);
	if(verbose==0) Rprintf("Step 4: lasso/elasticNet selection path.\n");
	for(ilambda = 0;ilambda<Nlambdas;ilambda++)
	{
		alpha_factor = alpha_factors[ilambda];
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, Nlambdas,lambda_factor);
		// call Weighted_LassoSf		Missing,
		lambda = Weighted_LassoSf_adaEN(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max,
					alpha_factor); 	// mueL not calculated
		if(verbose>0) Rprintf("\tlambda: %f\n", lambda);			
		// call QlambdaMiddle
		//QlambdaMiddle(Y,X, Q,B,f, mueL, sigma2,M,N);
		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];
	//----------------------------------------------------------R-package_stability selection		
		readPtr = &Bout[ilambda*MM];
		F77_CALL(dcopy)(&MM,B,&inci,readPtr,&incj);	
	//----------------------------------------------------------R-package_stability selection		
	}//ilambda; selection path
	//return B,f; 
	//correct positive
	//printMat(B,M,M);
	stat[0] = 0;// correct positive
	stat[2] = 0;//false positive
	stat[3] = 0;//positive detected
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index = j*M + i;
			if(Strue[index]==0 && B[index]!=0) stat[2] = stat[2] + 1;
			if(i!=j)
			{
				//stat[3]
				if(B[index]!=0) 
				{
					stat[3] = stat[3] + 1;
					//stat[0]
					if(Strue[index]!=0) stat[0] = stat[0] + 1;
				}
			}
		}
	}
	//power
	stat[4] = stat[0]/stat[1];
	stat[5] = stat[2]/stat[3];
	if(verbose==0) Rprintf("Step 5: Finish calculation; detection power in stat vector.\n");
	Free(Strue);
	Free(meanY);
	Free(meanX);
	Free(rho_factors);
	Free(Ycopy);
	Free(Xcopy);
	Free(W);
	Free(QIBinv);
	Free(Q);
	
	//------adaEN Apr. 2013		
	Free(ErrorEN);
	Free(lambdaEN);
	
		Free(ErrorEN_min);
	Free(steEN_min);
//------adaEN Apr. 2013
	//-------------------------------- some function changes Y, X permantly.
}




























































