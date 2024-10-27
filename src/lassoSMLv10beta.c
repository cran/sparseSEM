

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



//----------------------------------------------- TEST FUNCTION
int cv_gene_nets_support_Rdg(double *Y, double *X, int Kcv,double *lambda_factors, double *rho_factors, 
			int maxiter, int M, int N,int Nlambdas, int Nrho,int verbose,double *W, 			//double sigma2learnt,
			double *sigmaLasso)		//,double *IBinv
{	

	int MM = M*M;
	int MN = M*N;
	double *Q, *BL, *fL,*mueL; 
	int i,j,index;	
	Q = (double* ) R_Calloc(MM, double); //ptr Q
	BL = (double* ) R_Calloc(MM, double); //ptr BL
	fL = (double* ) R_Calloc(M, double); //ptr fL
	mueL = (double* ) R_Calloc(M, double);
	//
	double *BC, *fC;
	BC = (double* ) R_Calloc(MM, double);
	fC = (double* ) R_Calloc(M, double);
	
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
	
	Errs = (double* ) R_Calloc(i, double);
	Sigmas2 = (double* ) R_Calloc(Nlambdas, double);
	ErrorMean= (double* ) R_Calloc(Nlambdas, double);
	i = Nlambdas*Kcv;
	ErrorCV = (double* ) R_Calloc(i, double);
	
	//parameters inside the loop
	double *Ylearn, *Xlearn, *Ytest,*Xtest;
	int NlearnM = Nlearn*M;
	int NtestM = Ntest*M;
	Ylearn = (double* ) R_Calloc(NlearnM, double);
	Xlearn = (double* ) R_Calloc(NlearnM, double);
	Ytest = (double* ) R_Calloc(NtestM, double);
	Xtest = (double* ) R_Calloc(NtestM, double);

	//centerYX function
	double *Ylearn_centered, *Xlearn_centered, *meanY, *meanX;
	Ylearn_centered = (double* ) R_Calloc(NlearnM, double);
	Xlearn_centered = (double* ) R_Calloc(NlearnM, double);
	meanY = (double* ) R_Calloc(M, double);
	meanX = (double* ) R_Calloc(M, double);
	

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
	double *SL;				//zero or not-zero
	SL = (double* ) R_Calloc(MM, double);

	//convergence 
	
	if(verbose>1) Rprintf("\t\tEnter Function: cv_support. Nlambdas: %d; \t %d-fold cross validation.\n", Nlambdas,Kcv);
	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	int lda,ldb,ldc,ldk;
	double *Xptr, *XsubPtr; //submatrix
	double alpha, beta;


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
	mueR= &mueL[0];//
	

	double *testNoise,*ImB,*NOISE;
	NOISE =(double* ) R_Calloc(MN, double);

	testNoise = &NOISE[0];
	ImB = (double* ) R_Calloc(MM, double);
	char transa = 'N'; 
	char transb = 'N';
	lda = M;
	ldb = M;
	ldc = M;
	ldk = M;
	double testNorm;
	//result of CV 
	double rho_factor; // no pointers needed

	while(irho<Nrho && min_attained ==0)
	{

		i_rho_factor = rho_factors[irho];
		err_mean_prev = err_mean;
		err_mean = 0;
		//cross validation
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>6) Rprintf("\t\t\t\t\t\t\t crossValidation %d/Kcv\n\n",cv);
			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);

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
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc FCONE FCONE);
			for(i=0;i<M;i++)
			{
				Xptr = &Xtest[i];
				XsubPtr = &testNoise[i];
				//row i of noise
				alpha = -fR[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &ldk,XsubPtr, &M);
			}

			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				Xptr = &testNoise[i*M];
				F77_CALL(daxpy)(&M, &alpha,mueR, &inci,Xptr, &incj);
			}


			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);

			Errs[cv*Nrho+irho] = testNorm;
	
			
			
			err_mean = err_mean + testNorm;

		}

		
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
			
	}
	
	if(verbose>4) Rprintf("\t\t\t\t\tExit RidgeCV. sigma2R: %f\t",sigma2R);	
	
	if(irho == Nrho || irho == 1) //nonstop in while loop
	{	
		sigma2R = err_mean/(MN -1);
	}else
	{
		sigma2R = err_mean_prev/(MN -1);
	}

	rho_factor = rho_factors[irho-1]*N/(N-Ntest);
	if(verbose>4) Rprintf("sigma2learnt: %f\n",sigma2R);	

	if(verbose==0) Rprintf("Step 2: ridge CV; find rho : %f\n", rho_factor);
	sigma2learnt = constrained_ridge_cff(Y, X, rho_factor, M, N,BR,fR,mueR,verbose);	
	if(verbose==0) Rprintf("Step 3: ridge; calculate weights.\n");

	for(i=0;i<MM;i++) W[i] = 1/fabs(BL[i]+ 1e-10);
	

	sigmaLasso[0] = sigma2R;
	

	double *IBinv,*IBinvZero,*lambda_Max;
	double *IBinvPath,*IBinvPathZero;
	irho  = MM*Kcv;
	IBinvPath = (double* ) R_Calloc(irho, double);
	IBinvPathZero = (double* ) R_Calloc(irho, double);
	lambda_Max = (double* ) R_Calloc(Kcv, double);
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

		err_mean = 0;
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>3) Rprintf("\t\t %d/Kcv cross validation.\n", cv);

			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);

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
			

			
			//Ylearn_centered
			F77_CALL(dcopy)(&NlearnM,Xlearn,&inci,Xlearn_centered,&incj);
			F77_CALL(dcopy)(&NlearnM,Ylearn,&inci,Ylearn_centered,&incj);

			centerYX(Ylearn_centered,Xlearn_centered,meanY, meanX,M, Nlearn);

			
			if(ilambda == 0)
			{
				alpha = 1; 
				F77_CALL(dcopy)(&M,&alpha,&inc0,fL,&incj); // call dcopy(n, x, inci, y, incy)
				alpha = 0;
				F77_CALL(dcopy)(&MM,&alpha,&inc0,BL,&incj); 

				QlambdaStart(Ylearn_centered,Xlearn_centered, Q, sigma2learnt,M, Nlearn);

			
				lambda_Max[cv]	= lambdaMax(Ylearn_centered,Xlearn_centered,W,M, Nlearn);
			}
			lambda_max_cv = lambda_Max[cv];

			lambda_factor = lambda_factors[ilambda];

							
			IBinv = &IBinvPath[cv*MM];
			IBinvZero= &IBinvPathZero[cv*MM];
			lambda = Weighted_LassoSf_MLf(W, BL, fL, Ylearn,Xlearn, Q, lambda_factor, 
							lambda_factor_prev, sigma2learnt, maxiter,M, Nlearn, verbose,
							BC, fC, mueL,IBinv,IBinvZero,lambda_max_cv);


							
			if(verbose>3) Rprintf("\t\t\t step 1 SML lasso regression, lambda: %f.\n",lambda);

			QlambdaMiddleCenter(Ylearn_centered,Xlearn_centered, Q,BL,fL,sigma2learnt,M, Nlearn,IBinv);
			
			
			// constrained_MLf
			if(verbose>3) Rprintf("\t\t\t step 2 SML ZeroRegression.\n");
			

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
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc FCONE FCONE);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				XsubPtr = &NOISE[NtestM*cv+i];

				alpha = -fC[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &M,XsubPtr, &ldk);

			}

			alpha = -1;
			for(i=0;i<Ntest;i++)
			{

				Xptr = &NOISE[NtestM*cv+i*M];
				F77_CALL(daxpy)(&M, &alpha,mueL, &inci,Xptr, &incj);
			}

			
			
			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);

			index = cv*Nlambdas+ilambda;
			
			Errs[index] = testNorm;
	
	
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
		
		

		lambda_factor_prev = lambda_factor;
		
		if(verbose>2) Rprintf("\t\t\t %d/Nlambdas. %d fold cv; \t Err_Mean: %f; std:%f; \t sigma2learnt:%f.\n", ilambda,Kcv,err_mean,Sigmas2[ilambda],sigma2learnt);
		ilambda = ilambda + 1; 
	}
	
	double *errorMeanCpy;
	errorMeanCpy = (double* ) R_Calloc(Nlambdas, double);

	double minimumErr = -1e5 - MN;
	F77_CALL(dcopy)(&Nlambdas,&minimumErr,&inc0,errorMeanCpy,&inci);
	alpha = 1;
	F77_CALL(daxpy)(&Nlambdas,&alpha,ErrorMean,&inci,errorMeanCpy,&incj); // y = ax + y
	int ilambda_ms;	
	ilambda_ms = F77_CALL(idamax)(&Nlambdas, errorMeanCpy, &inci);//index of the max(errs_mean)<--min(mean)
	index = ilambda_ms - 1;


	//return index of lambda_factors
	if(verbose>1) Rprintf("\t\tExit Function: cv_support. optimal lambda index: %d.\n\n", ilambda_ms);

	R_Free(NOISE);
	R_Free(ImB);
	
	
	R_Free(Q);
	R_Free(BL);
	R_Free(fL);
	R_Free(mueL);
	R_Free(BC);
	R_Free(fC);
	
	R_Free(Errs);
	R_Free(Sigmas2);
	R_Free(ErrorMean);
	R_Free(ErrorCV);
	R_Free(errorMeanCpy);
	
	R_Free(Ylearn);
	R_Free(Xlearn);
	R_Free(Ytest);
	R_Free(Xtest);	


	R_Free(Ylearn_centered);
	R_Free(Xlearn_centered);	
	R_Free(meanY);
	R_Free(meanX);
	R_Free(SL);


	R_Free(IBinvPath);
	R_Free(IBinvPathZero);
	
	R_Free(lambda_Max);
	
	return ilambda_ms;
	
	
	
} //end of cv_gene_nets_support		

			
void mainSML(double *Y, double *X, int *m, int *n, int *Missing,double*B, double *f,double*stat,int*VB)
{

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
	Strue 		= (double* ) R_Calloc(MM, double);
	F77_CALL(dcopy)(&MM,B,&inci,Strue,&incj);
	stat[1] 	= 0;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index = j*M  +i;
			if(i!=j && B[index]!=0)
			{	stat[1] = stat[1] + 1;} //stat[1] total positive
		}
	}
	double alpha = 1;	
	F77_CALL(dcopy)(&M,&alpha,&inc0,f,&inci);

	alpha = 0;
	//F77_CALL(dscal)(&MM,&alpha,B,&inci);
	F77_CALL(dcopy)(&MM,&alpha,&inc0,B,&inci);

	for(i=0;i<MN;i++)
	{
		if(Missing[i] == 1) Y[i] = 0;
	}
	
	int maxiter 	= 500;
	int Kcv 		= 5;
	int L_lambda 	= 20; // number of lambdas in lambda_factors	stop at 0.001
	double *lambda_factors;
	lambda_factors 	= (double* ) R_Calloc(L_lambda, double);
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
	rho_factors 			= (double* ) R_Calloc(L, double);
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
	double *W;  //weight on diagonal?
	W = (double* ) R_Calloc(MM, double);

	double *QIBinv;
	QIBinv = (double* ) R_Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1
	
	ilambda_cv_ms = cv_gene_nets_support_Rdg(Y, X, Kcv,lambda_factors, rho_factors, 
			maxiter, M, N,Nlambdas, Nrho,verbose,W, &sigma2);
	
	
	if(verbose==0) Rprintf("Step 1: CV support; return number of lambda needed: %d\n", ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) R_Calloc(M, double);
	meanX = (double* ) R_Calloc(M, double);
	Ycopy = (double* ) R_Calloc(MN, double);
	Xcopy = (double* ) R_Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) R_Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda; 
	double lambda_max;
	lambda_max = lambdaMax(Ycopy,Xcopy,W,M, N);

	if(verbose==0) Rprintf("Step 4: lasso selection path.\n");

	for(ilambda = 0;ilambda<ilambda_cv_ms;ilambda++)
	{
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, ilambda_cv_ms,lambda_factor);

		lambda = Weighted_LassoSf(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max); 	// mueL not calculated
	
if(verbose>0) Rprintf("\tlambda: %f\n", lambda);						

		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];

	}//ilambda; selection path

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
	R_Free(Strue);
	R_Free(meanY);
	R_Free(meanX);
	R_Free(lambda_factors);
	R_Free(rho_factors);
	R_Free(Ycopy);
	R_Free(Xcopy);

	R_Free(W);
	R_Free(QIBinv);
	R_Free(Q);
}


double lambdaMax_adaEN(double *Y,double *X,double * Wori,int M, int N,double alpha_factor) 	//------adaEN 	
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
	int MN = M*N;
	int MM = M*M;
	//------adaEN 	
	double *W;
	W 				= (double* ) R_Calloc(MM, double);
	F77_CALL(dcopy)(&MM,Wori,&inci,W,&incj);
	F77_CALL(dscal)(&MM,&alpha_factor,W,&inci); 
	
	for(i=0;i<M;i++)
	{
		readPtr1 	= &X[i]; //ith row
		readPtr2 	= &Y[i];

		dxx[i] = F77_CALL(ddot)(&N,readPtr1,&lda,readPtr1,&M);
		rxy[i] 		= F77_CALL(ddot)(&N,readPtr1,&lda,readPtr2,&M);
		DxxRxy[i] 	= rxy[i]/dxx[i];		
	}
	

	double * XDxxRxy;

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

	R_Free(XDxxRxy);
	R_Free(YYXDR);
	//------adaEN 
	R_Free(W);

	return lambda_max;	
}


//
double Weighted_LassoSf_adaEN(double * Wori, double *B, double *f, double *Ycopy,double *Xcopy, 				//------adaEN 	
		double *Q, double lambda_factor, double lambda_factor_prev, double sigma2, int max_iter,
		int M, int N, int verbose,double *QIBinv,double lambda_max,			//double * mue,
		double alpha_factor)
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

	double *s, *S,*Wcopy, *W;
	S = (double* ) R_Calloc(MM, double);
	s = (double* ) R_Calloc(M, double);
	W 						= (double* ) R_Calloc(MM, double);
	F77_CALL(dcopy)(&MM,Wori,&inci,W,&incj);
	F77_CALL(dscal)(&MM,&alpha_factor,W,&inci); 
//------adaEN
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

			index = j*M  +i;
//------adaEN  
		
			if(fabs(Q[index] -(1-alpha_factor)*lambda*B[index])>= Wcopy[index] && i!= j)	
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
		//F77_CALL(dcopy)(&N,readPtr,&M,xi,&inci);
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
	//quadratic function
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
				if(verbose>6) Rprintf("\t\t\t\t\t updating Node %d \n",i);
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
//------adaEN  
r_ij = r_ij + (1 -alpha_factor)*lambda;

							beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
							
							if (fabs(m_ij)<1e-10) //go to the linear equation 
							{
								//
								if(verbose>7) Rprintf("\t\t\t\t\t\t\t Node %d \t interact with Node %d.\tLinear equation\n",i,j);
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
								if(verbose>7) Rprintf("\t\t\t\t\t\t\t Node %d \t interact with Node %d.\tQuadratic equation\n",i,j);

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
			}else//s[i]  no un-zero weight in this Node
			{
				readPtr = &B[i];
				F77_CALL(dcopy)(&M,&toyZero,&inc0,readPtr,&ldM);
				f[i] = f0[i];
			} // s[i]
		}//i= 1:M
		
		F77_CALL(dcopy)(&MM,B,&inci,BfNew,&incj);
		F1ptr = &BfNew[MM];
		F77_CALL(dcopy)(&M,f,&inci,F1ptr,&incj);
		
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
	
	//------adaEN 
	R_Free(W);
	return lambda;
}//weighted_LassoSf

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
	
//------adaEN
	double *s, *S,*Wcopy,*W; //global copy of W
	S = (double* ) R_Calloc(MM, double);
	s = (double* ) R_Calloc(M, double);
	W = (double* ) R_Calloc(MM, double);
	F77_CALL(dcopy)(&MM,Wori,&inci,W,&incj);
	F77_CALL(dscal)(&MM,&alpha_factor,W,&inci); //MAKE sure alpha_factor is bcasted

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

			index = j*M  +i;
//------adaEN		
// version 1_3: Qij			
			if(fabs(Q[index] -(1-alpha_factor)*lambda*B[index])>= Wcopy[index] && i!= j)	
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


		//dot product
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
	double *eiB;
	eiB = (double* ) R_Calloc(M, double);
	double *BiT;
	BiT = (double* ) R_Calloc(M, double);
	//quadratic function
	double d_ij, theta_ijp,k_ijp,q_ijp,Bijpp, Bijpm; //case (14)
	double q_ijm, theta_ijm, Bijmm, Bijmp,Lss,candsBij,LssCands;
	
	//converge of Node i
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
				if(verbose>6) Rprintf("\t\t\t\t\t updating Node %d \n",i);
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
//------adaEN 
r_ij = r_ij + (1 -alpha_factor)*lambda;

							beta_ij = F77_CALL(ddot)(&N, y_j, &inci,a_iT, &incj);
							
							if (fabs(m_ij)<1e-10) //go to the linear equation 
							{
								//
								if(verbose>7) Rprintf("\t\t\t\t\t\t\t Node %d \t interact with Node %d.\tLinear equation\n",i,j);
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
								if(verbose>7) Rprintf("\t\t\t\t\t\t\t Node %d \t interact with Node %d.\tQuadratic equation\n",i,j);
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
						}//if(j!= i)
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
			}else//s[i]  no un-zero weight in this Node
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
		// BLOCK COORDINATE ASCEND: Update IBinv
		UpdateIBinv(QIBinv, B,M);
		
		if(verbose>5) Rprintf("\t\t\t\t\t\tdelta_BF: %f\n",delta_BF);
		if(delta_BF<1e-3)		//break out
		{
			break;
		}

	}//while
	
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
	//quadratic function
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

		// inner loop
		for(i=0;i<M;i++)
		{
			if(s[i] >0)
			{
				//
				if(verbose>6) Rprintf("\t\t updating Node %d \n",i);
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
							if(verbose>7) Rprintf("\t\t\t Node %d \t interact with Node %d.\tLinear equation\n",i,j);
							B[j*M+i] = beta_ij/r_ij;

						}else //m_ij ~=0 go to the quadratic equation
						{
							//
							if(verbose>7) Rprintf("\t\t\t Node %d \t interact with Node %d.\tQuadratic equation\n",i,j);
							//					

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
	transa = 'N';
	alpha = 1;
	//beta = 1;
	ldk = M;
	F77_CALL(dgemv)(&transa, &M, &ldk,&alpha, IBinv, &ldM, meanY, &inci, &beta,mue, &incj FCONE);
	
	
	//----------------------------------------------------------------------------------END OF LASSO_ZERO

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
	//R_Free(zzi);
	R_Free(a_iT);
	
	R_Free(eiB);
	R_Free(BiT);
	R_Free(BfOld);
	R_Free(BfNew);
	R_Free(BfChange);
	
	//R_Free(ipiv);
	//------adaEN
	R_Free(W);

	return lambda;
}//weighted_LassoSf


int cv_gene_nets_support_adaEN(double *Y, double *X, int Kcv,double *lambda_factors, double *rho_factors, 
			int maxiter, int M, int N,int Nlambdas, int Nrho,int verbose,double *W, 			//double sigma2learnt,
			double *sigmaLasso,
			int i_alpha, double alpha_factor,double * ErrorEN,double *sigma2learnt_EN,
						double *ErrorEN_min, double* steEN_min)	
{	
	int MM = M*M;
	int MN = M*N;
	// to save for each lambda
	double *Q, *BL, *fL,*mueL; 
	int i,j,index;	
	Q = (double* ) R_Calloc(MM, double); //ptr Q
	BL = (double* ) R_Calloc(MM, double); //ptr BL
	fL = (double* ) R_Calloc(M, double); //ptr fL
	mueL = (double* ) R_Calloc(M, double);
	//
	double *BC, *fC;
	BC = (double* ) R_Calloc(MM, double);
	fC = (double* ) R_Calloc(M, double);
	
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
	
	Errs = (double* ) R_Calloc(i, double);
	Sigmas2 = (double* ) R_Calloc(Nlambdas, double);
	ErrorMean= (double* ) R_Calloc(Nlambdas, double);
	i = Nlambdas*Kcv;
	ErrorCV = (double* ) R_Calloc(i, double);
	
	//parameters inside the loop
	double *Ylearn, *Xlearn, *Ytest,*Xtest;
	int NlearnM = Nlearn*M;
	int NtestM = Ntest*M;
	Ylearn = (double* ) R_Calloc(NlearnM, double);
	Xlearn = (double* ) R_Calloc(NlearnM, double);
	Ytest = (double* ) R_Calloc(NtestM, double);
	Xtest = (double* ) R_Calloc(NtestM, double);

	//centerYX function
	double *Ylearn_centered, *Xlearn_centered, *meanY, *meanX;
	Ylearn_centered = (double* ) R_Calloc(NlearnM, double);
	Xlearn_centered = (double* ) R_Calloc(NlearnM, double);
	meanY = (double* ) R_Calloc(M, double);
	meanX = (double* ) R_Calloc(M, double);
	

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

	double *SL;				//zero or not-zero
	SL = (double* ) R_Calloc(MM, double);

	//convergence
	//------adaEN 
	if(verbose>1) Rprintf("\t\ti_alpha %d; \t Enter Function: cv_support. Nlambdas: %d; \t %d-fold cross validation.\n", i_alpha, Nlambdas,Kcv);
	if(verbose>4) Rprintf("\n\t\t\t\t\tEnter Function: ridge_cvf. %d-fold cross validation.\n", Kcv);

	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	int lda,ldb,ldc,ldk;
	double *Xptr, *XsubPtr; //submatrix
	double alpha, beta;

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
	
	double *testNoise,*ImB,*NOISE;
	NOISE =(double* ) R_Calloc(MN, double);
	testNoise = &NOISE[0];
	ImB = (double* ) R_Calloc(MM, double);
	char transa = 'N'; 
	char transb = 'N';
	lda = M;
	ldb = M;
	ldc = M;
	ldk = M;
	double testNorm;
	//result of CV 
	double rho_factor; // no pointers needed

	
//---------------------adaEN  
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
			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);

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
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc FCONE FCONE);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				XsubPtr = &testNoise[i];
				//row i of noise
				alpha = -fR[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &ldk,XsubPtr, &M);
			}//row i = 1:M

			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				Xptr = &testNoise[i*M];
				F77_CALL(daxpy)(&M, &alpha,mueR, &inci,Xptr, &incj);
			}

			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);

			Errs[cv*Nrho+irho] = testNorm;
	
			
			
			err_mean = err_mean + testNorm;
		}//cv = 0: Kcv
		
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

	if(verbose>4) Rprintf("\t\t\t\t\tExit RidgeCV. sigma2R: %f\t",sigma2R);	

	if(irho == Nrho || irho == 1) //nonstop in while loop
	{	
		sigma2R = err_mean/(MN -1);
	}else
	{
		sigma2R = err_mean_prev/(MN -1);
	}
	rho_factor = rho_factors[irho-1]*N/(N-Ntest);
	if(verbose>4) Rprintf("sigma2learnt: %f\n",sigma2R);	
	if(verbose==0) Rprintf("Step 1: ridge CV; find rho : %f\n", rho_factor);
	sigma2learnt = constrained_ridge_cff(Y, X, rho_factor, M, N,BR,fR,mueR,verbose);	
	if(verbose==0) Rprintf("Step 2: ridge; calculate weights.\n");
	for(i=0;i<MM;i++) W[i] = 1/fabs(BL[i]+ 1e-10);
	
	sigmaLasso[0] = sigma2R;
	sigma2learnt_EN[0] = sigma2learnt;
	}else	//i_alpha ==0
	{
		sigma2R 		= sigmaLasso[0];
		sigma2learnt 	= sigma2learnt_EN[0];
	}
//-----------------------------------------------adaEN  
	double *IBinv,*IBinvZero,*lambda_Max;
	double *IBinvPath,*IBinvPathZero;
	irho  = MM*Kcv;
	IBinvPath = (double* ) R_Calloc(irho, double);
	IBinvPathZero = (double* ) R_Calloc(irho, double);
	lambda_Max = (double* ) R_Calloc(Kcv, double);
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
		err_mean = 0;
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>3) Rprintf("\t\t %d/Kcv cross validation.\n", cv);
			// test, learn
			// start and end point
			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);
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
			
			//Ylearn_centered
			F77_CALL(dcopy)(&NlearnM,Xlearn,&inci,Xlearn_centered,&incj);
			F77_CALL(dcopy)(&NlearnM,Ylearn,&inci,Ylearn_centered,&incj);

			centerYX(Ylearn_centered,Xlearn_centered,meanY, meanX,M, Nlearn);
			//first ilambda
		
			if(ilambda == 0)
			{
				alpha = 1; 
				F77_CALL(dcopy)(&M,&alpha,&inc0,fL,&incj); // call dcopy(n, x, inci, y, incy)
				alpha = 0;
				F77_CALL(dcopy)(&MM,&alpha,&inc0,BL,&incj); 
				QlambdaStart(Ylearn_centered,Xlearn_centered, Q, sigma2learnt,M, Nlearn);
				lambda_Max[cv]	= lambdaMax_adaEN(Ylearn_centered,Xlearn_centered,W,M, Nlearn,alpha_factor); 	//------adaEN 	
			}//ilambda ==0
			lambda_max_cv = lambda_Max[cv];
			lambda_factor = lambda_factors[ilambda];
							
			IBinv = &IBinvPath[cv*MM];
			IBinvZero= &IBinvPathZero[cv*MM];
			lambda = Weighted_LassoSf_MLf_adaEN(W, BL, fL, Ylearn,Xlearn, Q, lambda_factor, 
							lambda_factor_prev, sigma2learnt, maxiter,M, Nlearn, verbose,
							BC, fC, mueL,IBinv,IBinvZero,lambda_max_cv,
							alpha_factor); 								//------adaEN  

							
			if(verbose>3) Rprintf("\t\t\t step 1 SML lasso regression, lambda: %f.\n",lambda);
			QlambdaMiddleCenter(Ylearn_centered,Xlearn_centered, Q,BL,fL,sigma2learnt,M, Nlearn,IBinv);
			
			// constrained_MLf
			if(verbose>3) Rprintf("\t\t\t step 2 SML ZeroRegression.\n");
			
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
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc FCONE FCONE);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				XsubPtr = &NOISE[NtestM*cv+i];
				alpha = -fC[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &M,XsubPtr, &ldk);
			}//row i = 1:M
			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				Xptr = &NOISE[NtestM*cv+i*M];
				F77_CALL(daxpy)(&M, &alpha,mueL, &inci,Xptr, &incj);
			}
			
			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);
			index = cv*Nlambdas+ilambda;
			
			Errs[index] = testNorm;

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
		
		lambda_factor_prev = lambda_factor;
		
		if(verbose>2) Rprintf("\t\t\t %d/Nlambdas. %d fold cv; \t Err_Mean: %f; std:%f; \t sigma2learnt:%f.\n", ilambda,Kcv,err_mean,Sigmas2[ilambda],sigma2learnt);
		ilambda = ilambda + 1; 
	}
	double *errorMeanCpy;
	errorMeanCpy = (double* ) R_Calloc(Nlambdas, double);
	//int inc0  = 0;
	double minimumErr = -1e5 - MN;
	F77_CALL(dcopy)(&Nlambdas,&minimumErr,&inc0,errorMeanCpy,&inci);
	alpha = 1;
	F77_CALL(daxpy)(&Nlambdas,&alpha,ErrorMean,&inci,errorMeanCpy,&incj); // y = ax + y
	int ilambda_ms;	
	ilambda_ms = F77_CALL(idamax)(&Nlambdas, errorMeanCpy, &inci);//index of the max(errs_mean)<--min(mean)
	index = ilambda_ms - 1;
	minimumErr = ErrorMean[index] + Sigmas2[index]; //actually max

	double distance;
	int tempIlambda_ms = index;
	if(index ==0)
	{
		tempIlambda_ms=1;
	}else
	{
		double minDist = fabs(ErrorMean[index -1] - minimumErr);
		for(i=index-1;i>0;i--) 
		{
			distance = fabs(ErrorMean[i] - minimumErr);
			if(distance <minDist)
			{
				minDist = distance;
				tempIlambda_ms = i + 1;
			}
		}
	}
	ilambda_ms = tempIlambda_ms;
	
	index = ilambda_ms - 1;
	ErrorEN[i_alpha] = ErrorMean[index];

//------adaEN	
	ErrorEN_min[i_alpha] = ErrorMean[index];
	steEN_min[i_alpha] = Sigmas2[index];

	if(verbose>1) Rprintf("\t\tExit Function: cv_support. optimal lambda index: %d.\n\n", ilambda_ms);

	R_Free(NOISE);
	R_Free(ImB);
	
	
	R_Free(Q);
	R_Free(BL);
	R_Free(fL);
	R_Free(mueL);
	R_Free(BC);
	R_Free(fC);
	
	R_Free(Errs);
	R_Free(Sigmas2);
	R_Free(ErrorMean);
	R_Free(ErrorCV);
	R_Free(errorMeanCpy);
	
	R_Free(Ylearn);
	R_Free(Xlearn);
	R_Free(Ytest);
	R_Free(Xtest);	

	R_Free(Ylearn_centered);
	R_Free(Xlearn_centered);	
	R_Free(meanY);
	R_Free(meanX);
	R_Free(SL);
	R_Free(IBinvPath);
	R_Free(IBinvPathZero);
	
	R_Free(lambda_Max);
	
	return ilambda_ms;

	
} //end of cv_gene_nets_support		


void mainSML_adaEN(double *Y, double *X, int *m, int *n, int *Missing,double*B, double *f,double*stat,
                   double*out_alpha,
                   double*out_lambda,
                   int*VB)
{
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
	Strue 		= (double* ) R_Calloc(MM, double);
	F77_CALL(dcopy)(&MM,B,&inci,Strue,&incj);
	stat[1] 	= 0;
	for(i=0;i<M;i++)
	{
		for(j=0;j<M;j++)
		{
			index = j*M  +i;
			if(i!=j && B[index]!=0)
			{	stat[1] = stat[1] + 1;} //stat[1] total positive
		}
	}
	double alpha = 1;	
	F77_CALL(dcopy)(&M,&alpha,&inc0,f,&inci);

	alpha = 0;
	F77_CALL(dcopy)(&MM,&alpha,&inc0,B,&inci);

	for(i=0;i<MN;i++)
	{
		if(Missing[i] == 1) Y[i] = 0;
	}
	
	//call cv_gene_nets_support ------------------------SYSTEM PARAMETERS
	int maxiter 	= 500;
	int Kcv 		= 5;
	int L_lambda 	= 20; // number of lambdas in lambda_factors	stop at 0.001
	double *lambda_factors;
	lambda_factors 	= (double* ) R_Calloc(L_lambda, double);
	double step 	= -0.2;
	for(i=0;i<L_lambda;i++) 
	{
		lambda_factors[i] 	= pow(10.0,step);
		step 				= step - 0.2;
	}
	
	step 					= -6;
	double * rho_factors;
	int L 					= 31; // number of rho_factors
	rho_factors 			= (double* ) R_Calloc(L, double);
	for(i=0;i<L;i++) 
	{
		rho_factors[i] 		= pow(10.0,step);
		step 				= step + 0.2;
	}
	//------adaEN
	double *alpha_factors, *ErrorEN,*lambdaEN, sigma2learnt;
	int L_alpha  		= 19; 
	alpha_factors 		= (double* ) R_Calloc(L_alpha, double); //variable in all ranks
	ErrorEN 			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	lambdaEN 			= (double* ) R_Calloc(L_alpha, double);

	double *ErrorEN_min, *steEN_min;
	ErrorEN_min			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	steEN_min 			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	
	step 					= 0.05;
	for(i=0;i<L_alpha;i++) 
	{
		alpha_factors[i]	= 0.95 - step*i;
	}

//------adaEN
	
	int Nlambdas,Nrho;
	Nlambdas 				= L_lambda;
	Nrho 					= L;
	double sigma2; //return value;

	double *W;  //weight on diagonal?
	W = (double* ) R_Calloc(MM, double);

	double *QIBinv;
	QIBinv = (double* ) R_Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1

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

	ilambda_cv_ms 			= lambdaEN[ind_minEN];
	alpha_factor  			= alpha_factors[ind_minEN];
	if(verbose==0) Rprintf("\tAdaptive_EN %d-fold CV, alpha: %f.\n", Kcv,alpha_factor);

	
	if(verbose==0) Rprintf("Step 3: CV support; alpha: %f, number of lambda needed: %d\n", alpha_factor,ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) R_Calloc(M, double);
	meanX = (double* ) R_Calloc(M, double);
	Ycopy = (double* ) R_Calloc(MN, double);
	Xcopy = (double* ) R_Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) R_Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda; 
	double lambda_max;
	
	lambda_max = lambdaMax_adaEN(Ycopy,Xcopy,W,M, N,alpha_factor);
	
	if(verbose==0) Rprintf("Step 4: lasso selection path.\n");
	
	for(ilambda = 0;ilambda<ilambda_cv_ms;ilambda++)
	{
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, ilambda_cv_ms,lambda_factor);
		// call Weighted_LassoSf		Missing,
		lambda = Weighted_LassoSf_adaEN(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max,
					alpha_factor); 	// mueL not calculated
			
if(verbose>0) Rprintf("\tlambda: %f\n", lambda);				

		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];
	}//ilambda: selection path
	
	out_alpha[0] = alpha_factor;
	out_lambda[0] = lambda_factor;
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
	R_Free(Strue);
	R_Free(meanY);
	R_Free(meanX);
	R_Free(lambda_factors);
	R_Free(rho_factors);
	R_Free(Ycopy);
	R_Free(Xcopy);
	R_Free(W);
	R_Free(QIBinv);
	R_Free(Q);
	
	//------adaEN	
	R_Free(alpha_factors);
	R_Free(ErrorEN);
	R_Free(lambdaEN);
	
	R_Free(ErrorEN_min);
	R_Free(steEN_min);

}


int cv_gene_nets_support_adaENcv(double *Y, double *X, int Kcv,double *lambda_factors, double *rho_factors, 
			int maxiter, int M, int N,int Nlambdas, int Nrho,int verbose,double *W, 			//double sigma2learnt,
			double *sigmaLasso,
			int i_alpha, double alpha_factor,double * ErrorEN,double *sigma2learnt_EN,//------adaEN
			double *mseStd, 		//
			double *ErrorEN_min, double* steEN_min)		
{
	int MM = M*M;
	int MN = M*N;
	// to save for each lambda
	double *Q, *BL, *fL,*mueL; 
	int i,j,index;	
	Q = (double* ) R_Calloc(MM, double); //ptr Q
	BL = (double* ) R_Calloc(MM, double); //ptr BL
	fL = (double* ) R_Calloc(M, double); //ptr fL
	mueL = (double* ) R_Calloc(M, double);
	//
	double *BC, *fC;
	BC = (double* ) R_Calloc(MM, double);
	fC = (double* ) R_Calloc(M, double);
	
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
	
	Errs = (double* ) R_Calloc(i, double);
	Sigmas2 = (double* ) R_Calloc(Nlambdas, double);
	ErrorMean= (double* ) R_Calloc(Nlambdas, double);
	i = Nlambdas*Kcv;
	ErrorCV = (double* ) R_Calloc(i, double);
	
	//parameters inside the loop
	double *Ylearn, *Xlearn, *Ytest,*Xtest;
	int NlearnM = Nlearn*M;
	int NtestM = Ntest*M;
	Ylearn = (double* ) R_Calloc(NlearnM, double);
	Xlearn = (double* ) R_Calloc(NlearnM, double);
	Ytest = (double* ) R_Calloc(NtestM, double);
	Xtest = (double* ) R_Calloc(NtestM, double);

	//centerYX function
	double *Ylearn_centered, *Xlearn_centered, *meanY, *meanX;
	Ylearn_centered = (double* ) R_Calloc(NlearnM, double);
	Xlearn_centered = (double* ) R_Calloc(NlearnM, double);
	meanY = (double* ) R_Calloc(M, double);
	meanX = (double* ) R_Calloc(M, double);
	

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

	double *SL;				//zero or not-zero
	SL = (double* ) R_Calloc(MM, double);

	//------adaEN
	if(verbose>1) Rprintf("\t\ti_alpha %d; \t Enter Function: cv_support. Nlambdas: %d; \t %d-fold cross validation.\n", i_alpha, Nlambdas,Kcv);
	if(verbose>4) Rprintf("\n\t\t\t\t\tEnter Function: ridge_cvf. %d-fold cross validation.\n", Kcv);
	int inci = 1;
	int incj = 1;
	int inc0 = 0;
	int lda,ldb,ldc,ldk;
	double *Xptr, *XsubPtr; //submatrix
	double alpha, beta;

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

	double *testNoise,*ImB,*NOISE;
	NOISE =(double* ) R_Calloc(MN, double);
	testNoise = &NOISE[0];
	ImB = (double* ) R_Calloc(MM, double);
	char transa = 'N'; 
	char transb = 'N';
	lda = M;
	ldb = M;
	ldc = M;
	ldk = M;
	double testNorm;
	//result of CV 
	double rho_factor; // no pointers needed

	
//---------------------adaEN 
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
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc FCONE FCONE);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				XsubPtr = &testNoise[i];
				//row i of noise
				alpha = -fR[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &ldk,XsubPtr, &M);
			}//row i = 1:M

			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				Xptr = &testNoise[i*M];
				F77_CALL(daxpy)(&M, &alpha,mueR, &inci,Xptr, &incj);
			}

			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);

			Errs[cv*Nrho+irho] = testNorm;
	
			
			
			err_mean = err_mean + testNorm;

		}//cv = 0: Kcv
		
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
	
	if(verbose>4) Rprintf("\t\t\t\t\tExit RidgeCV. sigma2R: %f\t",sigma2R);	

	if(irho == Nrho || irho == 1) //nonstop in while loop
	{	
		sigma2R = err_mean/(MN -1);
	}else
	{
		sigma2R = err_mean_prev/(MN -1);
	}

	rho_factor = rho_factors[irho-1]*N/(N-Ntest);
	if(verbose>4) Rprintf("sigma2learnt: %f\n",sigma2R);	

	if(verbose==0) Rprintf("Step 1: ridge CV; find rho : %f\n", rho_factor);
	sigma2learnt = constrained_ridge_cff(Y, X, rho_factor, M, N,BR,fR,mueR,verbose);	
	if(verbose==0) Rprintf("Step 2: ridge; calculate weights.\n");

	for(i=0;i<MM;i++) W[i] = 1/fabs(BL[i]+ 1e-10);

	sigmaLasso[0] = sigma2R;
	sigma2learnt_EN[0] = sigma2learnt;
	}else	//i_alpha ==0
	{
		sigma2R 		= sigmaLasso[0];
		sigma2learnt 	= sigma2learnt_EN[0];
	}
//-----------------------------------------------adaEN

	double *IBinv,*IBinvZero,*lambda_Max;
	double *IBinvPath,*IBinvPathZero;
	irho  = MM*Kcv;
	IBinvPath = (double* ) R_Calloc(irho, double);
	IBinvPathZero = (double* ) R_Calloc(irho, double);
	lambda_Max = (double* ) R_Calloc(Kcv, double);
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
		err_mean = 0;
		for(cv=0;cv<Kcv;cv++)
		{
			if(verbose>3) Rprintf("\t\t %d/Kcv cross validation.\n", cv);

			testStart = Ntest*cv + 1;
			testEnd = Ntest*(cv+1);
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
			
			F77_CALL(dcopy)(&NlearnM,Xlearn,&inci,Xlearn_centered,&incj);
			F77_CALL(dcopy)(&NlearnM,Ylearn,&inci,Ylearn_centered,&incj);

			centerYX(Ylearn_centered,Xlearn_centered,meanY, meanX,M, Nlearn);
			
			if(ilambda == 0)
			{
				alpha = 1; 
				F77_CALL(dcopy)(&M,&alpha,&inc0,fL,&incj); // call dcopy(n, x, inci, y, incy)
				alpha = 0;
				F77_CALL(dcopy)(&MM,&alpha,&inc0,BL,&incj); 

				QlambdaStart(Ylearn_centered,Xlearn_centered, Q, sigma2learnt,M, Nlearn);
				lambda_Max[cv]	= lambdaMax_adaEN(Ylearn_centered,Xlearn_centered,W,M, Nlearn,alpha_factor); 	//------adaEN	
			}//ilambda ==0
			lambda_max_cv = lambda_Max[cv];

			lambda_factor = lambda_factors[ilambda];

			IBinv = &IBinvPath[cv*MM];
			IBinvZero= &IBinvPathZero[cv*MM];
			lambda = Weighted_LassoSf_MLf_adaEN(W, BL, fL, Ylearn,Xlearn, Q, lambda_factor, 
							lambda_factor_prev, sigma2learnt, maxiter,M, Nlearn, verbose,
							BC, fC, mueL,IBinv,IBinvZero,lambda_max_cv,
							alpha_factor); 								//------adaEN

			if(verbose>3) Rprintf("\t\t\t step 1 SML lasso regression, lambda: %f.\n",lambda);

			QlambdaMiddleCenter(Ylearn_centered,Xlearn_centered, Q,BL,fL,sigma2learnt,M, Nlearn,IBinv);
			
			// constrained_MLf
			if(verbose>3) Rprintf("\t\t\t step 2 SML ZeroRegression.\n");
			
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
			F77_CALL(dgemm)(&transa, &transb,&M, &Ntest, &ldk,&alpha, ImB, &lda, Ytest, &ldb, &beta, testNoise, &ldc FCONE FCONE);
			for(i=0;i<M;i++)
			{
				// row i of X
				Xptr = &Xtest[i];
				XsubPtr = &NOISE[NtestM*cv+i];
				alpha = -fC[i];
				F77_CALL(daxpy)(&Ntest, &alpha,Xptr, &M,XsubPtr, &ldk);
			}//row i = 1:M
			alpha = -1;
			for(i=0;i<Ntest;i++)
			{
				Xptr = &NOISE[NtestM*cv+i*M];
				F77_CALL(daxpy)(&M, &alpha,mueL, &inci,Xptr, &incj);
			}

			testNorm = F77_CALL(ddot)(&NtestM, testNoise, &inci,testNoise, &incj);
			index = cv*Nlambdas+ilambda;
			
			Errs[index] = testNorm;

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
		
		lambda_factor_prev = lambda_factor;
		
		if(verbose>2) Rprintf("\t\t\t %d/Nlambdas. %d fold cv; \t Err_Mean: %f; std:%f; \t sigma2learnt:%f.\n", ilambda,Kcv,err_mean,Sigmas2[ilambda],sigma2learnt);
		ilambda = ilambda + 1; 
	}//while ilambda

	double *errorMeanCpy;
	errorMeanCpy = (double* ) R_Calloc(Nlambdas, double);
	//int inc0  = 0;
	double minimumErr = -1e5 - MN;
	F77_CALL(dcopy)(&Nlambdas,&minimumErr,&inc0,errorMeanCpy,&inci);
	alpha = 1;
	F77_CALL(daxpy)(&Nlambdas,&alpha,ErrorMean,&inci,errorMeanCpy,&incj); // y = ax + y
	int ilambda_ms;	
	ilambda_ms = F77_CALL(idamax)(&Nlambdas, errorMeanCpy, &inci);//index of the max(errs_mean)<--min(mean)
	index = ilambda_ms - 1;
	minimumErr = ErrorMean[index] + Sigmas2[index]; //actually max

//-------------------------------------------------------------adaEN	

	double *readPtr1;

	readPtr1 = &mseStd[0];
	F77_CALL(dcopy)(&Nlambdas,ErrorMean,&inci,readPtr1,&incj);
	readPtr1 = &mseStd[Nlambdas];
	F77_CALL(dcopy)(&Nlambdas,Sigmas2,&inci,readPtr1,&incj);

	double distance;
	int tempIlambda_ms = index;
	if(index ==0)
	{	
		tempIlambda_ms=1;
	}else
	{	
		double minDist = fabs(ErrorMean[index -1] - minimumErr);
		for(i=index-1;i>0;i--) 
		{
			distance = fabs(ErrorMean[i] - minimumErr);
			if(distance <minDist)
			{
				minDist = distance;
				tempIlambda_ms = i + 1;
			}
		}
	}
	ilambda_ms = tempIlambda_ms;

	index = ilambda_ms - 1;
	ErrorEN[i_alpha] = ErrorMean[index];

//------adaEN		
	ErrorEN_min[i_alpha] = ErrorMean[index];
	steEN_min[i_alpha] = Sigmas2[index];

	if(verbose>1) Rprintf("\t\tExit Function: cv_support. optimal lambda index: %d.\n\n", ilambda_ms);


	R_Free(NOISE);
	R_Free(ImB);
	
	
	R_Free(Q);
	R_Free(BL);
	R_Free(fL);
	R_Free(mueL);
	R_Free(BC);
	R_Free(fC);
	
	R_Free(Errs);
	R_Free(Sigmas2);
	R_Free(ErrorMean);
	R_Free(ErrorCV);
	R_Free(errorMeanCpy);
	
	R_Free(Ylearn);
	R_Free(Xlearn);
	R_Free(Ytest);
	R_Free(Xtest);	

	R_Free(Ylearn_centered);
	R_Free(Xlearn_centered);	
	R_Free(meanY);
	R_Free(meanX);
	R_Free(SL);


	R_Free(IBinvPath);
	R_Free(IBinvPathZero);
	
	R_Free(lambda_Max);
	
	return ilambda_ms;

	
} //end of cv_gene_nets_support		



void mainSML_adaENcv(double *Y, double *X, int *m, int *n, int *Missing, double*B, double *f,double*stat,
			double*alpha_factors,int *nAlpha, 	// must be scalar
			double *lambda_factors, int *nLambda,
			double*mse,double*mseSte,
			double *mseStd, int*kFold, 
			double*out_alpha,
			double*out_lambda,
			int*VB) 						// mseStd: nLmabda x 2 matrix, keep mse + std 
{

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
	Strue 					= (double* ) R_Calloc(MM, double);
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
	int Kcv 				= kFold[0];

	double step 	= -0.2;
	//rho factor
	step 					= -6;
	double * rho_factors;
	int L 					= 31; // number of rho_factors
	rho_factors 			= (double* ) R_Calloc(L, double);
	for(i=0;i<L;i++) 
	{
		rho_factors[i] 		= pow(10.0,step);
		step 				= step + 0.2;
	}
	//------adaEN 	
	double  *ErrorEN,*lambdaEN, sigma2learnt;	//*alpha_factors,
	int L_alpha  		= nAlpha[0]; 
	ErrorEN 			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	lambdaEN 			= (double* ) R_Calloc(L_alpha, double);

	double *ErrorEN_min, *steEN_min;
	ErrorEN_min			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	steEN_min 			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	
	step 					= 0.05;

	//---------------------------------------------END  SYSTEM PARAMETERS
	int Nlambdas,Nrho;
	//----------------------------------------------------------R-package
	Nlambdas 				= nLambda[0]; 			

	//----------------------------------------------------------R-package	
	Nrho 					= L;
	double sigma2; //return value;

	double *W;  //weight on diagonal?
	W = (double* ) R_Calloc(MM, double);
	double *QIBinv;
	QIBinv = (double* ) R_Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1

	int i_alpha;
	double alpha_factor;
	
	//--- ouput to R
	double*readPtr1, *readPtr2;
	for(i_alpha=0;i_alpha<L_alpha;i_alpha++)
	{	
	  alpha_factor 		= alpha_factors[i_alpha];		
		
		ilambda_cv_ms = cv_gene_nets_support_adaENcv(Y, X, Kcv,lambda_factors, rho_factors, 
			maxiter, M, N,Nlambdas, Nrho,verbose,W, &sigma2,
			i_alpha,alpha_factor,ErrorEN, &sigma2learnt,mseStd,
			ErrorEN_min,steEN_min);	 //version V1_1delta.c			
			
		lambdaEN[i_alpha] 	= ilambda_cv_ms;
		
		
		//--- output to R: mse, mseSte
		readPtr1 = &mse[i_alpha*Nlambdas];
		readPtr2 = &mseStd[0];
		F77_CALL(dcopy)(&Nlambdas,readPtr2,&inci,readPtr1,&incj);
		readPtr1 = &mseSte[i_alpha*Nlambdas];
		readPtr2 = &mseStd[Nlambdas];
		F77_CALL(dcopy)(&Nlambdas,readPtr2,&inci,readPtr1,&incj);
		
		
		
	}
	
	//find the min of ErrorEN;
	double minEN 			= ErrorEN[0];
	int ind_minEN 			= 0;
	
	if(L_alpha>1)
	{
		for(i_alpha = 1;i_alpha<L_alpha;i_alpha++)
		{
			if(minEN > ErrorEN[i_alpha])
			{
				minEN 			= ErrorEN[i_alpha];
				ind_minEN 		= i_alpha;
			}	
		}
		
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
		
	}
	ilambda_cv_ms 			= lambdaEN[ind_minEN];
	alpha_factor  			= alpha_factors[ind_minEN];
	if(verbose>=0) Rprintf("\tAdaptive_EN %d-fold CV, alpha: %f.\n", Kcv,alpha_factor);

	
	if(verbose==0) Rprintf("Step 3: CV support; alpha: %f, number of lambda needed: %d\n", alpha_factor,ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) R_Calloc(M, double);
	meanX = (double* ) R_Calloc(M, double);
	Ycopy = (double* ) R_Calloc(MN, double);
	Xcopy = (double* ) R_Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) R_Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda;
	double lambda_max;
	
	lambda_max = lambdaMax_adaEN(Ycopy,Xcopy,W,M, N,alpha_factor);
	
	if(verbose==0) Rprintf("Step 4: lasso selection path.\n");

	for(ilambda = 0;ilambda<ilambda_cv_ms;ilambda++)
	{
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, ilambda_cv_ms,lambda_factor);
		lambda = Weighted_LassoSf_adaEN(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max,
					alpha_factor); 	// mueL not calculated
					if(verbose>0) Rprintf("\tlambda: %f\n", lambda);
		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];

	}//ilambda; selection path
	out_alpha[0] = alpha_factor;
	out_lambda[0] = lambda_factor;
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
	

	
	
	R_Free(Strue);
	R_Free(meanY);
	R_Free(meanX);
	R_Free(rho_factors);
	R_Free(Ycopy);
	R_Free(Xcopy);
	R_Free(W);
	R_Free(QIBinv);
	R_Free(Q);
	
	//------adaEN 		
	R_Free(ErrorEN);
	R_Free(lambdaEN);
	
	R_Free(ErrorEN_min);
	R_Free(steEN_min);

}


//point estimate given 1 alpha, 1lambda, estimate B,F
//facilitate further computation such as:
//1.  	User decided lambda
//2.	Stability selection

void mainSML_adaENpointLambda(double *Y, double *X, int *m, int *n, int *Missing, double*B, double *f,double*stat,
			double*alpha_factors,//int *nAlpha, 	 must be scalar
			double *lambda_factors, int *nLambda, double *mseStd, int*VB) 						// mseStd: nLmabda x 2 matrix, keep mse + std 
{
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
	Strue 					= (double* ) R_Calloc(MM, double);
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
	rho_factors 			= (double* ) R_Calloc(L, double);
	for(i=0;i<L;i++) 
	{
		rho_factors[i] 		= pow(10.0,step);
		step 				= step + 0.2;
	}

	double  *ErrorEN,*lambdaEN, sigma2learnt;	//*alpha_factors,
	int L_alpha  		= 1; 
	ErrorEN 			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	lambdaEN 			= (double* ) R_Calloc(L_alpha, double);
	double *ErrorEN_min, *steEN_min;
	ErrorEN_min			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	steEN_min 			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.

	step 					= 0.05;

	//---------------------------------------------END  SYSTEM PARAMETERS
	int Nlambdas,Nrho;
	//----------------------------------------------------------R-package
	Nlambdas 				= 1; 			

	Nrho 					= L;
//call ridge_cvf
	double sigma2; //return value;
	double *W;  //weight on diagonal?
	W = (double* ) R_Calloc(MM, double);

	double *QIBinv;
	QIBinv = (double* ) R_Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1

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

	ilambda_cv_ms 			= lambdaEN[ind_minEN];
	alpha_factor  			= alpha_factors[ind_minEN];
	if(verbose>=0) Rprintf("\tAdaptive_EN %d-fold CV, alpha: %f.\n", Kcv,alpha_factor);
	//----------------------------------------------------------R-package
	Nlambdas 				= nLambda[0]; 			

	//----------------------------------------------------------R-package	
	if(verbose==0) Rprintf("Step 3: CV support; alpha: %f, number of lambda needed: %d\n", alpha_factor,ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) R_Calloc(M, double);
	meanX = (double* ) R_Calloc(M, double);
	Ycopy = (double* ) R_Calloc(MN, double);
	Xcopy = (double* ) R_Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) R_Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda; 
	double lambda_max;
	
	lambda_max = lambdaMax_adaEN(Ycopy,Xcopy,W,M, N,alpha_factor);

	if(verbose==0) Rprintf("Step 4: lasso selection path.\n");

	for(ilambda = 0;ilambda<Nlambdas;ilambda++)
	{
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, Nlambdas,lambda_factor);
		// call Weighted_LassoSf		Missing,
		lambda = Weighted_LassoSf_adaEN(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max,
					alpha_factor); 	// mueL not calculated
					if(verbose>0) Rprintf("\tlambda: %f\n", lambda);
		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];

	}//ilambda; selection path

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
	R_Free(Strue);
	R_Free(meanY);
	R_Free(meanX);
	R_Free(rho_factors);
	R_Free(Ycopy);
	R_Free(Xcopy);
	R_Free(W);
	R_Free(QIBinv);
	R_Free(Q);
	
	//------adaEN 		
	R_Free(ErrorEN);
	R_Free(lambdaEN);
		R_Free(ErrorEN_min);
	R_Free(steEN_min);
	
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
			double *Bout,
			int*kFold) 						// mseStd: nLmabda x 2 matrix, keep mse + std 
{
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
	Strue 					= (double* ) R_Calloc(MM, double);
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
	int Kcv 				= kFold[0];

	double step 	= -0.2;
	//rho factor
	step 					= -6;
	double * rho_factors;
	int L 					= 31; // number of rho_factors
	rho_factors 			= (double* ) R_Calloc(L, double);
	for(i=0;i<L;i++) 
	{
		rho_factors[i] 		= pow(10.0,step);
		step 				= step + 0.2;
	}

	double  *ErrorEN,*lambdaEN, sigma2learnt;	//*alpha_factors,
	int L_alpha  		= 1; 
	ErrorEN 			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	lambdaEN 			= (double* ) R_Calloc(L_alpha, double);

	double *ErrorEN_min, *steEN_min;
	ErrorEN_min			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	steEN_min 			= (double* ) R_Calloc(L_alpha, double);//rank0 variable actually.
	
	
	step 					= 0.05;
	
	
	//---------------------------------------------END  SYSTEM PARAMETERS
	int Nlambdas,Nrho;
	//----------------------------------------------------------R-package_1alpha point lambda
	Nlambdas 				= 1; 			

	//----------------------------------------------------------R-package_1alpha point lambda	
	Nrho 					= L;
	double sigma2; //return value;

	double *W;  //weight on diagonal?
	W = (double* ) R_Calloc(MM, double);

	double *QIBinv;
	QIBinv = (double* ) R_Calloc(MM, double);
	double beta = 0;
	F77_CALL(dcopy)(&MM,&beta,&inc0,QIBinv,&incj);
	for(i=0;i<M;i++) QIBinv[i*M+i] =1;
		
	//
	int ilambda_cv_ms; //index + 1

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

	ilambda_cv_ms 			= lambdaEN[ind_minEN];
	alpha_factor  			= alpha_factors[ind_minEN];
	if(verbose>=0) Rprintf("\tAdaptive_EN %d-fold CV for ridge-weight, alpha: %f.\n", Kcv,alpha_factor);
	//----------------------------------------------------------R-package
	Nlambdas 				= nLambda[0]; 			

	//----------------------------------------------------------R-package	
	if(verbose==0) Rprintf("Step 3: CV support; alpha: %f, number of lambda needed: %d\n", alpha_factor,ilambda_cv_ms);
	
	//call centerYX;
	double *meanY, *meanX, *Ycopy, *Xcopy;
	meanY = (double* ) R_Calloc(M, double);
	meanX = (double* ) R_Calloc(M, double);
	Ycopy = (double* ) R_Calloc(MN, double);
	Xcopy = (double* ) R_Calloc(MN, double);
	
	
	//copy Y,X
	F77_CALL(dcopy)(&MN,X,&inci,Xcopy,&incj);
	F77_CALL(dcopy)(&MN,Y,&inci,Ycopy,&incj);
	
	centerYX(Ycopy,Xcopy, meanY, meanX,M, N);
	// call Q_start
	double *Q;
	Q = (double* ) R_Calloc(MM, double);
	QlambdaStart(Ycopy,Xcopy, Q, sigma2, M, N);
	
	//selection path
	double lambda_factor_prev = 1.0;
	double lambda_factor;
	int ilambda;
	double lambda; lambda = 1.0;
	double lambda_max;
	//----------------------------------------------------------R-package_stability selection
	alpha_factor = 1;
	double *readPtr;
	//----------------------------------------------------------R-package_stability selection		
	lambda_max = lambdaMax_adaEN(Ycopy,Xcopy,W,M, N,alpha_factor);
	if(verbose==0) Rprintf("Step 4: lasso/elasticNet selection path.\n");
	for(ilambda = 0;ilambda<Nlambdas;ilambda++)
	{
		alpha_factor = alpha_factors[ilambda];
		lambda_factor = lambda_factors[ilambda];
		if(verbose>0) Rprintf("\t%d/%d lambdas. \tlambda_factor: %f", ilambda, Nlambdas,lambda_factor);
		lambda = Weighted_LassoSf_adaEN(W, B, f, Y,X, Q, lambda_factor,lambda_factor_prev, 
					sigma2, maxiter, M, N, verbose,QIBinv,lambda_max,
					alpha_factor); 	// mueL not calculated
		if(verbose>0) Rprintf("\tlambda: %f\n", lambda);			

		QlambdaMiddleCenter(Ycopy,Xcopy, Q,B,f,sigma2,M, N,QIBinv); //same set of Y,X 		<-- mueL not needed
		lambda_factor_prev = lambda_factors[ilambda];
	//----------------------------------------------------------R-package_stability selection		
		readPtr = &Bout[ilambda*MM];
		F77_CALL(dcopy)(&MM,B,&inci,readPtr,&incj);	
	//----------------------------------------------------------R-package_stability selection		
	}//ilambda; selection path

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
	R_Free(Strue);
	R_Free(meanY);
	R_Free(meanX);
	R_Free(rho_factors);
	R_Free(Ycopy);
	R_Free(Xcopy);
	R_Free(W);
	R_Free(QIBinv);
	R_Free(Q);
	
	//------adaEN 		
	R_Free(ErrorEN);
	R_Free(lambdaEN);
	
		R_Free(ErrorEN_min);
	R_Free(steEN_min);

}









