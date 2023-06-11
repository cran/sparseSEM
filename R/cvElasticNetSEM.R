#'cv function in lasso secenario: generally select (alpha, labmda ) that is within min(mse) + 1ste
#'
#'
elasticNetSEMcv <- function(Y,X,Missing = NULL,B=NULL,alpha_factors = seq(1,0.05, -0.05), lambda_factors =10^seq(-0.2,-4,-0.2), kFold = 5, verbose = 0){
	M 					= nrow(Y);
	N 					= ncol(Y);
	if (is.null(Missing)) Missing = matrix(0,M, N);
	if (is.null(B)) B = matrix(0,M,M);
	if(nrow(X) !=M){
	  if(verbose>=0) cat("error: sparseSEM currently support only the same dimension of X, Y.");
	  return( NULL);
	  
	}
	this.call=match.call()#returns a call in which all of the specified arguments are specified by their full names.
	if(verbose>=0) cat("\telastic net SML;",M, "Nodes, ", N , "samples; verbose: ", verbose, "\n\n")
	f 					= matrix(1,M,1);
	stat 				= rep(0,6);
#------------------------------------------------------R_package parameter
	nAlpha 	= length(alpha_factors);
	nLambda = length(lambda_factors);
	mse 	= rep(0,nLambda*nAlpha);
	mseSte = rep(0,nLambda*nAlpha);
	mseStd = rep(0, nLambda*2); 
	
	parameters = matrix(0,0,2);
	for (i in alpha_factors){
	  col1 = rep(i,nLambda);
	  col2 = lambda_factors;
	  para = cbind(col1,col2);
	  parameters = rbind(parameters,para)
	}
	
#------------------------------------------------------R_package parameter

	#dyn.load("elasticSMLv1.dll")
	tStart 				= proc.time();
	output<-.C("mainSML_adaENcv",
				Y 		= as.double(Y),
				X 		= as.double(X),
				M  		= as.integer(M),
				N 		= as.integer(N),			
				Missing 	= as.integer(Missing),
				B 		= as.double(B),
				f 		= as.double(f),
				stat 	= as.double(stat),
				alpha 	= as.double(alpha_factors),
				nAlpha 	= as.integer(nAlpha),
				lambda 	= as.double(lambda_factors),
				nLambda = as.integer(nLambda),
				mse 	= as.double(mse),
				mseSte 	= as.double(mseSte),
				mseStd 	= as.double(mseStd),
				kFold   = as.integer(kFold),
				para_alpha = as.double(0),
				para_lambda= as.double(0),
				verbose = as.integer(verbose),
				package = "sparseSEM"); 

	tEnd 				= proc.time();
	simTime 			= tEnd - tStart;
	#dyn.unload("elasticSMLv1.dll")
	if(verbose>=0) cat("\t computation time:", simTime[1], "sec\n");

	Bout = matrix(output$B,nrow= M, ncol = M, byrow = F);
	fout = matrix(output$f,nrow= M, ncol = 1, byrow = F);
	stat = matrix(output$stat,nrow = 6,ncol = 1, byrow = F);
	rownames(stat) <- c("correct_positive", "total_ground truth", "false_positive", "true_positive", "Power", "FDR")
#------------------------------------------------------R_package parameter
#	mseStd = matrix(output$mseStd,nrow= nLambda, ncol = 2, byrow = F);
	mse   = matrix(output$mse,nrow = nAlpha*nLambda, ncol =1, byrow= F)
	mseSte   = matrix(output$mseSte,nrow = nAlpha*nLambda, ncol =1, byrow= F)

	
	cvResults = cbind(parameters,mse,mseSte);
	colnames(cvResults)<-c("alpha","lambda","mean Error", "Ste");
#------------------------------------------------------R_package parameter
	hyperparameters= list(output$para_alpha, output$para_lambda);
	names(hyperparameters) = c("alpha", "lambda");
  fit = list(hyperparameters, Bout,fout,stat, simTime[1])
	
  names(fit) = c("hyperparameter", "weight","F","statistics","simTime")
  fit$call=this.call
  
  SMLresult 			<- list(cvResults, fit);
	names(SMLresult)	<-c("cv", "fit")
	return(SMLresult)
}

