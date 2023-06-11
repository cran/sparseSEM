
elasticNetSEMpoint <- function(Y,X,Missing = NULL,B = NULL,alpha_factor = 1, lambda_factor =0.01 ,verbose = 0){
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
	lambda_factor_cand 	= 10^seq(-0.2, -4, -0.2);
	index 				= which(lambda_factor_cand>=lambda_factor);
	lambda_factors 		= lambda_factor_cand[index]
	
	nAlpha 				= length(alpha_factor);
	nLambda 			= length(lambda_factors);
	mseStd 				= rep(0,nLambda*2);
#------------------------------------------------------R_package parameter

	#dyn.load("elasticSMLv1.dll")
	tStart 				= proc.time();
	output<-.C("mainSML_adaENpointLambda",
				Y 		= as.double(Y),
				X 		= as.double(X),
				M  		= as.integer(M),
				N 		= as.integer(N),			
				Missing 	= as.integer(Missing),
				B 		= as.double(B),
				f 		= as.double(f),
				stat 	= as.double(stat),
				alpha 	= as.double(alpha_factor),
				lambda 	= as.double(lambda_factors),
				nLambda = as.integer(nLambda),
				mseStd 	= as.double(mseStd),
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
	mseStd = matrix(output$mseStd,nrow= nLambda, ncol = 2, byrow = F);
	mseStd = cbind(lambda_factors,mseStd);
	colnames(mseStd)<-c("lambda","mean Error", "Std");
#------------------------------------------------------R_package parameter


	SMLresult 			<- list(Bout,fout,stat,simTime[1]);
	names(SMLresult)	<-c("weight","F","statistics","simTime")
	SMLresult$call = this.call;
	
	return(SMLresult)
}