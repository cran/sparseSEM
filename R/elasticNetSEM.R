elasticNetSEM<-
function(Y,X,Missing= NULL,B = NULL,verbose = 0){
	M = nrow(Y);
	N = ncol(Y);
	
	if (is.null(Missing)) Missing = matrix(0,M, N);
	if (is.null(B)) B = matrix(0,M,M);
	if(nrow(X) !=M){
	  if(verbose>=0) cat("error: sparseSEM currently support only the same dimension of X, Y.");
	  return( NULL);
	  
	}
	this.call=match.call()#returns a call in which all of the specified arguments are specified by their full names.
	if(verbose>=0) cat("\telastic net SML;",M, "Nodes, ", N , "samples; verbose: ", verbose, "\n\n")
	f = matrix(1,M,1);
	stat = rep(0,6);

	#dyn.load("elasticSMLv1.dll")
	tStart 	= proc.time();
	output<-.C("mainSML_adaEN",
				Y 	= as.double(Y),
				X 	= as.double(X),
				M  	= as.integer(M),
				N 		= as.integer(N),			
				Missing 	= as.integer(Missing),
				B 	= as.double(B),
				f = as.double(f),
				stat = as.double(stat),
				alpha = as.double(0),
				lambda = as.double(0),
				verbose = as.integer(verbose),
				package = "sparseSEM"); 

	tEnd = proc.time();
	simTime = tEnd - tStart;
	#dyn.unload("elasticSMLv1.dll")
	if(verbose>=0) cat("\t computation time:", simTime[1], "sec\n");

	Bout = matrix(output$B,nrow= M, ncol = M, byrow = F);
	fout = matrix(output$f,nrow= M, ncol = 1, byrow = F);
	stat = matrix(output$stat,nrow = 6,ncol = 1, byrow = F);
	rownames(stat) <- c("correct_positive", "total_ground truth", "false_positive", "true_positive", "Power", "FDR")
	
  hyperparameters = list(output$alpha, output$lambda);
	names(hyperparameters) = c("alpha", "lambda");
	SMLresult 			<- list(Bout,fout,stat,hyperparameters,simTime[1]);
	names(SMLresult)	<-c("weight","F","statistics","hyperparameters","runTime")
	
	SMLresult$call = this.call;
	return(SMLresult)
}
