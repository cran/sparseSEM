lassoSEM <-
function(Y,X,Missing = NULL,B= NULL,verbose = 5){
	M = nrow(Y);
	N = ncol(Y);
	if (is.null(Missing)) Missing = matrix(0,M, N);
	if (is.null(B)) B = matrix(0,M,M);
	if(nrow(X) !=M){
	  if(verbose>=0) cat("error: sparseSEM currently support only the same dimension of X, Y.");
	  return( NULL);
	  
	}
	
	if(verbose>=0) cat("\tLASSO SML version_1;",M, "Genes, ", N , "samples; verbose: ", verbose, "\n\n")
	f = matrix(1,M,1);
	stat = rep(0,6);

	#dyn.load("lassoSMLv10beta.dll")
	tStart 	= proc.time();
	output<-.C("mainSML",
			Y 	= as.double(Y),
			X 	= as.double(X),
			M  	= as.integer(M),
			N 		= as.integer(N),			
			Missing 	= as.integer(Missing),
			B 	= as.double(B),
			f = as.double(f),
			stat = as.double(stat),
			verbose = as.integer(verbose),
			package = "sparseSEM"); 

	tEnd = proc.time();
	simTime = tEnd - tStart;
	#dyn.unload("lassoSMLv9beta.dll")
	if(verbose>=0) cat("\t computation time:", simTime[1], "sec\n");

	Bout = matrix(output$B,nrow= M, ncol = M, byrow = F);
	fout = matrix(output$f,nrow= M, ncol = 1, byrow = F);
	stat = matrix(output$stat,nrow = 6,ncol = 1, byrow = F);
	stat

	SMLresult 			<- list(Bout,fout,stat,simTime[1]);
	names(SMLresult)	<-c("weight","F","statistics","simTime")
	return(SMLresult)
}
