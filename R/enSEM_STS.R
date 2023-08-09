#' internal function for enSEM_stability_selection
#' i is a parameter holds for parLapply
#'output a MM by nPara matrix with 0 denotes B[i*M+j] = 0 denotes edge (i,j) has no connection and 1 otherwise
#'
#'
#' 
enSEM_STS <- function(i=0,
                      Y,
                      X,
                      Missing = NULL,
                      B=NULL,
                      STS_para = NULL,
                      kFold = 5,
                      verbose = 0){
  
  
  M 					= nrow(Y);
  N 					= ncol(Y);
  if (is.null(Missing)) Missing = matrix(0,M, N);
  if (is.null(B)) B = matrix(0,M,M);
  if(nrow(X) !=M){
    if(verbose>=0) cat("error: sparseSEM currently support only the same dimension of X, Y.");
    return( NULL);
    
  }
  if (is.null(STS_para)){
    if(verbose>=0) cat("error: enSEM_STS is an internal function with STS_para setup as 2-column matrix of [alphas, lambdas].");
    return(NULL);
  }
  if(verbose>=1) {
    #cat("\tstability selection elastic net SEM;",M, "Nodes, ", N , "samples; verbose: ", verbose, "\n")
    if(i!=0) cat("\t bootstrapping: ", i, "\n\n") 
  }
  alpha_factors 			= STS_para[,1];
  lambda_factors 			= STS_para[,2];
  nAlpha 				= length(alpha_factors);
  nLambda 			= length(lambda_factors);
  #
  f 					= matrix(1,M,1);
  stat 				= rep(0,6);
  mseStd 				= rep(0,nLambda*2);
  #------------------------------------------------------R_package stability selection
  MM 					= M*M;
  Bout 				= matrix(0,MM,nAlpha);
  nKeep 				= floor(N/2);
  keepSample 			= sample(N,nKeep);
  exprY = Y
  genoX = X
  Y 					= exprY[,keepSample];
  X 					= genoX[,keepSample];
  Missing     = Missing[,keepSample]
  N 					= ncol(Y);	
  #Missing 			= matrix(0,N,M);
  #B 					= matrix(0,N,M);
  tStart 				= proc.time();
  output<-.C("mainSML_adaENstabilitySelection",
             Y 		= as.double(Y),
             X 		= as.double(X),
             M  		= as.integer(M),
             N 		= as.integer(N),			
             Missing 	= as.integer(Missing),
             B 		= as.double(B),
             f 		= as.double(f),
             stat 	= as.double(stat),
             alpha 	= as.double(alpha_factors),
             nAlpha = as.integer(nAlpha),				
             lambda 	= as.double(lambda_factors),
             nLambda = as.integer(nLambda),
             mseStd 	= as.double(mseStd),
             verbose = as.integer(verbose),
             Bout 	= as.double(Bout),
             kFold   = as.integer(kFold),
             package = "sparseSEM");
  
  tEnd 				= proc.time();
  simTime 			= tEnd - tStart;
  if(verbose>=1) cat("\t computation time:", simTime[1], "sec\n");
  
  temp 	= as.numeric(output$Bout != 0);
  Bout 	= matrix(temp,MM,nAlpha,byrow = F); #
  return(Bout)
  
}


