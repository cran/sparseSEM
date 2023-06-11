#' Stability Selection 
#' Ref1: Meinshausen N. and Buhlmann P, J. R. Statist. Soc. B (2010) 72
#' Ref2: Shah R. and Samworth R,J. R. Statist. Soc. B (2013) 75
#'output:
#'for a list of candidate threshold (denoted as pi in the reference paper):
#'col1: pi
#'col2: per-family error rate
#'col3: E(v)
#'col4: E(v)_ShahR
#'col5: nSTS - the number of stable selected edges
#'col6: FDR
#'col7: FDR_ShahR
#'
#'The optimal threshold from simulation is the one with the smallest FDR and largest true nSTS.
#' snow -> parallel
enSEM_stability_selection_parallel <- function(Y,
                                      X,
                                      Missing = NULL,
                                      B=NULL,
                                      alpha_factors = seq(1,0.05, -0.05), 
                                      lambda_factors =10^seq(-0.2,-4,-0.2), 
                                      kFold = 5,
                                      nBootstrap = 100,
                                      verbose = 0,
                                      clusters = NULL){
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Parallel STS requires parallel package", call. = FALSE)
  }
  requireNamespace("parallel")
  M = nrow(Y)
  N = ncol(Y)
  if (is.null(Missing)) Missing = matrix(0,M, N);
  if (is.null(B)) B = matrix(0,M,M);
  
  if(nrow(X) !=M){
    if(verbose>=0) cat("error: sparseSEM currently support only the same dimension of X, Y.");
    return( NULL);
    
  }
  if(is.null(clusters)){
    if(verbose>=0) cat("error: enSEM_stability_selection_parallel requires input of clusters. Please run: \n\t
                       libary(parallel) \n\t
                       clusters = makeCluster(16, type = 'sock')\n\t
                       clusterEvalQ(clusters, library(sparseSEM)) \n ");
    return(NULL)
  }
  
  this.call=match.call()#returns a call in which all of the specified arguments are specified by their full names.
  if(verbose>=0) {
    cat("\tstability selection elastic net SEM;",M, "Nodes, ", N , "samples; verbose: ", verbose, "\n")
    if(nBootstrap!=0) cat("\t bootstrapping: ", nBootstrap, "\n\n") 
  }
  
  f 					= matrix(1,M,1);
  stat 				= rep(0,6);
  tStart 				= proc.time();
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
  nPara= nrow(parameters)
  #------------------------------------------------------R_package parameter
  MM = M*M
  NtopEff 				= MM - M;
  qEffect 				= rep(0,nBootstrap);
  piThreshold 			= seq(0.6,0.95, by = 0.01);
  
  qsEffect 				= matrix(0,MM,nBootstrap);
  Bselection 				= vector("list",nBootstrap);
  #Stability selection
  Bselection				= parallel::parLapply(clusters,rep(1,nBootstrap),enSEM_STS,
                                   Y,X,Missing ,B,
                                   STS_para =parameters,
                                   kFold = kFold,
                                   verbose = verbose) ; # column wise
  
  nStep 				= nPara;
  nTopEff 			= MM - M;
  qsEffect 			= matrix(0,MM,nBootstrap);
  qABave 				= rep(0,nStep); #ave of (alpha,lambda) selected effects over nBootstrap
  qAB 				= matrix(0,MM,nStep); #each effect: how many times selected with this ab
  for(i_ab in 1:nStep)
  {
    for(j_repeat in 1:nBootstrap)
    {
      qsEffect[,j_repeat] 	= Bselection[[j_repeat]][,i_ab];
    }
    #average number of effect for this a,b
    qABave[i_ab]    = mean(colSums(qsEffect));
    #each effect: how many times selected with this ab
    qAB[,i_ab]      = rowSums(qsEffect);
  }
  
  
  #re-grid the set of shrinkage parameters
  
  #re-grid the set of shrinkage parameters
  
  piThreshold         = seq(0.6, 1, 0.01);
  nPi                 = length(piThreshold);
  stablySelected      = vector("list",nPi);
  preFDR              = rep(0,nPi);
  preFDR2             = rep(0,nPi);
  pCER                = rep(0,nPi);
  nSTS                = rep(0,nPi);
  
  nGrid               = nStep;
  qStable             = rep(0, MM);
  FDRsetS             = vector("list",nGrid);
  for(i_grid in 1:nGrid)
  {
    nStep                   = i_grid;
    qABave_igrid            = qABave[1:i_grid];#average number of effect for this a,b
    qAB_igrid               = qAB[,1:i_grid, drop = FALSE];  #each effect: how many times selected with this ab   
    for (i in 1:MM)
    {
      qStable[i]          = max(qAB_igrid[i,]);  #each effect: how many times selected with this ab-set
    }
    for(i in 1:nPi)
    {
      piThresholdCut      = piThreshold[i]*nBootstrap; 
      stablySelected[[i]] = which(qStable>= piThresholdCut);
      nSTS[i]             = length(stablySelected[[i]]);
      pCER[i]             = mean(qABave_igrid)^2/((2*piThreshold[i] -1)*nTopEff^2);
      preFDR[i]           = mean(qABave_igrid)^2/((2*piThreshold[i] -1)*nTopEff);
      if(piThreshold[i]<=0.75)
      {
        preFDR2[i]      = mean(qABave_igrid)^2/(2*(2*piThreshold[i] -1-1/(2*nBootstrap))*nTopEff);
      }else
      {
        preFDR2[i]      = 4*(1-piThreshold[i] + 1/(2*nBootstrap))*mean(qABave_igrid)^2/((1+1/(2*nBootstrap))*nTopEff);
      }
    }
    FDR 					= preFDR;
    FDR2 					= preFDR2;
    FDRset 					= cbind(piThreshold,pCER,preFDR,preFDR2,nSTS,FDR,FDR2);
    for(i in 1:nPi){
      piThresholdCut      = piThreshold[i]* nBootstrap;
      FDRset[i,6] 		=  FDRset[i,6]/length(which(qStable>=piThresholdCut));
      FDRset[i,7] 		=  FDRset[i,7] /length(which(qStable>=piThresholdCut));
    }
    FDRsetS[[i_grid]]         = FDRset;
  }
  
  #what is the optimal i_grid? 
  OUT <- elasticNetSEMcv(Y, X, Missing, B, alpha_factors,
                         lambda_factors, kFold = 5, verbose = -1);
  mse = OUT$cv[,3]
  index = which.min(mse)
  j = index
  while(j>0){
    if (OUT$cv[j,3] > (OUT$cv[index,3] + OUT$cv[index,4] ) ){
      break;
    }
    j=j-1;
    
  }
  
  index = j+1;
  
  
  alpha = OUT$cv[index,1]
  lambda =  OUT$cv[index,2]
  ind1 = which(parameters[,1] == alpha);
  ind2 = which(parameters[,2] == lambda)
  index = intersect(ind1,ind2);
  if (length(index)>1){
    index = index[1];
  }
  
  FDRset = FDRsetS[[index]]
  
  
  colnames(FDRset) <- c('threshold','pre-comparison error rate','E(v)','E(v)_ShahR','nSTS','FDR','FDR_ShahR');
  tEnd = proc.time()
  simTime = tEnd - tStart;
  
  j = which.min(FDRset[,7])
  #stableEffcts = stablySelected[j]
  stableEffcts_index = stablySelected[j][[1]]
  data = rep(0, MM)
  data[stableEffcts_index] = 1
  stableEffects = matrix(data, M, M, byrow = FALSE)
  
  stat = FDRset[j,]
  computeData=list(simTime, piThreshold, stablySelected, FDRset, nSTS, nBootstrap, qsEffect, qEffect, NtopEff,FDRsetS,Bselection,parameters,qABave,qAB);
  names(computeData)	<-c("simTime","piThreshold","stablySelected","FDRset","nSTS","nBootstrap","bootstrap_B","bootstrap_count","NtopEff","FDRsetS","selectionB","parameters","qABave","qAB");
  output = list(stableEffects, stat, computeData)
  names(output) = c("STS", "statistics", "STS data")
  
  output$call = this.call;
  return(output)
}

