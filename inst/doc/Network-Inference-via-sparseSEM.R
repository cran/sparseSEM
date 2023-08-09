## -----------------------------------------------------------------------------
rm(list = ls())
library(sparseSEM)

## -----------------------------------------------------------------------------
data(B);
data(Y);
data(X);
data(Missing);
cat("dimenstion of Y: ",dim(Y) )
cat("dimenstion of X: ",dim(X) )


## -----------------------------------------------------------------------------
set.seed(1)
output = elasticNetSEM(Y, X, Missing, B, verbose = 1); 
names(output)

## -----------------------------------------------------------------------------
fit_SEM = output$weight


## -----------------------------------------------------------------------------
library('plot.matrix')
par(mfrow = c(1, 2),mar=c(5.1, 4.1, 4.1, 4.1))
plot(B)
plot(fit_SEM)


## -----------------------------------------------------------------------------
set.seed(1)
cvfit = elasticNetSEMcv(Y, X, Missing, B, alpha_factors = c(0.75, 0.5, 0.25),
                      lambda_factors=c(0.1, 0.01, 0.001), kFold = 5, verbose  = 1);
names(cvfit)

## -----------------------------------------------------------------------------
head(cvfit$cv)

## -----------------------------------------------------------------------------
par(mfrow = c(1, 2),mar=c(5.1, 4.1, 4.1, 4.1))
plot(B)
plot(cvfit$fit$weight)

## -----------------------------------------------------------------------------
cvfit$fit$statistics


## -----------------------------------------------------------------------------
tStart = proc.time()
set.seed(0)
output = enSEM_stability_selection(Y,X, Missing,B,
                                 alpha_factors = seq(1,0.1, -0.1), 
                                 lambda_factors =10^seq(-0.2,-3,-0.2), 
                                 kFold = 5,
                                 nBootstrap = 100,
                                 verbose = -1)
tEnd = proc.time()
simTime = tEnd - tStart;
print(simTime)
names(output)
cat("nSTS = ", length(which(output$STS !=0)))

## -----------------------------------------------------------------------------
B[which(B!=0)] =1
par(mfrow = c(1, 2),mar=c(5.1, 4.1, 4.1, 4.1))
plot(B)
plot(output$STS)

## ---- eval=FALSE--------------------------------------------------------------
#  library(parallel)
#  cl<-makeCluster(4,type="SOCK")
#  clusterEvalQ(cl,{library(sparseSEM)})
#  output = enSEM_stability_selection_parallel(Y,X, Missing,B,
#                                        alpha_factors = seq(1,0.1, -0.1),
#                                        lambda_factors =10^seq(-0.2,-3,-0.2),
#                                        kFold = 3,
#                                        nBootstrap = 100,
#                                        verbose = -1,
#                                        clusters = cl)
#  stopCluster(cl)	

## ---- eval=FALSE--------------------------------------------------------------
#  rm(list = ls())
#  library(sparseSEM)
#  data(yeast)
#  output = elasticNetSEM(Y, X, verbose = 1)
#  # STS
#  STS = enSEM_stability_selection(Y,X,verbose = -1)

