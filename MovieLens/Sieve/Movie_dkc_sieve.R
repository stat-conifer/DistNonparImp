
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tmp <- proc.time()
source("Functions.R")

#########################################################################
load("MTransData.RData")
Feature <- MTransData[, -c(1, 2)]
MTransData <- cbind(MTransData[, c(1, 2)], Feature[, diag(var(Feature)) > 0.075])

MTransData[MTransData[, 2] == 0, 1] <- 0 
d <- length(MTransData[1, ]) - 2
N <- length(MTransData[, 1])

X <- MTransData[, -c(1,2)]
delta <- MTransData[, 2]
Y <- MTransData[, 1]
#########################################################################

kc_seq <- seq(0.1,2.1,0.4)

L <- c(100)
n <- floor(N/L)
ds <- c()
for (k in 1:length(kc_seq)) {
  kc <- kc_seq[k]
  ####
  KN <- ceiling(kc * N^{d/(1 + 2*d)})
  
  X_pl <- powers(X, KN)
  
  
  ## 
  
  mark <- 1:n
  
  Y_1 <- Y[mark]
  dlt_1 <- delta[mark]
  X_pl_1 <- X_pl[mark,]
  
  Y_1_obs <- Y_1[dlt_1==1]
  X_pl_1_obs <- X_pl_1[dlt_1==1, ]
  
  alpha <- 0.015 * log(KN, 10)^2 * sqrt(KN / n)
  hessian_ds <- t(X_pl_1_obs) %*% X_pl_1_obs / n
  inv_ds <- solve(hessian_ds + alpha * diag(dim(hessian_ds)[1]))
  rm(hessian_ds)
  beta_ds <- inv_ds %*% t(X_pl_1_obs) %*% Y_1_obs / n
  
  iter <- max(floor(30 * max(log(L, 10), 1) * (log(N, 10) + log(KN, 10)) / log(max(1 / alpha, log(KN, 10)), 10)),2)
  grad_fun <- function(l){
    mark <- ((l-1)*n+1):(l*n)
    
    Y_l <- Y[mark]
    dlt_l <- delta[mark]
    X_pl_l <- X_pl[mark,]
    
    Y_l_obs <- Y_l[dlt_l==1]
    X_pl_l_obs <- X_pl_l[dlt_l==1, ]
    
    return(-t(X_pl_l_obs) %*% (Y_l_obs - X_pl_l_obs %*% beta_ds))
    rm(X_pl_l, X_pl_l_obs)
  }
  
  for (t in 1:iter) {
    
    grad_seq <- matrix(unlist(lapply(1:L, grad_fun)),dim(X_pl_1_obs)[2])
    grad = rowSums(grad_seq)
    
    beta_ds <- beta_ds - inv_ds %*% grad / N
  }
  
  
  X_pl_mis_mean <- apply((1 - delta) * X_pl, 2, mean)
  ds[k] <- mean(delta * Y) + X_pl_mis_mean %*% beta_ds
  
}

print(ds)
print(proc.time()-tmp)



