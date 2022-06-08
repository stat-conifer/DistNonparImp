
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tmp <- proc.time()
source("Functions.R")

hok <- 20
v <- rep(0,hok / 2) # coefficients of Legendre kernel of order hok
for (j in 1:length(v)) {
  qt <- 2*(j - 1)
  v[j] <- (-1)^{qt/2}*factorial(qt) / (2^qt * factorial(qt/2)^2) * (2*qt + 1) #zero order term times inverse of norm^2
}
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

hc_seq <- seq(0.1,2.1,0.4)

L <- c(100)
n <- floor(N/L)


bws_fun <- function(l){
  mark <- ((l - 1) * n + 1) : (l * n)
  Y.mark <- Y[mark]
  dlt.mark <- delta[mark]
  X.mark <- as.matrix(X[mark,],n,d)
  
  ##
  
  mark <- rep(FALSE, n)
  mark[sample(n, round(n/2))] <- TRUE
  X_training <- X.mark[mark, ]
  Y_training <- Y.mark[mark]
  dlt_training <- dlt.mark[mark]
  X_training_obs <- X_training[dlt_training==1,]
  Y_training_obs <- Y_training[dlt_training==1]
  n_tr <- length(dlt_training)
  X_test <- X.mark[!mark, ]
  Y_test <- Y.mark[!mark]
  dlt_test <- dlt.mark[!mark]
  X_test_obs <- X_test[dlt_test==1,]
  Y_test_obs <- Y_test[dlt_test==1]
  X_test_mis <- X_test[dlt_test==0,]
  Y_test_mis <- Y_test[dlt_test==0]
  n_te_obs <- sum(dlt_test)
  # n_te_mis <- sum(1-dlt_test)
  # 
  # EX_test_mis <- cbind(rep(1,n_te_mis))
  
  X_M_te_mis <- moment(X_test_mis)
  mu_te <- colMeans(X_M_te_mis)
  X_M_te_obs <- moment(X_test_obs)
  w_te <- 1 - n_te_obs*X_M_te_obs%*%solve(t(X_M_te_obs)%*%X_M_te_obs)%*%(colMeans(X_M_te_obs)-mu_te)
  
  cv_dev_m <- rep(NA, length(hc_seq))
 
  for(k in 1:length(hc_seq)){
    hc <- hc_seq[k]
    hc_n_cv <- hc * ((1/(n_tr))^(1/(1 + 2 * d)))
    cv_dev_m[k] <- CV_m(Y_training_obs, Y_test_obs, X_training_obs, X_test_obs, w_te, hc_n_cv)
    #print(cv_dev_m[k,s])
  }
  return(cv_dev_m)
}

set.seed(d) 
cv_ms <- matrix(NA,length(hc_seq),L) 
ti <- proc.time()
t1 <- proc.time()
for (l in 1:L) {
  cv_ms[,l] = bws_fun(l)
}
t2 <- proc.time()

cvm <- rowMeans(cv_ms)
# 
####
hc_opt_m <- hc_seq[which.min(cvm)]
hN_opt <- hc_opt_m * ((log(L, 10) * (L/N))^(1/(1 + 2 * d)))

########################################################################    

ds <- rep(NA,L)
t3 <- proc.time()
for(l in 1:L){
  
  mark <- ((l-1)*n+1):(l*n)
  
  Y.l <- Y[mark]
  dlt.l <- delta[mark]
  X.l <- as.matrix(X[mark,],n,d)
  
  ds[l] <- ker_mean(Y.l,dlt.l,X.l,hN_opt)
  
}
t4 <- proc.time()

ds <- mean(ds)
tf = proc.time()

ds_time <- (tf-ti-(t4-t3+t2-t1)*((L-1)/L))[3]

print(ds)
print(cvm)
print(ds_time)
print(proc.time()-tmp)



