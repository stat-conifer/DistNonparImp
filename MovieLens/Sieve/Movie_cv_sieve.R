
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


bws_fun <- function(l){
  mark <- ((l - 1) * n + 1) : (l * n)
  Y.mark <- Y[mark]
  dlt.mark <- delta[mark]
  X.mark <- as.matrix(X[mark,],n,d)
  
  ##
  
  mark <- rep(FALSE, n)
  mark[sample(n, round(5*n/10))] <- TRUE
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
  
  kc_n_cv <- ceiling(kc * n_tr^{d/(1 + 2*d)})
  
  # t9 = proc.time()
  cv_val <- CV_m(Y_training_obs, Y_test_obs, X_training_obs, X_test_obs, w_te, kc_n_cv)
  # t10 = proc.time()-t9
  
  return(cv_val)
}

set.seed(d)
cvm <- c()

ti <- proc.time()[[3]]
t2 <- c()
for(j in 1:length(kc_seq)){
  kc <- kc_seq[j]
  t1 = proc.time()
  cv_ms <- c()
  for (l in 1:L) {
    cv_ms[l] <- bws_fun(l)
  }
  t2[j] = (proc.time()-t1)[[3]]
  cvm[j] <- mean(cv_ms)
  
  if(j>1){
    if((cvm[j-1]-cvm[j])/cvm[j-1]<0.01){
      break
    }
  }
  
}



####
kc_opt_m <- kc_seq[length(cvm)-1]
KN_opt <- ceiling(kc_opt_m * N^{d/(1 + 2*d)})

t1 <- proc.time()
X_pl <- powers(X, KN_opt)
t3 <- (proc.time()-t1)[[3]]

mark <- 1:n

Y_1 <- Y[mark]
dlt_1 <- delta[mark]
X_1 <- X[mark, ]
X_pl_1 <- powers(X_1, KN_opt)
Y_1_obs <- Y_1[dlt_1==1]
X_pl_1_obs <- X_pl_1[dlt_1==1, ]

alpha <- 0.015 * log(KN_opt, 10)^2 * sqrt(KN_opt / n)
hessian_ds <- t(X_pl_1_obs) %*% X_pl_1_obs / n
inv_ds <- solve(hessian_ds + alpha * diag(dim(hessian_ds)[1]))
rm(hessian_ds)
beta_ds <- inv_ds %*% t(X_pl_1_obs) %*% Y_1_obs / n

iter <- max(floor(30 * max(log(L, 10), 1) * (log(N, 10) + log(KN_opt, 10)) / log(max(1 / alpha, log(KN_opt, 10)), 10)),2)

t_par <- c()

for (t in 1:iter) {
  grad <- matrix(NA, length(X_pl_1[1, ]), L)
  grad[, 1] <- - t(X_pl_1_obs) %*% (Y_1_obs - X_pl_1_obs %*% beta_ds)
  
  t_tmp <- proc.time()
  
  for (l in 2:L) {
    
    mark <- ((l-1)*n+1):(l*n)
    
    Y_l <- Y[mark]
    dlt_l <- delta[mark]
    X_pl_l <- X_pl[mark,]
    
    Y_l_obs <- Y_l[dlt_l==1]
    X_pl_l_obs <- X_pl_l[dlt_l==1, ]
    
    grad[, l] <- - t(X_pl_l_obs) %*% (Y_l_obs - X_pl_l_obs %*% beta_ds)
    
    rm(X_pl_l, X_pl_l_obs)
  }
  
  t_par[t] <- (proc.time() - t_tmp)[[3]]
  grad <- apply(grad, 1, sum)
  beta_ds <- beta_ds - inv_ds %*% grad / N
}

X_pl_mis_mean <- apply((1 - delta) * X_pl, 2, mean)
ds <- mean(delta * Y) + X_pl_mis_mean %*% beta_ds

ds_time <- proc.time()[[3]]-ti-sum(t2)*((L-1)/L)-t3- sum(t_par)

print(ds)
print(cvm)
print(ds_time)
print(proc.time()-tmp)



