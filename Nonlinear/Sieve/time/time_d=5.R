
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tmp <- proc.time()
source("Functions.R")



mu <- 2
simul <- 20
N <- 200000
d <- 5
kc <- 0.5


L_seq <- c(10, 20, 50, 100, 200, 500)

sim_time <- function(sim) {
  set.seed(sim)
  X1 <- rnorm(N,0,1)
  X2 <- rnorm(N,0,1) 
  X3 <- rnorm(N,0,1)
  X4 <- 2 * runif(N,0,1) - 1
  X5 <- 2 * runif(N,0,1) - 1
  X <- cbind(X1, X2, X3, X4, X5)
  gamma <- c(1, 1, 1, 1, 1)
  p.miss <- 0.5 * exp(X %*% gamma)/(1 + exp(X %*% gamma)) + 0.5
  delta <- as.numeric(p.miss > runif(N,0,1))
  Y <- 2 + 1 * sin(X1 + X2 + X3)+2 * qnorm((X4 + 1) / 2)+2*qnorm((X5 + 1) / 2)+rnorm(N,0,1)
  
  ####
  t1 <- proc.time()
  
  KN <- ceiling(kc * N^{d/(1 + 2*d)})
  
  X_pl <- powers(X, KN)
  
  
  Y_obs <- Y[delta==1]
  X_pl_obs <- X_pl[delta==1, ]
  
  hessian <- t(X_pl_obs) %*% X_pl_obs / N
  inv <- solve(hessian + 1e-7 * diag(dim(hessian)[1])) 
  rm(hessian)
  beta <- inv %*% t(X_pl_obs) %*% Y_obs / N #compute estimate from a single machine
  rm(inv)
  
  X_pl_mis_mean <- apply((1 - delta) * X_pl[1:N, 1:(1 + KN)], 2, mean)
  all <- mean(delta * Y) + X_pl_mis_mean %*% beta
  
  all.time <- (proc.time()-t1)[[3]]
  
  ds.time <- rep(NA, length(L_seq))
  
  for(j in 1:length(L_seq)){
    L <- L_seq[j]
    n <- floor(N/L)
    
    t1 <- proc.time()
    
    mark <- 1:n
    
    Y_1 <- Y[mark]
    dlt_1 <- delta[mark]
    X_1 <- X[mark, ]
    X_pl_1 <- powers(X_1, KN)
    Y_1_obs <- Y_1[dlt_1==1]
    X_pl_1_obs <- X_pl_1[dlt_1==1, ]
    
    alpha <- 0.015 * log(KN, 10)^2 * sqrt(KN / n)
    hessian_ds <- t(X_pl_1_obs) %*% X_pl_1_obs / n
    inv_ds <- solve(hessian_ds + alpha * diag(dim(hessian_ds)[1]))
    rm(hessian_ds)
    beta_ds <- inv_ds %*% t(X_pl_1_obs) %*% Y_1_obs / n
    
    iter <- max(floor(30 * max(log(L, 10), 1) * (log(N, 10) + log(KN, 10)) / log(max(1 / alpha, log(KN, 10)), 10)),2)
    
    t_par <- c()
    grad <- matrix(NA, KN + 1, L)
    
    for (t in 1:iter) {
      
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
      grad_sum <- grad %*% rep(1, L)
      beta_ds <- beta_ds - inv_ds %*% grad_sum / N
    }
    
    Y_sum <- c()
    X_pl_mis_sum <- matrix(NA, KN + 1, L)
    Y_sum[1] <- sum(Y_1_obs)
    X_pl_mis_sum[, 1] <- t(X_pl_1[dlt_1 == 0, ]) %*% rep(1, n - sum(dlt_1))
    
    t_tmp <- proc.time()
    for (l in 2:L) {
      
      mark <- ((l-1)*n+1):(l*n)
      
      Y_l <- Y[mark]
      dlt_l <- delta[mark]
      X_pl_l <- X_pl[mark,]
      
      Y_l_obs <- Y_l[dlt_l==1]
      
      Y_sum[l] <- sum(Y_l_obs)
      X_pl_mis_sum[, l] <- t(X_pl_l[dlt_l == 0, ]) %*% rep(1, n - sum(dlt_l))
    }
    
    t_par2 <- (proc.time() - t_tmp)[3]
    
    X_pl_mis_sum <- X_pl_mis_sum %*% rep(1, L)
    ds <- (sum(Y_sum) + t(X_pl_mis_sum) %*% beta_ds) / N
    
    ds.time[j] <- (proc.time()-t1)[[3]] - sum(t_par) - t_par2
  }
  c(all.time, ds.time)
}

cl <- parallel::makeCluster(20, setup_strategy = "sequential")
parallel::clusterExport(cl, c("N", "d", "mu", "L_seq", "kc"))
parallel::clusterEvalQ(cl, c(library("MASS"), source('Functions.R')))
Res <- parallel::parSapply(cl, 1:simul, sim_time, simplify = "array")
parallel::stopCluster(cl)

print(apply(Res, 1, mean))
print(proc.time()-tmp)