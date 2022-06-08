

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tmp <- proc.time()
source("Functions.R")

###

sim_fun <- function(sim) {
  set.seed(sim)
  X1 <- rnorm(N,0,1)
  X2 <- rnorm(N,0,1) 
  X3 <- rnorm(N,0,1)
  X4 <- 2 * runif(N,0,1) - 1
  X5 <- 2 * runif(N,0,1) - 1
  X <- cbind(X1, X2, X3, X4, X5)
  gamma <- rep(1, 5)
  p.miss <- 0.5 * exp(X %*% gamma)/(1 + exp(X %*% gamma)) + 0.5
  delta <- as.numeric(p.miss > runif(N,0,1))
  Y <- 2 + 1 * sin(X1 + X2 + X3)+2 * qnorm((X4 + 1) / 2)+2*qnorm((X5 + 1) / 2)+rnorm(N,0,1)
  
  
  ####
  
  KN <- ceiling(kc * N^{d/(1 + 2*d)})
  
  X_pl <- powers(X, KN)
  
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
  
  ds.time <- (proc.time()-t1)[[3]] - sum(t_par) - t_par2
  
  res <- ds.time
  return(res)
}


mu <- 2
simul <- 20
N <- 20000
d <- 5


kc_seq <- seq(0.1,2.1,0.4)

L <- 100
n <- floor(N/L)

Res.time <- matrix(NA,length(kc_seq),simul)

for (k in 1:length(kc_seq)) {
  kc <- kc_seq[k]
  cl <- parallel::makeCluster(20, setup_strategy = "sequential")
  parallel::clusterExport(cl, c("N","L","n","d", "kc"))
  parallel::clusterEvalQ(cl, c(library("MASS"), source('Functions.R')))
  Res <- parallel::parSapply(cl, 1:simul, sim_fun)
  parallel::stopCluster(cl)
  Res.time[k,] <- Res
  ## save results
}



setting <- data.frame(d=d,
                      N=N,
                      L=L,
                      n=n,
                      simul=simul,
                      kc = kc)
print(setting)

L_d_Res_time <- paste("kc", as.character(kc),  "d", as.character(d), "Res_time", sep = "_", collapse = NULL)
assign(L_d_Res_time, Res.time)

## print results

print(rowMeans(Res.time))

print(proc.time()-tmp)




