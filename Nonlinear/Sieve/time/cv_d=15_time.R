

rm(list = ls())

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
  X6 <- rnorm(N,0,1)
  X7 <- rnorm(N,0,1) 
  X8 <- rnorm(N,0,1)
  X9 <- 2 * runif(N,0,1) - 1
  X10 <- 2 * runif(N,0,1) - 1
  X11 <- rnorm(N,0,1)
  X12 <- rnorm(N,0,1) 
  X13 <- rnorm(N,0,1)
  X14 <- 2 * runif(N,0,1) - 1
  X15 <- 2 * runif(N,0,1) - 1
  X <- cbind(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15)
  gamma <- c(rep(1, 5), rep(0, 10))
  p.miss <- 0.5 * exp(X %*% gamma)/(1 + exp(X %*% gamma)) + 0.5
  delta <- as.numeric(p.miss > runif(N,0,1))
  Y <- 2 + 1 * sin(X1 + X2 + X3)+2 * qnorm((X4 + 1) / 2)+2*qnorm((X5 + 1) / 2)+rnorm(N,0,1)
  
  ### sample
  ### find h_opt by CV
  
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
  grad <- matrix(NA, KN_opt + 1, L)
  
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
  X_pl_mis_sum <- matrix(NA, KN_opt + 1, L)
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
  
  with_w_time <- proc.time()[[3]] - ti - sum(t2)*((L-1)/L) - t3 - sum(t_par) - t_par2
  
  
  ##
  or <- mean(Y)
  cc <- mean(Y[delta==1])
  
  
  res <- list(or=or, cc=cc, ds = ds, cv = cvm, time = with_w_time)
  return(res)
  gc()
}


mu <- 2
simul <- 20
N <- 200000
d <- 15


kc_seq <- seq(0.1,2.1,0.4)

L_seq <- c(100)#c(10,20,50,100,200,500)

Res.bias.ds <- rep(NA, length(L_seq))
Res.sd.ds <- rep(NA, length(L_seq))


for(j in 1:length(L_seq)){
  L <- L_seq[j]
  n <- floor(N/L)
  
  cl <- parallel::makeCluster(20, setup_strategy = "sequential")
  parallel::clusterExport(cl, c("N","L","n","d", "kc_seq"))
  parallel::clusterEvalQ(cl, c(library("MASS"), source('Functions.R')))
  Res <- parallel::parLapply(cl, 1:simul, sim_fun)
  parallel::stopCluster(cl)
  
  Res.or <- matrix(NA,1,simul)
  Res.cc <- matrix(NA,1,simul)
  Res.ds <- matrix(NA,1,simul)
  Res.time <- matrix(NA,1,simul)
  
  Res.cv <- list()
  
  for (i in 1:simul) {
    Res.or[1,i] <- Res[[i]]$or
    Res.cc[1,i] <- Res[[i]]$cc
    Res.ds[1,i] <- Res[[i]]$ds
    Res.time[1,i] <- Res[[i]]$time
    Res.cv[[i]] <- Res[[i]]$cv
  }
  
  
  setting <- data.frame(d=d,
                        N=N,
                        L=L,
                        n=n,
                        simul=simul,
                        kc_Seq = kc_seq)
  print(setting)
  
  ## save results
  L_d_Res_cv <- paste("L", as.character(L),  "d", as.character(d), "Res_cv", sep = "_", collapse = NULL)
  assign(L_d_Res_cv, Res.cv)
  
  L_d_Res_time <- paste("L", as.character(L),  "d", as.character(d), "Res_time", sep = "_", collapse = NULL)
  assign(L_d_Res_time, Res.time)
  
  L_d_Res_ds <- paste("L", as.character(L),  "d", as.character(d), "Res_ds", sep = "_", collapse = NULL)
  assign(L_d_Res_ds, Res.ds)
  
  ## print results
  Res.bias.or <- mean(Res.or)-mu
  Res.sd.or <- sqrt(mean((Res.or-mean(Res.or))^2))
  Res.bias.cc <- mean(Res.cc)-mu
  Res.sd.cc <- sqrt(mean((Res.cc-mean(Res.cc))^2))
  
  Res.bias.ds[j] <- mean(abn_rmv(Res.ds))-mu
  Res.sd.ds[j] <- sqrt(mean((abn_rmv(Res.ds)-mean(abn_rmv(Res.ds)))^2))
  
  or_bias_sd <- list(OR.bias.sd = round(c(Res.bias.or,Res.sd.or),4))
  cc_bias_sd <- list(CC.bias.sd = round(c(Res.bias.cc,Res.sd.cc),4))
  ds_bias_sd <- list(DS.bias.sd = round(c(Res.bias.ds[j],Res.sd.ds[j]),4))
  
  print(or_bias_sd)
  print(cc_bias_sd)
  print(ds_bias_sd)
  print(mean(Res.time))
  
  
}

print(proc.time()-tmp)




