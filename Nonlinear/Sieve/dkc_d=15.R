

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
  ds <- mean(delta * Y) + X_pl_mis_mean %*% beta_ds
  
  ##
  or <- mean(Y)
  cc <- mean(Y[delta==1])
  
  
  res <- list(or=or, cc=cc, ds = ds)
  return(res)
  gc()
}


mu <- 2
simul <- 200
N <- 200000
d <- 15


kc_seq <- seq(0.1,2.1,0.4)

L_seq <- c(100)#c(10,20,50,100,200,500)

Res.bias.ds <- rep(NA, length(L_seq))
Res.sd.ds <- rep(NA, length(L_seq))


for(j in 1:length(L_seq)){
  L <- L_seq[j]
  n <- floor(N/L)
  for (k in 1:length(kc_seq)) {
    kc <- kc_seq[k]
    cl <- parallel::makeCluster(20, setup_strategy = "sequential")
    parallel::clusterExport(cl, c("N","L","n","d", "kc"))
    parallel::clusterEvalQ(cl, c(library("MASS"), source('Functions.R')))
    Res <- parallel::parLapply(cl, 1:simul, sim_fun)
    parallel::stopCluster(cl)
    
    Res.or <- matrix(NA,1,simul)
    Res.cc <- matrix(NA,1,simul)
    Res.ds <- matrix(NA,1,simul)
    
    for (i in 1:simul) {
      Res.or[1,i] <- Res[[i]]$or
      Res.cc[1,i] <- Res[[i]]$cc
      Res.ds[1,i] <- Res[[i]]$ds
    }
    
    
    setting <- data.frame(d=d,
                          N=N,
                          L=L,
                          n=n,
                          simul=simul,
                          kc = kc)
    print(setting)
    
    ## save results
    
    L_d_Res_ds <- paste("kc", as.character(kc),  "d", as.character(d), "Res_ds", sep = "_", collapse = NULL)
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
  }
}

print(proc.time()-tmp)




