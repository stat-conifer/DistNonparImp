rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('Functions.R')

tmp=proc.time()

sim_fun <- function(sim) {
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
  alpha <- c(0.5,0.5,0.5,0.5,1)
  Y <- 2 + 1 * X %*% alpha + rnorm(N,0,1)
  
  ##
  
  KN <- ceiling(kc * N^{d/(1 + 2*d)})
  
  Kn <- ceiling(kc * n^{d/(1 + 2*d)})
  
  hop <- 1 + sum(KN > ref) #high order of polynomial
  
  X_pl <- powers(X, KN)

  ##
   
  mark <- 1:n
  
  Y_1 <- Y[mark]
  dlt_1 <- delta[mark]
  X_pl_1 <- X_pl[mark,]
  
  Y_1_obs <- Y_1[dlt_1==1]
  X_pl_1_obs <- X_pl_1[dlt_1==1, ]
  X_pl_1_obs_sg <- X_pl_1[dlt_1==1, 1:(1+Kn)]
  
  hessian_sg <- t(X_pl_1_obs_sg) %*% X_pl_1_obs_sg / n
  inv_sg <- solve(hessian_sg + 1e-7 * diag(dim(hessian_sg)[1])) 
  rm(hessian_sg)
  beta_sg <- inv_sg %*% t(X_pl_1_obs_sg) %*% Y_1_obs / n #compute estimate from a single machine
  rm(inv_sg)

  X_pl_mis_mean_1 <- apply((1 - dlt_1) * X_pl[1:n, 1:(1 + Kn)], 2, mean)
  sg <- mean(dlt_1 * Y_1) + X_pl_mis_mean_1 %*% beta_sg
  
  
  ## 
  
  alpha <- 0.015 * log(KN, 10)^2 * sqrt(KN / n)
  hessian_ds <- t(X_pl_1_obs) %*% X_pl_1_obs / n
  inv_ds <- solve(hessian_ds + alpha * diag(dim(hessian_ds)[1]))
  rm(hessian_ds)
  beta_ds <- inv_ds %*% t(X_pl_1_obs) %*% Y_1_obs / n
  
  iter <- max(floor(30 * max(log(L, 10), 1) * (log(N, 10) + log(KN, 10)) / log(max(1 / alpha, log(KN, 10)), 10)),2)
  for (t in 1:iter) {
    grad <- - t(X_pl_1_obs) %*% (Y_1_obs - X_pl_1_obs %*% beta_ds)
    for (l in 2:L) {
      
      mark <- ((l-1)*n+1):(l*n)
      
      Y_l <- Y[mark]
      dlt_l <- delta[mark]
      X_pl_l <- X_pl[mark,]
      
      Y_l_obs <- Y_l[dlt_l==1]
      X_pl_l_obs <- X_pl_l[dlt_l==1, ]
      
      grad <- grad - t(X_pl_l_obs) %*% (Y_l_obs - X_pl_l_obs %*% beta_ds)
      
      rm(X_pl_l, X_pl_l_obs)
    }
    beta_ds <- beta_ds - inv_ds %*% grad / N
  }
  
  X_pl_mis_mean <- apply((1 - delta) * X_pl, 2, mean)
  ds <- mean(delta * Y) + X_pl_mis_mean %*% beta_ds
  
  ## 
  
  or <- mean(Y)
  cc <- mean(Y[delta==1])
  
  res <- list(or=or, cc=cc, sg = sg, ds = ds)
  
  return(res)

}

mu <- 2
simul <- 200
N <- 200000
d <- 5
kc <- 0.5


L_seq <- c(500,200,100,50,20,10)

Res.bias.sg <- rep(NA, length(L_seq))
Res.sd.sg <- rep(NA, length(L_seq))

Res.bias.ds <- rep(NA, length(L_seq))
Res.sd.ds <- rep(NA, length(L_seq))


for(j in 1:length(L_seq)){
  L <- L_seq[j]
  n <- floor(N/L)
  
  cl <- parallel::makeCluster(20, setup_strategy = "sequential")
  parallel::clusterExport(cl, c("N", "d", "n", "kc", "L"))
  parallel::clusterEvalQ(cl, c(library("MASS"), source('Functions.R')))
  Res <- parallel::parLapply(cl, 1:simul, sim_fun)
  parallel::stopCluster(cl)
  
  Res.or <- matrix(NA,1,simul)
  Res.cc <- matrix(NA,1,simul)
  Res.sg <- matrix(NA,1,simul)
  Res.ds <- matrix(NA,1,simul)
  
  for (i in 1:simul) {
    Res.or[1,i] <- Res[[i]]$or
    Res.cc[1,i] <- Res[[i]]$cc
    Res.sg[1,i] <- Res[[i]]$sg
    Res.ds[1,i] <- Res[[i]]$ds
  }
  
  
  setting <- data.frame(d=d,
                        # hop=hop,
                        kc = kc,
                        N=N,
                        L=L,
                        n=n,
                        simul=simul)
  print(setting)
  
  ## save results
  L_d_Res_sg <- paste("L", as.character(L),  "d", as.character(d), "Res_sg", sep = "_", collapse = NULL)
  L_d_Res_ds <- paste("L", as.character(L),  "d", as.character(d), "Res_ds", sep = "_", collapse = NULL)
  
  assign(L_d_Res_sg, Res.sg)
  assign(L_d_Res_ds, Res.ds)
  
  ## print results
  Res.bias.or <- mean(abn_rmv(Res.or))-mu
  Res.sd.or <- sqrt(mean((abn_rmv(Res.or)-mean(abn_rmv(Res.or)))^2))
  Res.bias.cc <- mean(abn_rmv(Res.cc))-mu
  Res.sd.cc <- sqrt(mean((abn_rmv(Res.cc)-mean(abn_rmv(Res.cc)))^2))
  
  Res.bias.sg[j] <- mean(abn_rmv(Res.sg))-mu
  Res.sd.sg[j] <- sqrt(mean((abn_rmv(Res.sg)-mean(abn_rmv(Res.sg)))^2))
  Res.bias.ds[j] <- mean(abn_rmv(Res.ds))-mu
  Res.sd.ds[j] <- sqrt(mean((abn_rmv(Res.ds)-mean(abn_rmv(Res.ds)))^2))
  
  or_bias_sd <- list(OR.bias.sd = round(c(Res.bias.or,Res.sd.or),4))
  cc_bias_sd <- list(CC.bias.sd = round(c(Res.bias.cc,Res.sd.cc),4))
  sg_bias_sd <- list(sg.bias.sd = round(c(Res.bias.sg[j],Res.sd.sg[j]),4))
  ds_bias_sd <- list(DS.bias.sd = round(c(Res.bias.ds[j],Res.sd.ds[j]),4))
  
  print(or_bias_sd)
  print(cc_bias_sd)
  print(sg_bias_sd)
  print(ds_bias_sd)
  
}

print(proc.time() - tmp)

save.image("C:/Users/smm/Desktop/DC-code/Setting2-Linear/Sieve/linear_main_d=5.RData")

