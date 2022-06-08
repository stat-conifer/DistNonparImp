

tmp <- proc.time()
source("Functions.R")

hok <- 20
v <- rep(0,hok / 2) # coefficients of Legendre kernel of order hok
for (j in 1:length(v)) {
  qt <- 2*(j - 1)
  v[j] <- (-1)^{qt/2}*factorial(qt) / (2^qt * factorial(qt/2)^2) * (2*qt + 1) #zero order term times inverse of norm^2
}

#########################################################################

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
  
  ############ sample
  ### find h_opt by CV
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
    
    w_te <- rep(1,n_te_obs)
    
    cv_dev_m <- rep(NA, length(hc_seq))
    
    for(k in 1:length(hc_seq)){
      hc <- hc_seq[k]
      hc_n_cv <- hc * ((1/(n_tr))^(1/(1 + 2 * d)))
      cv_dev_m[k] <- CV_m(Y_training_obs, Y_test_obs, X_training_obs, X_test_obs, w_te, hc_n_cv)
      #print(cv_dev_m[k,s])
    }
    return(cv_dev_m)
  }
  
  cv_ms <- matrix(NA, length(hc_seq),L) 
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
  
  no_w_time <- (tf-ti-(t4-t3+t2-t1)*((L-1)/L))[3]
  
  or <- mean(Y)
  cc <- mean(Y[delta==1])
  # ds <- cc * mean(delta)
  
  res <- list(or=or, cc=cc, ds = ds, cv = cvm, time = no_w_time)
  return(res)
  gc()
}


mu <- 2
simul <- 20
N <- 200000
d <- 5

hc_seq <- seq(0.1,2.1,0.4)

L_seq <- c(100)#c(10,20,50,100,200,500)


Res.bias.ds <- rep(NA, length(L_seq))
Res.sd.ds <- rep(NA, length(L_seq))


for(j in 1:length(L_seq)){
  L <- L_seq[j]
  n <- floor(N/L)
  
  cl <<- parallel::makeCluster(40, setup_strategy = "sequential")
  parallel::clusterExport(cl, c("N","v","L","n","hok","d","hc_seq"))
  parallel::clusterEvalQ(cl,source("Functions.R"))
  Res <- parallel::parLapply(cl, 1:simul, sim_fun)
  parallel::stopCluster(cl)
  
  Res.or <- matrix(NA,1,simul)
  Res.cc <- matrix(NA,1,simul)
  Res.ds <- matrix(NA,1,simul)
  Res.time <- matrix(NA,1,simul)
  Res.cv <- matrix(NA,length(hc_seq),simul)
  
  for (i in 1:simul) {
    Res.or[1,i] <- Res[[i]]$or
    Res.cc[1,i] <- Res[[i]]$cc
    Res.ds[1,i] <- Res[[i]]$ds
    Res.time[1,i] <- Res[[i]]$time
    Res.cv[,i] <- Res[[i]]$cv
  }
  
  
  setting <- data.frame(d=d,
                        hok=hok,
                        N=N,
                        L=L,
                        n=n,
                        simul=simul,
                        hc_Seq = hc_seq)
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



