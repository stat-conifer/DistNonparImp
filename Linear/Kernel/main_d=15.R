
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
  alpha <- c(0.5,0.5,0.5,0.5,1, rep(0, 10))
  Y <- 2 + 1 * X %*% alpha + rnorm(N,0,1)
  
  ####
  hc <- 1.7
  hN <- hc * ((log(L, 10) * (L/N))^(1/(1 + 2 * d)))
  hn <- hc * ((1/n)^(1/(1 + 2 * d)))
  ########################################################################    
  
  mark <- 1:n
  
  Y.1 <- Y[mark]
  dlt.1 <- delta[mark]
  X.1 <- as.matrix(X[mark,],n,d)
  
  sg <- ker_mean(Y.1,dlt.1,X.1,hN)
  
  ds <- rep(NA,L)

  for(j in 1:L){
    #set.seed(j)
    mark <- ((j-1)*n+1):(j*n)

    Y.j <- Y[mark]
    dlt.j <- delta[mark]
    X.j <- as.matrix(X[mark,],n,d)

    ds[j] <- ker_mean(Y.j,dlt.j,X.j,hN)

  }
  
  ds <- mean(ds)
  or <- mean(Y)
  cc <- mean(Y[delta==1])
  
  res <- list(or=or, cc=cc, sg = sg, ds = ds)
  return(res)
  gc()
}


mu <- 2
simul <- 200
N <- 200000
d <- 15

L_seq <- c(10,20,50,100,200,500)

Res.bias.sg <- rep(NA, length(L_seq))
Res.sd.sg <- rep(NA, length(L_seq))


Res.bias.ds <- rep(NA, length(L_seq))
Res.sd.ds <- rep(NA, length(L_seq))


for(j in 1:length(L_seq)){
  L <- L_seq[j]
  n <- floor(N/L)
  
  cl <<- parallel::makeCluster(20, setup_strategy = "sequential")
  parallel::clusterExport(cl, c("N","v","L","n","hok","d"))
  parallel::clusterEvalQ(cl,source("Functions.R"))
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
                       hok=hok,
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

print(proc.time()-tmp)

#save.image("C:/Users/smm/Desktop/vary_m_fix_h/res_d=5.RData")

