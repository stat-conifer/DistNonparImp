
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
  X <- cbind(X1, X2, X3, X4, X5)
  gamma <- c(1, 1, 1, 1, 1)
  p.miss <- 0.5 * exp(X %*% gamma)/(1 + exp(X %*% gamma)) + 0.5
  delta <- as.numeric(p.miss > runif(N,0,1))
  Y <- 2 + 1 * sin(X1 + X2 + X3)+2 * qnorm((X4 + 1) / 2)+2*qnorm((X5 + 1) / 2)+rnorm(N,0,1)
  
  ####
  l_hc <- function(l){
    mark <- ((l-1)*n+1):(l*n)
    Y.l <- Y[mark]
    dlt.l <- delta[mark]
    X.l <- as.matrix(X[mark,],n,d)
    ds_l <- c()
    for(j in 1:length(hc_seq)){
      hc <- hc_seq[j]
      hN <- hc * ((log(L, 10) * (L/N))^(1/(1 + 2 * d)))
      ds_l[j] <- ker_mean(Y.l,dlt.l,X.l,hN)
    }
    
    return(ds_l)
  }
  
  ds_hc <- matrix(unlist(lapply(1:L, l_hc)),length(hc_seq))
  ds <- rowMeans(ds_hc)
  or <- mean(Y)
  cc <- mean(Y[delta==1])
  # ds <- cc * mean(delta)
  
  res <- list(or=or, cc=cc, ds = ds)
  return(res)
  gc()
}


mu <- 2
simul <- 200
N <- 200000
d <- 5

hc_seq <- seq(0.1,2.1,0.4)

L <- 100
n <- floor(N/L)


cl <<- parallel::makeCluster(20, setup_strategy = "sequential")
parallel::clusterExport(cl, c("N","v","L","n","hok","d","hc_seq"))
parallel::clusterEvalQ(cl,source("Functions.R"))
Res <- parallel::parLapply(cl, 1:simul, sim_fun)
parallel::stopCluster(cl)

Res.or <- matrix(NA,1,simul)
Res.cc <- matrix(NA,1,simul)

Res.ds <- matrix(NA,length(hc_seq),simul)

for (i in 1:simul) {
  Res.or[1,i] <- Res[[i]]$or
  Res.cc[1,i] <- Res[[i]]$cc
  Res.ds[,i] <- Res[[i]]$ds
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
L_d_Res_ds <- paste("L", as.character(L),  "d", as.character(d), "Res_ds", sep = "_", collapse = NULL)
assign(L_d_Res_ds, Res.ds)

## print results
Res.bias.or <- mean(Res.or)-mu
Res.sd.or <- sqrt(mean((Res.or-mean(Res.or))^2))
Res.bias.cc <- mean(Res.cc)-mu
Res.sd.cc <- sqrt(mean((Res.cc-mean(Res.cc))^2))

Res.bias.ds <- c()
Res.sd.ds <- c()
for(j in 1:length(hc_seq)){
  Res.bias.ds[j] <- mean(abn_rmv(Res.ds[j,]))-mu
  Res.sd.ds[j] <- sqrt(mean((abn_rmv(Res.ds[j,])-mean(abn_rmv(Res.ds[j,])))^2))
}

or_bias_sd <- list(OR.bias.sd = round(c(Res.bias.or,Res.sd.or),4))
cc_bias_sd <- list(CC.bias.sd = round(c(Res.bias.cc,Res.sd.cc),4))
ds_bias_sd <- list(DS.bias = round(Res.bias.ds,4),
                   DS.sd = round(Res.sd.ds,4))

print(or_bias_sd)
print(cc_bias_sd)
print(ds_bias_sd)

print(proc.time()-tmp)



