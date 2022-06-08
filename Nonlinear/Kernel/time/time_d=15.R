
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

time_fun <- function(t){
  set.seed(t)
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
  hc <- 1.7
  hN_all <- hc * ((1/N)^(1/(1 + 2 * d)))

  ########################################################################    
  
  t1 <- proc.time()
  
  all <- ker_mean(Y,delta,X,hN_all)
  
  all.time <- (proc.time()-t1)[[3]]
  
  ds.time <- c()
  for(j in 1:length(L_seq)){
    L <- L_seq[j]
    n <- floor(N/L)
    
    hN_ds <- hc * ((log(L, 10) * (L/N))^(1/(1 + 2 * d)))
    
    
    ds <- rep(NA,L)
    
    ti <- proc.time()
    t1 <- proc.time()
    for(l in 1:L){
      #set.seed(l)
      
      mark <- ((l-1)*n+1):(l*n)
      
      Y.l <- Y[mark]
      dlt.l <- delta[mark]
      X.l <- as.matrix(X[mark,],n,d)
      
      ds[l] <- ker_mean(Y.l,dlt.l,X.l,hN_ds)
      
    }
    t2 <- proc.time()
    ds = mean(ds)
    tf = proc.time()
    
    ds.time[j] <- (tf-ti-(t2-t1)*((L-1)/L))[3]
  }
  
  
  res <- list(all.time = all.time,
              ds.time = ds.time)
  return(res)
}


mu <- 2

N <- 200000
d <- 15
simul <- 20
L_seq <- c(500)

cl <<- parallel::makeCluster(20, setup_strategy = "sequential")
parallel::clusterExport(cl, c("N","v","L_seq","hok","d"))
parallel::clusterEvalQ(cl,source("Functions.R"))
Res <- parallel::parLapply(cl, 1:simul, time_fun)
parallel::stopCluster(cl)

Res.time.all <- rep(NA, simul)
Res.time.ds <- matrix(NA, length(L_seq), simul)

for(i in  1:simul){
  Res.time.all[i] <- Res[[i]]$all.time
  Res.time.ds[,i] <- Res[[i]]$ds.time
}

# ## save results
# d_time_all <- paste("d", as.character(d), "time_all", sep = "_", collapse = NULL)
# assign(d_time_all, Res.time.all)
# 
# d_time_ds <- paste("d", as.character(d), "time_ds", sep = "_", collapse = NULL)
# assign(d_time_ds, Res.time.ds)


all.time <- mean(Res.time.all)
ds.time <- rowMeans(Res.time.ds)


print(all.time)

print(ds.time)

print(proc.time()-tmp)
