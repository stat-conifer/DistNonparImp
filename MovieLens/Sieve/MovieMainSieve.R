setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
load("MTransData.RData")
source("Functions.R")
Feature <- MTransData[, -c(1, 2)]
MTransData <- cbind(MTransData[, c(1, 2)], Feature[, diag(var(Feature)) > 0.075])


MTransData[MTransData[, 2] == 0, 1] <- 0 #missing data replaced by zero
L_seq <- c(10, 20, 50, 100, 200, 500) #numbers of machines
d <- length(MTransData[1, ]) - 2
N <- length(MTransData[, 1])



delta <- MTransData[, 2]
Y <- MTransData[, 1]

cc <- mean(Y[delta == 1])

t1 <- proc.time()

KN <- ceiling(0.9 * N^{d/(1 + 2*d)})

X_pl <- powers(MTransData[, -c(1,2)], KN)

Y_obs <- Y[delta==1]
X_pl_obs <- X_pl[delta==1, ]

hessian <- t(X_pl_obs) %*% X_pl_obs / N
inv <- solve(hessian + 1e-7 * diag(dim(hessian)[1])) 
rm(hessian)
beta <- inv %*% t(X_pl_obs) %*% Y_obs / N #compute estimate from a single machine

rm(inv)

X_pl_mis_mean <- apply((1 - delta) * X_pl[1:N, ], 2, mean)
all <- mean(delta * Y) + X_pl_mis_mean %*% beta

all.time <- (proc.time()-t1)[[3]]

ds.time <- rep(NA, length(L_seq))

ds <- c()
sg <- c()

for(j in 1:length(L_seq)){
  L <- L_seq[j]
  n <- floor(N/L)
  
  mark <- 1:n
  
  Kn <- ceiling(0.9 * n^{d/(1 + 2*d)})
  
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
  sg[j] <- mean(dlt_1 * Y_1) + X_pl_mis_mean_1 %*% beta_sg
  
  t1 <- proc.time()
  
  
  Y_1 <- Y[mark]
  dlt_1 <- delta[mark]
  X_pl_1 <- powers(MTransData[mark, -c(1,2)], KN)
  Y_1_obs <- Y_1[dlt_1==1]
  X_pl_1_obs <- X_pl_1[dlt_1==1, ]
  length(X_pl_1_obs[, 1])
  alpha <- 0.015 * log(KN, 10)^2 * sqrt(KN / n)
  hessian_ds <- t(X_pl_1_obs) %*% X_pl_1_obs / n
  inv_ds <- solve(hessian_ds + alpha * diag(dim(hessian_ds)[1]))
  rm(hessian_ds)
  beta_ds <- inv_ds %*% t(X_pl_1_obs) %*% Y_1_obs / n
  
  iter <- max(floor(30 * max(log(L, 10), 1) * (log(N, 10) + log(KN, 10)) / log(max(1 / alpha, log(KN, 10)), 10)),2)
  
  t_par <- c()

  for (t in 1:iter) {
    grad <- - t(X_pl_1_obs) %*% (Y_1_obs - X_pl_1_obs %*% beta_ds)
    
    t_tmp <- proc.time()
    
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
    
    t_par[t] <- (proc.time() - t_tmp)[[3]]
    
    beta_ds <- beta_ds - inv_ds %*% grad / N
  }
  
  X_pl_mis_mean <- apply((1 - delta) * X_pl, 2, mean)
  ds[j] <- mean(delta * Y) + X_pl_mis_mean %*% beta_ds
  
  ds.time[j] <- (proc.time()-t1)[[3]] - sum(t_par)
}




ResRD <- as.data.frame(t(c(all, cc, sg, ds)))
names(ResRD) <- c("all", "cc","sg10", "sg20", "sg50", "sg100", "sg200", "sg500", "ds10", "ds20", "ds50", "ds100", "ds200", "ds500")
cent <- ResRD$all
res <- as.vector(as.matrix(ResRD[, -1]))
z <- abs(unlist(res) - cent)
plot(rep(z[1], 6), type = "l", lty = 2, col = 2, ylim = c(0, 0.35), xaxt = "n",  ylab = "absolute difference", lwd = 1.8)
axis(1, at = 1:6, labels = L_seq, xlab = "L", lab.font = 3)
lines(z[2:7], lty = 3, col = 3, lwd = 1.8)
points(z[2:7], pch = 17, col = 3)
lines(z[8:13], lty = 3, col = 4, lwd = 1.8)
points(z[8:13], pch = 15, col = 4)

