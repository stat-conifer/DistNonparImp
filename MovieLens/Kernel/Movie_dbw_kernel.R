
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
load("MTransData.RData")
Feature <- MTransData[, -c(1, 2)]
MTransData <- cbind(MTransData[, c(1, 2)], Feature[, diag(var(Feature)) > 0.075])

MTransData[MTransData[, 2] == 0, 1] <- 0 
d <- length(MTransData[1, ]) - 2
N <- length(MTransData[, 1])

X <- MTransData[, -c(1,2)]
delta <- MTransData[, 2]
Y <- MTransData[, 1]
#########################################################################

hc_seq <- seq(0.1,2.1,0.4)

L <- c(100)
n <- floor(N/L)


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

print(ds)

print(proc.time()-tmp)



