
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
######################################################################

cc <- mean(Y[delta==1])
######################################################################

L_seq <- c(10,20,50,100,200,500)

ds_L <- rep(NA,length(L_seq))
sg_L <- rep(NA,length(L_seq))

time_ds_L <- rep(NA,length(L_seq))

for(j in 1:length(L_seq)){
  
  L <- L_seq[j]
  n <- floor(N/L)
  ######################################################################    
  hc <- 1.3
  hN <- hc * ((log(L, 10) * (L/N))^(1/(1 + 2 * d))) # ((max(log(L, 10), 1) * (L/N))^(1/(1 + 2 * d))) 
  hn <- hc * ((1/n)^(1/(1 + 2 * d)))

  ######SGM 
  
  mark <- 1:n
  Y_1 <-Y[mark]
  dlt_1 <- delta[mark]
  X_1 <- as.matrix(X[mark,],n,d)

  sg_L[j] <- ker_mean(Y_1,dlt_1,X_1,hn)
  #####KDI    
  
  
  ds <- rep(NA,L)
  time1 <- proc.time()
  for(l in 1:L){
    #set.seed(j)
    mark <- ((l-1)*n+1):(l*n)
    
    Y_l <- Y[mark]
    dlt_l <- delta[mark]
    X_l <- as.matrix(X[mark,],n,d)
    
    ds[l] <- ker_mean(Y_l,dlt_l,X_l,hN)
    
  }
  
  ds_L[j] <- mean(ds)
  
  time_ds_L[j] <- (proc.time()[[3]]-time1[[3]])/L
  
  print(time_ds_L[j])
  print(ds_L[j])
}

print(time_ds_L)

res <- list(cc = cc,
            sg = sg_L,
            ds = ds_L)
print(res)
######whole data
hc <- 1.3
hN <- hc * ((1/N)^(1/(1 + 2 * d)))
time1 <- proc.time()
whole <- 4#ker_mean(Y,delta,X,hN)
time_whole<-proc.time()[[3]]-time1[[3]]
print(time_whole)

library(latex2exp)

plot(1:length(L_seq),abs(rep(cc,length(L_seq))-whole), ylim = c(0,0.15),xlab = TeX(''), ylab=TeX(''),
     lty = 1, pch=1, type = 'o',cex = 0.5,cex.axis = 1.2, cex.lab = 1.4)
par(new = TRUE)
plot(1:length(L_seq),abs(res$sg-whole), ylim = c(0,0.15),xlab = '', ylab='',
     lty = 3,pch=4,type = 'o',cex = 0.5,cex.axis = 1.2, cex.lab = 1.4)
par(new = TRUE)
plot(1:length(L_seq),abs(res$ds-whole), ylim = c(0,0.15),xlab = '', ylab='',
     lty = 2,pch=2,type = 'o',cex = 0.5,cex.axis = 1.2, cex.lab = 1.4)
legend('topleft',c(TeX('cc'),TeX('sg'),TeX('ds')),lty=c(1,3,2),pch=c(1,4,2),cex=1)
#title(main = "curse of dimension")


print(proc.time()-tmp)