
powers <- function(x, k) {
  n <- length(x[, 1])
  d <- length(x[1, ])
  x <- apply(x, c(1, 2), truncov)
  x <- apply(x, 2, normcov) #if the return of apply function are vectors, the vectors will be arranged as a matrix
  tensor <- matrix(NA, n, k+1)
  for (i in 1:n) {
    v <- list()
    for (j in 1:d) {
      v[[j]] <- x[i, j]
    }
    a <- 1
    a <- c(a, unlist(v)[1:min(d, k)])
    l <- length(a)
    
      while (l < k + 1) {
       for (j in 1:d){
         l <- length(a)
         if (l >= k + 1) {
           break
         }
         b <- unlist(v[j:d])
         v[[j]] <- x[i, j] * b[1: min(length(b), k + 1 - l)]
         a <- c(a, unlist(v[[j]]))
       }
      }
    tensor[i, ] <- a
  }
  tensor
}

abn_rmv <- function(x){
  upper <- quantile(x, 0.75)
  lower <- quantile(x,0.25)
  x[(x <= upper + 1.5*(upper - lower)) & (x >= lower - 1.5*(upper - lower))]
}

truncov <- function(x)
{
  max(min(x, 0.4 * log(N)), -0.4 * log(N))
}

normcov <- function(x){(2*x - min(x) - max(x))/(max(x) - min(x))}


CV_m <- function(y_tr, y_te, x_tr, x_te, w_te, K){
  
  x_tr_pl <- powers(x_tr, K)
 
  hessian_tr <- t(x_tr_pl) %*% x_tr_pl / length(y_tr)
  inv_tr <- solve(hessian_tr)
  rm(hessian_tr)
  beta_tr <- inv_tr %*% t(x_tr_pl) %*% y_tr / length(y_tr)
  
  x_te_pl <- powers(x_te, K)
  
  im_y_te  <- x_te_pl %*% beta_tr
  
  # print(mean((abn_rmv((y_te-im_y_te)))^2))
  
  upper <- quantile(y_te-im_y_te, 0.75)
  lower <- quantile(y_te-im_y_te, 0.25)
  id_norm <- (y_te-im_y_te <= upper + 1.5*(upper - lower)) & (y_te-im_y_te >= lower - 1.5*(upper - lower))
  
  return(mean(w_te[id_norm] * (abn_rmv(y_te-im_y_te))^2,na.rm=TRUE))
  
}


moment <- function(x) {
  n <- length(x[, 1])
  d <- length(x[1, ])
  k <- choose(d + 2, 2) ##number of second order moment
  tensor <- matrix(NA, n, k+1)
  for (i in 1:n) {
    v <- list()
    for (j in 1:d) {
      v[[j]] <- x[i, j]
    }
    a <- 1
    a <- c(a, unlist(v)[1:min(d, k)])
    l <- length(a)
    
    while (l < k + 1) {
      for (j in 1:d){
        l <- length(a)
        if (l >= k + 1) {
          break
        }
        b <- unlist(v[j:d])
        v[[j]] <- x[i, j] * b[1: min(length(b), k + 1 - l)]
        a <- c(a, unlist(v[[j]]))
      }
    }
    tensor[i, ] <- a
  }
  tensor
}



