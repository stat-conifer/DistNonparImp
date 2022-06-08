
## kernel function

ker_fun <- function(x){
  P1 <- 1
  P2 <- x * (abs(x) <= 1)
  Kx <- v[1] * P1 
  for (j in 2:(hok/2)) {
    qt <- 2*(j-1)
    P1 <-  (2*qt - 1) / qt * x * P2 - (qt - 1) / qt * P1
    Kx <- Kx + v[j] * P1
    P2 <- (2*qt + 1) / (qt + 1) * x * P1 - qt / (qt + 1) * P2
  }
  Kx * (abs(x) <= 1) / 2
}


abn_rmv <- function(x){
  upper <- quantile(x, 0.75)
  lower <- quantile(x,0.25)
  x[(x <= upper + 1.5*(upper - lower)) & (x >= lower - 1.5*(upper - lower))]
}

ker_mean <- function(y,dlt,x,bw){
  
  x_mis = x[dlt==0,]
  x_obs = x[dlt==1,]
  y_obs <- y[dlt == 1]
  
  
  im_fun <- function(z){
    
    zx = t((z-t(x_obs))/(resample::colStdevs(x)*bw))
    
    
    ker_zx = ker_fun(zx)
    
    I <- rowSums(ker_zx==0)==0
    
    if (sum(I) == 0) {
      return(0)
    }
    
    M_tmp <- matrix(ker_zx[I, ], sum(I))
    ker_zx_prod <- (1 - 2 * rowSums(M_tmp < 0) %%2) * exp(rowSums(log(abs(M_tmp))))
    
    
    deno =  sum((1/(length(y)*bw^(ncol(x))))*ker_zx_prod)#rep(1,nt)
    nume = sum((1/(length(y)*bw^(ncol(x))))*ker_zx_prod*y_obs[I])
    
    
    sign = abs(deno) < 1e-120
    
    deno[sign] = 1e-120
    
    return(nume/deno)
  }
  
  # tmp = proc.time()
  im = apply(x_mis,1,im_fun)
  # print(proc.time()-tmp)
  
  mu_est = mean(y[dlt==1]) * mean(dlt) + mean(abn_rmv(im)) * mean(1 - dlt)
  return(mu_est)
}

# ker_fun_lower <- function(x){
#   (abs(x) <= 1) / 2
# }
# 
# ker_pi <- function(dlt, x, bw){
#   
#   im_fun <- function(z){
#     
#     zx = t((z-t(x))/(resample::colStdevs(x)*bw))
#     
#     
#     ker_zx = ker_fun_lower(zx)
#     
#     
#     # kernel.t = apply(vec.kernel.t,1,prod)
#     
#     I <- rowSums(ker_zx==0)==0
#     
#     if (sum(I) == 0) {
#       return(0)
#     }
#     
#     M_tmp <- matrix(ker_zx[I, ], sum(I))
#     ker_zx_prod<- (1 - 2 * rowSums(M_tmp < 0) %%2) * exp(rowSums(log(abs(M_tmp))))
#     
#     
#     deno =  sum((1/(length(dlt)*bw^(ncol(x))))*ker_zx_prod)#rep(1,nt)
#     nume = sum((1/(length(dlt)*bw^(ncol(x))))*ker_zx_prod*dlt[I])
#     
#     
#     sign = abs(deno) < 1e-120
#     
#     deno[sign] = 1e-120
#     
#     return(nume/deno)
#   }
#   
#   # tmp = proc.time()
#   im = apply(x,1,im_fun)
#   # print(proc.time()-tmp)
#   
#   return(im)
# }
# 
# 
# CV_pi <- function(y_tr, y_te, x_tr, x_te, bw){
#   
#   im_fun <- function(z){
#     
#     zx = t((z-t(x_tr))/(resample::colStdevs(x_tr)*bw))
#     
#     ker_zx = ker_fun_lower(zx)
#     
#     # kernel.t = apply(vec.kernel.t,1,prod)
#     
#     I <- rowSums(ker_zx==0)==0
#     
#     if (sum(I) == 0) {
#       return(0)
#     }
#     
#     M_tmp <- matrix(ker_zx[I, ], sum(I))
#     ker_zx_prod <- (1 - 2 * rowSums(M_tmp < 0) %%2) * exp(rowSums(log(abs(M_tmp))))
#     
#     deno =  sum((1/(length(y_tr)*bw^(ncol(x_tr))))*ker_zx_prod)#rep(1,nt)
#     nume = sum((1/(length(y_tr)*bw^(ncol(x_tr))))*ker_zx_prod*y_tr[I])
#     
#     
#     sign = abs(deno)<1e-120
#     
#     deno[sign] = 1e-120
#     
#     return(nume/deno)
#   }
#   
#   im_y_te = apply(x_te,1,im_fun)
#   
#   
#   return(mean(((y_te-im_y_te))^2))
#   
# }
# 
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

CV_m <- function(y_tr_obs, y_te_obs, x_tr_obs, x_te_obs, w_te, bw){
  
  im_fun <- function(z){
    
    zx = t((z-t(x_tr_obs))/(resample::colStdevs(x_tr_obs)*bw))
    
    ker_zx = ker_fun(zx)
    
    # kernel.t = apply(vec.kernel.t,1,prod)
    
    I <- rowSums(ker_zx==0)==0
    
    if (sum(I) == 0) {
      return(0)
    }
    
    M_tmp <- matrix(ker_zx[I, ], sum(I))
    ker_zx_prod <- (1 - 2 * rowSums(M_tmp < 0) %%2) * exp(rowSums(log(abs(M_tmp))))
    
    deno =  sum((1/(length(y_tr_obs)*bw^(ncol(x_tr_obs))))*ker_zx_prod)#rep(1,nt)
    nume = sum((1/(length(y_tr_obs)*bw^(ncol(x_tr_obs))))*ker_zx_prod*y_tr_obs[I])
    
    
    sign = abs(deno)<1e-120
    
    deno[sign] = 1e-120
    
    return(nume/deno)
  }
  
  im_y_te = apply(x_te_obs,1,im_fun)
  
  upper <- quantile(y_te_obs-im_y_te, 0.75)
  lower <- quantile(y_te_obs-im_y_te, 0.25)
  id_norm <- (y_te_obs-im_y_te <= upper + 1.5*(upper - lower)) & (y_te_obs-im_y_te >= lower - 1.5*(upper - lower))
  
  return(mean(w_te[id_norm] * (abn_rmv(y_te_obs-im_y_te))^2,na.rm=TRUE))
  
}
