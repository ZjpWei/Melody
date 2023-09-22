GDM.reg <- function(Y, X, max.iter){
  rm.ind = which(colSums(Y) == 0)
  if(length(rm.ind) > 0)
    Y = Y[,-rm.ind,drop=FALSE]
  keep.ind = which(rowSums(Y) != 0)
  Y = Y[keep.ind,,drop=FALSE]
  X = X[keep.ind,,drop=FALSE]
  avg.rel = colMeans(Y/rowSums(Y))
  Y = Y[,order(avg.rel, decreasing = TRUE)]
  K = ncol(Y) - 1
  est.para.lst = lapply(1:K, function(x){
    Y.tmp = cbind(Y[,x], rowSums(Y[,(x+1):(K+1),drop=FALSE]))
    # id.sel = rowSums(Y.tmp) != 0
    # Y.tmp = Y.tmp[id.sel,,drop=FALSE]
    # X.tmp = X[id.sel,,drop=FALSE]
    alpha0 = matrix(0.001, nrow = ncol(X), ncol = 1); phi0 = 1
    .GDM_EM(Y.tmp, X, alpha0, phi0, max.iter = max.iter)
  })
  est.para = do.call(Map, c(f = cbind, est.para.lst))
  return(est.para)
}

r.GDM <- function(X, mod, SeqDepth = 1e3){
  alpha = mod$alpha.est
  phi = mod$phi.est
  N = nrow(X); K = ncol(alpha)
  if(length(SeqDepth) == 1){
    SeqDepth = rep(SeqDepth, N)
  }else{
    if(length(SeqDepth) != N) stop("The length of Sequencing depth vector does not match the length of covariates!")
  }
  mu.mat = 1/(1+exp(-(X %*% alpha)))
  phi.mat = matrix(phi, N, K, byrow = TRUE)
  a.mat = mu.mat * phi.mat
  b.mat = (1 - mu.mat) * phi.mat
  
  Z = matrix(rbeta(N*K, a.mat, b.mat), N, K)
  P = matrix(NA, N, K+1)
  for (j in 1:K) {
    P[,j] = Z[,j]*matrixStats::rowProds(1 - Z[,0:(j-1), drop = FALSE])
  }
  P[, ncol(P)] =  1 - rowSums(P[, -ncol(P), drop = FALSE])
  P[P<0] = 0
  # Y = c()
  # for (i in 1:N) {
  #   Y = rbind(Y, t(rmultinom(1, SeqDepth[i], P[i,])))
  # }
  # return(Y)
  return(P)
}

.GDM_EM <- function(Y, X, alpha0, phi0, tol = 0.0001, max.iter = 1000){
  CONV = 0
  CONV.iter = max.iter
  n = nrow(Y)
  K = ncol(Y) - 1
  d = ncol(X)
  
  alpha.last = alpha0; phi.last = phi0
  alpha.now = alpha0; phi.now = phi0
  
  AB.R = matrix(NA, nrow = n, ncol = 2*K)
  for (l in 1:max.iter) {
    # E-step
    # print(paste("====== ", l, "th ======", sep=""))
    
    for(i in 1:n){
      tmp = exp(X[i,] %*% alpha.last)
      mv = as.numeric(tmp/(1+tmp))
      sv = phi.last
      
      av = mv * sv
      bv = (1 - mv) * sv
      
      par.post = .eZ(av, bv, Y[i,])
      tmp = par.post$av.post + par.post$bv.post
      AB.R[i,] = c(digamma(par.post$av.post) - digamma(tmp), digamma(par.post$bv.post) - digamma(tmp))
    }
    A.R = AB.R[, 1:K, drop = FALSE]
    B.R = AB.R[, K+1:K, drop = FALSE]
    
    # M-step
    for (j in 1:K) {
      tmp = .BetaOptim(A.R[,j], B.R[,j], X, alpha.last[,j], phi.last[j])
      alpha.now[,j ] = tmp[1:d]; phi.now[j] = tmp[-(1:d)]
    }
    
    diff = 0
    diff = diff + sum(abs(alpha.now - alpha.last))
    diff = diff + sum(abs(phi.now - phi.last))
    if(diff < tol){
      CONV = 1; CONV.iter = l
      break
    }else{
      alpha.last = alpha.now; phi.last = phi.now
    }
  }
  return(list(alpha.est = alpha.now, phi.est = phi.now, 
              CONV = CONV, CONV.iter = CONV.iter))
}

# .eZ: expectation of Z in the GD
.eZ <- function(av, bv, Y){
  K = length(Y)-1
  N = sum(Y)
  
  av.prim = av + Y[1:K]  
  bv.prim = bv
  for(j in 1:K){
    bv.prim[j] = bv.prim[j] + (N - sum(Y[1:j]))
  }
  return(list(av.post = av.prim, bv.post = bv.prim))
}
.BetaOptim <- function(A, B, X, alpha.ini, phi.ini){
  # A = A.R[,j]; B = B.R[,j]; alpha.ini = alpha.last[,j]; phi.ini = phi.last[j]
  Beta.par.ini = c(alpha.ini, phi.ini)
  Beta.data = list(A = A, B = B, X = X)
  return(optim(par = Beta.par.ini, fn = .BetaNegLoglik, gr = .BetaNegScore, data = Beta.data, method = "BFGS")$par)
}
# negative log pseudo-likelihood function of the GDM model
# minimizing logL_gdm is equivalent to maximizing the pseudo-likelihood
.BetaNegLoglik <- function(par, data){
  # par = Beta.par.ini; data = Beta.data
  d = ncol(data$X)
  alpha = par[1:d]
  phi = par[-(1:d)]
  
  tmp = as.numeric(exp(data$X %*% alpha))
  mu.tmp = as.numeric( tmp/(1+tmp) )
  
  a = phi * mu.tmp
  b = phi * (1 - mu.tmp)
  a[a < 0] = 0
  b[b < 0] = 0
  
  return(-sum( -lbeta(a, b) + data$A * (a - 1) + data$B * (b - 1)))
}

# score function, obtained by differentiating the log-likelihood function w.r.t unknown parameters
.BetaNegScore <- function(par, data){
  d = ncol(data$X)
  alpha = par[1:d]
  phi = par[-(1:d)]
  
  tmp = as.numeric(exp(data$X %*% alpha))
  mu.tmp = as.numeric( tmp/(1+tmp) )
  
  a = phi * mu.tmp
  b = phi * (1 - mu.tmp)
  a[a < 0] = 0
  b[b < 0] = 0
  
  zstar = data$A - data$B
  mustar = digamma(a) - digamma(b)
  
  return(-c(colSums(phi * mu.tmp * (1/(1+tmp)) * (zstar - mustar) * data$X),
            sum(mu.tmp * (zstar - mustar) + data$B - digamma(b) + digamma(phi))))
}
