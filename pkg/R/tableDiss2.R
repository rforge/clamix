# -------------------------------------------------------------------------- #

check <- function(x,penalty=1e6){
  x[which(is.nan(x))] <- 0 # (0/0 becomes 0)
  x[which(is.infinite(x))] <- penalty # X/0 becomes big number
  x
}

# distSO and the functions for different dissimilarities

distSO <- function(X,Y,nVar,alpha,diss){
  d <- numeric(nVar)
  for(i in 1:nVar) {
    ln <- length(X[[i]]);
    nX <- as.double(X[[i]][[ln]]); nY <- as.double(Y[[i]][[ln]])
    d[[i]] <- diss(X[[i]][-ln],Y[[i]][-ln],nX,nY)
  }
  dis <- as.numeric(d %*% alpha)
  if (is.na(dis)) dis <- Inf
  return(dis)
}

deltas <- list(
  d1 = function(x,y,nx,ny){
    nx*sum((x-y)^2)  
  },
  d1w = function(x,y,nx,ny){
    nx*sum((x/nx-y/ny)^2)
  },
  d2 = function(x,y,nx,ny){
    temp = ((x-y)/y)^2
    temp = check(temp)
    nx*sum(temp)   
  },  
  d3 = function(x,y,nx,ny){
    temp = (x-y)^2/y
    temp = check(temp)
    nx*sum(temp)   
  },
  d4 = function(x,y,nx,ny){
    temp = ((x-y)/x)^2
    temp = check(temp)
    nx*sum(temp)   
  },
  d5 = function(x,y,nx,ny){
    temp = (x-y)^2/x
    temp = check(temp)
    nx*sum(temp)   
  },
  d6 = function(x,y,nx,ny){
    temp = (x-y)^2/(x*y)
    temp = check(temp)
    nx*sum(temp)   
  }  
)

computePQHG <- function(X,ivar){
  ln = length(X[[1]][[ivar]])
  P = Q = H = G = rep(0,(ln-1)); 
  w = 0
  for(j in 1:length(X)){
    P = P + X[[j]][[ivar]][-ln]*X[[j]][[ivar]][ln]
    Q = Q + X[[j]][[ivar]][-ln]^2*X[[j]][[ivar]][ln]
    H = H + check(X[[j]][[ivar]][ln]/X[[j]][[ivar]][-ln])
    G = G + check(X[[j]][[ivar]][ln]/X[[j]][[ivar]][-ln]^2)
    w = w + X[[j]][[ivar]][ln]
  }
  return(list(P=P,Q=Q,H=H,G=G,w=w))
}

# -------------------------------------------------------------------------- #
# computes distC and functions for different dissimilarities

distC <- function(X,Y,Z,nVar,alpha,SOx,SOy,diss){
  d <- numeric(nVar)
  for(i in 1:nVar) {
    ln <- length(X[[i]]);
    nX <- as.double(X[[i]][[ln]]); nY <- as.double(Y[[i]][[ln]])
    d[[i]] <- diss(X[[i]][-ln],Y[[i]][-ln],Z[[i]][-ln],SOx,SOy,i,nX,nY)
  }
  return(as.numeric(d %*% alpha))
}

Ds <- list(
  d1 = function(x,y,z,X,Y,ivar,nx,ny){
    nx*ny/(nx+ny)*sum((x-y)^2)
  },
  d1w = function(x,y,z,X,Y,ivar,nx,ny){
    nx*ny/(nx+ny)*sum((x/nx-y/ny)^2)
  },
  d2 = function(x,y,z,X,Y,ivar,nx,ny){
    temp1 <- ((x - z)/z)^2; temp1 = check(temp1)
    temp2 <- ((y - z)/z)^2; temp2 = check(temp2)
    u = computePQHG(X,ivar); v = computePQHG(Y,ivar);
    sum(check(u$P/x)*temp1) + sum(check(v$P/y)*temp2)
  },
  d3 = function(x,y,z,X,Y,ivar,nx,ny){
    temp1 <- (x - z)^2/z; temp1 = check(temp1)
    temp2 <- (y - z)^2/z; temp2 = check(temp2)
    nx*sum(temp1) + ny*sum(temp2)
  },
  d4 = function(x,y,z,X,Y,ivar,nx,ny){
    temp1 <- (x - z)^2
    temp2 <- (y - z)^2
    u = computePQHG(X,ivar); v = computePQHG(Y,ivar);
    sum(u$G*temp1) + sum(v$G*temp2)
  },
  d5 = function(x,y,z,X,Y,ivar,nx,ny){
    temp1 <- (x - z)^2/x; temp1 = check(temp1)
    temp2 <- (y - z)^2/y; temp2 = check(temp2)
    nx*sum(temp1) + ny*sum(temp2)
  },
  d6 = function(x,y,z,X,Y,ivar,nx,ny){
    temp1 <- ((x - z)/(x*z))^2; temp1 = check(temp1)
    temp2 <- ((y - z)/(y*z))^2; temp2 = check(temp2)
    u = computePQHG(X,ivar); v = computePQHG(Y,ivar);
    sum(check(u$P/x)*temp1) + sum(check(v$P/y)*temp2)
  }
)

# -------------------------------------------------------------------------- #
# computes t (leaders for leaderSO) and functions for different dissimilarities

tlead <- function(SOs,numSO,nVar,maxL,so,clust,ldr){
  L <- vector("list", maxL)
  for (k in 1:maxL) {
    L[[k]] <- so
    names(L)[[k]] <- paste("L", k, sep = "")
    X = SOs[k==clust]
    L[[k]] = ldr(L[[k]],X,nVar)
  }
  return(L)
}

ts <- list(
  d1 = function(L,X,nVar){
    for(i in 1:nVar){
      temp = computePQHG(X,i)
      L[[i]] = c(temp$P/temp$w,temp$w)}
    L
  },
  d2 = function(L,X,nVar){
    for(i in 1:nVar){
      temp = computePQHG(X,i)
      L[[i]] = c(temp$Q/temp$P,temp$w)}
    L
  },
  d1w = function(L,X,nVar){
    for(i in 1:nVar){
      for(j in 1:length(X))
        L[[i]] = L[[i]] + X[[j]][[i]]}
    L
  },
  d3 = function(L,X,nVar){
    for(i in 1:nVar){
      temp = computePQHG(X,i)
      L[[i]] = c(sqrt(temp$Q/temp$w),temp$w)}
    L
  },
  d4 = function(L,X,nVar){
    for(i in 1:nVar){
      temp = computePQHG(X,i)
      L[[i]] = c(check(temp$H/temp$G),temp$w)}
    L
  },
  d5 = function(L,X,nVar){
    for(i in 1:nVar){
      temp = computePQHG(X,i)
      L[[i]] = c(check(temp$w/temp$H),temp$w)}
    L
  },
  d6 = function(L,X,nVar){
    for(i in 1:nVar){
      temp = computePQHG(X,i)
      L[[i]] = c(sqrt(check(temp$P/temp$H)),temp$w)}
    L
  }
)

# -------------------------------------------------------------------------- #
# computes z (leader for hclustSO) and functions for different dissimilarities

zlead <- function(Lu,Lv,nVar,so,SOu,SOv,ldr){
  L <- so
  for (i in 1:nVar) {
    ln <- length(Lu[[i]])
    L[[i]] <- ldr(Lu[[i]][-ln],Lv[[i]][-ln],Lu[[i]][ln],Lv[[i]][ln],SOu,SOv,i)
  }
  return(L)
}

zs <- list(
  d1 = function(u,v,nu,nv,U,V,ivar){
    c((nu*u + nv*v)/(nu + nv),nu+nv)
  },
  d1w = function(u,v,nu,nv,U,V,ivar){
    c(u+v,nu+nv)
  },
  d2 = function(u,v,nu,nv,U,V,ivar){
    Pu = computePQHG(U,ivar)$P; Pv = computePQHG(V,ivar)$P;
    c(check((u*Pu + v*Pv)/(Pu + Pv)),nu+nv)
  },
  d3 = function(u,v,nu,nv,U,V,ivar){
    c(sqrt((nu*u^2 +nv*v^2)/(nu+nv)),nu+nv)
  },
  d4 = function(u,v,nu,nv,U,V,ivar){
    Hu = computePQHG(U,ivar)$H; Hv = computePQHG(V,ivar)$H;
    c(check((Hu + Hv)/(check(Hu/u) + check(Hv/v))),nu+nv)
  },
  d5 = function(u,v,nu,nv,U,V,ivar){
    Hu = computePQHG(U,ivar)$H; Hv = computePQHG(V,ivar)$H;
    c(check((nu + nv)/(Hu + Hv)),nu+nv)
  },
  d6 = function(u,v,nu,nv,U,V,ivar){
    Pu = computePQHG(U,ivar)$P; Pv = computePQHG(V,ivar)$P;
    c(sqrt(check((Pu + Pv)/(check(Pu/u^2) + check(Pv/v^2)))),nu+nv)
  }
)

