# Clamix
# clustering symbolic objects described by discrete distributions
# V. Batagelj, S. Korenjak-Èerne, N. Kejžar
# University of Ljubljana
# july 2010
# -------------------------------------------------------------------
# parameters
#   nVar - number of variables in SO description
#   w    - weights in dissimilarity
#   so   - empty SO
# are provided as global variables

# http://www.math.tau.ac.il/~yekutiel/CDA/Thompson manual.pdf
# page 10 & 11
# http://books.google.com/books?id=hpEzw4T0sPUC&printsec=frontcover&dq=analysis+of+categorical+data
# page 20
# V. Batagelj / 28.7.2010
# ---------------------------------------------------------------------
acceptbin<-function(x,n,p){
# computes the acceptability of p when x is observed and X is Bin(n,p)
# adapted from Blaker (2000)
  p1 <- 1-pbinom(x-1,n,p)
  p2 <- pbinom(x,n,p)
  a1 <- p1 + pbinom(qbinom(p1,n,p)-1,n,p)
  a2 <- p2 + 1 - pbinom(qbinom(1-p2,n,p),n,p)
  min(a1,a2)
}
acceptinterval<-function(x, n, level, tolerance=1e-04){
# computes acceptability interval for p at 1-alpha equal to level
# (in (0,1)) when x is an observed value of X which is Bin(n,p)
# adapted from Blaker (2000)
  lower <- 0; upper <- 1
  if(x) {lower <- qbeta((1-level)/2,x,n-x+1)
    while(acceptbin(x,n,lower) < (1-level)) lower <- lower+tolerance
  }
  if(x!=n) {upper <- qbeta(1-(1-level)/2,x+1,n-x)
    while(acceptbin(x,n,upper) < (1-level)) upper <- upper-tolerance
  }
  c(lower=lower,upper=upper)
}
# acceptinterval(x=0,n=25,level=.95)
# lower upper
# 0 0.1275852

# -----------------------------------------------------------------------

# creates an empty SO of type "symObject"
empty.symObject <- function(nCats){
  nVar <- length(nCats)
  s <- vector("list", nVar)
  for(i in 1:nVar) s[[i]] <- double(nCats[i]+1)
  class(s) <- "symObject"
  return(s)
}

# computes weighted squared Euclidean dissimilarity between SOs
distSO <- function(X,Y,nVar,alpha){
  d <- numeric(nVar)
  for(i in 1:nVar) {
    ln <- length(X[[i]]);
    nX <- as.double(X[[i]][[ln]]); nY <- as.double(Y[[i]][[ln]])
    d[[i]] <- nX*sum((X[[i]][-ln]/nX-Y[[i]][-ln]/nY)**2)
  }
  dis <- as.numeric(d %*% alpha)/2
  if (is.na(dis)) dis <- Inf
  return(dis)
}

# computes dissimilarity between clusters
distCl <- function(X,Y,nVar,alpha){
  d <- numeric(nVar)
  for(i in 1:nVar) {
    ln <- length(X[[i]]);
    nX <- as.double(X[[i]][[ln]]); nY <- as.double(Y[[i]][[ln]])
    d[[i]] <- nX*nY/(nX+nY)*sum((X[[i]][-ln]/nX-Y[[i]][-ln]/nY)**2)
  }
  return(as.numeric(d %*% alpha)/2)
}

# encode numerical vector using given encoding
encodeSO <- function(x,encoding,codeNA){
  if(is.na(x)) return(codeNA)
  for(i in 1:length(encoding)) if(encoding[[i]](x)) return(i)
}

# adapted leaders method
#---------------------------------------------
# VB, 16. julij 2010
#   dodaj omejitev na najmanjše število enot v skupini
#   omejeni polmer

leaderSO <- function(dataset,maxL){
  #attach(dataset)
  so <- dataset$so; alpha <- dataset$alpha
  SOs <- dataset$SOs; nVar <- length(so)
  numSO <- length(SOs)
  L <- vector("list",maxL); Ro <- numeric(maxL)
# random partition into maxL clusters
  clust <- sample(1:maxL,numSO,replace=TRUE)
  tim <- 1; step <- 0
  repeat {
    step <- step+1
  # new leaders - determine the leaders of clusters in current partition
    for(k in 1:maxL){L[[k]] <- so; names(L)[[k]] <- paste("L",k,sep="")}
    for(i in 1:numSO){j <- clust[[i]];
      for(k in 1:nVar) L[[j]][[k]] <- L[[j]][[k]] + SOs[[i]][[k]] }
  # new partition - assign each unit to the nearest new leader
    clust <- integer(numSO)
    R <- numeric(maxL); p <- double(maxL)
    for(i in 1:numSO){d <- double(maxL)
      for(k in 1:maxL){d[[k]] <- distSO(SOs[[i]],L[[k]],nVar,alpha)}
      r <- min(d); j <- which(d==r)
      if(length(j)==0){
        cat("unit",i,"\n",d,"\n"); flush.console(); print(SOs[[i]]); flush.console()
        u <- which(is.na(d))[[1]]; cat("leader",u,"\n"); print(L[[u]])
        stop()}
      j <- which(d==r)[[1]];
      clust[[i]] <- j
      p[[j]] <- p[[j]] + r; if(R[[j]]<r) R[[j]]<- r
    }
  # report
    cat("\nStep",step,"\n")
    print(table(clust)); print(R); print(Ro-R); Ro <- R;
    print(p); print(sum(p)); flush.console()
    tim <- tim-1
    if(tim<1){
      ans <- readline("Times repeat = ")
      tim <- as.integer(ans); if (tim < 1)  break
    }
  # in the case of empty cluster put in the most distant SO
    repeat{
      t <- table(clust); em <- setdiff(1:maxL,as.integer(names(t)))
      if(length(em)==0) break
      j <- em[[1]]; rmax <- 0; imax <- 0
      for(i in 1:numSO){d <- double(maxL)
        for(k in 1:maxL){d[[k]] <- distSO(SOs[[i]],L[[k]],nVar,alpha)}
        r <- max(d);  if(rmax<r) {rmax <- r; imax <- i}
      }
      clust[[imax]] <- j; L[[j]] <- SOs[[imax]]
      cat("*** empty cluster",j,"- SO",imax,"transfered, rmax =",rmax,"\n")
      flush.console()
    }
  }
  #detach(dataset)
  leaders <- dataset
  leaders$SOs <- L
  return (list(clust=clust,leaders_symData=leaders,R=R,p=p))
}

# adapted Ward's hierarchical clustering
#---------------------------------------------
# VB, 16. julij 2010
#   imena notranjih toèk v drevesu

hclustSO <- function(dataset){
# compute order of dendrogram
  orDendro <- function(i){if(i<0) return(-i)
    return(c(orDendro(m[i,1]),orDendro(m[i,2])))}
  #attach(dataset)
  so <- dataset$so; alpha <- dataset$alpha
  L <- dataset$SOs; nVar <- length(so)
  numL <- length(L); numLm <- numL-1
# each unit is a cluster; compute inter-cluster dissimilarity matrix
  D <- matrix(nrow=numL,ncol=numL); diag(D) <- Inf
  for(i in 1:numLm) for(j in (i+1):numL) {
      D[i,j] <- distCl(L[[i]],L[[j]],nVar,alpha); D[j,i] <- D[i,j]
  }
  active <- 1:numL; m <- matrix(nrow=numLm,ncol=2)
  node <- rep(0,numL); h <- numeric(numLm); LL <- vector("list",numLm)
  for(k in 1:numLm) LL[[k]] <- so
  for(k in 1:numLm){
  # determine the closest pair of clusters (p,q)
    ind <- active[sapply(active,function(i) which.min(D[i,active]))]
    dd <- sapply(active,function(i) min(D[i,active]))
    pq <- which.min(dd)
  # join the closest pair of clusters
    p<-active[pq]; q <- ind[pq]; h[k] <- D[p,q]
    if(node[p]==0){m[k,1] <- -p; Lp <- L[[p]]
      } else {m[k,1] <- node[p]; Lp <- LL[[node[p]]]}
    if(node[q]==0){m[k,2] <- -q; Lq <- L[[q]]
      } else {m[k,2] <- node[q]; Lq <- LL[[node[q]]]}
    for(t in 1:nVar) LL[[k]][[t]] <- Lp[[t]] + Lq[[t]]
    active <- setdiff(active,p)
    Lpq <- LL[[k]]
  # determine dissimilarities to the new cluster
    for(s in setdiff(active,q)){
      if(node[s]==0){Ls <- L[[s]]} else {Ls <- LL[[node[s]]]}
      D[q,s] <- distCl(Lpq,Ls,nVar,alpha); D[s,q] <- D[q,s]
    }
    node[[q]] <- k
  }
  hc <- list(merge=m,height=h,order=orDendro(numLm),labels=names(L),
    method="adapted ward",call=match.call(),dist.method="squared euclidean",leaders=LL)
  class(hc) <- "hclust"
  #detach(dataset)
  return(hc)
}

# prints values that are significantely different from the average
testSOvar <- function(SO,var,total,a){
  ln <- length(total[[var]]); lnM <- ln-1
  pt <- total[[var]][-ln]/total[[var]][[ln]]
  pj <- SO[[var]][-ln]/SO[[var]][[ln]]
  for(k in 1:lnM) {
    ai <- acceptinterval(x=total[[var]][[k]], n=total[[var]][[ln]], level=1-a)
    if(pj[[k]] > ai[[2]]) cat(format(names(total[[var]])[[k]],width=15,justify="right"),
      format(c(pj[[k]],ai[[2]],pt[[k]],ai[[1]],(pj[[k]]-pt[[k]])/(ai[[2]]-pt[[k]])),
      width=12,justify="left",digits=5,nsmall=7),"\n")
    if(pj[[k]] < ai[[1]]) cat(format(names(total[[var]])[[k]],width=15,justify="right"),
      format(c(pj[[k]],ai[[2]],pt[[k]],ai[[1]],(pt[[k]]-pj[[k]])/(pt[[k]]-ai[[1]])),
      width=12,justify="left",digits=5,nsmall=7),"\n")
  }
}

# prints values that are significantely larger than the average
testSOvarP <- function(SO,var,total,a){
  ln <- length(total[[var]]); lnM <- ln-1
  pt <- total[[var]][-ln]/total[[var]][[ln]]
  pj <- SO[[var]][-ln]/SO[[var]][[ln]]
  for(k in 1:lnM) {
    ai <- acceptinterval(x=total[[var]][[k]], n=total[[var]][[ln]], level=1-a)
    if(pj[[k]] > ai[[2]]) cat(format(names(total[[var]])[[k]],width=15,justify="right"),
      format(c(pj[[k]],ai[[2]],pt[[k]],ai[[1]],(pj[[k]]-pt[[k]])/(ai[[2]]-pt[[k]])),
      width=12,justify="left",digits=5,nsmall=7),"\n")
  }
}

# prints values that are significantely different from the average
testSOvarQ <- function(SO,var,total,a,qmin){
  ln <- length(total[[var]]); lnM <- ln-1
  pt <- total[[var]][-ln]/total[[var]][[ln]]
  pj <- SO[[var]][-ln]/SO[[var]][[ln]]
  for(k in 1:lnM) {
    ai <- acceptinterval(x=total[[var]][[k]], n=total[[var]][[ln]], level=1-a)
    if(pj[[k]] > ai[[2]]) {q <- (pj[[k]]-pt[[k]])/(ai[[2]]-pt[[k]])
      if(q>qmin) cat(format(names(total[[var]])[[k]],width=15,justify="right"),
      format(c(pj[[k]],ai[[2]],pt[[k]],ai[[1]],q),
      width=12,justify="left",digits=5,nsmall=7),"\n")}
    if(pj[[k]] < ai[[1]]) {q <- (pt[[k]]-pj[[k]])/(pt[[k]]-ai[[1]])
      if(q>qmin) cat(format(names(total[[var]])[[k]],width=15,justify="right"),
      format(c(pj[[k]],ai[[2]],pt[[k]],ai[[1]],q),
      width=12,justify="left",digits=5,nsmall=7),"\n")}
  }
}

# prints "equiprobabilistic" encodings for numerical variables
makeEnc <- function(var,name,k,file=""){
  n <- length(var)-length(which(is.na(var)))
  v <- sort(var)
  steps <- seq.int(0,n,n %/% (k-1))
  cuts <- c(0,v[steps[1:(k-1)]],max(v,na.rm=TRUE))
  cat("enc",name," <- list(\n",sep="",file=file,append=TRUE)
  cat('  "[0]" = function(x) x<=0,\n',file=file,append=TRUE)
  for(j in 2:k) cat('  "(',cuts[j-1],',',cuts[j],']" = function(x) x<=',cuts[j]+1e-7,',\n',sep='',file=file,append=TRUE)
  cat('  "NA" = function(x) TRUE )\n',file=file,append=TRUE)
}

# computes totals
computeTotal <- function(dataset){
  #attach(dataset)
  so <- dataset$so; namedSO <- dataset$namedSO
  L <- dataset$SOs; nVar <- length(so)
  numL <- length(L); total <- namedSO 
  for(i in 1:numL) for(j in 1:nVar) total[[j]] <- total[[j]] + L[[i]][[j]]
  #detach(dataset)
  return(total)
}

print.symObject <- function(x,...){
  # print symObject
  cat("symObject",":","\n",sep=" ")
  cat("Number of variables:", length(x) ,"\n",sep=" ")
  for (i in 1:length(x)){
    cat(names(x)[i],"\n",sep=" ")
    print(x[[i]])
  }
}

summary.symData <- function(x,...){
  # summary of symData - vector of symObjects
  cat("summary symData",":","\n\n",sep=" ")
  cat("Dimension (units x variables):", length(x$SOs),"x",length(x$so),"\n",sep=" ")
  cat("Type of distributions:", x$type,"\n", sep=" ")
  cat("Outlook of symObject: ","\n")
  print(x$namedSO)
}

create.symData <- function(datalist,type="gDist",alpha=NULL){
  # TYPE - PDIST - SUM!?!?
  # creates object symData
  # from list of dataframes
  nVar <- length(datalist)
  if (is.null(alpha))
    alpha <- rep(1/nVar,nVar)
  n_units <- nrow(datalist[[1]])
  nVar_names <- names(datalist)
  unit_names <- rownames(datalist[[1]])
 for (i in 1:nVar){
    temp <- nrow(datalist[[i]])
    if (temp != n_units)
      stop("number of units must match!")
  }
  dataset <- vector(mode="list",length=nVar)
  nCats <- vector(mode="numeric",length=nVar)
  for (i in 1:nVar){
    dataset[[i]] <- cbind(datalist[[i]],0,apply(datalist[[i]],1,sum))
    coldata <- ncol(dataset[[i]])
    nCats[i] <- coldata-1
    colnames(dataset[[i]])[[coldata-1]] <- "NA"; colnames(dataset[[i]])[[coldata]] <- "num"
  }
  SOs <- vector(mode="list",length=n_units)
  for (i in 1:n_units){
    s <- vector("list", nVar)
    for (j in 1:nVar){
      s[[j]] <- as.double(dataset[[j]][i,])
    }
    names(s) <- nVar_names
    class(s) <- "symObject"
    SOs[[i]] <- s
 }
  names(SOs) <- unit_names
  so <- empty.symObject(nCats)
  names(so) <- nVar_names
  namedSO <- so
  for (i in 1:nVar)
    names(namedSO[[i]]) <- colnames(dataset[[i]])

  # make an object
  object <- list(SOs=SOs,so=so,namedSO=namedSO,alpha=alpha,type=type)
  class(object) <- "symData"
  return(object)
}

