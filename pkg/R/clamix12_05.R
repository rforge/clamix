# Clamix
# clustering symbolic objects described by discrete distributions
# V. Batagelj, S. Korenjak-Cerne, N. Kejzar
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

## clamix12_01 - tableDiss2.R klicemo sami -- S 10. 7. 2016, dodane vse razlicnosti d1 - d6
## clamix12_02 - v "create.symData" dodana moznost "rDist" za samo relativne porazdelitve, vsota frekvenc = 1
## clamix12_04 - "create.symData" ni NA vrednosti (te naj se dodajo kot svoja kategorija v data.frame-u)
## clamix12_05 - "leaderSO" dodana 2 if stavka za "report", 
##             - "create.symData" - dodana rDist.KT in pDist.KT (Krichevsky-Trofimov estimator za relativne frekvence)


# -----------------------------------------------------------------------

#source("tableDiss2.R")   ### klicemo sami

# creates an empty SO of type "symObject"
empty.symObject <- function(nCats){
  nVar <- length(nCats)
  s <- vector("list", nVar)
  for(i in 1:nVar) s[[i]] <- double(nCats[i]+1)
  class(s) <- "symObject"
  return(s)
}

# encode numerical vector using given encoding
encodeSO <- function(x,encoding,codeNA){
  if(is.na(x)) return(codeNA)
  for(i in 1:length(encoding)) if(encoding[[i]](x)) return(i)
}

# adapted leaders method
#---------------------------------------------
# VB, 16. julij 2010
#   dodana omejitev na najmanj?e ?tevilo enot v skupini (1.7.2016)
    # iz katere jemljemo najbolj oddaljen element v primeru prazne skupine
#   omejeni polmer

leaderSO <- function(dataset,maxL,initial=NULL,stabil=1e-6,report=TRUE,interact=TRUE,type="d1"){
  #attach(dataset)
  # initial - initial clustering into maxL clusters 
  so <- dataset$so; alpha <- dataset$alpha
  SOs <- dataset$SOs; nVar <- length(so)
  numSO <- length(SOs)
  L <- vector("list",maxL); Ro <- numeric(maxL)
  if(is.null(initial)){
  # random partition into maxL clusters
    clust <- sample(1:maxL,numSO,replace=TRUE)
  }else{clust <- initial}
  tim <- 1; step <- 0
  pOld <- double(maxL) ### important when interact = FALSE
  repeat {
    step <- step+1
  # new leaders - determine the leaders of clusters in current partition
    L <- tlead(SOs,numSO,nVar,maxL,so,clust,ts[[type]])
  # new partition - assign each unit to the nearest new leader
    clust <- integer(numSO)
    R <- numeric(maxL); p <- double(maxL)
    for(i in 1:numSO){d <- double(maxL)
      for(k in 1:maxL){
        d[[k]] <- distSO(SOs[[i]],L[[k]],nVar,alpha,deltas[[type]])}
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
    if(report){
      cat("\nStep",step,"\n")
      print(table(clust)); print(R); print(Ro-R); Ro <- R;
      print(p); print(sum(p)); flush.console()
    }
    tim <- tim-1
    if(tim<1 && interact){   ## added for possibility interact = FALSE
      ans <- readline("Times repeat = ")
      tim <- as.integer(ans); if (tim < 1)  break
    }else{                  ## added for possibility interact = FALSE
      if(!interact){
        if(tim!=0){ ## at least one run of leaders already run and stabilized
          if(report){
            print(pOld);print(p);flush.console()}
          if(max(pOld-p)<stabil) break}
        pOld <- p}}
  # in the case of empty cluster put in the most distant SO (from cluster with at least 2 units)
    repeat{
      t <- table(clust); em <- setdiff(1:maxL,as.integer(names(t)))
      if(length(em)==0) break
      j <- em[[1]]; rmax <- 0; imax <- 0
      ## only eligible (units from cluster with at least 2 units)
      eligibleClu = as.integer(names(which(t!=1))); eligibleSOind = which(clust %in% eligibleClu)
      for(i in eligibleSOind){d <- double(maxL)
        for(k in 1:maxL){
          d[[k]] <- distSO(SOs[[i]],L[[k]],nVar,alpha,deltas[[type]])}
        r <- max(d);  if(rmax<r) {rmax <- r; imax <- i}
      }
      clust[[imax]] <- j; L[[j]] <- SOs[[imax]]
      if(report){
        cat("*** empty cluster",j,"- SO",imax,"transfered, rmax =",rmax,"\n");flush.console()}
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
#   imena notranjih to?k v drevesu
## 4.7.2016 - spremeni klice funkcij za izracun razlicnosti
## dodaj razlicnosti d2, d4, d5, d6

hclustSO <- function(dataset,type="d1"){
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
      Z <- zlead(L[[i]],L[[j]],nVar,so,L[i],L[j],zs[[type]])
      D[i,j] <- distC(L[[i]],L[[j]],Z,nVar,alpha,L[i],L[j],Ds[[type]])
      D[j,i] <- D[i,j]
  }
  active <- 1:numL; m <- matrix(nrow=numLm,ncol=2)
  node <- rep(0,numL); h <- numeric(numLm); LL <- vector("list",numLm)
  temp <- vector("list",numL)
  for(k in 1:numLm) {LL[[k]] <- so; temp[[k]] <- k} ## elements in cluster k
  temp[[numL]] <- numL # the last element too
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
    LL[[k]] <- zlead(Lq,Lp,nVar,so,L[temp[[q]]],L[temp[[p]]],zs[[type]])
    active <- setdiff(active,p)
    Lpq <- LL[[k]]
    temp[[p]] <- temp[[q]] <- c(temp[[p]],temp[[q]])
  # determine dissimilarities to the new cluster
    for(s in setdiff(active,q)){
      if(node[s]==0){Ls <- L[[s]]} else {Ls <- LL[[node[s]]]}
      Z <- zlead(Lpq,Ls,nVar,so,L[temp[[q]]],L[temp[[s]]],zs[[type]])
      D[q,s] <- distC(Lpq,Ls,Z,nVar,alpha,L[temp[[q]]],L[temp[[s]]],Ds[[type]])
      D[s,q] <- D[q,s]
    }
    node[[q]] <- k
  }
  hc <- list(merge=m,height=h,order=orDendro(numLm),labels=names(L),
    method="adapted ward",call=match.call(),dist.method="squared euclidean",leaders=LL)
  class(hc) <- "hclust"
  #detach(dataset)
  return(hc)
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

summary.symData <- function(object,...){
  # summary of symData - vector of symObjects
  cat("summary symData",":","\n\n",sep=" ")
  cat("Dimension (units x variables):", length(object$SOs),"x",length(object$so),"\n",sep=" ")
  cat("Type of distributions:", object$type,"\n", sep=" ")
  cat("Outlook of symObject: ","\n")
  print(object$namedSO)
}

# create.symData <- function(datalist,type="gDist",alpha=NULL){
#   # TYPE - PDIST - SUM!?!? uses relative distribution, but save absolute total=n
#   # type rDist uses only relative distribution, total=1
#   # creates object symData
#   # from list of dataframes
#   nVar <- length(datalist)
#   if (is.null(alpha))
#     alpha <- rep(1/nVar,nVar)
#   n_units <- nrow(datalist[[1]])
#   nVar_names <- names(datalist)
#   unit_names <- rownames(datalist[[1]])
#  for (i in 1:nVar){
#     temp <- nrow(datalist[[i]])
#     if (temp != n_units)
#       stop("number of units must match!")
#   }
#   dataset <- vector(mode="list",length=nVar)
#   nCats <- vector(mode="numeric",length=nVar)
#   for (i in 1:nVar){
#     ### changed for probabilities - if clause ###
#     n <- apply(datalist[[i]],1,sum)
#     if(type=="pDist")   ## use relative distribution, but save absolute total=n
#       dataset[[i]] <-
#         cbind(datalist[[i]]/n,0,n)
#     else
# 	    if(type=="rDist")   ## use only relative distribution, total=1
#            dataset[[i]] <-
#            cbind(datalist[[i]]/n,0,1)
#         else
#            dataset[[i]] <- cbind(datalist[[i]],0,n)
#     coldata <- ncol(dataset[[i]])
#     nCats[i] <- coldata-1
#     colnames(dataset[[i]])[[coldata-1]] <- "NA"; colnames(dataset[[i]])[[coldata]] <- "num"
#   }
#   SOs <- vector(mode="list",length=n_units)
#   for (i in 1:n_units){
#     s <- vector("list", nVar)
#     for (j in 1:nVar){
#       s[[j]] <- as.double(dataset[[j]][i,])
#     }
#     names(s) <- nVar_names
#     class(s) <- "symObject"
#     SOs[[i]] <- s
#  }
#   names(SOs) <- unit_names
#   so <- empty.symObject(nCats)
#   names(so) <- nVar_names
#   namedSO <- so
#   for (i in 1:nVar)
#     names(namedSO[[i]]) <- colnames(dataset[[i]])
# 
#   # make an object
#   object <- list(SOs=SOs,so=so,namedSO=namedSO,alpha=alpha,type=type)
#   class(object) <- "symData"
#   return(object)
# }

###### no NA values !!! ######

create.symData <- function(datalist,type="pDist",alpha=NULL){
  # type pDist uses relative distribution, but save absolute total=n
  # type rDist uses only relative distribution, total=1
  # type pDist.KT, rDist.KT - Krichevsky-Trofimov estimator for relative dist. (total=n, total=1)
  # KT relative estimates for p_i = (x_i + 1/nCat)/(sum(x_i) + 1)
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
    ### changed for probabilities - if clause ###
    n <- apply(datalist[[i]],1,sum)
    nCats[i] <- ncol(datalist[[i]])
    if(type %in% c("pDist","rDist","pDist.KT","rDist.KT")){
     if(type=="pDist")   ## use relative distribution, but save absolute total=n
      dataset[[i]] <- cbind(datalist[[i]]/n,n)
     if(type=="rDist")   ## use only relative distribution, total=1
      dataset[[i]] <- cbind(datalist[[i]]/n,1)
     if(type=="pDist.KT")   ## use KT relative distribution, but save absolute total=n
      dataset[[i]] <- cbind((datalist[[i]]+1/nCats[i])/(n+1),n)
     if(type=="rDist.KT")   ## use KT relative distribution, total=1
      dataset[[i]] <- cbind((datalist[[i]]+1/nCats[i])/(n+1),1)
    }else
      dataset[[i]] <- cbind(datalist[[i]],n)
    colnames(dataset[[i]])[[nCats[i]+1]] <- "num"
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

