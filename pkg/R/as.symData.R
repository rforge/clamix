as.symData <- function(L,nvnames= NULL){
  # checks if L is list of symObjects
  # makes symData if is true
  # or returns FALSE
  if (!is.list(L) && class(L)!="symObject")
    return(FALSE)
  else{
    if (class(L)=="symObject")
      L1 <- L
    else
      L1 <- L[[1]]
    nVar <- length(L1)
    nCats <- NULL
    nnames <- vector(mode="list",length(L1))
    for (i in (1:length(L1))){
      nCats <- c(nCats,length(L1[[i]])-1)
      if (is.null(nvnames)){
        nnames[[i]] <- names(L1[[i]])}
    }
    so <- empty.symObject(nCats)
    namedSO <- so
    names(namedSO) <- names(L1)
    if (is.null(nvnames)) nvnames <- nnames
    for (i in 1:nVar)
      names(namedSO[[i]]) <- nvnames[[i]]
    alpha <- rep(1,length(L1))/length(L1)
    type <- "gDist"
  }
  newdata <- list(SOs=L,so=so,namedSO=namedSO,alpha=alpha,type=type)
  class(newdata) <- "symData"
  return(newdata)
}
