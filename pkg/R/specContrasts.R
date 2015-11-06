specContrasts <- function(rs,rc,j){
    ## rs = leader of whole sample (as symObject)
    ## rc = leader of one cluster (as symObject)
  mj = length(rs[[j]])  #components + 1
  sizers = rs[[j]][mj]
  sizerc = rc[[j]][mj]
  prs = rs[[j]][-c((mj-1),mj)]/sizers # as probability, last two are NA and N
  prc = rc[[j]][-c((mj-1),mj)]/sizerc # as probability, last two are NA and N
  spec = sum(abs(prs-prc))/2
  contr = prc/prs
  contr = sapply(contr,FUN = function(x)return(ifelse(!is.na(x)&x<1,-x^(-1),x)))
  return(list(specificity = spec, constrasts = contr))
}
  
