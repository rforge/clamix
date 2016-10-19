# change because NA not in the creating of the SDA data set

specContrasts <- function(rs,rc,j){
    ## rs = leader of whole sample (as symObject)
    ## rc = leader of one cluster (as symObject)
  mj = length(rs[[j]])  #components + 1
  sizers = rs[[j]][mj]
  sizerc = rc[[j]][mj]
  prs = rs[[j]][-mj]/sizers # as probability, last is N
  prc = rc[[j]][-mj]/sizerc # as probability, last is N
  spec = sum(abs(prs-prc))/2
  contr = prc/prs
  contr = sapply(contr,FUN = function(x)return(ifelse(!is.na(x)&x<1,-x^(-1),x)))
  return(list(specificity = spec, constrasts = contr))
}
  
