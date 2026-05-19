logit <- function(p){
  return(log(p/(1-p)))
}

logistic <- function(q){
  return(1/(exp(-q)+1))
}