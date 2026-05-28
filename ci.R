ci <- function(mu,sigma2,level=0.95){
  alpha <- 1-level
  plusminus <- sqrt(sigma2)*qnorm(1-alpha/2)
  return(cbind(mu-plusminus,mu,mu+plusminus))
}