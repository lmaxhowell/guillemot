trans <- function(param){
  # param is a matrix of time and age dependent
  # parameters, so could be or r
  # transfroms FROM matrices TO vectors
  theta <- unique(c(param))
  return(theta)
}

untrans <- function(theta,age,tme){
  # age is the age structure of what to keep constant
  # tme is the time structure of what to keep constant
  # both should be a list of indices
  # so list(1,2,3:Time) for age means that age is specific
  # to first years, second years and then constant after
  # and the same applies to first years
  Time <- length(unlist(tme))
  Age <- length(unlist(age))
  param <- array(dim=c(Time,Age))
  for(a in 1:length(age)){
    for(t in 1:length(tme)){
      param[tme[[t]],age[[a]]] <- theta[1]
      if(length(theta)>1){
        theta <- theta[2:length(theta)]
      }
    }
  }
  return(param)
}