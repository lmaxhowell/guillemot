trans <- function(param){
  # param is a matrix of state and time dependent
  # parameters, so could be or r
  # transfroms FROM matrices TO vectors
  theta <- unique(c(param))
  return(theta)
}

untrans <- function(theta,tme,state){
  # tme is the Time structure of what to keep constant
  # state is the state structure of what to keep constant
  # both should be a list of indices
  # so list(1,2,3:Time) for tme means that tme is specific
  # to first years, second years and then constant after
  # and the same applies to first years
  State <- length(unlist(state))
  Time <- length(unlist(tme))
  param <- array(dim=c(State,Time))
  for(a in 1:length(tme)){
    for(t in 1:length(state)){
      param[state[[t]],tme[[a]]] <- theta[1]
      if(length(theta)>1){
        theta <- theta[2:length(theta)]
      }
    }
  }
  return(param)
}