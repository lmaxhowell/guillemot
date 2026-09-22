source("suff.stat.R")

il.alt <- function(phi,psi,mv){
  Indicator <- function(prob){
    # a function that checks if the probability pu into it will
    # "cause problems" in the likelihood - i.e. if we log
    # this probability will it be infinite and cause the whole ll
    # to be infinite as a result
    # or has a small amount of numerical error crept in
    # and made a probability negative?
    # in which case just set it to zero
    # will add a warning in in this case
    # this function just means we can call log(prob)
    # without worrying about the errors this can cause
    if(length(prob)==0){
      print(c(parent.frame()$i,parent.frame()$t,parent.frame()$Time))
      print(c("r",parent.frame()$current_state_index))
      print(c("s",parent.frame()$next_state_index))
      print(c("t",parent.frame()$current_time))
      print(c("a",parent.frame()$current_age))
      print(parent.frame()$ch)
      View(parent.frame()$transitions)
    }
    if(!exists("prob") | is.na(prob)){
      print(parent.frame()$ch)
      print(parent.frame()$transitions)
    }
    if(prob<0){
      return(0)
    }else{
      return(ifelse(is.finite(log(prob)),log(prob),0)) 
    }
  }
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  m <- mv$m
  v <- mv$v
  Time <- dim(mv$m)[3]
  Age <- Time
  
  if(is.null(unlist(mv$violations))==FALSE){
    warning("some transitions not allowed")
  }
  
  probs.m <- array(0,dim=dim(m)) # r,s,t,a
  probs.v <- array(0,dim=dim(v)) # r,t,a
  
  skip <- which(row.names(psi)=="S")
  for(r in 1:length(states)){
    for(t in 1:Time){
      for(a in 1:Age){
        for(s in 1:length(states)){
          if(m[r,s,t,a]>0){
            probs.m[r,s,t,a] <- log(phi[t,a,r]*psi[r,s,t,a])
          }
        }
        if(r %in% 1:2){
          if(v[r,t,a]>0){
            # print(c(r,t,a))
            probs.v[r,t,a] <- log(Chi(r,t,a,phi,psi))
          }
        }else if(r %in% 3:7){
          if(t<(Time-1)){
            if(v[r,t,a]>0){
              # print(c(r,t,a))
              probs.v[r,t,a] <- log(phi[t,a,r]*psi[r,skip,t,a]*(1-phi[t+1,a+1,skip]) + (1-phi[t,a,r]))
            }
          }else if(t==(Time-1)){
            # if we are in a breeding state at Time-1 and then a zero at Time
            # then they dont have to have died they could still be in the skipping state
            # with no information about them being alive
            if(v[r,t,a]>0){
              # print(c(r,t,a))
              probs.v[r,t,a] <- log(phi[t,a,r]*psi[r,skip,t,a] + (1-phi[t,a,r]))
            }
          }
        }
      }
    }
  }
  ll <- probs.m%*%mv$m + probs.v%*%mv$v
  return(as.numeric(ll))
}

ll.il.alt <- function(theta,phi.ind,delt.ind,kap.ind,rho.ind,gam.ind,eps.ind,p.ind,struc,mv){
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  
  phi <- untrans(logistic(theta[phi.ind]),struc$phi$age,struc$phi$time,struc$phi$state)
  delt <- untrans(logistic(theta[delt.ind]),struc$delt$age,struc$delt$time,struc$delt$state)
  kap <- untrans(logistic(theta[kap.ind]),struc$kap$age,struc$kap$time,struc$kap$state)
  rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
  gam <- untrans(logistic(theta[gam.ind]),struc$gam$age,struc$gam$time,struc$gam$state)
  eps <- untrans(logistic(theta[eps.ind]),struc$eps$age,struc$eps$time,struc$eps$state)
  p <- untrans(logistic(theta[p.ind]),struc$p$age,struc$p$time,struc$p$state)
  
  psi <- make.psi(delt,kap,rho,gam,eps,p)
  # print(psi)
  
  ll <- il.alt(phi,psi,mv)
  
  return(ll)
}

ll.il.alt2 <- function(theta,ageclasses,timeclasses,beta.struc,mv){
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  
  phi.trans <- function(alpha,beta,beta.struc){
    # inputs should be vectors
    alpha <- c(0,alpha) # to constrain alpha_0=0
    Ages <- length(alpha)
    Time_1 <- length(unlist(beta.struc)) # cause this is Time-1
    phi <- array(0,dim=c(Time_1,Time_1+1,length(states)))
    for(a in 1:(Time_1+1)){
      for(t in 1:Time_1){
        # print(c(t,a,min(a,Ages)))
        for(i in 1:length(beta.struc)){
          if(t %in% beta.struc[[i]]){
            break
          }
        }
        phi[t,a,] <- alpha[min(a,Ages)]+beta[i]
      }
    }
    return(phi)
  }
  
  rho.trans <- function(rho.input){
    rho <- array(0,dim=c(Time-1,Time,length(states)))
    # rho[,1,] <- 0 # sets as default anyway # age 0
    rho[,2:3,] <- rho.input[1] # age one and two
    rho[,4,] <- rho.input[2] # age 3
    rho[,5:Time,] <- rho.input[3] # age 4+
    return(rho)
  }
  
  alpha <- theta[1:(ageclasses-1)]
  beta <- theta[ageclasses:(timeclasses+ageclasses-1)]
  rho.input <- theta[(timeclasses+ageclasses):(timeclasses+ageclasses+2)]
  phi <- logistic(phi.trans(alpha,beta,beta.struc))
  
  rho <- logistic(rho.trans(rho.input))
  
  Time <- dim(mv$m)[3]
  
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  
  theta <- theta[-c(1:(timeclasses+ageclasses+2))]
  
  delt <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
  kap <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(logistic(theta[3]),struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(theta[4]),struc.const$age,struc.const$time,struc.const$state)
  p <- untrans(logistic(theta[5]),struc.const$age,struc.const$time,struc.const$state)

  psi <- make.psi(delt,kap,rho,gam,eps,p)

  ll <- il.alt(phi,psi,mv)
  
  return(ll)
}