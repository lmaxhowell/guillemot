l <- as.numeric(commandArgs(trailingOnly=TRUE)[1]) #needed to queue jobs on storm
# file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
# file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
# file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
# lapply(file.names, source) # source all functions needed
# load("guillemot.RData")
# save(list=ls(),file="gf.storm.RData")

if(l==1){
  load("gf.storm.RData")
  library(parallel)
  phi.trans <- function(alpha,beta){
    # inputs should be vectors
    alpha <- c(0,alpha) # to constrain alpha_1=0
    Ages <- length(alpha)
    Time <- length(beta)
    phi <- array(0,dim=c(Time,Time+1,length(states)))
    for(a in 1:(Time+1)){
      for(t in 1:Time){
        # print(c(t,a,min(a,Ages)))
        phi[t,a,] <- alpha[min(a,Ages)]+beta[t]
      }
    }
    return(phi)
  }
  
  rho.trans <- function(rho.input){
    rho <- array(0,dim=c(Time-1,Time,length(states)))
    # rho[,1:2,] <- 0 # sets as default anyway # ages 0 and 1
    rho[,3,] <- rho.input[1] # age two
    rho[,4,] <- rho.input[2] # age 3
    rho[,5:Time,] <- rho.input[3] # age 4+
    return(rho)
  }
  
  theta.s <- function(n,state.depen=FALSE){ # want n to be desired number of age classes
    x <- ifelse(state.depen==FALSE,21,30)
    return(logit(rep(0.55,x+n)))
  }
  
  ll.il.ms <- function(theta,ageclasses,ch){
    alpha <- theta[1:(ageclasses-1)]
    beta <- theta[ageclasses:(Time+ageclasses-2)]
    rho.input <- theta[(Time+ageclasses-1):(Time+ageclasses+1)]
    phi <- logistic(phi.trans(alpha,beta))
    # print(length(alpha))
    # print(length(beta))
    # print(length(rho.input))
    # print(dim(phi))
    rho <- logistic(rho.trans(rho.input))
    # rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
    
    # Ages <- length(alpha)
    Time <- ncol(ch)-1
    
    # struc.time <- list("age"=list(1:Time),"time"=as.list(1:(Time-1)),"state"=list(1:length(states)))
    struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
    # struc.state <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
    
    theta <- theta[-c(1:(Time+ageclasses+1))]
    
    delt <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
    kap <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
    gam <- untrans(logistic(theta[3]),struc.const$age,struc.const$time,struc.const$state)
    eps <- untrans(logistic(theta[4]),struc.const$age,struc.const$time,struc.const$state)
    # print(eps)
    # print(theta)
    
    # print(dim(delt))
    # print(dim(kap))
    # print(dim(rho))
    # print(dim(gam))
    # print(dim(eps))
    psi <- make.psi(delt,kap,rho,gam,eps)
    # print(dim(psi))
    # print(phi[,1:ageclasses,1])
    
    ll <- il(ch,phi,psi)
    return(ll)
  }
  
  ll.il.ms.s <- function(theta,ageclasses,ch){
    alpha <- theta[1:(ageclasses-1)]
    beta <- theta[ageclasses:(Time+ageclasses-2)]
    rho.input <- theta[(Time+ageclasses-1):(Time+ageclasses+1)]
    phi <- logistic(phi.trans(alpha,beta))
    rho <- logistic(rho.trans(rho.input))
    # rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
    
    # Ages <- length(alpha)
    Time <- ncol(ch)-1
    
    struc.time <- list("age"=list(1:Time),"time"=as.list(1:(Time-1)),"state"=list(1:length(states)))
    struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
    struc.state <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
    struc.kap <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(c(1,2,8),3,4,5,6,7))
    struc.delt <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(c(1,2),3,4,5,6,7,8))
    
    # so kappa needs 5 parameters and delta needs 6
    # add a zero at the beginning for the rows where that parameter doesn't actually show up
    delt <- untrans(logistic(c(0,theta[(Time+ageclasses+2):(Time+ageclasses+2+5)])),struc.delt$age,struc.delt$time,struc.delt$state)
    kap <- untrans(logistic(c(0,theta[(Time+ageclasses+2+5+1):(Time+ageclasses+2+5+1+4)])),struc.kap$age,struc.kap$time,struc.kap$state)
    
    theta <- theta[-c(1:(Time+ageclasses+2+5+1+4))]
    
    gam <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
    eps <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
    
    psi <- make.psi(delt,kap,rho,gam,eps)
    # print(psi)
    
    ll <- il(ch,phi,psi)
    return(ll)
  }
  
  
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  
  AgeClasses <- 2:8
  cores <- length(AgeClasses)
  op <- mclapply(AgeClasses,function(x) optim(theta.s(x,FALSE),ll.il.ms,ageclasses=x,ch=ch,
                                              control=list(fnscale=-1,maxit=25000),method="Nelder-Mead",hessian=TRUE),mc.cores=cores)
  op.s <- mclapply(AgeClasses,function(x) optim(theta.s(x,TRUE),ll.il.ms.s,ageclasses=x,ch=ch,
                                                control=list(fnscale=-1,maxit=25000),method="Nelder-Mead",hessian=TRUE),mc.cores=cores)
  aics1 <- sapply(AgeClasses,function(x) 2*length(theta.s(x))-2*op[[x-1]]$value)
  aics2 <- sapply(AgeClasses,function(x) 2*length(theta.s(x))-2*op.s[[x-1]]$value)
  
  df.aic <- data.frame("aic"=c(aics1,aics2),
                       "AgeClasses"=rep(AgeClasses,2),
                       "StateDependence"=rep(c(FALSE,TRUE),each=length(AgeClasses)),
                       "Convergence"=c(sapply(AgeClasses,function(x) op[[x-1]]$convergence),
                                       sapply(AgeClasses,function(x) op.s[[x-1]]$convergence)))
  
  par.name.fun <- function(ageclass,state.depen=FALSE){
    alphas <- paste0("alpha",(1:(ageclass-1))+1)
    betas <- paste0("beta",1:(Time-1))
    rhos <- paste0("rho",3:5)
    if(state.depen==TRUE){
      deltas <- paste0("delta",3:8)
      kappas <- paste0("kappa",3:7)
    }else{
      deltas <- "delta"
      kappas <- "kappa"
    }
    return(c(alphas,betas,rhos,deltas,kappas,"gamma","epsilon"))
  }
  
  names1 <- lapply(AgeClasses,par.name.fun,state.depen=FALSE)
  names2 <- lapply(AgeClasses,par.name.fun,state.depen=TRUE)
  
  save(df.aic,op,op.s,file="ModelSelectionStorm.RData")
  
  par.fun <- function(op.par,ageclass){
    output <- c(op.par[1:(Time+ageclass)],logistic(op.par[(Time+ageclass+1):length(op.par)]))
    return(output)
  }
  
  df.par <- data.frame("par"=c(unlist(names1),
                               unlist(names2)),
                       "MLE"=c(unlist(lapply(AgeClasses,function(x) par.fun(op[[x-1]]$par,x))),
                               unlist(lapply(AgeClasses,function(x) par.fun(op.s[[x-1]]$par,x)))),
                       "AgeClasses"=c(rep(AgeClasses,sapply(names1,length)),
                                      rep(AgeClasses,sapply(names2,length))),
                       "StateDependance"=rep(c(FALSE,TRUE),c(length(unlist(names1)),length(unlist(names2)))),
                       "AIC"=c(rep(aics1,sapply(names1,length)),
                               rep(aics2,sapply(names2,length))),
                       "convergence"=c(rep(sapply(AgeClasses,function(x) op[[x-1]]$convergence),sapply(names1,length)),
                                       rep(sapply(AgeClasses,function(x) op.s[[x-1]]$convergence),sapply(names2,length))))
  
  save(df.par,df.aic,op,op.s,file="ModelSelectionStorm.RData")
  
}