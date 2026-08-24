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
  
  theta.s <- function(n,state.depen=FALSE){ # want n to be desired number of age classes
    x <- ifelse(state.depen==FALSE,8,19) # based on two timeclasses
    return(logit(runif(x+n,0.1,0.9)))
  }
  
  ll.il.ms <- function(theta,ageclasses,timeclasses,beta.struc,ch){
    # beta.struc should be a list corresponding to
    # where the breaks in beta go e.g. list(1:10,11:15)
    # would be a change in beta at the 11 time point
    alpha <- theta[1:(ageclasses-1)]
    beta <- theta[ageclasses:(timeclasses+ageclasses-1)]
    rho.input <- theta[(timeclasses+ageclasses):(timeclasses+ageclasses+2)]
    phi <- logistic(phi.trans(alpha,beta,beta.struc))
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
    
    theta <- theta[-c(1:(timeclasses+ageclasses+2))]
    
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
    return(-ll)
  }
  
  # ll.il.ms.s <- function(theta,ageclasses,ch){
  #   alpha <- theta[1:(ageclasses-1)]
  #   beta <- theta[ageclasses:(Time+ageclasses-2)]
  #   rho.input <- theta[(Time+ageclasses-1):(Time+ageclasses+1)]
  #   phi <- logistic(phi.trans(alpha,beta))
  #   rho <- logistic(rho.trans(rho.input))
  #   # rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
  #   
  #   # Ages <- length(alpha)
  #   Time <- ncol(ch)-1
  #   
  #   struc.time <- list("age"=list(1:Time),"time"=as.list(1:(Time-1)),"state"=list(1:length(states)))
  #   struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  #   struc.state <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
  #   struc.kap <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(c(1,2,8),3,4,5,6,7))
  #   struc.delt <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(c(1,2),3,4,5,6,7,8))
  #   
  #   # so kappa needs 5 parameters and delta needs 6
  #   # add a zero at the beginning for the rows where that parameter doesn't actually show up
  #   delt <- untrans(logistic(c(0,theta[(Time+ageclasses+2):(Time+ageclasses+2+5)])),struc.delt$age,struc.delt$time,struc.delt$state)
  #   kap <- untrans(logistic(c(0,theta[(Time+ageclasses+2+5+1):(Time+ageclasses+2+5+1+4)])),struc.kap$age,struc.kap$time,struc.kap$state)
  #   
  #   theta <- theta[-c(1:(Time+ageclasses+2+5+1+4))]
  #   
  #   gam <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
  #   eps <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
  #   
  #   psi <- make.psi(delt,kap,rho,gam,eps)
  #   # print(psi)
  #   
  #   ll <- il(ch,phi,psi)
  #   return(ll)
  # }
  
  TimeClasses <- 2:15
  cores <- length(TimeClasses)
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  op <- mclapply(TimeClasses, function(x){
          timeclasses <- 2
          ageclasses <- 4
          beta.struc <- list(1:(x-1),x:15)
          starter <- theta.s(4)
          op <- nlm(f=ll.il.ms, # likelihood function
                    p=starter, # initial values
                    ageclasses=ageclasses,
                    timeclasses=timeclasses,
                    beta.struc=beta.struc,
                    ch=ch,
                    hessian=TRUE)
          return(c(op,"start"=list(starter)))
        },mc.cores=cores)
  
  AICc <- function(ll,n,ac=4){
    p <- length(theta.s(ac))
    aic <- 2*ll + 2*p # not minus ll as assuming its coming out of nlm already its negative
    aicc <- aic + 2*p*(p+1)/(n-p-1)
    return(aicc)
  }
  n1 <- nrow(ch) # number of individuals
  n2 <- nrow(uch) # number of unique chs
  # n3 # number of transition pairs
  aics1 <- sapply(TimeClasses,function(x) AICc(op[[x-1]]$minimum,n1))
  aics2 <- sapply(TimeClasses,function(x) AICc(op[[x-1]]$minimum,n2))
  
  df.aic1 <- data.frame("AICc"=aics1,
                        "TimeChange"=TimeClasses,
                        "convergence"=sapply(TimeClasses,function(x) op[[x-1]]$code==1))
  df.aic2 <- data.frame("AICc"=aics2,
                        "TimeChange"=TimeClasses,
                        "convergence"=sapply(TimeClasses,function(x) op[[x-1]]$code==1))
  par.fun <- function(mle){
    c(mle[1:5],logistic(mle[6:length(mle)]))
  }
  
  save(df.aic1,df.aic2,file="MS.time.RData")
  
  nam <- c(paste0("alpha",1:3),paste0("beta",1:2),paste0("rho",c(23,4,5)),"delta","kappa","gamma","epsilon")
  
  df.par <- data.frame("par"=rep(nam,length(TimeClasses)),
                       "mle"=c(unlist(sapply(TimeClasses,function(x) par.fun(op[[x-1]]$estimate)))),
                       "TimeChange"=rep(TimeClasses,each=length(theta.s(4))),
                       "convergence"=rep(sapply(TimeClasses,function(x) op[[x-1]]$code==1),each=length(theta.s(4))),
                       "AICc"=rep(aics1,each=length(theta.s(4))),
                       "AICc_un"=rep(aics2,each=length(theta.s(4))),
                       "start"=c(unlist(sapply(TimeClasses,function(x) op[[x-1]]$start))))
  
  save(df.aic1,df.aic2,df.par,file="MS.time.RData")
}