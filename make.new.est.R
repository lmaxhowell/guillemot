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
  
  ll.il.c <- function(theta,ageclasses,ch){
    Time <- ncol(ch)-1
    struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
    
    phi <- untrans(logistic(theta[1:ageclasses]),c(as.list(1:ageclasses),list((ageclasses+1):Time)),struc.const$time,struc.const$state)
    
    theta <- theta[-c(1:(ageclasses))]
    rho <- logistic(rho.trans(theta[1:3]))
    delt <- untrans(logistic(theta[4]),struc.const$age,struc.const$time,struc.const$state)
    kap <- untrans(logistic(theta[5]),struc.const$age,struc.const$time,struc.const$state)
    gam <- untrans(logistic(theta[6]),struc.const$age,struc.const$time,struc.const$state)
    eps <- untrans(logistic(theta[7]),struc.const$age,struc.const$time,struc.const$state)
    
    psi <- make.psi(delt,kap,rho,gam,eps)
    ll <- il(ch,phi,psi)
    
    return(-ll)
  }
  cores <- 30
  n <- 30
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  timeclasses <- 2
  ageclasses <- 4
  beta.struc <- list(1:(11-1),11:15)
  tb <- table(unlist(ch3[,1:16]))
  
  Bs <- tb["LB"] + tb["L_B"]
  B_s <- tb["LB_"] + tb["L_B_"]
  Ls <- tb["LB"] + tb["LB_"]
  L_s <- tb["L_B"] + tb["L_B_"]
  
  gamm <- unname(Bs/(Bs+B_s))
  del <- unname(L_s/(Ls+L_s))
  
  par.name <- c(paste0("alpha",1:3),paste0("beta",10:11),paste0("rho",c(23,4,5)),"delta","kappa","gamma","epsilon")
  
  theta.s <- function(n=4){ # want n to be desired number of age classes
    x <- n-1+2+3 # based on two timeclasses
    rtn <- c(logit(runif(x,0.1,0.9)),logit(del),logit(runif(1,0.1,0.9)),logit(gamm),logit(runif(1,0.1,0.9))) # so need 12 values
    return(rtn)
  }
  
  op <- mclapply(1:n,function(x){
                          start <- theta.s()
                          opl <- optim(start,
                                       ll.il.ms,
                                       ageclasses=ageclasses,
                                       timeclasses=timeclasses,
                                       beta.struc=beta.struc,
                                       ch=ch3,
                                       control=list(maxit=10000))
                          return(c(opl,list("start"=start)))},mc.cores=cores)
  
  par.fun <- function(mle){
    return(c(mle[1:5],logistic(mle[6:length(mle)])))
  }
  
  df <- data.frame("par"=rep(par.name,n),
                   "mle"=c(sapply(1:n,function(x) par.fun(op[[x]]$par))),
                   "value"=rep(sapply(1:n,function(x) op[[x]]$value),each=length(par.name)),
                   "convergence"=rep(sapply(1:n,function(x) op[[x]]$convergence),each=length(par.name)),
                   "start"=c(sapply(1:n,function(x) par.fun(op[[x]]$start))))
  
  theta.s2 <- function(n=4){ # want n to be desired number of age classes
    x <- n+3 
    rtn <- c(logit(runif(x,0.1,0.9)),logit(del),logit(runif(1,0.1,0.9)),logit(gamm),logit(runif(1,0.1,0.9))) # so need 12 values
    return(rtn)
  }
  
  par.name2 <- c(paste0("phi",1:4),paste0("rho",c(23,4,5)),"delta","kappa","gamma","epsilon")
  
  op2 <- mclapply(1:n,function(x){
                          start <- theta.s2()
                          opl <- optim(start,
                                       ll.il.c,
                                       ageclasses=ageclasses,
                                       ch=ch3,
                                       control=list(maxit=10000))
                          return(c(opl,list("start"=start)))},mc.cores=cores)
  
  df2 <- data.frame("par"=rep(par.name2,n),
                    "mle"=c(sapply(1:n,function(x) logistic(op2[[x]]$par))),
                    "value"=rep(sapply(1:n,function(x) op2[[x]]$value),each=length(par.name2)),
                    "convergence"=rep(sapply(1:n,function(x) op2[[x]]$convergence),each=length(par.name2)),
                    "start"=c(sapply(1:n,function(x) logistic(op2[[x]]$start))))
  
  
  save(df,df2,file="estimate.RData")
  
}