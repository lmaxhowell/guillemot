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
    # rho[,1:2,] <- 0 # sets as default anyway
    rho[,3,] <- rho.input[1]
    rho[,4,] <- rho.input[2]
    rho[,5:Time,] <- rho.input[3]
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
  
  
  ll.il.ms.pl <- function(theta,fix,at,ageclasses,ch){
    theta.full <- append(theta,fix,after=(at-1))
    ll <- ll.il.ms(theta.full,ageclasses,ch)
    return(ll)
  }
  
  len <- 60
  cores <- len
  maxits <- 25000
  pl.log <- seq(-3,3,length.out=len)
  pl <- logit(seq(0.05,0.95,length.out=len))
  
  for(i in 1:18){
    assign(paste0("pl",i),mclapply(1:len, function(x) optim(theta.s(4)[1:24],ll.il.ms.pl,
                                                 fix=pl.log[x],at=i,
                                                 ageclasses=4,ch=ch,
                                                 control=list(fnscale=-1,maxit=maxits),
                                                 method="Nelder-Mead"),mc.cores=cores))
    print(i)
  }
  for(i in 19:25){
    assign(paste0("pl",i),mclapply(1:len, function(x) optim(theta.s(4)[1:24],ll.il.ms.pl,
                                                 fix=pl[x],at=i,
                                                 ageclasses=4,ch=ch,
                                                 control=list(fnscale=-1,maxit=maxits),
                                                 method="Nelder-Mead"),mc.cores=cores))
    print(i)
  }
  
  save(list=paste0("pl",1:25),file="AgeClass4PLs.RData")
  
  op <- optim(theta.s(4),ll.il.ms,ageclasses=4,ch=ch,
              control=list(fnscale=-1,maxit=maxits*2),method="Nelder-Mead",hessian=TRUE)
  
  df.pl <- data.frame("value"=unlist(lapply(1:25,function(i) sapply(1:len,function(x) eval(as.name(paste0("pl",i)))[[x]]$value))),
                      "fix"=c(rep(pl.log,18),rep(seq(0.05,0.95,length.out=len),25-18)),
                      "convergence"=rep(unlist(lapply(1:25,function(i) sapply(1:len,function(x) eval(as.name(paste0("pl",i)))[[x]]$convergence))),each=len),
                      "par"=rep(c(paste0("alpha",1:3),paste0("beta",1:15),paste0("rho",1:3),"delta","kappa","gamma","epsilon"),each=len))
  
  df.op <- data.frame("value"=rep(op$value,25),
                      "fix"=c(op$par[1:18],logistic(op$par[19:25])),
                      "convergence"=rep(op$convergence,25),
                      "par"=c(paste0("alpha",1:3),paste0("beta",1:15),paste0("rho",1:3),"delta","kappa","gamma","epsilon"))
  
  save(df.pl,df.op,op,file="AgeClass4PLdf.RData")
}