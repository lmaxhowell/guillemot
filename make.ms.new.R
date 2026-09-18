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
    
    psi <- make.psi(delt,kap,rho,gam,eps)
    
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
  cores <- 89
  n <- 10
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  timeclasses <- 2
  ageclasses <- 2:6
  beta.struc <- lapply(1:14,function(x) list(1:(x),(x+1):15))
  tb <- table(unlist(ch3[,1:16]))
  
  Bs <- tb["LB"] + tb["L_B"]
  B_s <- tb["LB_"] + tb["L_B_"]
  Ls <- tb["LB"] + tb["LB_"]
  L_s <- tb["L_B"] + tb["L_B_"]
  
  gamm <- unname(Bs/(Bs+B_s))
  del <- unname(L_s/(Ls+L_s))
  
  theta.s <- function(n=4,timee=FALSE){ # want n to be desired number of age classes
    x <- ifelse(timee,n-1+2+3,n+3) # alphas+betas+rho vs phis+rho
    rtn <- c(logit(runif(x,0.1,0.9)),logit(del),logit(runif(1,0.1,0.9)),logit(gamm),logit(runif(1,0.1,0.9)))
    return(rtn)
  }
  
  par.fun <- function(mle,alphabet){
    return(c(mle[1:alphabet],logistic(mle[(alphabet+1):length(mle)])))
  }
  
  loopn <- rep(1:10,length(beta.struc))
  loopb <- rep(beta.struc,each=n)
  
  for(a in ageclasses){
    assign(paste0("op",a),mclapply(1:(n*length(beta.struc)),function(x){
                              start <- theta.s(a,timee=TRUE)
                              opl <- optim(start,
                                           ll.il.ms,
                                           ageclasses=a,
                                           timeclasses=timeclasses,
                                           beta.struc=loopb[[x]],
                                           ch=ch3,
                                           control=list(maxit=10000))
                              return(c(opl,list("start"=start)))},mc.cores=cores))
    
    seq1 <- seq(1,length(beta.struc)*n,n)
    for(i in 1:length(beta.struc)){
      selector <- sapply(seq1[i]:(seq1[i]+n-1),function(x) eval(as.name(paste0("op",a)))[[x]]$value )
      par.name <- c(paste0("alpha",1:(a-1)),paste0("beta",c(tail(beta.struc[[i]][[1]],1),beta.struc[[i]][[2]][1])),paste0("rho",c(23,4,5)),"delta","kappa","gamma","epsilon")
      selected <- eval(as.name(paste0("op",a)))[[(seq1[i]:(seq1[i]+n-1))[which.min(selector)]]]
      
      print(par.name)
      print(selected$par)
      assign(paste0("dff",i), data.frame("par"=par.name,
                                         "mle"=par.fun(selected$par,a-1+i),
                                         "value"=rep(selected$value,length(par.name)),
                                         "aic"=rep(selected$value+2*length(par.name) ,length(par.name)),
                                         "convergence"=rep(selected$convergence,length(par.name)),
                                         "ageclasses"=rep(a,length(par.name)),
                                         "betatimesplit"=rep(i,length(par.name)),
                                         "start"=par.fun(selected$start,a-1+i)))
    }
    
    assign(paste0("df",a), do.call("rbind", lapply( paste0("dff", 1:length(beta.struc)) , get) ))  
    print(a)
  }
  
  df <- do.call("rbind", lapply( paste0("df", ageclasses) , get) )
  save(df,file="ms1.RData")
}

if(l==2){
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
    
    psi <- make.psi(delt,kap,rho,gam,eps)
    
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
  cores <- 89
  n <- 10
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  timeclasses <- 2
  ageclasses <- 2:6
  beta.struc <- lapply(1:14,function(x) list(1:(x),(x+1):15))
  tb <- table(unlist(ch3[,1:16]))
  
  Bs <- tb["LB"] + tb["L_B"]
  B_s <- tb["LB_"] + tb["L_B_"]
  Ls <- tb["LB"] + tb["LB_"]
  L_s <- tb["L_B"] + tb["L_B_"]
  
  gamm <- unname(Bs/(Bs+B_s))
  del <- unname(L_s/(Ls+L_s))
  
  theta.s <- function(n=4,timee=FALSE){ # want n to be desired number of age classes
    x <- ifelse(timee,n-1+2+3,n+3) # alphas+betas+rho vs phis+rho
    rtn <- c(logit(runif(x,0.1,0.9)),logit(del),logit(runif(1,0.1,0.9)),logit(gamm),logit(runif(1,0.1,0.9)))
    return(rtn)
  }
  
  par.fun <- function(mle){
    return(logistic(mle))
  }
  
  for(a in ageclasses){
    assign(paste0("op",a),mclapply(1:(n),function(x){
                              start <- theta.s(a,timee=FALSE)
                              opl <- optim(start,
                                           ll.il.c,
                                           ageclasses=a,
                                           ch=ch3,
                                           control=list(maxit=10000))
                              return(c(opl,list("start"=start)))},mc.cores=cores))
    
    selector <- sapply(1:n,function(x) eval(as.name(paste0("op",a)))[[x]]$value )
    par.name <- c(paste0("phi",1:a),paste0("rho",c(23,4,5)),"delta","kappa","gamma","epsilon")
    selected <- eval(as.name(paste0("op",a)))[[which.min(selector)]]
      
    assign(paste0("df",a), data.frame("par"=par.name,
                                      "mle"=par.fun(selected$par),
                                      "value"=rep(selected$value,length(par.name)),
                                      "aic"=rep(selected$value+2*length(par.name) ,length(par.name)),
                                      "convergence"=rep(selected$convergence,length(par.name)),
                                      "ageclasses"=rep(a,length(par.name)),
                                      "betatimesplit"=rep(NA,length(par.name)),
                                      "start"=par.fun(selected$start)))
    print(a)
  }
  
  df <- do.call("rbind", lapply( paste0("df", ageclasses) , get) )
  save(df,file="ms2.RData")
}