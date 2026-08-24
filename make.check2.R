l <- as.numeric(commandArgs(trailingOnly=TRUE)[1]) #needed to queue jobs on storm
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
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  timeclasses <- 2
  ageclasses <- 4
  beta.struc <- list(1:(11-1),11:15)
  reps <- 89
  cores <- reps
  
  par.fun <- function(mle){
    c(mle[1:5],logistic(mle[6:length(mle)]))
  }
  par.name <- c(paste0("alpha",1:3),paste0("beta",10:11),paste0("rho",c(23,4,5)),"delta","kappa","gamma","epsilon")
  
  nf <- function(n){
    starter <- theta.s(4) # just zero starts for now
    op <- nlm(f=ll.il.ms, # likelihood function
              p=starter, # initial values
              ageclasses=ageclasses,
              timeclasses=timeclasses,
              beta.struc=beta.struc,
              ch=ch,
              hessian=TRUE)
    op0 <- c(op,"start"=list(starter))
    
  
    time.start <- Sys.time()
    for(i in 1:n){
      assign(paste0("op",i),c(nlm(f=ll.il.ms, # likelihood function
                                  p=eval(as.name(paste0("op",i-1)))$estimate, # initial values
                                  ageclasses=ageclasses,
                                  timeclasses=timeclasses,
                                  beta.struc=beta.struc,
                                  ch=ch,
                                  hessian=TRUE),"start"=list(eval(as.name(paste0("op",i-1)))$estimate)))
      cat("\r",i)
    }
    time.end <- Sys.time()
    
    df <- data.frame("par"=rep(par.name,n+1),
                     "n"=rep(0:n,each=length(par.name)),
                     "mle"=c(sapply(0:n,function(x) par.fun(eval(as.name(paste0("op",x)))$estimate))),
                     "value"=c(sapply(0:n,function(x) eval(as.name(paste0("op",x)))$minimum)),
                     "aic"=c(sapply(0:n,function(x) eval(as.name(paste0("op",x)))$minimum+24)),
                     "convergence"=rep(sapply(0:n,function(x) eval(as.name(paste0("op",x)))$code),each=length(par.name)),
                     "start"=c(sapply(0:n,function(x) eval(as.name(paste0("op",x)))$start)))
    df2 <- data.frame("phi"=logistic(c(phi.trans(op5$estimate[1:3],op5$estimate[4:5],beta.struc)[,1:4,1])),
                      "ageclass"=as.factor(rep(1:4,each=15)),
                      "time"=rep(1:15,4))
    return(list("df.mle"=df,"df.mle.phi"=df2,"time"=time.end-time.start))
  }
  
  n <- 5
  biglist <- mclapply(rep(n,reps),nf,mc.cores=cores)
  
  save(biglist,file="code.check2op.RData")
  
  bldf <- lapply(1:reps,function(x) biglist[[x]]$df.mle)
  
  bigdf <- do.call("rbind",bldf)
  
  save(bigdf,file="code.check2.RData")
  
  bigdf <- cbind(bigdf,"rep"=rep(1:reps,each=length(par.name)*(n+1)))
  
  save(bigdf,file="code.check2.RData")
}