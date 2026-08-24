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
    return(logit(runif(64+n,0.1,0.9)))
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
    
    struc.time <- list("age"=list(1:Time),"time"=as.list(1:(Time-1)),"state"=list(1:length(states)))
    struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
    # struc.state <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
    
    theta <- theta[-c(1:(timeclasses+ageclasses+2))]
    
    delt <- untrans(logistic(theta[1:15]),struc.time$age,struc.time$time,struc.time$state)
    kap <- untrans(logistic(theta[16:30]),struc.time$age,struc.time$time,struc.time$state)
    gam <- untrans(logistic(theta[31:45]),struc.time$age,struc.time$time,struc.time$state)
    eps <- untrans(logistic(theta[46:60]),struc.time$age,struc.time$time,struc.time$state)
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
  n <- 200
  cores <- 89
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  timeclasses <- 2
  ageclasses <- 4
  alph <- c(0.4,0.8,1.2)
  bet <- c(-0.4,-0.9)
  rho <- c(0.1,0.3,0.6)
  del <- c(0.5, 0.2, 0.6, 0.15, 0.65, 0.8, 0.4, 0.35, 0.55, 0.45, 0.85, 0.7, 0.3, 0.1, 0.2)
  kap <- c(0.15, 0.3, 0.55, 0.8, 0.7, 0.65, 0.9, 0.25, 0.35, 0.45, 0.2, 0.1, 0.75, 0.5, 0.65)
  gam <- c(0.65, 0.45, 0.4, 0.3, 0.35, 0.25, 0.55, 0.7, 0.8, 0.6, 0.5, 0.2, 0.85, 0.15, 0.55)
  eps <- c(0.4, 0.1, 0.9, 0.15, 0.5, 0.25, 0.6, 0.7, 0.55, 0.75, 0.65, 0.8, 0.35, 0.85, 0.15)
  otherparm <- c(del,kap,gam,eps)
  theta <- c(alph,bet,logit(rho),logit(otherparm))
  thetat <- c(alph,bet,rho,otherparm) # non transformed
  beta.struc <- list(1:(11-1),11:15)
  
  dat.sim.wrap2 <- function(theta,ageclasses,timeclasses,beta.struc,seed){
    alpha <- theta[1:(ageclasses-1)]
    beta <- theta[ageclasses:(timeclasses+ageclasses-1)]
    rho.input <- theta[(timeclasses+ageclasses):(timeclasses+ageclasses+2)]
    phi <- logistic(phi.trans(alpha,beta,beta.struc))
    rho <- logistic(rho.trans(rho.input))
    
    struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
    struc.time <- list("age"=list(1:Time),"time"=as.list(1:(Time-1)),"state"=list(1:length(states)))
    
    theta <- theta[-c(1:(timeclasses+ageclasses+2))]
    
    delt <- untrans(logistic(theta[1:15]),struc.time$age,struc.time$time,struc.time$state)
    kap <- untrans(logistic(theta[16:30]),struc.time$age,struc.time$time,struc.time$state)
    gam <- untrans(logistic(theta[31:45]),struc.time$age,struc.time$time,struc.time$state)
    eps <- untrans(logistic(theta[46:50]),struc.time$age,struc.time$time,struc.time$state)
    df <- dat.simulate2(phi,delt,kap,rho,gam,eps,ni2,seed)
    return(df)
  }
  set.seed(2083304330)
  seeds <- sample(1:.Machine$integer.max,n)
  ni2 <- lapply(1:15,function(x) ni2[[x]])
  dat <- mclapply(1:n, function(x) dat.sim.wrap2(theta,ageclasses,
                                                 timeclasses,beta.struc,
                                                 seeds[x]),mc.cores=TRUE)
  op <- mclapply(1:n, function(x){
    starter <- rep(0,68) # just zero starts for now
    op <- nlm(f=ll.il.ms, # likelihood function
              p=starter, # initial values
              ageclasses=ageclasses,
              timeclasses=timeclasses,
              beta.struc=beta.struc,
              ch=dat[[x]],
              hessian=TRUE)
    return(c(op,"start"=list(starter)))
  },mc.cores=cores)
  
  par.name <- c(paste0("alpha",1:3),paste0("beta",10:11),paste0("rho",c(23,4,5)),paste0("delta",1:15),paste0("kappa",1:15),paste0("gamma",1:15),paste0("epsilon",1:15))
  par.group <- c(rep(c("alpha","beta","rho","delta","kappa","gamma","epsilon"),c(3,2,3,15,15,15,15)))
  par.fun <- function(mle){
    c(mle[1:5],logistic(mle[6:length(mle)]))
  }
  df.mle <- data.frame("par"=rep(par.name,n),
                       "mle"=c(sapply(1:n,function(x) par.fun(op[[x]]$estimate))),
                       "convergence"=rep(sapply(1:n,function(x) op[[x]]$code),each=68),
                       "rep"=rep(1:n,each=68),
                       "seed"=rep(seeds,each=68),
                       "group"=rep(par.group,n))
  df.true <- data.frame("par"=par.name,
                        "value"=thetat,
                        "group"=par.group)
  save(df.mle,df.true,file="sim.check2.RData")
}