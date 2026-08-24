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
  
  timestruc <- function(np){
    if(np==1){
      bas <- list(1:15)
    }
    else{
      bas <- unname(split(1:15, cut(1:15, breaks = np, labels = FALSE)))
    }
    return(bas)
  }
  
  ll.il.ms <- function(theta,np,ch){
    ageclasses <- 4
    timeclasses <- 2
    beta.struc <- list(1:(11-1),11:15)
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
    
    struc.time <- list("age"=list(1:Time),"time"=timestruc(np),"state"=list(1:length(states)))
    struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
    # struc.state <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
    
    theta <- theta[-c(1:(timeclasses+ageclasses+2))]
    
    delt <- untrans(logistic(theta[1:np]),struc.time$age,struc.time$time,struc.time$state)
    kap <- untrans(logistic(theta[(np+1):(2*np)]),struc.time$age,struc.time$time,struc.time$state)
    gam <- untrans(logistic(theta[(2*np+1):(3*np)]),struc.time$age,struc.time$time,struc.time$state)
    eps <- untrans(logistic(theta[(3*np+1):(4*np)]),struc.time$age,struc.time$time,struc.time$state)
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
  n <- 89*3 # 267
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
  
  dat.sim.wrap2 <- function(theta,np,seed){
    ageclasses <- 4
    timeclasses <- 2
    beta.struc <- list(1:(11-1),11:15)
    alpha <- theta[1:(ageclasses-1)]
    beta <- theta[ageclasses:(timeclasses+ageclasses-1)]
    rho.input <- theta[(timeclasses+ageclasses):(timeclasses+ageclasses+2)]
    phi <- logistic(phi.trans(alpha,beta,beta.struc))
    rho <- logistic(rho.trans(rho.input))
    
    struc.const <- list("age"=list(1:Time),"time"=timestruc(np),"state"=list(1:length(states)))
    struc.time <- list("age"=list(1:Time),"time"=as.list(1:(Time-1)),"state"=list(1:length(states)))
    
    theta <- theta[-c(1:(timeclasses+ageclasses+2))]
    
    delt <- untrans(logistic(theta[1:np]),struc.time$age,struc.time$time,struc.time$state)
    kap <- untrans(logistic(theta[(np+1):(2*np)]),struc.time$age,struc.time$time,struc.time$state)
    gam <- untrans(logistic(theta[(2*np+1):(3*np)]),struc.time$age,struc.time$time,struc.time$state)
    eps <- untrans(logistic(theta[(3*np+1):(4*np)]),struc.time$age,struc.time$time,struc.time$state)
    df <- dat.simulate2(phi,delt,kap,rho,gam,eps,ni2,seed)
    return(df)
  }
  set.seed(2083304330)
  seeds <- sample(1:.Machine$integer.max,n)
  ni2 <- lapply(1:15,function(x) ni2[[x]])
  
  for(np in 1:15){
    otherparm <- c(del[1:np],kap[1:np],gam[1:np],eps[1:np])
    theta <- c(alph,bet,logit(rho),logit(otherparm))
    thetat <- c(alph,bet,rho,otherparm) # non transformed
    dat <- mclapply(1:n, function(x) dat.sim.wrap2(theta,np,
                                                   seeds[x]),mc.cores=TRUE)
    
    op <- mclapply(1:n, function(x){
                    starter <- rep(0,length(theta)) # just zero starts for now
                    op <- nlm(f=ll.il.ms, # likelihood function
                              p=starter, # initial values
                              np=np,
                              ch=dat[[x]],
                              hessian=TRUE)
                    return(c(op,"start"=list(starter)))
                  },mc.cores=cores)
    par.name <- c(paste0("alpha",1:3),paste0("beta",10:11),paste0("rho",c(23,4,5)),paste0("delta",1:np),paste0("kappa",1:np),paste0("gamma",1:np),paste0("epsilon",1:np))
    par.group <- c(rep(c("alpha","beta","rho","delta","kappa","gamma","epsilon"),c(3,2,3,np,np,np,np)))
    par.fun <- function(mle){
      c(mle[1:5],logistic(mle[6:length(mle)]))
    }
    
    assign(paste0("df.mle",np),data.frame("par"=rep(par.name,n),
                                          "mle"=c(sapply(1:n,function(x) par.fun(op[[x]]$estimate))),
                                          "convergence"=rep(sapply(1:n,function(x) op[[x]]$code),each=length(theta)),
                                          "rep"=rep(1:n,each=length(theta)),
                                          "seed"=rep(seeds,each=length(theta)),
                                          "group"=rep(par.group,n),
                                          "numtest"=rep(np,length(theta)*n)))
    assign(paste0("df.true",np),data.frame("par"=par.name,
                                           "value"=thetat,
                                           "group"=par.group,
                                           "numtest"=rep(np,length(thetat))))
    print(np)
  }
  
  df.mle <- do.call(rbind, lapply( paste0("df.mle", 1:15) , get) )
  df.true <- do.call(rbind, lapply( paste0("df.true", 1:15) , get) )
  
  save(df.mle,df.true,file="sim.check3.RData")
}