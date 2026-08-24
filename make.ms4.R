# trying ms on the real data BUT allowing two parameters in delta and kappa
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
  library(matrixcalc,include.only="is.singular.matrix")
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
  
  theta.s <- function(){ # want n to be desired number of age classes
    return(logit(runif(14,0.1,0.9)))
  }
  
  ll.il.ms <- function(theta,struc,ch){
    beta.struc <- struc[[1]]
    delta.struc <- struc[[2]]
    kappa.struc <- struc[[3]]
    ageclasses <- 4
    timeclasses <- 2
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
    
    # struc.time <- list("age"=list(1:Time),"time"=timestruc(np),"state"=list(1:length(states)))
    struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
    # struc.state <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
    
    theta <- theta[-c(1:(timeclasses+ageclasses+2))]
    
    delt <- untrans(logistic(theta[1:2]),struc.const$age,delta.struc,struc.const$state)
    kap <- untrans(logistic(theta[3:4]),struc.const$age,kappa.struc,struc.const$state)
    gam <- untrans(logistic(theta[5]),struc.const$age,struc.const$time,struc.const$state)
    eps <- untrans(logistic(theta[6]),struc.const$age,struc.const$time,struc.const$state)
    # print(eps)
    # print(theta)
    # print(delt[,,1])
    # print(kap[,,1])
    
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
  strucs <- lapply(2:15,function(x)list(1:(x-1),x:(Time-1)) )
  cores <- length(strucs)
  timeclasses <- 2
  ageclasses <- 4
  tl <- 14 # length of theta/how many parameters
  par.name <- c(paste0("alpha",1:3),paste0("beta",10:11),paste0("rho",c(23,4,5)),paste0("delta",1:2),paste0("kappa",1:2),paste0("gamma"),paste0("epsilon"))
  par.group <- c(rep(c("alpha","beta","rho","delta","kappa","gamma","epsilon"),c(3,2,3,2,2,1,1)))
  par.fun <- function(mle){
    c(mle[1:5],logistic(mle[6:length(mle)]))
  }
  boundry.flag <- function(mle){
    mle1 <- mle[1:5]
    mle2 <- mle[6:length(mle)]
    rtn <- c(rep(NA,5),ifelse(mle2>0.95 | mle2<0.05,"YES","NO"))
    return(rtn)
  }
  hess.fun <- function(mat){
    if(is.singular.matrix(mat)){
      rtn <- rep(NA,tl)
    }else{
      rtn <- diag(solve(mat))
    }
    return(rtn)
  }
  for(b in 1:length(strucs)){ # b for beta
    for(d in 1:length(strucs)){ # d for delta
      # for(k in 1:length(strucs)){ # k for kappa
      # can use this loop if decide to parallelise over random starts
        starter <- rep(0,tl) # just zero starts for now
        op <- mclapply(1:length(strucs), function(k){
                    struc <- list(strucs[[b]],strucs[[d]],strucs[[k]])
                    ops <- nlm(f=ll.il.ms, # likelihood function
                               p=starter, # initial values
                               ch=ch,
                               struc=struc,
                               hessian=TRUE)
                    return(c(ops,"start"=list(starter)))
                  },mc.cores=cores)
      # }
      assign(paste0("df.mle",b,d),data.frame("par"=rep(par.name,length(strucs)),
                                             "value"=c(sapply(1:length(strucs),function(k) par.fun(op[[k]]$estimate) )),
                                             "convergence"=sapply(1:length(strucs),function(k) op[[k]]$code ),
                                             "beta.cut"=rep(b+1,each=tl*length(strucs)),
                                             "delta.cut"=rep(d+1,each=tl*length(strucs)),
                                             "kappa.cut"=rep((1:length(strucs))+1,each=tl),
                                             "boundry"=c(sapply(1:length(strucs),function(k) boundry.flag(par.fun(op[[k]]$estimate)) )),
                                             "var_hess"=c(sapply(1:length(strucs),function(k) hess.fun(op[[k]]$hessian)))))
      assign(paste0("df.aic",b,d),data.frame("aic"=c(sapply(1:length(strucs),function(k) 2*op[[k]]$minimum + 2*tl)),
                                             "beta.cut"=rep(b+1,length(strucs)),
                                             "delta.cut"=rep(d+1,length(strucs)),
                                             "kappa.cut"=(1:length(strucs))+1,
                                             "convergence"=c(sapply(1:length(strucs),function(k) op[[k]]$code)),
                                             "boundry"=sapply(1:length(strucs),function(k) sum(boundry.flag(par.fun(op[[k]]$estimate))=="YES",na.rm=TRUE)>0 )))
      print(c(b,d))
    }
  }
  
  for(i in 1:length(strucs)){
    assign(paste0("df.mle",i),do.call(rbind, lapply( paste0("df.mle",i, 1:length(strucs)) , get) ))
    assign(paste0("df.aic",i),do.call(rbind, lapply( paste0("df.aic",i, 1:length(strucs)) , get) ))
  }
  
  df.mle <- do.call(rbind, lapply( paste0("df.mle", 1:length(strucs)) , get) )
  df.aic <- do.call(rbind, lapply( paste0("df.aic", 1:length(strucs)) , get) )
  
  save(df.mle,df.aic,file="ms4.RData")
}