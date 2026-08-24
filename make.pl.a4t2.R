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
  tb <- table(unlist(ch[,1:16]))
  
  Bs <- tb["LB"] + tb["L_B"]
  B_s <- tb["LB_"] + tb["L_B_"]
  Ls <- tb["LB"] + tb["LB_"]
  L_s <- tb["L_B"] + tb["L_B_"]
  
  gamm <- unname(Bs/(Bs+B_s))
  del <- unname(Ls/(Ls+L_s))
  
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
  
  theta.s <- function(n){ # want n to be desired number of age classes
    x <- 6
    return(logit(runif(x+n,0.1,0.9)))
  }
  
  ll.il.ms <- function(theta,beta.struc,ch){
    # beta.struc should be a list corresponding to
    # where the breaks in beta go e.g. list(1:10,11:15)
    # would be a change in beta at the 11 time point
    alpha <- theta[1:(4-1)]
    beta <- theta[4:(2+4-1)]
    rho.input <- theta[(2+4):(2+4+2)]
    phi <- logistic(phi.trans(alpha,beta,beta.struc))
    rho <- logistic(rho.trans(rho.input))
    
    Time <- ncol(ch)-1
    
    struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
    
    theta <- theta[-c(1:(2+4+2))]
    
    delt <- untrans(del,struc.const$age,struc.const$time,struc.const$state)
    kap <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
    gam <- untrans(gamm,struc.const$age,struc.const$time,struc.const$state)
    eps <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)

    psi <- make.psi(delt,kap,rho,gam,eps)

    ll <- il(ch,phi,psi)
    return(-ll)
  }
  
  ll.il.ms.fix <- function(theta,fix,at,beta.struc,ch){
    # beta.struc should be a list corresponding to
    # where the breaks in beta go e.g. list(1:10,11:15)
    # would be a change in beta at the 11 time point
    # alpha <- theta[1:(4-1)]
    alpha <- append(theta[1:2],fix,at-1)
    beta <- theta[3:4]
    rho.input <- theta[5:7]
    phi <- logistic(phi.trans(alpha,beta,beta.struc))
    rho <- logistic(rho.trans(rho.input))
    
    Time <- ncol(ch)-1
    
    struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
    
    theta <- theta[-c(1:(7))]
    
    delt <- untrans(del,struc.const$age,struc.const$time,struc.const$state)
    kap <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
    gam <- untrans(gamm,struc.const$age,struc.const$time,struc.const$state)
    eps <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
    
    psi <- make.psi(delt,kap,rho,gam,eps)
    
    ll <- il(ch,phi,psi)
    return(-ll)
  }
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  # ageclasses <- 4
  # timeclasses <- 2
  beta.struc <- list(1:(11-1),11:15)
  
  starter <- rep(0,10)
  op <- nlm(f=ll.il.ms, # likelihood function
            p=starter, # initial values
            beta.struc=beta.struc,
            ch=ch,
            hessian=TRUE)
  par.fun <- function(mle){
    c(mle[1:5],logistic(mle[6:length(mle)]))
  }
  
  lo <- 89
  cores <- 89
  over <- seq(-4,4,length.out=lo)
  
  op.a1 <- mclapply(1:lo,function(x) nlm(f=ll.il.ms.fix, # likelihood function
                                         p=starter[1:9], # initial values
                                         beta.struc=beta.struc,
                                         fix=over[x],
                                         at=1,
                                         ch=ch,
                                         hessian=TRUE) )
  print(1)
  op.a2 <- mclapply(1:lo,function(x) nlm(f=ll.il.ms.fix, # likelihood function
                                         p=starter, # initial values
                                         beta.struc=beta.struc,
                                         fix=over[x],
                                         at=2,
                                         ch=ch,
                                         hessian=TRUE) )
  print(2)
  op.a3 <- mclapply(1:lo,function(x) nlm(f=ll.il.ms.fix, # likelihood function
                                         p=starter, # initial values
                                         beta.struc=beta.struc,
                                         fix=over[x],
                                         at=3,
                                         ch=ch,
                                         hessian=TRUE) )
  print(3)
  
  pl.a1 <- sapply(1:lo,function(x) op.a1[[x]]$minimum )
  pl.a2 <- sapply(1:lo,function(x) op.a2[[x]]$minimum )
  pl.a3 <- sapply(1:lo,function(x) op.a3[[x]]$minimum )
  
  con.a1 <- sapply(1:lo,function(x) op.a1[[x]]$code )
  con.a2 <- sapply(1:lo,function(x) op.a2[[x]]$code )
  con.a3 <- sapply(1:lo,function(x) op.a3[[x]]$code )
  
  df <- data.frame("over"=rep(over,3),
                   "value"=c(pl.a1,pl.a2,pl.a3),
                   "par"=rep(c("alpha1","alpha2","alpha3"),each=lo),
                   "convergence"=c(con.a1,con.a2,con.a3))
  df.true <- data.frame("par"=c(paste0("alpha",1:3),paste0("beta",1:2),paste0("rho",c(23,4,5)),"kappa","epsilon"),
                        "mle"=par.fun(op$estimate),
                        "convergence"=rep(op$code,10),
                        "value"=rep(op$minimum,10))
  save(df,df.true,file="pl.phi.RData")
  
}