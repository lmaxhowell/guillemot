# trying to see if there is a step change in delta
file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed
load("guillemot.RData")

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
  x <- ifelse(state.depen==FALSE,22,31) # plus one compared to other scripts bc two delta parameters
  return(logit(rep(0.55,x+n)))
}

ll.il.delt <- function(theta,step.change,ageclasses=4,ch=ch){
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
  struc.age <- list("age"=list(1:step.change,(step.change+1):Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  # chosen so 1 is the smallest step.change can be and can only go up to (maximum(age)==Time)-1
  # struc.state <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
  
  theta <- theta[-c(1:(Time+ageclasses+1))]
  
  delt <- untrans(logistic(theta[1:2]),struc.age$age,struc.age$time,struc.age$state)
  kap <- untrans(logistic(theta[3]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(logistic(theta[4]),struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(theta[5]),struc.const$age,struc.const$time,struc.const$state)
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

ll.il.null <- function(theta,ageclasses=4,ch=ch){
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

Ages <- 1:14
timer(op.delt <- mclapply(Ages, function(x) optim(theta.s(4),ll.il.delt,
                                             step.change=x,
                                             ageclasses=4,ch=ch,
                                             control=list(fnscale=-1,maxit=1000),
                                             method="Nelder-Mead"),mc.cores=cores))
# Time difference of 52.37804 mins

# three in age classes EVEN THOUGH its four age classes as we want one less parameter to account for the one delta parameter
timer(op.null <- optim(theta.s(3),ll.il.null,
                 ageclasses=4,ch=ch,
                 control=list(fnscale=-1,maxit=1000),
                 method="Nelder-Mead"))

aics <- c(2*length(theta.s(3))-2*op.null$value,
          sapply(Ages,function(x) 2*length(theta.s(4))-2*op.delt[[x]]$value))
df.aic.delt <- data.frame("AIC"=aics,
                          "StepAge"=c(0,Ages),
                          "convergence"=c(op.null$convergence,
                                          sapply(Ages,function(x) op.delt[[x]]$convergence)),
                          "DeltaAIC"=aics-min(aics))
selected <- which.min(aics)
if(selected>1){ # its one with a step change in delta
  op.sel <- op.delt[[selected-1]] # minus one to reset which index of the list it is
}else if(selected==1){ # it selected the null model
  op.sel <- op.null
}
par.fun <- function(op.par){
  return(c(op.par[1:18],logistic(op.par[19:length(op.par)])))
}
par.name.fun <- function(op.par){
  if(length(op.par)==25){
    return(c(paste0("alpha",1:3),paste0("beta",1:15),paste0("rho",1:3),"delta","kappa","gamma","epsilon"))
  }else if(length(op.par)==26){
    return(c(paste0("alpha",1:3),paste0("beta",1:15),paste0("rho",1:3),"delta1","delta2","kappa","gamma","epsilon"))
  }
}
df.par <- data.frame("value"=par.fun(op.sel$par),
                     "par"=par.name.fun(op.sel$par),
                     "convergence"=rep(op.sel$convergence,length(op.sel$par)),
                     "StepAge"=rep(selected-1,length(op.sel$par)))

df.phi.delt <- data.frame("Time"=rep(2010:2024,4),"phi"=c(logistic(phi.trans(op.sel$par[1:3],op.sel$par[4:18]))[,1:4,1]),
                          "AgeClass"=as.factor(rep(1:4,each=(Time-1))))
ggplot(df.phi.delt,aes(Time,phi,col=AgeClass)) + geom_line() + geom_point() + theme_bw()


