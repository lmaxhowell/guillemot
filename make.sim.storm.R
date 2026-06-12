# file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
# file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
# file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
# lapply(file.names, source) # source all functions needed
# save(list=ls(),file="guill.func.RData")

l <- as.numeric(commandArgs(trailingOnly=TRUE)[1]) #needed to queue jobs on storm
if(l==1){
  load("guill.func.RData")
  library(plyr, include.only = c("count"))
  library(parallel)
  states <- c("N","N_","B1","LB","L_B","LB_","L_B_","S")
  ni <- c(48, 30, 31, 65, 42, 49, 64, 48, 61, 75, 65, 87, 62, 73, 77)
  Time <- length(ni)+1
  ni2 <- list()
  for(i in 1:(length(ni)-1)){
    ni2[[i]] <- rep(1,ni[i])
  }
  
  struc <- list("phi"=list("age"=list(1,2:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                 "delt"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                 "kap"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                 "rho"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                 "gam"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                 "eps"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
  struc2 <- list("phi"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                "delt"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                "kap"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                "rho"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                "gam"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                "eps"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
  
  
  theta <- logit(c(0.3,0.7,0.3,0.2,0.5,0.8,0.6))
  theta2 <- logit(c(0.7,0.3,0.2,0.5,0.8,0.6))
  # phi3 <- untrans2(logistic(theta[1:2]),struc3$phi$age,struc3$phi$time,struc3$phi$state)
  # delt3 <- untrans2(logistic(theta[3]),struc3$delt$age,struc3$delt$time,struc3$delt$state)
  # kap3 <- untrans2(logistic(theta[4]),struc3$kap$age,struc3$kap$time,struc3$kap$state)
  # rho3 <- untrans2(logistic(theta[5]),struc3$rho$age,struc3$rho$time,struc3$rho$state)
  # gam3 <- untrans2(logistic(theta[6]),struc3$gam$age,struc3$gam$time,struc3$gam$state)
  # eps3 <- untrans2(logistic(theta[7]),struc3$eps$age,struc3$eps$time,struc3$eps$state)
  
  
  
  n <- 200
  set.seed(722461813)
  seeds <- sample(1:.Machine$integer.max,n)
  cores <- 50
  sim.dat <- mclapply(1:n,function(x) dat.sim.wrap(theta,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
                                                    eps.ind=7,struc,ni2,seeds[x]),mc.cores=cores)
  sim.dat2 <- mclapply(1:n,function(x) dat.sim.wrap(theta2,phi.ind=1,delt.ind=2,kap.ind=3,rho.ind=4,gam.ind=5,
                                                    eps.ind=6,struc2,ni2,seeds[x]),mc.cores=cores)
  
  timer(op.n <- mclapply(1:n,function(x) optim(theta,ll.il,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
                                                eps.ind=7,struc=struc,ch=sim.dat[[x]],control=list(fnscale=-1),
                                                method="BFGS"),mc.cores=cores))
  timer(op.n2 <- mclapply(1:n,function(x) optim(theta2,ll.il,phi.ind=1,delt.ind=2,kap.ind=3,rho.ind=4,gam.ind=5,
                                                eps.ind=6,struc=struc2,ch=sim.dat2[[x]],control=list(fnscale=-1),
                                                method="BFGS"),mc.cores=cores))
  
  df.sim <- data.frame("par"=rep(c("phi1","phi2","delta","kappa","rho","gamma","epsilon"),n),
                        "MLE"=c(sapply(1:n,function(x) logistic(op.n[[x]]$par))),
                        "convergence"=rep(sapply(1:n,function(x) op.n[[x]]$convergence),each=7))
  df.sim2 <- data.frame("par"=rep(c("phi","delta","kappa","rho","gamma","epsilon"),n),
                        "MLE"=c(sapply(1:n,function(x) logistic(op.n2[[x]]$par))),
                        "convergence"=rep(sapply(1:n,function(x) op.n2[[x]]$convergence),each=6))
  
  save(df.sim,df.sim2,file="guill.sim.RData")
}