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
  
  cutoff <- 7
  struclist <- lapply(1:cutoff,function(x) c(as.list(1:x),list((x+1):Time)))
  struc.s <- lapply(1:cutoff,function(x){
    list("phi"=list("age"=struclist[[x]],"time"=list(1:Time),"state"=list(1:length(states))),
                 "delt"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                 "kap"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                 "rho"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                 "gam"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
                 "eps"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
  })
  
  theta <- lapply(1:cutoff,function(x) logit(sample(seq(0.1,0.9,0.05),5+length(struclist[[x]]),replace = TRUE)))
  
  n <- 200
  set.seed(722461813)
  seeds <- sample(1:.Machine$integer.max,n)
  cores <- 89
  
  sim.dat <- lapply(1:cutoff,function(x) mclapply(1:n,function(y) dat.sim.wrap(theta[[x]],
                                                    phi.ind=1:length(struclist[[x]]),
                                                    delt.ind=(length(struclist[[x]])+1),
                                                    kap.ind=(length(struclist[[x]])+2),
                                                    rho.ind=(length(struclist[[x]])+3),
                                                    gam.ind=(length(struclist[[x]])+4),
                                                    eps.ind=(length(struclist[[x]])+5),
                                                    struc.s[[x]],ni2,seeds[y]),mc.cores=cores))
  
  op.n <- lapply(1:cutoff,function(x) mclapply(1:n,function(y) optim(theta[[x]],ll.il,
                                                    phi.ind=1:length(struclist[[x]]),
                                                    delt.ind=(length(struclist[[x]])+1),
                                                    kap.ind=(length(struclist[[x]])+2),
                                                    rho.ind=(length(struclist[[x]])+3),
                                                    gam.ind=(length(struclist[[x]])+4),
                                                    eps.ind=(length(struclist[[x]])+5),
                                                    struc=struc.s[[x]],ch=sim.dat[[x]][[y]],control=list(fnscale=-1,maxit=1000),
                                                method="BFGS"),mc.cores=cores))
  
  df.sim <- data.frame("par"=unlist(sapply(1:cutoff,function(x) rep(c(paste0("phi",1:(x+1)),"delta","kappa","rho","gamma","epsilon"),n))),
                       "MLE"=unlist(sapply(1:cutoff,function(x) c(sapply(1:n,function(y) logistic(op.n[[x]][[y]]$par))))),
                       "convergence"=unlist(sapply(1:cutoff,function(x) rep(sapply(1:n,function(y) op.n[[x]][[y]]$convergence),each=length(c(paste0("phi",1:(x+1)),"delta","kappa","rho","gamma","epsilon"))))),
                       "ageclasses"=unlist(sapply(1:cutoff,function(x) rep(x+1,n*(length(struclist[[x]])+5)))))

  df.true <- data.frame("par"=unlist(sapply(1:cutoff,function(x) c(paste0("phi",1:(x+1)),"delta","kappa","rho","gamma","epsilon"))),
                        "MLE"=unlist(lapply(1:cutoff,function(x) logistic(theta[[x]]))),
                        "ageclasses"=unlist(sapply(1:cutoff,function(x) rep(x+1,(length(struclist[[x]])+5)))))
  
  save(df.sim,df.true,file="phi.sim.RData")
}