file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed

load("matrix_AUK.RData")
load("age_mat.RData")

states <- c("N", "B1","LB","L_B","LB_","L_B_","S")
states_code <- c(1,2,3,5,4,6) # no skipped code needed so this vector is one shorter length than states
ch <- matrix_AUK
for(i in 1:length(states_code)){
  ch[which(ch==states_code[i],arr.ind=TRUE)] <- states[i]
}


# parameters are phi,delta,kappa,rho,gamma

# Time <- 7
# Ages <- Time
# Sexs <- 2

# phi <- array(dim=c(Time,Ages,Sexs))
# delta <- array(dim=c(Time,Ages,Sexs))
# kap <- array(dim=c(Time,Ages,Sexs))
# rho <- array(dim=c(Time,Ages,Sexs))
# gam <- array(dim=c(Time,Ages,Sexs))

phi <- array(rep(c(0.3,0.7),c(1,6)),dim=c(length(states)))
delt <- array(0.3,dim=c(length(states)))
kap <- array(0.2,dim=c(length(states)))
rho <- array(0.5,dim=c(length(states)))
gam <- array(0.8,dim=c(length(states)))

test <- rbind(c("N","N","N","B1","LB","LB","0","L_B_","0"),
              c("N","N","N","B1","LB","LB","0","L_B_","0"))
il(test,phi,psi)

# calculate raw likelihood to test against
2*(log(phi[1]*(1-rho[1])) + # first N->N
  log(phi[1]*(1-rho[1])) + # second N->N
  log(phi[1]*rho[1]) + # third N->B1
  log(phi[2]*(1-kap[2])*(1-delt[2])*gam[2]) + # fourth B1->LB
  log(phi[3]*(1-kap[3])*(1-delt[3])*gam[3]) + # fifth LB->LB
  log(phi[3]*kap[3]) + # sixth LB->S
  log(phi[7]*delt[7]*(1-gam[7])) + # seventh S->L_B_
  log(phi[6]*kap[6]*(1-phi[7])+(1-phi[6]))) # eighth L_B_->0 could have skipped or died?
il(test,phi,psi) # slightly over so its double counting

test_il <- c(log(phi[1]*(1-rho[1])), # first N->N
              log(phi[1]*(1-rho[1])), # second N->N
              log(phi[1]*rho[1]), # third N->B1
              log(phi[2]*(1-kap[2])*(1-delt[2])*gam[2]), # fourth B1->LB
              log(phi[3]*(1-kap[3])*(1-delt[3])*gam[3]), # fifth LB->LB
              log(phi[3]*kap[3]), # sixth LB->S
              log(phi[7]*delt[7]*(1-gam[7])), # seventh S->L_B_
              log(phi[6]*kap[6]*(1-phi[7])+(1-phi[6]))) # eighth L_B_->0 could have skipped or died?
sum(test_il)
log(Pr_r0(6,phi,psi))
test_il[8]


psi <- make.psi(delt,kap,rho,gam)

il(ch,phi,psi)

theta <- logit(c(0.3,0.7,0.3,0.2,0.5,0.8))
struc <- list("phi"=list(1,2:length(states)),
              "delt"=list(1:length(states)),
              "kap"=list(1:length(states)),
              "rho"=list(1:length(states)),
              "gam"=list(1:length(states)))
ll.il(theta,1:2,3,4,5,6,struc,ch)

timer(op <- optim(theta,ll.il,phi.ind=1:2,delt.ind=3,kap.ind=4,
            rho.ind=5,gam.ind=6,struc=struc,ch=ch,
            control=list(fnscale=-1),hessian=TRUE)) # Time difference of 3.178106 mins, convergence 1
timer(op2 <- optim(theta,ll.il,phi.ind=1:2,delt.ind=3,kap.ind=4,
                  rho.ind=5,gam.ind=6,struc=struc,ch=ch,
                  control=list(fnscale=-1),method="BFGS",hessian=TRUE)) # Time difference of 2.249783 mins, convergence 0
op$convergence
c("phi1","phiA","delta","kappa","rho","gamma")
logistic(op$par)
op$hessian
diag(solve(-op$hessian))

# adding time dependance to function
Time <- ncol(ch)
struc <- list("phi"=list("state"=list(1,2:length(states)),"time"=list(1:Time)),
              "delt"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "kap"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "rho"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "gam"=list("state"=list(1:length(states)),"time"=list(1:Time)))
ll.il(theta,1:2,3,4,5,6,struc,ch)
theta2 <- logit(c(0.3,0.7,0.3,0.2,0.7,0.5,0.8))
struc2 <- list("phi"=list("state"=list(1,2:length(states)),"time"=list(1:Time)),
              "delt"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "kap"=list("state"=list(1:length(states)),"time"=list(1:8,9:Time)),
              "rho"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "gam"=list("state"=list(1:length(states)),"time"=list(1:Time)))
ll.il(theta2,1:2,3,4:5,6,7,struc2,ch)
timer(op3 <- optim(theta2,ll.il,phi.ind=1:2,delt.ind=3,kap.ind=4:5,
                   rho.ind=6,gam.ind=7,struc=struc2,ch=ch,
                   control=list(fnscale=-1),method="BFGS",hessian=TRUE)) # Time difference of 1.921875 mins, convergence 0
op3$convergence
logistic(op3$par)
diag(solve(-op3$hessian))
ci(op3$par,diag(solve(-op3$hessian)))

struc3 <- list("phi"=list("state"=list(1,2:length(states)),"time"=list(1:Time)),
               "delt"=list("state"=list(1:length(states)),"time"=list(1:Time)),
               "kap"=list("state"=list(1:length(states)),"time"=list(1:Time)),
               "rho"=list("state"=list(1:length(states)),"time"=list(1:Time)),
               "gam"=list("state"=list(1:length(states)),"time"=list(1:Time)))
timer(op4 <- optim(theta,ll.il,phi.ind=1:2,delt.ind=3,kap.ind=4,
                   rho.ind=5,gam.ind=6,struc=struc3,ch=ch,
                   control=list(fnscale=-1),method="BFGS",hessian=TRUE)) # Time difference of 1.506441 mins, convergence 0
op4$convergence
c("phi1","phiA","delta","kappa","rho","gamma")
logistic(op4$par)

#################################### New il to add in Nbar state
rm(list=ls())
file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
# file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed
states <- c("N","N_","B1","LB","L_B","LB_","L_B_","S")
states_code <- c(1,13,2,3,5,4,6)
# no skipped code needed so this vector is one shorter length than states
# the 13 is so no element of nch will match and get N_ and keep the indicies correct
load("matrix_AUK.RData")
nch <- matrix_AUK
Time <- ncol(nch)
for(i in 1:length(states_code)){
  nch[which(nch==states_code[i],arr.ind=TRUE)] <- states[i]
}
for(i in 1:nrow(nch)){
  # want to add Nbar as a state between any Ns or until N becomes B1
  # and want to add S as a state between the breeding states
  wn <- which(nch[i,]=="N")
  wb1 <- which(nch[i,]=="B1")
  wnlss <- wn[wn<wb1] # which non breeding states are before the first breeding state (as opposed to the non breeding due to skipping)
  if(length(wnlss)==0){
    if(length(wn)>=2){
      for(j in 1:(Time-1)){
        if(j>min(wn) & j<max(wn) & nch[i,j]=="0"){
          # newch[j] <- "N_"
          nch[i,j] <- "N_"
        }
      }
    }
  }else{
    if(length(wnlss)>0){
      for(j in 1:(Time-1)){
        if(j>min(wn) & j<wb1 & nch[i,j]=="0"){
          # newch[j] <- "N_"
          nch[i,j] <- "N_"
        }
      }
    }
  }
  wbb <- which(nch[i,] %in% c("LB","L_B","LB_","L_B_"))
  if(length(wbb)>1){
    for(j in 1:(length(wbb)-1)){
      if((wbb[j+1]-wbb[j])>1){
        nch[i,(wbb[j]+1):(wbb[j+1]-1)] <- "S"
      }
    }
  }
}

struc <- list("phi"=list("state"=list(1,2:length(states)),"time"=list(1:Time)),
               "delt"=list("state"=list(1:length(states)),"time"=list(1:Time)),
               "kap"=list("state"=list(1:length(states)),"time"=list(1:Time)),
               "rho"=list("state"=list(1:length(states)),"time"=list(1:Time)),
               "gam"=list("state"=list(1:length(states)),"time"=list(1:Time)),
               "ups"=list("state"=list(1:length(states)),"time"=list(1:Time)))
delt <- untrans(0.3,struc$delt$time,struc$delt$state)
kap <- untrans(0.2,struc$kap$time,struc$kap$state)
rho <- untrans(0.5,struc$rho$time,struc$rho$state)
gam <- untrans(0.8,struc$gam$time,struc$gam$state)
ups <- untrans(0.6,struc$ups$time,struc$ups$state)
make.psi(delt,kap,rho,gam,ups)


theta <- logit(c(0.3,0.7,0.3,0.2,0.5,0.8,0.6))
ll.il(theta,1:2,3,4,5,6,7,struc,nch)

test <- rbind(c("N","N","N","B1","LB","LB","S","L_B_","0"),
              c("N","N","N","B1","LB","LB","S","L_B_","0"))

# calculate raw likelihood to test against
test_il <- c(log(phi[1,1]*(1-rho[1,1])*ups[1,1]), # first N->N
             log(phi[1,2]*(1-rho[1,2])*ups[1,2]), # second N->N
             log(phi[1,3]*rho[1,3]), # third N->B1
             log(phi[2,4]*(1-kap[2,4])*(1-delt[2,4])*gam[2,4]), # fourth B1->LB
             log(phi[3,5]*(1-kap[3,5])*(1-delt[3,5])*gam[3,5]), # fifth LB->LB
             log(phi[3,6]*kap[3,6]), # sixth LB->S
             log(phi[7,7]*delt[7,7]*(1-gam[7,7])), # seventh S->L_B_
             log(phi[6,8]*kap[6,8]*(1-phi[7,8])+(1-phi[6,8]))) # eighth L_B_->0 could have skipped or died?
2*sum(test_il)
il(test,phi,psi)


library(plyr, include.only = c("count"))
uch <- count(nch)
# want to see how many individuals are ringed at each time
ni <- rep(0,ncol(uch)-1) # number of individuals
for(i in 1:nrow(uch)){
  ni[which(uch[i,1:16]!=0)[1]] <- ni[which(uch[i,1:16]!=0)[1]] + uch[i,17]
}
mean(ni) # 54.8 -> 55 individuals roughly each time

phi <- untrans(c(0.3,0.7),struc$phi$time,struc$phi$state)
sim <- dat.simulate(phi,delt,kap,rho,gam,ups,ni[1:15],722461813)
ll.il(theta,1:2,3,4,5,6,7,struc,sim)

op <- optim(theta,ll.il,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
            ups.ind=7,struc=struc,ch=sim,control=list(fnscale=-1),method="BFGS")
op$convergence
logistic(op$par)
op$value

library(parallel)
n <- 20
set.seed(722461813)
seeds <- sample(1:.Machine$integer.max,n)
cores <- 2
sim.dat <- mclapply(1:n,function(x) dat.simulate(phi,delt,kap,rho,gam,ups,ni[1:15],seeds[x]),mc.cores=cores)

timer(op.n <- mclapply(1:n,function(x) optim(theta,ll.il,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
                                       ups.ind=7,struc=struc,ch=sim.dat[[x]],control=list(fnscale=-1),
                                       method="BFGS"),mc.cores=2))

df.sim <- data.frame("par"=rep(c("phi1","phi2","delta","kappa","rho","gamma","upsilon"),n),
                     "MLE"=c(sapply(1:n,function(x) logistic(op.n[[x]]$par))),
                     "convergence"=rep(sapply(1:n,function(x) op.n[[x]]$convergence),each=7))
library(ggplot2)
ggplot(df.sim,aes(par,MLE)) + geom_boxplot()
df.true <- data.frame("par"=c("phi1","phi2","delta","kappa","rho","gamma","upsilon"),
                      "MLE"=logistic(theta))
ggplot(df.sim,aes(par,MLE,col=par)) + geom_boxplot() + geom_point(data=df.true,aes(par,MLE),shape=5,col="black")

####################
# trying again but with constant phi
####################
struc2 <- list("phi"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "delt"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "kap"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "rho"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "gam"=list("state"=list(1:length(states)),"time"=list(1:Time)),
              "ups"=list("state"=list(1:length(states)),"time"=list(1:Time)))
phi2 <- untrans(c(0.6),struc2$phi$time,struc2$phi$state)
theta2 <- logit(c(0.6,0.3,0.2,0.5,0.8,0.6))

sim.dat2 <- mclapply(1:n,function(x) dat.simulate(phi2,delt,kap,rho,gam,ups,ni[1:15],seeds[x]),mc.cores=cores)

timer(op.n2 <- mclapply(1:n,function(x) optim(theta2,ll.il,phi.ind=1,delt.ind=2,kap.ind=3,rho.ind=4,gam.ind=5,
                                             ups.ind=6,struc=struc2,ch=sim.dat2[[x]],control=list(fnscale=-1),
                                             method="BFGS"),mc.cores=2))

df.sim2 <- data.frame("par"=rep(c("phi","delta","kappa","rho","gamma","upsilon"),n),
                     "MLE"=c(sapply(1:n,function(x) logistic(op.n2[[x]]$par))),
                     "convergence"=rep(sapply(1:n,function(x) op.n2[[x]]$convergence),each=6))
ggplot(df.sim2,aes(par,MLE)) + geom_boxplot()
df.true2 <- data.frame("par"=c("phi","delta","kappa","rho","gamma","upsilon"),
                      "MLE"=logistic(theta2))
ggplot(df.sim2,aes(par,MLE,col=par)) + geom_boxplot() + geom_point(data=df.true2,aes(par,MLE),shape=5,col="black")

####################################
# Adding in age dependance
####################################
struc3 <- list("phi"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1,2:length(states))),
               "delt"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "kap"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "rho"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "gam"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "ups"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
phi3 <- untrans2(logistic(theta[1:2]),struc3$phi$age,struc3$phi$time,struc3$phi$state)
delt3 <- untrans2(logistic(theta[3]),struc3$delt$age,struc3$delt$time,struc3$delt$state)
kap3 <- untrans2(logistic(theta[4]),struc3$kap$age,struc3$kap$time,struc3$kap$state)
rho3 <- untrans2(logistic(theta[5]),struc3$rho$age,struc3$rho$time,struc3$rho$state)
gam3 <- untrans2(logistic(theta[6]),struc3$gam$age,struc3$gam$time,struc3$gam$state)
ups3 <- untrans2(logistic(theta[7]),struc3$ups$age,struc3$ups$time,struc3$ups$state)

psi3 <- make.psi2(delt3,kap3,rho3,gam3,ups3)
nch3 <- cbind(nch,"age"=rep(1,nrow(nch)))

ll.il(theta,1:2,3,4,5,6,7,struc,nch)
ll.il2(theta,1:2,3,4,5,6,7,struc3,nch3)
# they produce the same number so adding in age worked!

ni2 <- list()
for(i in 1:(length(ni)-1)){
  ni2[[i]] <- rep(1,ni[i])
}
sim2 <- dat.simulate2(phi3,delt3,kap3,rho3,gam3,ups3,ni2,722461813)
ll.il2(theta,1:2,3,4,5,6,7,struc3,sim2)
op2 <- optim(theta,ll.il2,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
            ups.ind=7,struc=struc3,ch=sim2,control=list(fnscale=-1),method="BFGS")
op2$convergence
logistic(op2$par)


n <- 20
set.seed(722461813)
seeds <- sample(1:.Machine$integer.max,n)
cores <- 2
sim.dat2 <- mclapply(1:n,function(x) dat.simulate2(phi3,delt3,kap3,rho3,gam3,ups3,ni2,seeds[x]),mc.cores=cores)

timer(op.n2 <- mclapply(1:n,function(x) optim(theta,ll.il2,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
                                             ups.ind=7,struc=struc3,ch=sim.dat2[[x]],control=list(fnscale=-1),
                                             method="BFGS"),mc.cores=cores))

df.sim2 <- data.frame("par"=rep(c("phi1","phi2","delta","kappa","rho","gamma","upsilon"),n),
                     "MLE"=c(sapply(1:n,function(x) logistic(op.n2[[x]]$par))),
                     "convergence"=rep(sapply(1:n,function(x) op.n2[[x]]$convergence),each=7))
library(ggplot2)
ggplot(df.sim2,aes(par,MLE)) + geom_boxplot()
df.true <- data.frame("par"=c("phi1","phi2","delta","kappa","rho","gamma","upsilon"),
                      "MLE"=logistic(theta))
ggplot(df.sim2,aes(par,MLE,col=par)) + geom_boxplot() + geom_point(data=df.true,aes(par,MLE),shape=5,col="black")


####### quick storm sim
load("guill.sim.RData")
ggplot(df.sim3,aes(par,MLE,col=par)) + geom_boxplot() + geom_point(data=df.true,aes(par,MLE),shape=5,col="black")

###### trying other structures
struc4 <- list("phi"=list("age"=list(1,2:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "delt"=list("age"=list(1:Time),"time"=list(1:8,9:Time),"state"=list(1:length(states))),
               "kap"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "rho"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "gam"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "ups"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
theta4 <- logit(c(0.3,0.7,0.2,0.6,0.2,0.7,0.5,0.8))
sim.dat4 <- dat.sim.wrap(theta4,1:2,3:4,5,6,7,8,struc4,ni2,1330080468)
op4 <- optim(theta4,ll.il2,phi.ind=1:2,delt.ind=3:4,kap.ind=5,rho.ind=6,gam.ind=7,
             ups.ind=8,struc=struc4,ch=sim.dat4,control=list(fnscale=-1),method="BFGS")
op4$convergence
logistic(op4$par)


############## trying more direct calculations of the likelihood
phi4 <- untrans2(logistic(theta4[1:2]),struc4$phi$age,struc4$phi$time,struc4$phi$state)
delt4 <- untrans2(logistic(theta4[3:4]),struc4$delt$age,struc4$delt$time,struc4$delt$state)
kap4 <- untrans2(logistic(theta4[5]),struc4$kap$age,struc4$kap$time,struc4$kap$state)
rho4 <- untrans2(logistic(theta4[6]),struc4$rho$age,struc4$rho$time,struc4$rho$state)
gam4 <- untrans2(logistic(theta4[7]),struc4$gam$age,struc4$gam$time,struc4$gam$state)
ups4 <- untrans2(logistic(theta4[8]),struc4$ups$age,struc4$ups$time,struc4$ups$state)

sim.dat4[418,]
# if state comes from column
test41 <- c(phi4[10,1,1]*(1-rho4[10,1,1])*ups4[10,1,1], # N -> N
            phi4[11,2,3]*rho4[11,2,3], # N -> B1
            phi4[12,3,7]*(1-kap4[12,3,7])*delt4[12,3,7]*(1-gam4[12,3,7]), # B1 -> L_B_
            phi4[13,4,8]*kap4[13,4,8], # L_B_ -> S
            phi4[14,5,6]*(1-delt4[14,5,6])*(1-gam4[14,5,6]), # S -> LB_
            (1-phi4[15,6,6])+phi4[15,6,6]*kap4[15,6,8]*(1-phi4[16,7,6])) # LB_ -> 0
# if state comes from row
test41 <- c(phi4[10,1,1]*(1-rho4[10,1,1])*ups4[10,1,1], # N -> N
            phi4[11,2,1]*rho4[11,2,1], # N -> B1
            phi4[12,3,3]*(1-kap4[12,3,3])*delt4[12,3,3]*(1-gam4[12,3,3]), # B1 -> L_B_
            phi4[13,4,7]*kap4[13,4,7], # L_B_ -> S
            phi4[14,5,8]*(1-delt4[14,5,8])*(1-gam4[14,5,8]), # S -> LB_
            (1-phi4[15,6,6])+phi4[15,6,6]*kap4[15,6,8]*(1-phi4[16,7,8])) # LB_ -> 0

sim.dat4[455,]
# if state comes from column
test42 <- c(phi4[10,1,2]*(1-rho4[10,1,2])*(1-ups4[10,1,2]), # N -> N_ state is 2 as its the column/the state its going to
            phi4[11,2,3]*rho4[11,2,3], # N_ -> B1
            phi4[12,3,8]*kap4[12,3,8], # B1 -> S
            phi4[13,4,7]*delt4[13,4,7]*(1-gam4[13,4,7]), # S -> L_B_
            phi4[14,5,4]*(1-kap4[14,5,4])*(1-delt4[14,5,4])*gam4[14,5,4], # L_B_ -> LB
            phi4[15,6,4]*(1-kap4[15,6,4])*(1-delt4[15,6,4])*gam4[15,6,4]) # LB -> LB
# if state comes from row
test42 <- c(phi4[10,1,1]*(1-rho4[10,1,1])*(1-ups4[10,1,1]), # N -> N_
            phi4[11,2,2]*rho4[11,2,2], # N_ -> B1
            phi4[12,3,3]*kap4[12,3,3], # B1 -> S
            phi4[13,4,8]*delt4[13,4,8]*(1-gam4[13,4,8]), # S -> L_B_
            phi4[14,5,7]*(1-kap4[14,5,7])*(1-delt4[14,5,7])*gam4[14,5,7], # L_B_ -> LB
            phi4[15,6,4]*(1-kap4[15,6,4])*(1-delt4[15,6,4])*gam4[15,6,4]) # LB -> LB

sum(log(test41),log(test42))

ll.il2(theta4,1:2,3:4,5,6,7,8,struc4,sim.dat4[c(418,455),])

sim.dat4[761,]
test43 <- sum(c((1-phi4[14,1,1]), # it can die
                phi4[14,1,1]*(1-rho4[14,1,1])*(1-ups4[14,1,1])*(1-phi4[15,2,2]), # it can survive, not recruit and die
                phi4[14,1,1]*(1-rho4[14,1,1])*(1-ups4[14,1,1])*phi4[15,2,2]*(1-rho4[15,2,2])*(1-ups4[15,2,2])))  # it can survive, not recruit and then survive and not recruit

ll.il2(theta4,1:2,3:4,5,6,7,8,struc4,sim.dat4[c(761,761),])/2

sum(log(test41),log(test42),log(test43))
ll.il2(theta4,1:2,3:4,5,6,7,8,struc4,sim.dat4[c(418,455,761),])


####################################
# Changing states!
####################################
rm(list=ls())
file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed
states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
states_code <- c(1,13,2,3,5,4,6)
# no skipped code needed so this vector is one shorter length than states
# the 13 is so no element of nch will match and get E and keep the indicies correct
load("matrix_AUK.RData")
nch <- matrix_AUK
Time <- ncol(nch)
for(i in 1:length(states_code)){
  nch[which(nch==states_code[i],arr.ind=TRUE)] <- states[i]
}
for(i in 1:nrow(nch)){
  # want to add Nbar as a state between any Ns or until N becomes B1
  # and want to add S as a state between the breeding states
  wn <- which(nch[i,]=="N")
  wb1 <- which(nch[i,]=="B1")
  wnlss <- wn[wn<wb1] # which non breeding states are before the first breeding state (as opposed to the non breeding due to skipping)
  if(length(wnlss)==0){
    if(length(wn)>=2){
      for(j in 1:(Time-1)){
        if(j>min(wn) & j<max(wn) & nch[i,j]=="0"){
          nch[i,j] <- "E"
        }
      }
    }
  }else{
    if(length(wnlss)>0){
      for(j in 1:(Time-1)){
        if(j>min(wn) & j<wb1 & nch[i,j]=="0"){
          nch[i,j] <- "E"
        }
      }
    }
  }
  wbb <- which(nch[i,] %in% c("B1","LB","L_B","LB_","L_B_"))
  if(length(wbb)>1){
    for(j in 1:(length(wbb)-1)){
      if((wbb[j+1]-wbb[j])>1){
        nch[i,(wbb[j]+1):(wbb[j+1]-1)] <- "S"
      }
    }
  }
}

theta <- logit(c(0.3,0.7,0.3,0.2,0.5,0.8,0.6))
struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1,2:length(states))),
               "delt"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "kap"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "rho"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "gam"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "eps"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))
phi <- untrans(logistic(theta[1:2]),struc$phi$age,struc$phi$time,struc$phi$state)
delt <- untrans(logistic(theta[3]),struc$delt$age,struc$delt$time,struc$delt$state)
kap <- untrans(logistic(theta[4]),struc$kap$age,struc$kap$time,struc$kap$state)
rho <- untrans(logistic(theta[5]),struc$rho$age,struc$rho$time,struc$rho$state)
gam <- untrans(logistic(theta[6]),struc$gam$age,struc$gam$time,struc$gam$state)
eps <- untrans(logistic(theta[7]),struc$eps$age,struc$eps$time,struc$eps$state)

psi <- make.psi(delt,kap,rho,gam,eps)

funs <-  sapply(14:1,function(x){
              time.start <- Sys.time()
              res <- Pr_r0(1,x,1,phi,psi)
              time.end <- Sys.time()
              tm <- difftime(time.end,time.start,units="secs")
              print(c(x,tm))
              return(c(res,tm))})

qdf <- as.data.frame(matrix(c(14.000000000,0.001687527,
                       13.000000000,0.004917145,
                       12.0000000,0.0585835,
                       11.0000000,0.2531517,
                       10.00000,1.71282,
                       9.00000,12.48036,
                       8.000000,101.4336,
                       7.0000,798.5432),ncol=2,byrow = TRUE))
colnames(qdf) <- c("x","time")
ggplot(qdf,aes(x,time)) + geom_point()
qdf <- cbind(as.data.frame(t(funs)),"V3"=14:1)
colnames(qdf) <- c("value","time","x")
ggplot(qdf,aes(x,time)) + geom_point()

uch <- count(nch)
View(uch)


ni <- rep(0,ncol(uch)-1) # number of individuals
for(i in 1:nrow(uch)){
  ni[which(uch[i,1:16]!=0)[1]] <- ni[which(uch[i,1:16]!=0)[1]] + uch[i,17]
}
ni2 <- list()
for(i in 1:(length(ni)-1)){
  ni2[[i]] <- rep(1,ni[i])
}
sim.dat <- dat.sim.wrap(theta,1:2,3,4,5,6,7,struc,ni2,1330080468)

ll.il(theta,1:2,3,4,5,6,7,struc,sim.dat)

mch <- cbind(nch,"age"=rep(1,nrow(nch)))
ll.il(theta,1:2,3,4,5,6,7,struc,mch)

op <- optim(theta,ll.il,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
            eps.ind=7,struc=struc,ch=sim.dat,control=list(fnscale=-1),method="BFGS")
op$convergence
logistic(op$par)

library(parallel)
n <- 20
set.seed(722461813)
seeds <- sample(1:.Machine$integer.max,n)
cores <- 2
sim.dat.n <- mclapply(1:n,function(x) dat.simulate2(phi,delt,kap,rho,gam,eps,ni2,seeds[x]),mc.cores=cores)

timer(op.n <- mclapply(1:n,function(x) optim(theta,ll.il,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
                                             eps.ind=7,struc=struc,ch=sim.dat.n[[x]],control=list(fnscale=-1),
                                             method="BFGS"),mc.cores=cores))

df.sim <- data.frame("par"=rep(c("phi1","phi2","delta","kappa","rho","gamma","epsilon"),n),
                     "MLE"=c(sapply(1:n,function(x) logistic(op.n[[x]]$par))),
                     "convergence"=rep(sapply(1:n,function(x) op.n[[x]]$convergence),each=7))
# ggplot(df.sim,aes(par,MLE)) + geom_boxplot()
df.true <- data.frame("par"=c("phi1","phi2","delta","kappa","rho","gamma","epsilon"),
                      "MLE"=logistic(theta))
ggplot(df.sim,aes(par,MLE,col=par)) + geom_boxplot() + geom_point(data=df.true,aes(par,MLE),shape=5,col="black")

nch[7,]
Chi(1,14,1,phi,psi)
sum((1-phi[14,1,1]),
    phi[14,1,1]*psi[1,2,14,1]*(1-phi[15,2,2]),
    phi[14,1,1]*psi[1,2,14,1]*phi[15,2,2]*psi[2,2,15,2])
uch[153,]
Chi(1,13,1,phi,psi)
sum((1-phi[13,1,1]),
    phi[13,1,1]*psi[1,2,13,1]*(1-phi[14,2,2]),
    phi[13,1,1]*psi[1,2,13,1]*phi[14,2,2]*psi[2,2,14,2]*(1-phi[15,3,2]),
    phi[13,1,1]*psi[1,2,13,1]*phi[14,2,2]*psi[2,2,14,2]*phi[15,3,2]*psi[2,2,15,3])



theta2 <- logit(c(0.3,0.7,0.3,0.2,0.5,0.8,0.6,0.3))
struc2 <- list("phi"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1,2:length(states))),
               "delt"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "kap"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "rho"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "gam"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "eps"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
               "iot"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))

ll.il2(theta2,1:2,3,4,5,6,7,8,struc2,mch)
ll.il2(theta2,1:2,3,4,5,6,7,8,struc2,sim.dat)

op2 <- optim(theta2,ll.il2,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
            eps.ind=7,iot.ind=8,struc=struc2,ch=sim.dat,control=list(fnscale=-1),method="BFGS")
op2$convergence
logistic(op2$par)
op3 <- optim(theta2,ll.il2,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
             eps.ind=7,iot.ind=7,struc=struc2,ch=sim.dat,control=list(fnscale=-1),method="BFGS")
op3$convergence
logistic(op3$par)
logistic(op$par)==logistic(op3$par[1:7])

timer(op.n2 <- mclapply(1:n,function(x) optim(theta2,ll.il2,phi.ind=1:2,delt.ind=3,kap.ind=4,rho.ind=5,gam.ind=6,
                                             eps.ind=7,iot.ind=8,struc=struc2,ch=sim.dat.n[[x]],control=list(fnscale=-1),
                                             method="BFGS"),mc.cores=cores))
df.sim2 <- data.frame("par"=rep(c("phi1","phi2","delta","kappa","rho","gamma","epsilon","iota"),n),
                     "MLE"=c(sapply(1:n,function(x) logistic(op.n2[[x]]$par))),
                     "convergence"=rep(sapply(1:n,function(x) op.n2[[x]]$convergence),each=8))
# ggplot(df.sim,aes(par,MLE)) + geom_boxplot()
theta2 <- logit(c(0.3,0.7,0.3,0.2,0.5,0.8,0.6,0.6))
df.true2 <- data.frame("par"=c("phi1","phi2","delta","kappa","rho","gamma","epsilon","iota"),
                      "MLE"=logistic(theta2))
ggplot(df.sim2,aes(par,MLE,col=par)) + geom_boxplot() + geom_point(data=df.true2,aes(par,MLE),shape=5,col="black")


###############################
# Calculating some manual likelihoods to check
###############################
uch[153,]
Chi(1,13,6,phi,psi)
sum((1-phi[13,6,1]),
    phi[13,6,1]*psi[1,2,13,6]*(1-phi[14,2,7]),
    phi[13,6,1]*psi[1,2,13,6]*phi[14,2,7]*psi[2,2,14,7]*(1-phi[15,8,2]),
    phi[13,6,1]*psi[1,2,13,6]*phi[14,2,7]*psi[2,2,14,7]*phi[15,8,2]*psi[2,2,15,8])

test.l <- c(phi[8,1,1]*psi[1,2,8,1], # N -> E
            phi[9,2,2]*psi[2,2,9,2], # E -> E
            phi[10,3,2]*psi[2,1,10,3], #  E -> N
            phi[11,4,1]*psi[1,1,11,4], # N -> N
            phi[12,5,1]*psi[1,1,12,5], # N -> N
            sum((1-phi[13,6,1]),
                phi[13,6,1]*psi[1,2,13,6]*(1-phi[14,2,7]),
                phi[13,6,1]*psi[1,2,13,6]*phi[14,2,7]*psi[2,2,14,7]*(1-phi[15,8,2]),
                phi[13,6,1]*psi[1,2,13,6]*phi[14,2,7]*psi[2,2,14,7]*phi[15,8,2]*psi[2,2,15,8])) # N -> 0
sum(log(test.l))*2
ll.il(theta,1:2,3,4,5,6,7,struc,cbind(nch[c(325,325),],"age"=rep(1,2)))


nch[377,]
test.l2 <- c(phi[10,1,1]*psi[1,2,10,1], # N -> E
             phi[11,2,2]*psi[2,2,11,2], # E -> E
             phi[12,3,2]*psi[2,1,12,3], # E -> N
             phi[13,4,1]*psi[1,3,13,4], # N -> B1
             phi[14,5,3]*psi[3,8,14,5], # B1 -> S
             phi[15,6,8]*psi[8,7,15,6]) # S -> L_B_
sum(log(test.l2))*2
ll.il(theta,1:2,3,4,5,6,7,struc,cbind(nch[c(377,377),],"age"=rep(1,2)))

sum(log(test.l))+sum(log(test.l2))
ll.il(theta,1:2,3,4,5,6,7,struc,cbind(nch[c(325,377),],"age"=rep(1,2)))

load("guill.sim.RData")
df.true <- data.frame("par"=c("phi1","phi2","delta","kappa","rho","gamma","epsilon"),
                      "MLE"=logistic(theta))
ggplot(df.sim,aes(par,MLE,col=par)) + geom_boxplot() + geom_point(data=df.true,aes(par,MLE),shape=5,col="black")


