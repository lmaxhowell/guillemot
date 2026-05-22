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



# want to see how many individuals are ringed at each time
ni <- rep(0,ncol(uch)-1) # number of individuals
for(i in 1:nrow(uch)){
  ni[which(uch[i,1:16]!=0)[1]] <- ni[which(uch[i,1:16]!=0)[1]] + uch[i,17]
}
mean(ni) # 54.8 -> 55 individuals roughly each time




