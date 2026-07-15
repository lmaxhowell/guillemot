file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed

phi.trans <- function(alpha,beta){
  # inputs should be vectors
  alpha <- c(0,alpha) # to constrain alpha_1=0
  Ages <- length(alpha)
  Time <- length(beta)
  phi <- array(0,dim=c(Time,Time,length(states)))
  for(a in 1:Time){
    for(t in 1:Time){
      phi[t,a,] <- alpha[min(a,Ages)]+beta[t]
    }
  }
  return(phi)
}

rho.trans <- function(rho.input){
  rho <- array(0,dim=c(Time,Time,length(states)))
  # rho[,1:2,] <- 0 # sets as default anyway
  rho[,3,] <- rho.input[1]
  rho[,4,] <- rho.input[2]
  rho[,5:Time,] <- rho.input[3]
  return(rho)
}

ll.il.wrap <- function(theta,alp.ind,bet.ind,rho.ind,ch){
  alpha <- theta[alp.ind]
  beta <- theta[bet.ind]
  rho.input <- theta[rho.ind]
  phi <- logistic(phi.trans(alpha,beta))
  rho <- logistic(rho.trans(rho.input))
  # rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
  
  # Ages <- length(alpha)
  Time <- ncol(ch)-1
  
  struc.time <- list("age"=list(1:Time),"time"=as.list(1:Time),"state"=list(1:length(states)))
  struc.const <- list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states)))
  struc.state <- list("age"=list(1:Time),"time"=list(1:Time),"state"=as.list(1:length(states)))
  
  theta <- theta[-c(alp.ind,bet.ind,rho.ind)]
  
  delt <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
  kap <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(logistic(theta[3]),struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(theta[4]),struc.const$age,struc.const$time,struc.const$state)
  
  psi <- make.psi(delt,kap,rho,gam,eps)
  # print(psi)
  
  ll <- il(ch,phi,psi)
  return(ll)
}

# three survival age classes with alpha1=0
# means that alpha is of length 2 
# then beta is of length Time
# then we have rho input being of length 3
# then need four more parameters for the other matrices

load("guillemot.RData")
Time <- ncol(ch)-1
theta <- logit(rep(0.55,Time+2+3+4))

states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
ll.il.wrap(theta,1:2,3:(Time+2),(Time+3):(Time+5),ch)
par.fun <- function(op.par,n){
  output <- c(op.par[1:(n-1)],logistic(op.par[n:length(op.par)]))
  return(output)
}

op1 <- optim(theta,ll.il.wrap,alp.ind=1:2,bet.ind=3:(Time+2),
             rho.ind=(Time+3):(Time+5),ch=ch,
             control=list(fnscale=-1,maxit=1000),method="BFGS",hessian=TRUE)
op1$convergence
par.fun(op1$par,19)

op2 <- optim(theta,ll.il.wrap,alp.ind=1:2,bet.ind=3:(Time+2),
             rho.ind=(Time+3):(Time+5),ch=ch,
             control=list(fnscale=-1,maxit=1000),method="Nelder-Mead",hessian=TRUE)
op2$convergence
c(logistic(op2$par[1:2]),op2$par[3:18],logistic(op2$par[19:25]))
par.fun(op2$par,19)

df.phi <- data.frame("Time"=2010:2025,"phi"=c(logistic(phi.trans(op2$par[1:2],op2$par[3:18]))[,1:3,1]),
                     "AgeClass"=as.factor(rep(1:3,each=Time)))

ggplot(df.phi,aes(Time,phi,linetype=AgeClass,col=AgeClass)) + geom_line() + geom_point()

data.frame("par"=c(paste0("alpha",2:3),paste0("beta",1:Time),paste0("rho",3:5),"delta","kappa","gamma","epsilon"),
           "MLE"=par.fun(op2$par,19))


################################################
# Now want to make the function do additional age classes
################################################
theta2 <- logit(rep(0.55,Time+2+3+4+1)) # plus one compared to theta as one more survival age class
op3 <- optim(theta2,ll.il.wrap,alp.ind=1:3,bet.ind=4:(Time+3),
             rho.ind=(Time+4):(Time+6),ch=ch,
             control=list(fnscale=-1,maxit=1000),method="Nelder-Mead",hessian=TRUE)
op3$convergence
data.frame("par"=c(paste0("alpha",2:4),paste0("beta",1:Time),paste0("rho",3:5),"delta","kappa","gamma","epsilon"),
           "MLE"=par.fun(op3$par,20))
data.frame("par"=c(paste0("alpha",2:4),paste0("beta",1:Time),paste0("rho",3:5),"delta","kappa","gamma","epsilon"),
           "MLE_4AC"=par.fun(op3$par,20),
           "MLE_3AC"=c(op2$par[1:2],NA,op2$par[3:18],logistic(op2$par[19:25])))
df.phi2 <- data.frame("Time"=2010:2025,"phi"=c(logistic(phi.trans(op3$par[1:3],op3$par[4:19]))[,1:4,1]),
                      "AgeClass"=as.factor(rep(1:4,each=Time)))
ggplot(df.phi2,aes(Time,phi,linetype=AgeClass,col=AgeClass)) + geom_line() + geom_point()

################################################
# aaaaaaaand 5 age classes
################################################
theta3 <- logit(rep(0.55,Time+2+3+4+2)) # plus two compared to theta as two more survival age classs
op4 <- optim(theta3,ll.il.wrap,alp.ind=1:4,bet.ind=5:(Time+4),
             rho.ind=(Time+5):(Time+7),ch=ch,
             control=list(fnscale=-1,maxit=1000),method="Nelder-Mead",hessian=TRUE)
op4$convergence
data.frame("par"=c(paste0("alpha",2:5),paste0("beta",1:Time),paste0("rho",3:5),"delta","kappa","gamma","epsilon"),
           "MLE"=par.fun(op4$par,21))
data.frame("par"=c(paste0("alpha",2:5),paste0("beta",1:Time),paste0("rho",3:5),"delta","kappa","gamma","epsilon"),
           "MLE_5AC"=par.fun(op4$par,21),
           "MLE_4AC"=c(op3$par[1:3],NA,op3$par[4:19],logistic(op3$par[20:26])),
           "MLE_3AC"=c(op2$par[1:2],NA,NA,op2$par[3:18],logistic(op2$par[19:25])))
df.phi3 <- data.frame("Time"=2010:2025,"phi"=c(logistic(phi.trans(op4$par[1:4],op4$par[5:20]))[,1:5,1]),
                      "AgeClass"=as.factor(rep(1:5,each=Time)))
ggplot(df.phi3,aes(Time,phi,linetype=AgeClass,col=AgeClass)) + geom_line() + geom_point()

df.phi4 <- cbind(rbind(df.phi,df.phi2,df.phi3),"NumOfAgeClasses"=rep(c(3,4,5),c(nrow(df.phi),nrow(df.phi2),nrow(df.phi3))))
ggplot(df.phi4,aes(Time,phi,linetype=AgeClass,col=AgeClass)) + geom_line() + geom_point() + facet_wrap(~NumOfAgeClasses,ncol=1)

################################################
# Adding state dependence into delta and kappa
################################################
ll.il.wrap2 <- function(theta,alp.ind,bet.ind,rho.ind,delt.ind,kap.ind,ch){
  alpha <- theta[alp.ind]
  beta <- theta[bet.ind]
  rho.input <- theta[rho.ind]
  phi <- logistic(phi.trans(alpha,beta))
  rho <- logistic(rho.trans(rho.input))
  # rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
  
  # Ages <- length(alpha)
  Time <- ncol(ch)-1
  
  struc.time <- list("age"=list(1:Time),"time"=as.list(1:Time),"state"=list(1:length(states)))
  struc.const <- list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states)))
  struc.kap <- list("age"=list(1:Time),"time"=list(1:Time),"state"=list(c(1,2,8),3,4,5,6,7))
  struc.delt <- list("age"=list(1:Time),"time"=list(1:Time),"state"=list(c(1,2),3,4,5,6,7,8))
  
  # so kappa needs 5 parameters and delta needs 6
  # add a zero at the beginning for the rows where that parameter doesn't actually show up
  delt <- untrans(logistic(c(0,theta[delt.ind])),struc.delt$age,struc.delt$time,struc.delt$state)
  kap <- untrans(logistic(c(0,theta[kap.ind])),struc.kap$age,struc.kap$time,struc.kap$state)
  
  
  theta <- theta[-c(alp.ind,bet.ind,rho.ind,delt.ind,kap.ind)]
  gam <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
  
  psi <- make.psi(delt,kap,rho,gam,eps)
  # print(psi)
  
  ll <- il(ch,phi,psi)
  return(ll)
}

# three age classes
theta.s <- function(n){ # want n to be desired number of age classes
  return(logit(rep(0.55,31+n)))
}
op.s <- optim(theta.s(3),ll.il.wrap2,alp.ind=1:2,bet.ind=3:(Time+2),
              rho.ind=(Time+3):(Time+5),delt.ind=(Time+6):(Time+10),
              kap.ind=(Time+11):(Time+16),ch=ch,
              control=list(fnscale=-1,maxit=1000),method="Nelder-Mead",hessian=TRUE)
op.s$convergence
par.fun(op.s$par,19)
df.phi.s <- data.frame("Time"=2010:2025,"phi"=c(logistic(phi.trans(op.s$par[1:2],op.s$par[3:18]))[,1:3,1]),
                     "AgeClass"=as.factor(rep(1:3,each=Time)))
ggplot(df.phi,aes(Time,phi,linetype=AgeClass,col=AgeClass)) + geom_line() + geom_point()

# four age classes
op.s2 <- optim(theta.s(4),ll.il.wrap2,alp.ind=1:3,bet.ind=4:(Time+3),
              rho.ind=(Time+4):(Time+6),delt.ind=(Time+7):(Time+11),
              kap.ind=(Time+12):(Time+17),ch=ch,
              control=list(fnscale=-1,maxit=1000),method="Nelder-Mead",hessian=TRUE)
op.s2$convergence
par.fun(op.s2$par,20)
df.phi.s2 <- data.frame("Time"=2010:2025,"phi"=c(logistic(phi.trans(op.s2$par[1:3],op.s2$par[4:19]))[,1:4,1]),
                      "AgeClass"=as.factor(rep(1:4,each=Time)))

# five age classes
op.s3 <- optim(theta.s(5),ll.il.wrap2,alp.ind=1:4,bet.ind=5:(Time+4),
               rho.ind=(Time+5):(Time+7),delt.ind=(Time+8):(Time+12),
               kap.ind=(Time+13):(Time+18),ch=ch,
               control=list(fnscale=-1,maxit=1000),method="Nelder-Mead",hessian=TRUE)
op.s3$convergence
par.fun(op.s3$par,21)
df.phi.s3 <- data.frame("Time"=2010:2025,"phi"=c(logistic(phi.trans(op.s3$par[1:4],op.s3$par[5:20]))[,1:5,1]),
                      "AgeClass"=as.factor(rep(1:5,each=Time)))

# combined df for survival
df.phi.s4 <- cbind(rbind(df.phi.s,df.phi.s2,df.phi.s3),"NumOfAgeClasses"=rep(c(3,4,5),c(nrow(df.phi.s),nrow(df.phi.s2),nrow(df.phi.s3))))
ggplot(df.phi.s4,aes(Time,phi,linetype=AgeClass,col=AgeClass)) + geom_line() + geom_point() + facet_wrap(~NumOfAgeClasses,ncol=1)

df.phi.all <- cbind(rbind(df.phi4,df.phi.s4),"StateDependence"=rep(c(FALSE,TRUE),each=nrow(df.phi4)))
ggplot(df.phi.all,aes(Time,phi,linetype=AgeClass,col=AgeClass)) + geom_line() + geom_point() + facet_grid(NumOfAgeClasses~StateDependence)

df.delta <- data.frame("MLE"=c(par.fun(op.s$par,19)[22:26],
                               par.fun(op.s2$par,20)[23:27],
                               par.fun(op.s3$par,21)[24:28]),
                       "AgeClasses"=as.factor(rep(3:5,each=5)),
                       "par"=rep(paste0("delta",1:5),3))
df.kappa <- data.frame("MLE"=c(par.fun(op.s$par,19)[27:32],
                               par.fun(op.s2$par,20)[28:33],
                               par.fun(op.s3$par,21)[29:34]),
                       "AgeClasses"=as.factor(rep(3:5,each=6)),
                       "par"=rep(paste0("kappa",1:6),3))
df.states <- cbind(rbind(df.delta,df.kappa))
ggplot(df.states,aes(par,MLE,col=AgeClasses)) + geom_point()
ggplot(df.states,aes(par,MLE,col=AgeClasses)) + geom_point() + facet_wrap(~AgeClasses,ncol=3)

aics <- c(2*length(theta.s(3))-2*op2$value,
          2*length(theta.s(4))-2*op3$value,
          2*length(theta.s(5))-2*op4$value,
          2*length(theta.s(3))-2*op.s$value,
          2*length(theta.s(4))-2*op.s2$value,
          2*length(theta.s(5))-2*op.s3$value)
df.aic <- data.frame("aic"=aics,
                     "AgeClasses"=rep(3:5,2),
                     "StateDependence"=rep(c(FALSE,TRUE),each=3))

# selects no state dependence and four age classes which is op3
View(solve(-op3$hessian))
data.frame("par"=c(paste0("alpha",2:4),paste0("beta",1:Time),paste0("rho",3:5),"delta","kappa","gamma","epsilon"),
           "MLE"=par.fun(op3$par,20),
           "var"=diag(solve(-op3$hessian)),
           "sd"=sqrt(diag(solve(-op3$hessian))))

timer(op.fin <- optim(theta2,ll.il.wrap,alp.ind=1:3,bet.ind=4:(Time+3),
                    rho.ind=(Time+4):(Time+6),ch=ch,
                    control=list(fnscale=-1,maxit=10000),method="Nelder-Mead",hessian=TRUE))
# Time difference of 1.630507 hours
save(op.fin,file="op.fin.RData")
par.fun(op.fin$par,20)
df.fin <- data.frame("Time"=2010:2025,"phi"=c(logistic(phi.trans(op.fin$par[1:3],op.fin$par[4:19]))[,1:4,1]),
                     "AgeClass"=as.factor(rep(1:4,each=Time)))
ggplot(df.fin,aes(Time,phi,linetype=AgeClass,col=AgeClass)) + geom_line() + geom_point()
timer(op.fin2 <- optim(theta2,ll.il.wrap,alp.ind=1:3,bet.ind=4:(Time+3),
                       rho.ind=(Time+4):(Time+6),ch=ch,
                       control=list(fnscale=-1,maxit=50000),method="Nelder-Mead",hessian=TRUE))
# Time difference of 2.492975 hours
timer(op.fin3 <- optim(theta2,ll.il.wrap,alp.ind=1:3,bet.ind=4:(Time+3),
                       rho.ind=(Time+4):(Time+6),ch=ch,
                       control=list(fnscale=-1,maxit=10000),method="BFGS",hessian=TRUE))
# Time difference of 55.78593 mins

###########################################
# Seeing if delta is dependant on age
###########################################
ll.il.wrap3 <- function(theta,AgeClasses,ch){
  alpha <- theta[1:3] # choosing four age classes -> three alphas
  beta <- theta[4:19] # Time betas
  rho.input <- theta[20:22] # then three rhos
  phi <- logistic(phi.trans(alpha,beta))
  rho <- logistic(rho.trans(rho.input))
  # rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
  
  # Ages <- length(alpha)
  Time <- ncol(ch)-1
  
  struc.const <- list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states)))
  if(AgeClasses==1){
    struc.age <- list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states)))
  }else if(AgeClasses>1){
    struc.age <- list("age"=c(as.list(1:(AgeClasses-1)),list(AgeClasses:Time)),"time"=list(1:Time),"state"=list(1:length(states)))
  }
  
  
  theta <- theta[-c(1:22)]
  
  delt <- untrans(logistic(theta[1:AgeClasses]),struc.const$age,struc.const$time,struc.const$state)
  
  theta <- theta[-c(1:AgeClasses)]
  kap <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(theta[3]),struc.const$age,struc.const$time,struc.const$state)
  
  psi <- make.psi(delt,kap,rho,gam,eps)
  # print(psi)
  
  ll <- il(ch,phi,psi)
  return(ll)
}

DeltaAges <- 1:16
theta.s <- function(n){ # want n to be desired number of age classes
  return(logit(rep(0.55,25+n)))
}
DeltaAge <- mclapply(DeltaAges,function(x) optim(theta.s(x),ll.il.wrap3,AgeClasses=x,ch=ch,
                                   control=list(fnscale=-1,maxit=1000),method="Nelder-Mead",hessian=TRUE),mc.cores=2)
DeltaAIC <- data.frame("AgeClass"=DeltaAges,
                       "AIC"=sapply(DeltaAges,function(x) 2*length(theta.s(x))-2*DeltaAge[[x]]$value ))
which.min(DeltaAIC$AIC) # 4
# sapply(DeltaAges,function(x) par.fun(DeltaAge[[x]]$par,16+x))
par.fun(DeltaAge[[4]]$par,20)
df.delta4 <- data.frame("Time"=2010:2025,"phi"=c(logistic(phi.trans(par.fun(DeltaAge[[4]]$par,20)[1:3],par.fun(DeltaAge[[4]]$par,20)[4:19]))[,1:4,1]),
                        "AgeClass"=as.factor(rep(1:4,each=Time)))
ggplot(df.delta4,aes(Time,phi,col=AgeClass)) + geom_line() + geom_point()



load("ModelSelection.RData")
ggplot(df.par,aes(par,MLE,col=StateDependance)) + geom_point() + facet_wrap(~as.factor(AgeClasses))

df.selected <- data.frame("Time"=2010:2025,"phi"=c(logistic(phi.trans(df.par$MLE[df.par$AgeClasses==6 & df.par$StateDependance==TRUE][1:5],op.s3$par[6:21]))[,1:6,1]),
                          "AgeClass"=as.factor(rep(1:6,each=Time)))
ggplot(df.selected,aes(Time,phi,col=AgeClass)) + geom_line() + geom_point()

