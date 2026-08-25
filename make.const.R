load("gf.storm.RData")

Time <- ncol(ch)-1
states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "eps"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))))

start <- rep(0,6)
ll.il(start,phi.ind=1,delt.ind=2,kap.ind=3,rho.ind=4,gam.ind=5,eps.ind=6,struc=struc,ch=ch)
ll.il.m <- function(theta,phi.ind,delt.ind,kap.ind,rho.ind,gam.ind,eps.ind,struc,ch){
  ll <- ll.il(theta=theta,phi.ind=phi.ind,delt.ind=delt.ind,kap.ind=kap.ind,rho.ind=rho.ind,gam.ind=gam.ind,eps.ind=eps.ind,struc=struc,ch=ch)
  return(-ll)
}

op <- nlm(f=ll.il.m,
          p=start,
          phi.ind=1,
          delt.ind=2,
          kap.ind=3,
          rho.ind=4,
          gam.ind=5,
          eps.ind=6,
          struc=struc,
          ch=ch)
op$code
logistic(op$estimate)

op2 <- optim(start,
             ll.il.m,
             phi.ind=1,
             delt.ind=2,
             kap.ind=3,
             rho.ind=4,
             gam.ind=5,
             eps.ind=6,
             struc=struc,
             ch=ch)
op2$convergence
logistic(op2$par)
op2$value

op3 <- optim(start,
             ll.il.m,
             phi.ind=1,
             delt.ind=2,
             kap.ind=3,
             rho.ind=4,
             gam.ind=5,
             eps.ind=6,
             struc=struc,
             ch=ch,
             control=list(maxit=1000))
op3$convergence
logistic(op3$par)
op3$value

timer(op4 <- optim(start,
                   ll.il.m,
                   phi.ind=1,
                   delt.ind=2,
                   kap.ind=3,
                   rho.ind=4,
                   gam.ind=5,
                   eps.ind=6,
                   struc=struc,
                   ch=ch,
                   control=list(maxit=5000)))
op4$convergence
logistic(op4$par)
op4$value
df <- data.frame("par"=c("phi","delta","kappa","rho","gamma","epsilon"),
                 "mle"=logistic(op4$par),
                 "convergence"=rep(op4$convergence,6))


op5 <- optim(start,
             ll.il.m,
             phi.ind=1,
             delt.ind=2,
             kap.ind=3,
             rho.ind=4,
             gam.ind=5,
             eps.ind=6,
             struc=struc,
             ch=ch,
             method="BFGS")
op5$convergence
logistic(op5$par)
op5$value


struc2 <- list("phi"=list("age"=list(1,2,3,4,5,6:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
               "eps"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))))

start2 <- rep(0,6+5)
timer(op6 <- optim(start2,
                   ll.il.m,
                   phi.ind=1:6,
                   delt.ind=7,
                   kap.ind=8,
                   rho.ind=9,
                   gam.ind=10,
                   eps.ind=11,
                   struc=struc2,
                   ch=ch,
                   control=list(maxit=5000)))
# Time difference of 26.54476 mins
op6$convergence
logistic(op6$par)
op6$value
df2 <- data.frame("par"=c(paste0("phi",1:6),"delta","kappa","rho","gamma","epsilon"),
                  "mle"=logistic(op6$par),
                  "convergence"=rep(op6$convergence,11))
save(list=ls(),file="it.works.backup.RData")

# high estimate for phi2 so trying with lower maxit and seeing
# if, even if it doesnt "converge" it is slightly better
timer(op7 <- optim(start2,
                   ll.il.m,
                   phi.ind=1:6,
                   delt.ind=7,
                   kap.ind=8,
                   rho.ind=9,
                   gam.ind=10,
                   eps.ind=11,
                   struc=struc2,
                   ch=ch,
                   control=list(maxit=1000)))
op7$convergence
logistic(op7$par)
op7$value

tb <- table(unlist(ch[,1:16]))

Bs <- tb["LB"] + tb["L_B"]
B_s <- tb["LB_"] + tb["L_B_"]
Ls <- tb["LB"] + tb["LB_"]
L_s <- tb["L_B"] + tb["L_B_"]

gamm <- unname(Bs/(Bs+B_s))
del <- unname(L_s/(Ls+L_s))

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
  
  # struc.time <- list("age"=list(1:Time),"time"=as.list(1:(Time-1)),"state"=list(1:length(states)))
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  # struc.state <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
  
  theta <- theta[-c(1:(timeclasses+ageclasses+2))]
  
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
  return(-ll)
}
timeclasses <- 2
ageclasses <- 4
beta.struc <- list(1:(11-1),11:15)

par.name <- c(paste0("alpha",1:3),paste0("beta",11:12),paste0("rho",c(23,4,5)),"delta","kappa","gamma","epsilon")

start3 <- logit(c(rep(0.5,8),del,0.5,gamm,0.5))

timer(op8 <- optim(start3,
                   ll.il.ms,
                   ageclasses=ageclasses,
                   timeclasses=timeclasses,
                   beta.struc=beta.struc,
                   ch=ch,
                   control=list(maxit=10000)))

df3 <- data.frame("par"=par.name,
                  "mle"=logistic(op8$par),
                  "convergence"=rep(op8$convergence,length(par.name)))
logistic(phi.trans(op8$par[1:3],op8$par[4:5],beta.struc)[10:11,1:4,1])


load("guillemot2.RData")
timer(op9 <- optim(start3,
                   ll.il.ms,
                   ageclasses=ageclasses,
                   timeclasses=timeclasses,
                   beta.struc=beta.struc,
                   ch=ch3,
                   control=list(maxit=10000)))

df3 <- data.frame("par"=par.name,
                  "mle"=logistic(op8$par),
                  "convergence"=rep(op8$convergence,length(par.name)))
logistic(phi.trans(op8$par[1:3],op8$par[4:5],beta.struc)[10:11,1:4,1])


