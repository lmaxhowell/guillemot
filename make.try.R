library("gf.storm.RData")
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
  
  delt <- untrans(del,struc.const$age,struc.const$time,struc.const$state)
  kap <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(gamm,struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
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
states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
Time <- ncol(ch)-1
timeclasses <- 2
ageclasses <- 4
beta.struc <- list(1:(11-1),11:15)

starter <- theta.s(4)
op <- nlm(f=ll.il.ms, # likelihood function
          p=starter, # initial values
          ageclasses=ageclasses,
          timeclasses=timeclasses,
          beta.struc=beta.struc,
          ch=ch,
          hessian=TRUE)
par.fun <- function(mle){
  c(mle[1:5],logistic(mle[6:length(mle)]))
}

par.fun(op$estimate)

ll.il.ms2 <- function(theta,ageclasses,timeclasses,beta.struc,ch){
  # beta.struc should be a list corresponding to
  # where the breaks in beta go e.g. list(1:10,11:15)
  # would be a change in beta at the 11 time point
  alpha <- theta[1:(ageclasses-1)]
  beta <- theta[ageclasses:(timeclasses+ageclasses-1)]
  phi <- logistic(phi.trans(alpha,beta,beta.struc))
  # print(length(alpha))
  # print(length(beta))
  # print(length(rho.input))
  # print(dim(phi))
  rho <- logistic(rho.trans(op$estimate[6:8]))

  Time <- ncol(ch)-1
  
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))

  
  delt <- untrans(del,struc.const$age,struc.const$time,struc.const$state)
  kap <- untrans(logistic(op$estimate[9]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(gamm,struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(op$estimate[10]),struc.const$age,struc.const$time,struc.const$state)

  psi <- make.psi(delt,kap,rho,gam,eps)
  ll <- il(ch,phi,psi)
  return(-ll)
}

starter2 <- theta.s(4)[1:5]
op2 <- nlm(f=ll.il.ms2, # likelihood function
           p=rep(0,5), # initial values
           ageclasses=ageclasses,
           timeclasses=timeclasses,
           beta.struc=beta.struc,
           ch=ch,
           hessian=TRUE,
           print.level=2)

phi.trans <- function(theta){ # designed to be length(theta)=4
  phi <- untrans(theta,list(1:3,4:16),list(1:10,11:15),list(1:8))
  return(phi)
}

ll.il.ms3 <- function(theta,ch){
  phi <- logistic(phi.trans(theta))
  rho <- logistic(rho.trans(op$estimate[6:8]))
  
  Time <- ncol(ch)-1
  
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  
  
  delt <- untrans(del,struc.const$age,struc.const$time,struc.const$state)
  kap <- untrans(logistic(op$estimate[9]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(gamm,struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(op$estimate[10]),struc.const$age,struc.const$time,struc.const$state)
  
  psi <- make.psi(delt,kap,rho,gam,eps)
  ll <- il(ch,phi,psi)
  return(-ll)
}

op3 <- nlm(f=ll.il.ms3, # likelihood function
           p=rep(0,4), # initial values
           ch=ch,
           hessian=TRUE)

phi.trans <- function(theta){ # designed to be length(theta)=1
  phi <- untrans(theta,list(1:16),list(1:15),list(1:8))
  return(phi)
}

op4 <- nlm(f=ll.il.ms3, # likelihood function
           p=0, # initial values
           ch=ch,
           hessian=TRUE,
           steptol=1e-5)

