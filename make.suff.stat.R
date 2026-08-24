load("gf.storm.RData")
load("matrix_AUK.RData")
states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")

suff.stat <- function(ch){
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  Age <- Time
  
  Agesi <- ch[,ncol(ch)]
  ch <- ch[,-ncol(ch)]
  
  # m is r->s
  # v is r->0
  m <- array(0,dim=c(length(states),length(states),Time,Age)) # r,s,t,a
  v <- array(0,dim=c(length(states),Time,Age)) # r,t,a
  
  violations <- list(vector(mode = "list", length = 8),
                     vector(mode = "list", length = 8),
                     vector(mode = "list", length = 8),
                     vector(mode = "list", length = 8),
                     vector(mode = "list", length = 8),
                     vector(mode = "list", length = 8),
                     vector(mode = "list", length = 8),
                     vector(mode = "list", length = 8))
  isitavio <- function(r,s){
    rtn <- FALSE
    if((r %in% 3:7) & (s %in% 1:2)){ # if its in a breeding state it cant return to nonbreeding
      rtn <- TRUE
    }else if(r==8 & s==8){ # if it skips multiple times
      rtn <- TRUE
    }
    return(rtn)
  }
  
  for(i in 1:nrow(ch)){
    wh <- which(ch[i,]!=0)
    if(length(wh)==1){
      r <- which(states==ch[i,wh])
      v[r,wh,1] <- v[r,wh,1] + 1 
    }else{
      for(j in 1:(length(wh)-1)){
        r <- which(states==ch[i,wh[j]])
        s <- which(states==ch[i,wh[j+1]])
        m[r,s,wh[j],j+Agesi[i]-1] <- m[r,s,wh[j],j+Agesi[i]-1] + 1
        if(isitavio(r,s)){
          print(c(r,s))
          violations[[r]][[s]] <- append(violations[[r]][[s]],i)
        }
      }
      r <- which(states==ch[i,wh[length(wh)]])
      v[r,wh[length(wh)],1] <- v[r,wh[length(wh)],1] + 1
    }
  }
  m2 <- rowSums(m,dims = 2)
  dimnames(m2) <- list(states,states)
  
  return(list("m"=m,"v"=v,"m2"=m2,"violations"=violations))
}

mv <- suff.stat(ch2)
# [1] 4 1
# [1] 6 1
# [1] 8 8
# [1] 7 1
# [1] 3 1
# [1] 3 1
# [1] 3 1
# [1] 6 1
# [1] 3 1
# [1] 6 1
# [1] 5 1
# [1] 3 1
vs <- rowSums(mv$v,dims = 1)
ms <- rowSums(mv$m,dims = 2)
dimnames(ms) <- list(states,states)
ms

codes <- data.frame("state"=c("N","E","B1","LB","L_B","LB_","L_B_","S"),
               "code"=c(1,13,2,3,5,4,6,13))

ch[mv$violations[[8]][[8]],]
matrix_AUK[mv$violations[[8]][[8]],]


vios <- do.call("rbind",lapply(unlist(mv$violations),function(i) cbind(rbind(ch[i,1:16],matrix_AUK[i,]),"row"=c(i,i))))

# double skips
c(129)
# observed not breeding after first breeding
c(253,293,374,398,611)
# observed not breeding after breeding (not first breeding)
c(69,492,125,393,422,170)

chn <- ch[-unlist(mv$violations),]
mv2 <- suff.stat(chn)
vs2 <- rowSums(mv2$v,dims = 1)
ms2 <- rowSums(mv2$m,dims = 2)
dimnames(ms2) <- list(states,states)
ms2


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
  x <- 8
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
          ch=chn,
          hessian=TRUE)
par.fun <- function(mle){
  c(mle[1:5],logistic(mle[6:length(mle)]))
}
par.fun(op$estimate)
logistic(phi.trans(op$estimate[1:3],op$estimate[4:5],beta.struc))[,1:4,1]


opr <- lapply(1:20, function(x){ 
           starter <- theta.s(4)
           return(c(nlm(f=ll.il.ms, # likelihood function
           p=starter, # initial values
           ageclasses=ageclasses,
           timeclasses=timeclasses,
           beta.struc=beta.struc,
           ch=chn,
           hessian=TRUE),list("start"=starter)))})

sapply(1:20,function(x) opr[[x]]$code)
sapply(1:20,function(x) par.fun(opr[[x]]$estimate))

##########
# compare with simulated dataset
dat.sim.wrap2 <- function(theta,ageclasses,timeclasses,beta.struc,seed){
  alpha <- theta[1:(ageclasses-1)]
  beta <- theta[ageclasses:(timeclasses+ageclasses-1)]
  rho.input <- theta[(timeclasses+ageclasses):(timeclasses+ageclasses+2)]
  phi <- logistic(phi.trans(alpha,beta,beta.struc))
  rho <- logistic(rho.trans(rho.input))
  
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  
  theta <- theta[-c(1:(timeclasses+ageclasses+2))]
  
  delt <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
  kap <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(logistic(theta[3]),struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(theta[4]),struc.const$age,struc.const$time,struc.const$state)
  df <- dat.simulate2(phi,delt,kap,rho,gam,eps,ni2[1:15],seed)
  return(df)
}

alph <- c(0.8,1.1,1.4)
bet <- c(0.2,-0.2)
rho <- c(0.1,0.3,0.6)
otherparm <- c(0.4,0.2,0.6,0.5) # delta,kappa,gamma,epsilon
theta <- c(alph,bet,logit(rho),logit(otherparm))
dat <- dat.sim.wrap2(theta,4,2,beta.struc,1612049052)

mv.s <- suff.stat(dat)
ms3 <- rowSums(mv.s$m,dims = 2)
dimnames(ms3) <- list(states,states)
ms3

op2 <- nlm(f=ll.il.ms, # likelihood function
           p=theta, # initial values
           ageclasses=ageclasses,
           timeclasses=timeclasses,
           beta.struc=beta.struc,
           ch=dat,
           hessian=TRUE)
par.fun(op2$estimate)
logistic(phi.trans(op2$estimate[1:3],op2$estimate[4:5],beta.struc))[10:11,1:4,1]
dat1.5 <- dat.sim.wrap2(op2$estimate,4,2,beta.struc,1612049052)
mv.s1.5 <- suff.stat(dat1.5)
mv.s1.5$m2


set.seed(1612049052)
seeds <- sample(1:.Machine$integer.max,20)
dat2 <- lapply(1:20,function(x) dat.sim.wrap2(theta,4,2,beta.struc,seeds[x]))
ms4 <- lapply(1:20,function(x){
              mv.s <- suff.stat(dat2[[x]])
              ms <- rowSums(mv.s$m,dims = 2)
              dimnames(ms) <- list(states,states)
              return(ms)
            })

dat3 <- dat.sim.wrap2(rowMeans(sapply(1:20,function(x) par.fun(opr[[x]]$estimate))),4,2,beta.struc,1612049052)

mv3 <- suff.stat(dat3)
ms5 <- rowSums(mv3$m,dims = 2)
dimnames(ms5) <- list(states,states)
ms5


ll.il.simp <- function(theta,ch){
  timeclasses <- 2
  ageclasses <- 4
  beta.struc <- list(1:(11-1),11:15)
  alpha <- theta[1:(ageclasses-1)]
  beta <- theta[ageclasses:(timeclasses+ageclasses-1)]
  phi <- logistic(phi.trans(alpha,beta,beta.struc))
  
  # Ages <- length(alpha)
  Time <- ncol(ch)-1
  
  # struc.time <- list("age"=list(1:Time),"time"=as.list(1:(Time-1)),"state"=list(1:length(states)))
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  # struc.state <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=as.list(1:length(states)))
  
  theta <- theta[-c(1:(5))]
  
  rho <- untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state)
  delt <- untrans(logistic(theta[2]),struc.const$age,struc.const$time,struc.const$state)
  kap <- untrans(logistic(theta[3]),struc.const$age,struc.const$time,struc.const$state)
  gam <- untrans(logistic(theta[4]),struc.const$age,struc.const$time,struc.const$state)
  eps <- untrans(logistic(theta[5]),struc.const$age,struc.const$time,struc.const$state)
  psi <- make.psi(delt,kap,rho,gam,eps)
  
  ll <- il(ch,phi,psi)
  return(-ll)
}

op3 <- nlm(f=ll.il.simp, # likelihood function
           p=rep(0,10), # initial values
           ch=chn,
           hessian=TRUE)
par.fun(op3$estimate)

op4 <- optim(rep(0,10),
             ll.il.simp, # likelihood function
             ch=chn,
             hessian=TRUE)
op4$convergence
par.fun(op4$par)
logistic(phi.trans(op4$par[1:3],op4$par[4:5],beta.struc))[10:11,1:4,1]

opr2 <- lapply(1:20, function(x){ 
             starter <- theta.s(4)
             return(c(optim(starter,
                          ll.il.ms, # likelihood function
                          ageclasses=ageclasses,
                          timeclasses=timeclasses,
                          beta.struc=beta.struc,
                          ch=chn,
                          hessian=TRUE),list("start"=starter)))})
sapply(1:20,function(x) opr2[[x]]$convergence )
sapply(1:20,function(x) par.fun(opr2[[x]]$par) )
rowMeans(sapply(1:20,function(x) par.fun(opr2[[x]]$par) ))
which.min(sapply(1:20,function(x) opr2[[x]]$value ))
par.fun(opr2[[2]]$par)
dat4 <- dat.sim.wrap2(opr2[[2]]$par,4,2,beta.struc,1612049052)
suff.stat(dat4)$m2
