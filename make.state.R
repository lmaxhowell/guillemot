load("gf.storm.RData")
library(parallel)
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

theta.s <- function(state.depen=c(F,F,F,F)){
  n <- 8
  x <- sum(ifelse(state.depen,c(6,5,6,2),c(1,1,1,1)))
  return(logit(runif(x+n,0.1,0.9)))
}

ll.il.ms <- function(theta,state.depen=c(F,F,F,F),ch){
  # state4 is about the four final parameters and whether they should be state
  # dependant or not (false for no and true for yes)
  timeclasses <- 2
  ageclasses <- 4
  beta.struc <- list(1:(11-1),11:15)
  alpha <- theta[1:(ageclasses-1)]
  beta <- theta[ageclasses:(timeclasses+ageclasses-1)]
  rho.input <- theta[(timeclasses+ageclasses):(timeclasses+ageclasses+2)]
  phi <- logistic(phi.trans(alpha,beta,beta.struc))
  rho <- logistic(rho.trans(rho.input))

  Time <- ncol(ch)-1
  
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  struc.state <- function(n){ # n should be 1,2,3 or 4
    if(n==1){ # delta can only vary over states 3:8
      lsss <- c(list(1:3),as.list(4:8))
    }else if(n==2){ # kappa can only vary over states 3:7
      lsss <- c(list(1:3),as.list(4:6),list(7:8))
    }else if(n==3){ # gamma can only vary over states 3:8
      lsss <- c(list(1:3),as.list(4:8))
    }else if(n==4){ # epsilon can only vary over states 1:2
      lsss <- c(list(1),list(2:8))
    }
    return(lsss)  
  }
  
  theta <- theta[-c(1:(timeclasses+ageclasses+2))]
  
  nams <- c("delt","kap","gam","eps")
  for(i in 1:4){
    if(state.depen[i]){
      assign(nams[i],
             untrans(logistic(theta[1:length(struc.state(i))]),struc.const$age,struc.const$time,
                     struc.state(i)))
      theta <- theta[-(1:length(struc.state(i)))]
    }else{
      assign(nams[i],untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state))
      theta <- theta[-1]
    }
  }
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
par.fun <- function(mle){
  c(mle[1:5],logistic(mle[6:length(mle)]))
}

delt.dep <- c(T,F,F,F)
starter <- theta.s(delt.dep) # trying just state dependence on delta
timer(op <- nlm(f=ll.il.ms, # likelihood function
            p=starter, # initial values
            ch=ch,
            state.depen=delt.dep,
            hessian=TRUE))
par.fun(op$estimate)
logistic(phi.trans(op$estimate[1:3],op$estimate[4:5],list(1:(11-1),11:15)))[10:11,1:4,1]
cbind(par.nam(delt.dep),par.fun(op$estimate))

par.nam <- function(state.depen){
  al <- paste0("alpha",1:3)
  be <- paste0("beta",10:11)
  ro <- paste0("rho",c(23,4,5))
  allrep <- lapply(1:4,function(i,x,y) paste0(x[i],y[[i]]),x=c("delta","kappa","gamma","epsilon"),y=list(3:8,3:7,3:8,1:2))
  just <- c("delta","kappa","gamma","epsilon")
  return(c(al,be,ro,unlist(ifelse(state.depen,allrep,just))))
}

all.dep <- c(T,T,T,F)
starter2 <- theta.s(all.dep) # trying just state dependence on delta, kappa and gamma
timer(op2 <- nlm(f=ll.il.ms, # likelihood function
             p=starter2, # initial values
             ch=ch,
             state.depen=all.dep,
             hessian=TRUE))
# Time difference of 2.426483 mins
cbind(par.nam(all.dep),par.fun(op2$estimate))

# removing alpha/beta model and trying age classes

theta.s2 <- function(state.depen=c(F,F,F,F)){ # removing beta and increasing age classes
  n <- 9
  x <- sum(ifelse(state.depen,c(6,5,6,2),c(1,1,1,1)))
  return(logit(runif(x+n,0.1,0.9)))
}

par.nam2 <- function(state.depen){
  ph <- paste0("phi",1:6)
  ro <- paste0("rho",c(23,4,5))
  allrep <- lapply(1:4,function(i,x,y) paste0(x[i],y[[i]]),x=c("delta","kappa","gamma","epsilon"),y=list(3:8,3:7,3:8,1:2))
  just <- c("delta","kappa","gamma","epsilon")
  return(c(ph,ro,unlist(ifelse(state.depen,allrep,just))))
}

ll.il.ms2 <- function(theta,state.depen=c(F,F,F,F),ch){
  # state4 is about the four final parameters and whether they should be state
  # dependant or not (false for no and true for yes)
  struc.age <- list("age"=list(1,2,3,4,5,6:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  phi <- untrans(logistic(theta[1:6]),struc.age$age,struc.age$time,struc.age$state)
  rho.input <- theta[7:9]
  rho <- logistic(rho.trans(rho.input))
  
  Time <- ncol(ch)-1
  
  struc.state <- function(n){ # n should be 1,2,3 or 4
    if(n==1){ # delta can only vary over states 3:8
      lsss <- c(list(1:3),as.list(4:8))
    }else if(n==2){ # kappa can only vary over states 3:7
      lsss <- c(list(1:3),as.list(4:6),list(7:8))
    }else if(n==3){ # gamma can only vary over states 3:8
      lsss <- c(list(1:3),as.list(4:8))
    }else if(n==4){ # epsilon can only vary over states 1:2
      lsss <- c(list(1),list(2:8))
    }
    return(lsss)  
  }
  
  
  theta <- theta[-c(1:9)]
  
  nams <- c("delt","kap","gam","eps")
  for(i in 1:4){
    if(state.depen[i]){
      assign(nams[i],
             untrans(logistic(theta[1:length(struc.state(i))]),struc.const$age,struc.const$time,
                     struc.state(i)))
      theta <- theta[-(1:length(struc.state(i)))]
    }else{
      assign(nams[i],untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state))
      theta <- theta[-1]
    }
  }
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

two.dep <- c(T,T,F,F)
starter3 <- theta.s2(two.dep) # trying just state dependence on delta and kappa
timer(op3 <- nlm(f=ll.il.ms2, # likelihood function
                 p=starter3, # initial values
                 ch=ch,
                 state.depen=two.dep,
                 hessian=TRUE))
# Time difference of 9.764191 mins
cbind(par.nam2(two.dep),logistic(op3$estimate))

# trying one age class
theta.s3 <- function(state.depen=c(F,F,F,F)){ # removing beta and increasing age classes
  n <- 4
  x <- sum(ifelse(state.depen,c(6,5,6,2),c(1,1,1,1)))
  return(logit(runif(x+n,0.1,0.9)))
}

par.nam3 <- function(state.depen){
  ph <- paste0("phi")
  ro <- paste0("rho",c(23,4,5))
  allrep <- lapply(1:4,function(i,x,y) paste0(x[i],y[[i]]),x=c("delta","kappa","gamma","epsilon"),y=list(3:8,3:7,3:8,1:2))
  just <- c("delta","kappa","gamma","epsilon")
  return(c(ph,ro,unlist(ifelse(state.depen,allrep,just))))
}

ll.il.ms3 <- function(theta,state.depen=c(F,F,F,F),ch){
  # state4 is about the four final parameters and whether they should be state
  # dependant or not (false for no and true for yes)
  struc.age <- list("age"=list(1,2,3,4,5,6:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  phi <- untrans(logistic(theta[1]),struc.age$age,struc.age$time,struc.age$state)
  rho.input <- theta[2:4]
  rho <- logistic(rho.trans(rho.input))
  
  Time <- ncol(ch)-1
  
  struc.state <- function(n){ # n should be 1,2,3 or 4
    if(n==1){ # delta can only vary over states 3:8
      lsss <- c(list(1:3),as.list(4:8))
    }else if(n==2){ # kappa can only vary over states 3:7
      lsss <- c(list(1:3),as.list(4:6),list(7:8))
    }else if(n==3){ # gamma can only vary over states 3:8
      lsss <- c(list(1:3),as.list(4:8))
    }else if(n==4){ # epsilon can only vary over states 1:2
      lsss <- c(list(1),list(2:8))
    }
    return(lsss)  
  }
  
  
  theta <- theta[-c(1:4)]
  
  nams <- c("delt","kap","gam","eps")
  for(i in 1:4){
    if(state.depen[i]){
      assign(nams[i],
             untrans(logistic(theta[1:length(struc.state(i))]),struc.const$age,struc.const$time,
                     struc.state(i)))
      theta <- theta[-(1:length(struc.state(i)))]
    }else{
      assign(nams[i],untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state))
      theta <- theta[-1]
    }
  }
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

two.dep <- c(T,T,F,F)
starter4 <- theta.s3(two.dep) # trying just state dependence on delta and kappa
timer(op4 <- nlm(f=ll.il.ms3, # likelihood function
                 p=starter4, # initial values
                 ch=ch,
                 state.depen=two.dep,
                 hessian=TRUE))
# Time difference of 3.930135 mins
cbind(par.nam3(two.dep),logistic(op4$estimate))

# trying loads of age classes
theta.s4 <- function(state.depen=c(F,F,F,F)){ # removing beta and increasing age classes
  n <- 3 + 13
  x <- sum(ifelse(state.depen,c(6,5,6,2),c(1,1,1,1)))
  return(logit(runif(x+n,0.1,0.9)))
}

par.nam4 <- function(state.depen){
  ph <- paste0("phi",1:13)
  ro <- paste0("rho",c(23,4,5))
  allrep <- lapply(1:4,function(i,x,y) paste0(x[i],y[[i]]),x=c("delta","kappa","gamma","epsilon"),y=list(3:8,3:7,3:8,1:2))
  just <- c("delta","kappa","gamma","epsilon")
  return(c(ph,ro,unlist(ifelse(state.depen,allrep,just))))
}

ll.il.ms4 <- function(theta,state.depen=c(F,F,F,F),ch){
  # state.depen is about the four final parameters and whether they should be state
  # dependant or not (false for no and true for yes)
  struc.age <- list("age"=c(as.list(1:12),list(13:Time)),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  struc.const <- list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states)))
  phi <- untrans(logistic(theta[1:13]),struc.age$age,struc.age$time,struc.age$state)
  rho.input <- theta[14:16]
  rho <- logistic(rho.trans(rho.input))
  
  Time <- ncol(ch)-1
  
  struc.state <- function(n){ # n should be 1,2,3 or 4
    if(n==1){ # delta can only vary over states 3:8
      lsss <- c(list(1:3),as.list(4:8))
    }else if(n==2){ # kappa can only vary over states 3:7
      lsss <- c(list(1:3),as.list(4:6),list(7:8))
    }else if(n==3){ # gamma can only vary over states 3:8
      lsss <- c(list(1:3),as.list(4:8))
    }else if(n==4){ # epsilon can only vary over states 1:2
      lsss <- c(list(1),list(2:8))
    }
    return(lsss)  
  }
  
  
  theta <- theta[-c(1:16)]
  
  nams <- c("delt","kap","gam","eps")
  for(i in 1:4){
    if(state.depen[i]){
      assign(nams[i],
             untrans(logistic(theta[1:length(struc.state(i))]),struc.const$age,struc.const$time,
                     struc.state(i)))
      theta <- theta[-(1:length(struc.state(i)))]
    }else{
      assign(nams[i],untrans(logistic(theta[1]),struc.const$age,struc.const$time,struc.const$state))
      theta <- theta[-1]
    }
  }
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

two.dep <- c(T,T,F,F)
starter5 <- theta.s4(two.dep) # trying just state dependence on delta and kappa
timer(op5 <- nlm(f=ll.il.ms4, # likelihood function
                 p=starter5, # initial values
                 ch=ch,
                 state.depen=two.dep,
                 hessian=TRUE))
# Time difference of 11.37321 mins
cbind(par.nam4(two.dep),logistic(op5$estimate))

# okay now lets do some loops

# defining n as the number of unique chs as thats the smallest number
n <- nrow(count(ch))
# can use number of transitions which is largest possible n
# ss <- suff.stat(ch) # requires suff.stat from make.suff.stat.R
# n <- sum(sapply(1:16,function(x) colSums(ss$v,dims=1)[x,1]*(16-x) )) + sum(ss$m) 
# or could use number of individuals
# n <- nrow(ch)
dep.as.str <- sapply(1:16,function(x) paste((as.integer(intToBits(x)[1:4])), collapse=""))
for(i in 1:16){
  dep <- as.logical(as.numeric(unlist(strsplit(dep.as.str[i],""))))
  k <- length(par.nam(dep))
  starter <- theta.s(dep)
  ts <- Sys.time()
  assign(paste0("opl",i),nlm(f=ll.il.ms, # likelihood function
                             p=starter, # initial values
                             ch=ch,
                             state.depen=dep,
                             hessian=TRUE))
  te <- Sys.time()
  print(difftime(te,ts,units = "mins"))
  assign(paste0("df",i),data.frame("par"=par.nam(dep),
                                   "mle"=par.fun(eval(as.name(paste0("opl",i)))$estimate),
                                   "state.dep"=dep.as.str[i],
                                   "convergence"=eval(as.name(paste0("opl",i)))$code,
                                   "value"=eval(as.name(paste0("opl",i)))$minimum,
                                   "aic"=(eval(as.name(paste0("opl",i)))$minimum + 2*k),
                                   "aicc"=(eval(as.name(paste0("opl",i)))$minimum + 2*k) + ((2*k^2+2*k)/(n-k-1)),
                                   "time"=as.numeric(difftime(te,ts,units = "mins"))))
  print(i)
}

df <- do.call("rbind",lapply(1:16,function(x) eval(as.name(paste0("df",x))) ))
ggplot(df,aes(par,mle)) + geom_point() + facet_wrap(~state.dep) + lims(y=c(0,1))

for(i in 1:16){
  dep <- as.logical(as.numeric(unlist(strsplit(dep.as.str[i],""))))
  k <- length(par.nam(dep))
  starter <- theta.s(dep)
  ts <- Sys.time()
  assign(paste0("opl2",i),optim(starter,
                                ll.il.ms, # likelihood function
                                ch=ch,
                                state.depen=dep,
                                hessian=TRUE))
  te <- Sys.time()
  print(difftime(te,ts,units = "mins"))
  assign(paste0("df2",i),data.frame("par"=par.nam(dep),
                                    "mle"=par.fun(eval(as.name(paste0("opl2",i)))$par),
                                    "state.dep"=dep.as.str[i],
                                    "convergence"=eval(as.name(paste0("opl2",i)))$convergence,
                                    "value"=eval(as.name(paste0("opl2",i)))$value,
                                    "aic"=(eval(as.name(paste0("opl2",i)))$value + 2*k),
                                    "aicc"=(eval(as.name(paste0("opl2",i)))$value + 2*k) + ((2*k^2+2*k)/(n-k-1)),
                                    "time"=as.numeric(difftime(te,ts,units = "mins"))))
  print(i)
}

df.2 <- do.call("rbind",lapply(1:6,function(x) eval(as.name(paste0("df2",x))) ))

