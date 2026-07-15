ch2 <- read.csv("sim.data.csv")
ch2 <- ch2[,-c(1,18)]

states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
ch <- as.data.frame(array(0,dim=dim(ch2)))
for(j in 1:nrow(ch2)){
  for(i in 1:length(states)){
    ch[j,which(ch2[j,]==states[i])] <- as.numeric(i)
  }
  ch[j,which(ch2[j,]=="0")] <-as.numeric(0)
}

nobs <- nrow(ch)
Time <- ncol(ch)
m <- array(0,dim=c(8,8))
v <- array(0,dim=c(8,Time))

for (i in 1:nobs){
  for (t in 1:(Time-1)){
    r <- ch[i,t]
    s <- ch[i,t+1]
    if (r>0 & s>0){
      m[r,s] <- m[r,s]+1
    }
    if (r>0 && s==0){
      v[r,t] <- v[r,t]+1
    }
  }
}

MS_model <- function(theta,ch,m,v){
  nobs <-nrow(ch)
  Time<-ncol(ch)
  
  phi<-1/(1+exp(-theta[1]))
  delta<-1/(1+exp(-theta[2]))
  kappa<-1/(1+exp(-theta[3]))
  rho<-1/(1+exp(-theta[4]))
  gam<-1/(1+exp(-theta[5]))
  epsilon<-1/(1+exp(-theta[6]))
  
  psi<-array(0,dim=c(8,8))
  for (r in 1:2){
    psi[r,1]<-(1-epsilon)*(1-rho)
    psi[r,2]<-epsilon
    psi[r,3]<-(1-epsilon)*rho
  }
  for (r in 3:7){
    psi[r,4]<-(1-kappa)*(1-delta)*gam
    psi[r,5]<-(1-kappa)*delta*gam
    psi[r,6]<-(1-kappa)*(1-delta)*(1-gam)
    psi[r,7]<-(1-kappa)*delta*(1-gam)
    psi[r,8]<-kappa
  }
  psi[8,4]<-(1-delta)*gam
  psi[8,5]<-delta*gam
  psi[8,6]<-(1-delta)*(1-gam)
  psi[8,7]<-delta*(1-gam)
  
  prob_m <- array(0,dim=c(8,8))
  
  for (r in 1:8){
    for (s in 1:8){
      prob_m[r,s]<-phi*psi[r,s]
    }
  }
  
  chi <- array(1,dim=c(8,Time))
  for (t in (Time-1):1){
    for (r in 3:7){
      chi[r,t] <- (1-phi)+phi*psi[r,8]*(1-phi)
    }
    chi[1,t]<- (1-phi)+phi*epsilon*chi[1,t+1]
  }
  
  L1<-0
  L2<-0
  
  for (r in 1:8){
    for (s in 1:8){
      probm<-prob_m[r,s]
      if (probm>0){
        L1 <- L1 + m[r,s]*log(probm)
      }
    }
  }
  
  for (r in 1:8){
    for (t in 1:Time){
      probchi<-chi[r,t]
      if (probchi>0){
        L2 <- L2 + v[r,t]*log(probchi)
      }
    }
  }
  L <- -(L1+L2)
  return(L)
  # return(c(L,L1,L2))
}

L <-  MS_model(c(0.847297860387203, -0.847297860387204, -1.38629436111989, 0, 
                 1.38629436111989, 0.405465108108164),ch=ch,m=m,v=v)
print(L)

MSout<-optim(par = c(0,0,0,0,0,0), fn = MS_model, ch = ch, m=m,v=v,control = list(maxit = 1000000),method="BFGS")

MSout$convergence
1/(1+exp(-MSout$par))
MSout$value

obs_delta <- 1-(sum(ch==4)+sum(ch==6))/(sum(ch==4)+sum(ch==5)+sum(ch==6)+sum(ch==7))
obs_gamma <- (sum(ch==4)+sum(ch==5))/(sum(ch==4)+sum(ch==5)+sum(ch==6)+sum(ch==7))

obs_delta
obs_gamma