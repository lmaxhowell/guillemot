rm(list=ls())
file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed

# note that we gotta make sure that ll.il is returning a vector of the m 
# related likelihood and v likelihood
ch <- read.csv("sim.data.csv")
ch <- ch[,-1]

ch2 <- read.csv("sim.data.csv")
ch2 <- ch2[,-c(1,18)]

states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
ch3 <- as.data.frame(array(0,dim=dim(ch2)))
for(j in 1:nrow(ch2)){
  for(i in 1:length(states)){
    ch3[j,which(ch2[j,]==states[i])] <- as.numeric(i)
  }
  ch3[j,which(ch2[j,]=="0")] <-as.numeric(0)
}

nobs <- nrow(ch3)
Time <- ncol(ch3)
m <- array(0,dim=c(8,8))
v <- array(0,dim=c(8,Time))

for (i in 1:nobs){
  for (t in 1:(Time-1)){
    r <- ch3[i,t]
    s <- ch3[i,t+1]
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
  L <- (L1+L2)
  # return(L)
  return(c(L,L1,L2))
}

theta <- logit(c(0.7,0.3,0.2,0.5,0.8,0.6))
struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))),
              "eps"=list("age"=list(1:Time),"time"=list(1:Time),"state"=list(1:length(states))))

# total likelihoods
lrt <- MS_model(theta,ch=ch3,m=m,v=v)
lmt <- ll.il(theta,1,2,3,4,5,6,struc,ch)

df <- as.data.frame(array(0,dim=c(nrow(ch),7)))
colnames(df) <- c("rach.l","rach.lm","rach.lv","me.l","me.lm","me.lv","diff")
for(i in 1:nrow(ch)){
  # first construct rachels suff stats
  m2 <- array(0,dim=c(8,8))
  v2 <- array(0,dim=c(8,Time))
  for (t in 1:(Time-1)){
    r <- ch3[i,t]
    s <- ch3[i,t+1]
    if (r>0 & s>0){
      m2[r,s] <- m2[r,s]+1
    }
    if (r>0 && s==0){
      v2[r,t] <- v2[r,t]+1
    }
  }
  ch.i <- ch[c(i,i),]
  ch.i2 <- ch3[i,]
  lr <- MS_model(theta,ch=ch.i2,m=m2,v=v2)
  lm <- ll.il(theta,1,2,3,4,5,6,struc,ch.i)/2
  # if(lr[1]!=lm[1]){
  #   print(c(lr[1],lm[1]))
  #   print(ch.i)
  # }
  df[i,] <- c(lr,lm[1:3],lr[1]-lm[1])
  cat("\ritr",i," of ",nrow(ch))
}

df2 <- cbind(df,"zero"=(round(df$diff,6)==0))
table(df2$zero)
# 137 rows not the same
df2[!df2$zero,]
View(ch[!df2$zero,])


df3 <- as.data.frame(array(0,dim=c(nrow(ch),7)))
colnames(df3) <- c("rach.l","rach.lm","rach.lv","me.l","me.lm","me.lv","diff")
for(i in 1:nrow(ch)){
  # first construct rachels suff stats
  m2 <- array(0,dim=c(8,8))
  v2 <- array(0,dim=c(8,Time))
  for (t in 1:(Time-1)){
    r <- ch3[i,t]
    s <- ch3[i,t+1]
    if (r>0 & s>0){
      m2[r,s] <- m2[r,s]+1
    }
    if (r>0 && s==0){
      v2[r,t] <- v2[r,t]+1
    }
  }
  ch.i <- ch[c(i,i),]
  ch.i2 <- ch3[i,]
  lr <- MS_model(theta,ch=ch.i2,m=m2,v=v2)
  lm <- ll.il(theta,1,2,3,4,5,6,struc,ch.i)/2
  # if(lr[1]!=lm[1]){
  #   print(c(lr[1],lm[1]))
  #   print(ch.i)
  # }
  df3[i,] <- c(lr,lm[1:3],lr[1]-lm[1])
  cat("\ritr",i," of ",nrow(ch))
}

df4 <- cbind(df3,"zero"=(round(df3$diff,6)==0))
table(df4$zero)
df4[!df4$zero,]
View(ch[!df4$zero,])



