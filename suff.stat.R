suff.stat <- function(ch){
  states <- c("N","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  B <- states[-1]
  
  Agesi <- ch[,ncol(ch)]
  ch <- ch[,-ncol(ch)]
  
  Age <- Time + max(Agesi) - 1
  
  # m is r->s
  # v is r->0
  mB <- array(0,dim=c(length(states),length(states),Time,Age)) # r,s,t,a
  mN <- array(0,dim=c(2,Time,Age,Time)) # s,t,a,t' # 2 rather than length states because those are the only states N can go to
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
        t <- wh[j]
        a <- Agesi[i]+t-wh[1]
        if(s %in% 1:2){
          tprime <- wh[j+1]
          mN[s,t,a,tprime] <- mN[s,t,a,tprime] + 1
        }else{
          mB[r,s,t,a] <- mB[r,s,t,a] + 1
        }

        if(isitavio(r,s)){
          print(c(r,s))
          violations[[r]][[s]] <- append(violations[[r]][[s]],i)
        }
      }
      r <- which(states==ch[i,wh[length(wh)]])
      v[r,wh[length(wh)],wh[length(wh)]+Agesi[i]-1] <- v[r,wh[length(wh)],wh[length(wh)]+Agesi[i]-1] + 1
    }
  }
  mB2 <- rowSums(mB,dims = 2)
  dimnames(m2) <- list(states,states)
  v2 <- rowSums(v,dims = 1)
  names(v2) <- states
  
  mN2 <- array(0,dim=c(2,Time)) # s,t'
  for(i in 1:dim(mN)[2]){
    for(j in 1:dim(mN)[3]){
      mN2 <- mN2 + mN[,i,j,]
    }
  }
  
  return(list("mB"=mB,"mN"=mN,"v"=v,"mB2"=mB2,"mN2"=mN2,"v2"=v2,"violations"=violations))
}