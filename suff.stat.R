suff.stat <- function(ch){
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(ch)-1
  
  Agesi <- ch[,ncol(ch)]
  ch <- ch[,-ncol(ch)]
  
  Age <- Time + max(Agesi) - 1
  
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
        # if((j+Agesi[i]-1)==0){
        #   print(c(i,j))
        # }
        m[r,s,wh[j],j+Agesi[i]-1] <- m[r,s,wh[j],j+Agesi[i]-1] + 1
        if(isitavio(r,s)){
          print(c(r,s))
          violations[[r]][[s]] <- append(violations[[r]][[s]],i)
        }
      }
      r <- which(states==ch[i,wh[length(wh)]])
      v[r,wh[length(wh)],length(wh)+Agesi[i]-1] <- v[r,wh[length(wh)],length(wh)+Agesi[i]-1] + 1
    }
  }
  m2 <- rowSums(m,dims = 2)
  dimnames(m2) <- list(states,states)
  v2 <- rowSums(v,dims = 1)
  names(v2) <- states
  
  return(list("m"=m,"v"=v,"m2"=m2,"v2"=v2,"violations"=violations))
}