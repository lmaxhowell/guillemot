load("new_data.RData")
states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
Time <- ncol(new_data)-1

ch3 <- new_data[,1:16]
colnames(ch3) <- c(2010:2025)
states_code <- c(1,13,"B1","LSB","DSB","LUB","DUB",13)

for(s in 1:length(states_code)){
  ch3[which(ch3==states_code[s],arr.ind=TRUE)] <- states[s]
}

for(i in 1:nrow(ch3)){
  # want to add Nbar as a state between any Ns or until N becomes B1
  # and want to add S as a state between the breeding states
  wn <- which(ch3[i,]=="N")
  wb1 <- which(ch3[i,]=="B1")
  wnlss <- wn[wn<wb1] # which non breeding states are before the first breeding state (as opposed to the non breeding due to skipping)
  if(length(wnlss)==0){
    if(length(wn)>=2){
      for(j in 1:(Time-1)){
        if(j>min(wn) & j<max(wn) & ch3[i,j]=="0"){
          ch3[i,j] <- "E"
        }
      }
    }
  }else{
    if(length(wnlss)>0){
      for(j in 1:(Time-1)){
        if(j>min(wn) & j<wb1 & ch3[i,j]=="0"){
          ch3[i,j] <- "E"
        }
      }
    }
  }
  wbb <- which(ch3[i,] %in% c("B1","LB","L_B","LB_","L_B_"))
  if(length(wbb)>1){
    for(j in 1:(length(wbb)-1)){
      if((wbb[j+1]-wbb[j])>1){
        ch3[i,(wbb[j]+1):(wbb[j+1]-1)] <- "S"
      }
    }
  }
}

ch3 <- cbind(ch3,"age"=rep(1,nrow(ch3)))

ni3 <- rep(0,16)
for(i in 1:nrow(ch3)){
  ni3[which(ch3[i,1:16]!=0)[1]] <- ni3[which(ch3[i,1:16]!=0)[1]]+1
}

ni4 <- list()
for(i in 1:length(ni3)){
  ni4[[i]] <- rep(1,ni3[i])
}

save(ch3,ni3,ni4,file="guillemot2.RData")