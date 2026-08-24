states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
states_code <- c(1,13,2,3,5,4,6)
# no skipped code needed so this vector is one shorter length than states
# the 13 is so no element of ch will match and get E and keep the indicies correct
load("matrix_AUK.RData")
ch <- matrix_AUK
Time <- ncol(ch)
for(i in 1:length(states_code)){
  ch[which(ch==states_code[i],arr.ind=TRUE)] <- states[i]
}
for(i in 1:nrow(ch)){
  # want to add Nbar as a state between any Ns or until N becomes B1
  # and want to add S as a state between the breeding states
  wn <- which(ch[i,]=="N")
  wb1 <- which(ch[i,]=="B1")
  wnlss <- wn[wn<wb1] # which non breeding states are before the first breeding state (as opposed to the non breeding due to skipping)
  if(length(wnlss)==0){
    if(length(wn)>=2){
      for(j in 1:(Time-1)){
        if(j>min(wn) & j<max(wn) & ch[i,j]=="0"){
          ch[i,j] <- "E"
        }
      }
    }
  }else{
    if(length(wnlss)>0){
      for(j in 1:(Time-1)){
        if(j>min(wn) & j<wb1 & ch[i,j]=="0"){
          ch[i,j] <- "E"
        }
      }
    }
  }
  wbb <- which(ch[i,] %in% c("B1","LB","L_B","LB_","L_B_"))
  if(length(wbb)>1){
    for(j in 1:(length(wbb)-1)){
      if((wbb[j+1]-wbb[j])>1){
        ch[i,(wbb[j]+1):(wbb[j+1]-1)] <- "S"
      }
    }
  }
}

ni <- rep(0,16)
for(i in 1:nrow(ch)){
  ni[which(ch[i,1:16]!=0)[1]] <- ni[which(ch[i,1:16]!=0)[1]]+1
}
ni2 <- list()
for(i in 1:length(ni)){
  ni2[[i]] <- rep(1,ni[i])
}

ch <- as.data.frame(ch)
ch <- cbind(ch,"age"=rep(1,nrow(ch)))
ch2 <- ch

ch <- ch[-c(253, 293, 374, 398, 611, 69, 492, 125, 393, 422, 170, 129),]

save(ch,ch2,ni,ni2,file="guillemot.RData")