Pr_rs <- function(r,s,phi,psi){
  prob <- phi[r]*psi[r,s]
  return(prob)
}

Pr_r0s <- function(r,s,phi,psi){
  skip <- which(row.names(psi)=="S")
  prob <- phi[r]*psi[r,skip]*phi[skip]*psi[skip,s]
  return(prob)
}

Pr_r0 <- function(r,phi,psi){
  skip <- which(row.names(psi)=="S")
  # prob <- phi[r]*psi[r,skip]*(1-phi[skip]) + (1-phi[skip]) # this is the OG from the document
  # prob <- phi[r]*psi[r,skip] + (1-phi[skip]) # this is what I thought
  prob <- phi[r]*psi[r,skip]*(1-phi[skip]) + (1-phi[r]) # this is third one/blue equation
  return(prob)
}

make.psi <- function(delta,kap,rho,gam){
  psi <- array(0,dim=c(length(states),length(states)))
  psi[1,] <- c(1-rho[1],rho[2],0,0,0,0,0)
  row.names(psi) <- states
  colnames(psi) <- states
  for(i in 3:7){
    psi[i-1,] <- c(0,0,
                   (1-kap[i])*(1-delta[i])*gam[i],
                   (1-kap[i])*delta[i]*gam[i],
                   (1-kap[i])*(1-delta[i])*(1-gam[i]),
                   (1-kap[i])*delta[i]*(1-gam[i]),
                   kap[i])
  }
  psi[7,] <- c(0,0,
               (1-delta[7])*gam[7],
               delta[7]*gam[7],
               (1-delta[7])*(1-gam[7]),
               delta[7]*(1-gam[7]),
               0)
  return(psi)
}

# previous il that tries to go though things according to time
# il <- function(ch,phi,psi){ # il is "individual likelihood"
#   Indicator <- function(prob){
#     # a function that checks if the probability pu into it will
#     # "cause problems" in the likelihood - i.e. if we log
#     # this probability will it be infinite and cause the whole ll
#     # to be infinite as a result
#     # or has a small amount of numerical error crept in
#     # and made a probability negative?
#     # in which case just set it to zero
#     # will add a warning in in this case
#     # this function just means we can call log(prob)
#     # without worrying about the errors this can cause
#     if(prob<0){
#       return(0)
#     }else{
#      return(ifelse(is.finite(log(prob)),log(prob),0)) 
#     }
#   }
#   library(plyr, include.only = c("count"))
#   # just need the count function from this package
#   # to figure out how many unique chs there are
#   uch <- count(ch) # unique number of capture histories
#   # what are the states
#   states <- c("N", "B1","LB","L_B","LB_","L_B_","S")
#   # number of columns of the chs is Time
#   Time <- ncol(ch)
#   # how many unique capture histories are there
#   Individuals <- nrow(uch)
#   ll_i <- rep(0,Individuals) # want to create an individual likelihood for each capture history
#   for(i in 1:Individuals){
#     ch <- uch[i,1:Time] # the current capture history we are looking at
#     # need to deal with all the capture histories that
#     # only have one state in them first
#     if(sum(ch!="0")==1){
#       current_state <- ch[ch!="0"]
#       if(which(ch!="0")==Time){ # if this occurs at the final time...do nothing
#         next
#       }else{ # else we are in a Pr_r0 situation
#         current_state_index <- which(states==current_state) # what number is the current state
#         ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,phi,psi))
#       }
#     }else{ # therefore ch has multiple observations
#       for(t in 1:(Time-1)){# for each time, compare with the next time point 
#         current_state <- unlist(ch[t]) # what state are we in right now?
#         next_state <- unlist(ch[t+1]) # what is the next state?
#         current_state_index <- which(states==current_state) # what number is the current state
#         next_state_index <- which(states==next_state) # what number is the next state
#         # print(c(current_state,next_state,current_state_index,next_state_index))
#         if(current_state!="0"){ # only want to evaluate if we are in an observed state
#           if(next_state=="0"){ # if next state is a zero then check
#             # if we are in the final time as this is the only time
#             # that Pr_r0 can happen, so we shall deal with that first
#             if((t+1)==Time){
#               ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,phi,psi))
#             }else{ # otherwise its a Pr_r0s
#               print(ch)
#               plus <- 2
#               next_state <- ch[t+plus] # what is the NEXT next state?
#               while(next_state=="0"){
#                 plus <- plus + 1
#                 print(plus)
#                 next_state <- ch[t+plus] # what is the NEXT next state?
#               }
#               next_state_index <- which(states==next_state) # what number is the next state
#               # next_state_index <- ifelse(next_state=="0",stop(print(ch)),which(states==next_state)) # what number is the next state
#               ll_i[i] <- ll_i[i] + Indicator(Pr_r0s(current_state_index,next_state_index,phi,psi))
#               next # might need a double next
#             } # end else t+1==Time
#           }else if(next_state!="0"){ # if it isnt a zero then we have an r->s transition
#             ll_i[i] <- ll_i[i] + Indicator(Pr_rs(current_state_index,next_state_index,phi,psi))
#           } # end else next state not equal zero
#         } # end if current state not zero
#       } # end t
#     } # end i
#   } # end else ch has multiple observations
#   ll <- sum(uch$freq*ll_i) # multiply the log-likelihood for each unique capture history
#   # by how many times that capture history occurs and sum it up
#   # to make the multinomial loglikelihood
#   return(ll)
# }

# find all the state transitions in the relevant capture history
find.transitions <- function(ch){
  Time <- length(ch)
  where <- which(ch!="0")
  df <- as.data.frame(array(NA,dim=c(max(length(where)-1,1),4))) # want it to be a character df but full of nothing
  colnames(df) <- c("r","s","t_r","t_s") # state r, state s, time at state r, time at states
  if(nrow(df)==1){ # if only one observation
    if(length(where)==1){ # if only one state observed
      df[1,] <- list(ch[ch!="0"],"0",where,where+1)
      return(df)
    }else if(length(where)==2){ # the case where there is a state at the final time
      df[1,] <- list(ch[where[1]],ch[where[2]],where[1],where[2])
      return(df)
    }
  }else{ # more than one observation
    for(t in 1:(length(where)-1)){ # go through all the states
      df[t,] <- list(ch[where[t]],ch[where[t+1]],where[t],where[t+1])
    } # end for t in 1:length(where)-1
    lw <- where[length(where)] # the last element of where
    if(lw!=Time){ # if the last place that a state is observed ISNT the last time
      df <- rbind(df,
                  list(ch[lw],"0",lw,lw+1))
    }
    return(df)
  } # end else more than one observation
}

# new il that tries to go though things according to length of observed states
il <- function(ch,phi,psi){ # il is "individual likelihood"
  Indicator <- function(prob){
    # a function that checks if the probability pu into it will
    # "cause problems" in the likelihood - i.e. if we log
    # this probability will it be infinite and cause the whole ll
    # to be infinite as a result
    # or has a small amount of numerical error crept in
    # and made a probability negative?
    # in which case just set it to zero
    # will add a warning in in this case
    # this function just means we can call log(prob)
    # without worrying about the errors this can cause
    if(length(prob)==0){
      print(c(parent.frame()$i,parent.frame()$t,parent.frame()$Time))
      View(parent.frame()$transitions)
    }
    if(prob<0){
      return(0)
    }else{
      return(ifelse(is.finite(log(prob)),log(prob),0)) 
    }
  }
  library(plyr, include.only = c("count"))
  # just need the count function from this package
  # to figure out how many unique chs there are
  uch <- count(ch) # unique number of capture histories
  # what are the states
  states <- c("N", "B1","LB","L_B","LB_","L_B_","S")
  # number of columns of the chs is Time
  Time <- ncol(ch)
  # how many unique capture histories are there
  Individuals <- nrow(uch)
  ll_i <- rep(0,Individuals) # want to create an individual likelihood for each capture history
  for(i in 1:Individuals){
    ch <- uch[i,1:Time] # the current capture history we are looking at
    transitions <- find.transitions(ch)
    # need to deal with all the capture histories that
    # only have one state in them first
    if(nrow(transitions)==1){
      current_state <- transitions[1,1]
      if(transitions[1,4]==(Time+1)){ # if this occurs at the final time...do nothing
        next
      }else{ # else we are in a Pr_r0 situation
        current_state_index <- which(states==current_state) # what number is the current state
        ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,phi,psi))
      }
    }else{ # therefore ch has multiple observations
      for(t in 1:nrow(transitions)){
        if(transitions[t,4]==(Time+1)){ # in the situation where the final time point is an observed state and not a zero
          next # then we just skip it
        }else if((transitions[t,4]-transitions[t,3])==1){ # aka are the states next to each other and not a state in the final time
          current_state_index <- which(states==transitions[t,1]) # what number is the current state
          # need to deal with if the transition is the final one, a state to zero
          if(transitions[t,2]=="0"){ # skips shouldnt end up in this loop because they wont have sequential times
            ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,phi,psi))
          }else{
            next_state_index <- which(states==transitions[t,2]) # what number is the next state
            ll_i[i] <- ll_i[i] + Indicator(Pr_rs(current_state_index,next_state_index,phi,psi))
          } # end else
        }else if((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]!="N"){ # we have a skipped a time and it was previously breeding
          current_state_index <- which(states==transitions[t,1]) # what number is the current state
          next_state_index <- which(states==transitions[t,2]) # what number is the next state
          ll_i[i] <- ll_i[i] + Indicator(Pr_r0s(current_state_index,next_state_index,phi,psi))
        }else if(((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]=="N") | 
                 ((transitions[t,4]-transitions[t,3])>2)){
          # we have a skipped a time and it was NOT previously breeding OR
          # we have skipped multiple times (and so theoretically should not be breeding) OR
          current_state_index <- which(states==transitions[t,1]) # what number is the current state
          ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,phi,psi))
        }
        # print(ll_i)
      } # end for t in transitions
    } # end else ch has multiple observations
  } # end i
  ll <- sum((uch$freq)*(ll_i)) # multiply the log-likelihood for each unique capture history
  # by how many times that capture history occurs and sum it up
  # to make the multinomial loglikelihood
  return(ll)
}

ll.il <- function(theta,phi.ind,delt.ind,kap.ind,rho.ind,gam.ind,struc,ch){
  states <- c("N", "B1","LB","L_B","LB_","L_B_","S")
  
  phi <- untrans(logistic(theta[phi.ind]),list(1),struc$phi)
  delt <- untrans(logistic(theta[delt.ind]),list(1),struc$delt)
  kap <- untrans(logistic(theta[kap.ind]),list(1),struc$kap)
  rho <- untrans(logistic(theta[rho.ind]),list(1),struc$rho)
  gam <- untrans(logistic(theta[gam.ind]),list(1),struc$gam)
  
  psi <- make.psi(delt,kap,rho,gam)
  
  ll <- il(ch,phi,psi)
  
  return(ll)
}

Pr_rs <- function(r,s,t,phi,psi){
  prob <- phi[r,t]*psi[r,s,t]
  return(prob)
}

Pr_r0s <- function(r,s,t,phi,psi){
  skip <- which(row.names(psi)=="S")
  prob <- phi[r,t]*psi[r,skip,t]*phi[skip,t]*psi[skip,s,t]
  return(prob)
}

Pr_r0 <- function(r,t,phi,psi){
  skip <- which(row.names(psi)=="S")
  # prob <- phi[r]*psi[r,skip]*(1-phi[skip]) + (1-phi[skip]) # this is the OG from the document
  # prob <- phi[r]*psi[r,skip] + (1-phi[skip]) # this is what I thought
  prob <- phi[r,t]*psi[r,skip,t]*(1-phi[skip,t]) + (1-phi[r,t]) # this is third one/blue equation
  return(prob)
}

make.psi <- function(delta,kap,rho,gam){
  Time <- ncol(delta)
  if(ncol(kap)!=Time | ncol(rho)!=Time | ncol(gam)!=Time){
    stop("delta,kappa,rho and gamma must all have the same number of columns, corresponding to time")
  }
  psi <- array(0,dim=c(length(states),length(states),Time))
  for(t in 1:Time){
    psi[1,,t] <- c(1-rho[1,t],rho[2,t],0,0,0,0,0)
    row.names(psi) <- states
    colnames(psi) <- states
    for(i in 3:7){
      psi[i-1,,t] <- c(0,0,
                     (1-kap[i,t])*(1-delta[i,t])*gam[i,t],
                     (1-kap[i,t])*delta[i,t]*gam[i,t],
                     (1-kap[i,t])*(1-delta[i,t])*(1-gam[i,t]),
                     (1-kap[i,t])*delta[i,t]*(1-gam[i,t]),
                     kap[i,t])
    }
    psi[7,,t] <- c(0,0,
                 (1-delta[7,t])*gam[7,t],
                 delta[7,t]*gam[7,t],
                 (1-delta[7,t])*(1-gam[7,t]),
                 delta[7,t]*(1-gam[7,t]),
                 0)
  }
  return(psi)
}

# want to do the same as previous il but make it time dependant
il <- function(ch,phi,psi){ # il is "individual likelihood"
  Indicator <- function(prob){
    # a function that checks if the probability pu into it will
    # "cause problems" in the likelihood - i.e. if we log
    # this probability will it be infinite and cause the whole ll
    # to be infinite as a result
    # or has a small amount of numerical error crept in
    # and made a probability negative?
    # in which case just set it to zero
    # will add a warning in in this case
    # this function just means we can call log(prob)
    # without worrying about the errors this can cause
    if(length(prob)==0){
      print(c(parent.frame()$i,parent.frame()$t,parent.frame()$Time))
      View(parent.frame()$transitions)
    }
    if(prob<0){
      return(0)
    }else{
      return(ifelse(is.finite(log(prob)),log(prob),0)) 
    }
  }
  library(plyr, include.only = c("count"))
  # just need the count function from this package
  # to figure out how many unique chs there are
  uch <- count(ch) # unique number of capture histories
  # what are the states
  states <- c("N", "B1","LB","L_B","LB_","L_B_","S")
  # number of columns of the chs is Time
  Time <- ncol(ch)
  # how many unique capture histories are there
  Individuals <- nrow(uch)
  ll_i <- rep(0,Individuals) # want to create an individual likelihood for each capture history
  for(i in 1:Individuals){
    ch <- uch[i,1:Time] # the current capture history we are looking at
    transitions <- find.transitions(ch)
    # need to deal with all the capture histories that
    # only have one state in them first
    if(nrow(transitions)==1){
      current_state <- transitions[1,1]
      current_time <- transitions[1,3]
      if(transitions[1,4]==(Time+1)){ # if this occurs at the final time...do nothing
        next
      }else{ # else we are in a Pr_r0 situation
        current_state_index <- which(states==current_state) # what number is the current state
        ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,phi,psi))
      }
    }else{ # therefore ch has multiple observations
      for(t in 1:nrow(transitions)){
        if(transitions[t,4]==(Time+1)){ # in the situation where the final time point is an observed state and not a zero
          next # then we just skip it
        }else if((transitions[t,4]-transitions[t,3])==1){ # aka are the states next to each other and not a state in the final time
          current_state_index <- which(states==transitions[t,1]) # what number is the current state
          current_time <- transitions[t,3]
          # need to deal with if the transition is the final one, a state to zero
          if(transitions[t,2]=="0"){ # skips shouldnt end up in this loop because they wont have sequential times
            ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,phi,psi))
          }else{
            next_state_index <- which(states==transitions[t,2]) # what number is the next state
            ll_i[i] <- ll_i[i] + Indicator(Pr_rs(current_state_index,next_state_index,current_time,phi,psi))
          } # end else
        }else if((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]!="N"){ # we have a skipped a time and it was previously breeding
          current_state_index <- which(states==transitions[t,1]) # what number is the current state
          next_state_index <- which(states==transitions[t,2]) # what number is the next state
          current_time <- transitions[t,3]
          ll_i[i] <- ll_i[i] + Indicator(Pr_r0s(current_state_index,next_state_index,current_time,phi,psi))
        }else if(((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]=="N") | 
                 ((transitions[t,4]-transitions[t,3])>2)){
          # we have a skipped a time and it was NOT previously breeding OR
          # we have skipped multiple times (and so theoretically should not be breeding) OR
          current_state_index <- which(states==transitions[t,1]) # what number is the current state
          current_time <- transitions[t,3]
          ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,phi,psi))
        }
        # print(ll_i)
      } # end for t in transitions
    } # end else ch has multiple observations
  } # end i
  ll <- sum((uch$freq)*(ll_i)) # multiply the log-likelihood for each unique capture history
  # by how many times that capture history occurs and sum it up
  # to make the multinomial loglikelihood
  return(ll)
}

ll.il <- function(theta,phi.ind,delt.ind,kap.ind,rho.ind,gam.ind,struc,ch){
  states <- c("N", "B1","LB","L_B","LB_","L_B_","S")
  
  phi <- untrans(logistic(theta[phi.ind]),struc$phi$time,struc$phi$state)
  delt <- untrans(logistic(theta[delt.ind]),struc$delt$time,struc$delt$state)
  kap <- untrans(logistic(theta[kap.ind]),struc$kap$time,struc$kap$state)
  rho <- untrans(logistic(theta[rho.ind]),struc$rho$time,struc$rho$state)
  gam <- untrans(logistic(theta[gam.ind]),struc$gam$time,struc$gam$state)
  
  psi <- make.psi(delt,kap,rho,gam)
  # print(psi)
  
  ll <- il(ch,phi,psi)
  
  return(ll)
}

#################################################
# want to add in Nbar state for the non breeders
#################################################
make.psi <- function(delta,kap,rho,gam,upsilon){
  states <- c("N","N_","B1","LB","L_B","LB_","L_B_","S")
  Time <- ncol(delta)
  if(ncol(kap)!=Time | ncol(rho)!=Time | ncol(gam)!=Time | ncol(upsilon)!=Time){
    stop("delta,kappa,rho, gamma and upsilon must all have the same number of columns, corresponding to time")
  }
  psi <- array(0,dim=c(length(states),length(states),Time))
  for(t in 1:Time){
    psi[1,,t] <- c((1-rho[1,t])*upsilon[1,t],(1-rho[2,t])*(1-upsilon[2,t]),rho[3,t],0,0,0,0,0)
    psi[2,,t] <- c((1-rho[1,t])*upsilon[1,t],(1-rho[2,t])*(1-upsilon[2,t]),rho[3,t],0,0,0,0,0)
    row.names(psi) <- states
    colnames(psi) <- states
    # this is if these rows depend on the ROWS theyre in
    # for(i in 4:8){
    #   psi[i-1,,t] <- c(0,0,0,
    #                    (1-kap[i,t])*(1-delta[i,t])*gam[i,t],
    #                    (1-kap[i,t])*delta[i,t]*gam[i,t],
    #                    (1-kap[i,t])*(1-delta[i,t])*(1-gam[i,t]),
    #                    (1-kap[i,t])*delta[i,t]*(1-gam[i,t]),
    #                    kap[i,t])
    # }
    # psi[8,,t] <- c(0,0,0,
    #                (1-delta[8,t])*gam[8,t],
    #                delta[8,t]*gam[8,t],
    #                (1-delta[8,t])*(1-gam[8,t]),
    #                delta[8,t]*(1-gam[8,t]),
    #                0)
    # this is if these rows depend on the columns theyre in
    for(i in 4:8){
      psi[i-1,,t] <- c(0,0,0,
                       (1-kap[4,t])*(1-delta[4,t])*gam[4,t],
                       (1-kap[5,t])*delta[5,t]*gam[5,t],
                       (1-kap[6,t])*(1-delta[6,t])*(1-gam[6,t]),
                       (1-kap[7,t])*delta[7,t]*(1-gam[7,t]),
                       kap[8,t])
    }
    psi[8,,t] <- c(0,0,0,
                   (1-delta[4,t])*gam[4,t],
                   delta[5,t]*gam[5,t],
                   (1-delta[6,t])*(1-gam[6,t]),
                   delta[7,t]*(1-gam[7,t]),
                   0)
  }
  return(psi)
}

Pr_rs <- function(r,s,t,phi,psi){
  prob <- phi[r,t]*psi[r,s,t]
  return(prob)
}

Pr_r0 <- function(r,t,phi,psi){
  skip <- which(row.names(psi)=="S")
  prob <- phi[r,t]*psi[r,skip,t]*(1-phi[skip,t]) + (1-phi[r,t])
  return(prob)
}

# want to add in Nbar state
il <- function(ch,phi,psi){ # il is "individual likelihood"
  Indicator <- function(prob){
    # a function that checks if the probability pu into it will
    # "cause problems" in the likelihood - i.e. if we log
    # this probability will it be infinite and cause the whole ll
    # to be infinite as a result
    # or has a small amount of numerical error crept in
    # and made a probability negative?
    # in which case just set it to zero
    # will add a warning in in this case
    # this function just means we can call log(prob)
    # without worrying about the errors this can cause
    if(length(prob)==0){
      print(c(parent.frame()$i,parent.frame()$t,parent.frame()$Time))
      View(parent.frame()$transitions)
    }
    if(prob<0){
      return(0)
    }else{
      return(ifelse(is.finite(log(prob)),log(prob),0)) 
    }
  }
  library(plyr, include.only = c("count"))
  # just need the count function from this package
  # to figure out how many unique chs there are
  uch <- count(ch) # unique number of capture histories
  # what are the states
  states <- c("N","N_","B1","LB","L_B","LB_","L_B_","S")
  # number of columns of the chs is Time
  Time <- ncol(ch)
  # how many unique capture histories are there
  Individuals <- nrow(uch)
  ll_i <- rep(0,Individuals) # want to create an individual likelihood for each capture history
  for(i in 1:Individuals){
    ch <- uch[i,1:Time] # the current capture history we are looking at
    transitions <- find.transitions(ch)
    # print(transitions)
    # need to deal with all the capture histories that
    # only have one state in them first
    if(nrow(transitions)==1){
      current_state <- transitions[1,1]
      current_time <- transitions[1,3]
      if(transitions[1,4]==(Time+1)){ # if this occurs at the final time...do nothing
        next
      }else{ # else we are in a Pr_r0 situation
        current_state_index <- which(states==current_state) # what number is the current state
        ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,phi,psi))
      }
    }else{ # therefore ch has multiple observations
      for(t in 1:nrow(transitions)){
        if(transitions[t,4]==(Time+1)){ # in the situation where the final time point is an observed state and not a zero
          next # then we just skip it
        }else if((transitions[t,4]-transitions[t,3])==1){ # aka are the states next to each other and not a state in the final time
          current_state_index <- which(states==transitions[t,1]) # what number is the current state
          current_time <- transitions[t,3]
          # need to deal with if the transition is the final one, a state to zero
          if(transitions[t,2]=="0"){
            ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,phi,psi))
          }else{
            next_state_index <- which(states==transitions[t,2]) # what number is the next state
            ll_i[i] <- ll_i[i] + Indicator(Pr_rs(current_state_index,next_state_index,current_time,phi,psi))
          } # end else
        } # all the states should now be next to each other so commented out the next bit
        # }else if((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]!="N"){ # we have a skipped a time and it was previously breeding
        #   current_state_index <- which(states==transitions[t,1]) # what number is the current state
        #   next_state_index <- which(states==transitions[t,2]) # what number is the next state
        #   current_time <- transitions[t,3]
        #   ll_i[i] <- ll_i[i] + Indicator(Pr_r0s(current_state_index,next_state_index,current_time,phi,psi))
        # }else if(((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]=="N") | 
        #          ((transitions[t,4]-transitions[t,3])>2)){
        #   # we have a skipped a time and it was NOT previously breeding OR
        #   # we have skipped multiple times (and so theoretically should not be breeding) OR
        #   current_state_index <- which(states==transitions[t,1]) # what number is the current state
        #   current_time <- transitions[t,3]
        #   ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,phi,psi))
        # }
        # print(c(t,ll_i))
      } # end for t in transitions
    } # end else ch has multiple observations
  } # end i
  ll <- sum((uch$freq)*(ll_i)) # multiply the log-likelihood for each unique capture history
  # by how many times that capture history occurs and sum it up
  # to make the multinomial loglikelihood
  return(ll)
}

ll.il <- function(theta,phi.ind,delt.ind,kap.ind,rho.ind,gam.ind,ups.ind,struc,ch){
  states <- c("N","N_","B1","LB","L_B","LB_","L_B_","S")
  
  phi <- untrans(logistic(theta[phi.ind]),struc$phi$time,struc$phi$state)
  delt <- untrans(logistic(theta[delt.ind]),struc$delt$time,struc$delt$state)
  kap <- untrans(logistic(theta[kap.ind]),struc$kap$time,struc$kap$state)
  rho <- untrans(logistic(theta[rho.ind]),struc$rho$time,struc$rho$state)
  gam <- untrans(logistic(theta[gam.ind]),struc$gam$time,struc$gam$state)
  ups <- untrans(logistic(theta[ups.ind]),struc$ups$time,struc$ups$state)
  
  psi <- make.psi(delt,kap,rho,gam,ups)
  # print(psi)
  
  ll <- il(ch,phi,psi)
  
  return(ll)
}

####################
# want to add in age dependance
####################
make.psi2 <- function(delta,kap,rho,gam,upsilon){
  states <- c("N","N_","B1","LB","L_B","LB_","L_B_","S")
  # all five input arrays should be indexed by [time,age,state]
  Time <- dim(delta)[1]
  Ages <- dim(delta)[2]
  if(ncol(kap)!=Time | ncol(rho)!=Time | ncol(gam)!=Time | ncol(upsilon)!=Time){
    stop("delta,kappa,rho, gamma and upsilon must all have the same number of columns, corresponding to time")
  }
  psi <- array(0,dim=c(length(states),length(states),Time,Ages))
  for(t in 1:Time){
    for(a in 1:Ages){
      psi[1,,t,a] <- c((1-rho[t,a,1])*upsilon[t,a,1],(1-rho[t,a,2])*(1-upsilon[t,a,2]),rho[t,a,3],0,0,0,0,0)
      psi[2,,t,a] <- c((1-rho[t,a,1])*upsilon[t,a,1],(1-rho[t,a,2])*(1-upsilon[t,a,2]),rho[t,a,3],0,0,0,0,0)
      row.names(psi) <- states
      colnames(psi) <- states
      # this is if these rows depend on the ROWS theyre in
      # for(i in 4:8){
      #   psi[i-1,,t] <- c(0,0,0,
      #                    (1-kap[i,t])*(1-delta[i,t])*gam[i,t],
      #                    (1-kap[i,t])*delta[i,t]*gam[i,t],
      #                    (1-kap[i,t])*(1-delta[i,t])*(1-gam[i,t]),
      #                    (1-kap[i,t])*delta[i,t]*(1-gam[i,t]),
      #                    kap[i,t])
      # }
      # psi[8,,t] <- c(0,0,0,
      #                (1-delta[8,t])*gam[8,t],
      #                delta[8,t]*gam[8,t],
      #                (1-delta[8,t])*(1-gam[8,t]),
      #                delta[8,t]*(1-gam[8,t]),
      #                0)
      # this is if these rows depend on the columns theyre in
      for(i in 4:8){
        psi[i-1,,t,a] <- c(0,0,0,
                         (1-kap[t,a,4])*(1-delta[t,a,4])*gam[t,a,4],
                         (1-kap[t,a,5])*delta[t,a,5]*gam[t,a,5],
                         (1-kap[t,a,6])*(1-delta[t,a,6])*(1-gam[t,a,6]),
                         (1-kap[t,a,7])*delta[t,a,7]*(1-gam[t,a,7]),
                         kap[t,a,8])
      }
      psi[8,,t,a] <- c(0,0,0,
                     (1-delta[t,a,4])*gam[t,a,4],
                     delta[t,a,5]*gam[t,a,5],
                     (1-delta[t,a,6])*(1-gam[t,a,6]),
                     delta[t,a,7]*(1-gam[t,a,7]),
                     0)
    }
  }
  return(psi)
}

find.transitions2 <- function(ch){ # need to add in age as component
  age <- as.numeric(ch[length(ch)])
  ch <- ch[1:(length(ch)-1)]
  Time <- length(ch)
  where <- which(ch!="0")
  df <- as.data.frame(array(NA,dim=c(max(length(where)-1,1),5))) # want it to be a character df but full of nothing
  colnames(df) <- c("r","s","t_r","t_s","age") # state r, state s, time at state r, time at states
  if(nrow(df)==1){ # if only one observation
    if(length(where)==1){ # if only one state observed
      df[1,] <- list(ch[ch!="0"],"0",where,where+1,age)
      return(df)
    }else if(length(where)==2){ # the case where there is a state at the final time
      df[1,] <- list(ch[where[1]],ch[where[2]],where[1],where[2],age)
      return(df)
    }
  }else{ # more than one observation
    for(t in 1:(length(where)-1)){ # go through all the states
      df[t,] <- list(ch[where[t]],ch[where[t+1]],where[t],where[t+1],age+t-1)
      # the age+t-1 ONLY WORKS if the capture history has NO gaps once the bird is seen
    } # end for t in 1:length(where)-1
    lw <- where[length(where)] # the last element of where
    if(lw!=Time){ # if the last place that a state is observed ISNT the last time
      df <- rbind(df,
                  list(ch[lw],"0",lw,lw+1,age+length(where)-1))
    }
    return(df)
  } # end else more than one observation
}

Pr_rs2 <- function(r,s,t,a,phi,psi){
  prob <- phi[t,a,r]*psi[r,s,t,a]
  return(prob)
}

Pr_r02 <- function(r,t,a,phi,psi){
  skip <- which(row.names(psi)=="S")
  prob <- phi[t,a,r]*psi[r,skip,t,a]*(1-phi[t,a,skip]) + (1-phi[t,a,r])
  return(prob)
}

il2 <- function(ch,phi,psi){ # il is "individual likelihood"
  # phi is indexed by [time,age,state]
  # psi is indexed by [state_r,state_s,time,age]
  Indicator <- function(prob){
    # a function that checks if the probability pu into it will
    # "cause problems" in the likelihood - i.e. if we log
    # this probability will it be infinite and cause the whole ll
    # to be infinite as a result
    # or has a small amount of numerical error crept in
    # and made a probability negative?
    # in which case just set it to zero
    # will add a warning in in this case
    # this function just means we can call log(prob)
    # without worrying about the errors this can cause
    if(length(prob)==0){
      print(c(parent.frame()$i,parent.frame()$t,parent.frame()$Time))
      print(c("s",parent.frame()$current_state_index))
      print(c("t",parent.frame()$current_time))
      print(c("a",parent.frame()$current_age))
      print(parent.frame()$ch)
      View(parent.frame()$transitions)
    }
    if(prob<0){
      return(0)
    }else{
      return(ifelse(is.finite(log(prob)),log(prob),0)) 
    }
  }
  library(plyr, include.only = c("count"))
  # just need the count function from this package
  # to figure out how many unique chs there are
  uch <- count(ch) # unique number of capture histories
  # what are the states
  states <- c("N","N_","B1","LB","L_B","LB_","L_B_","S")
  # number of columns of the chs minus is Time to account for the age column
  Time <- ncol(ch)-1
  # how many unique capture histories are there
  Individuals <- nrow(uch)
  ll_i <- rep(0,Individuals) # want to create an individual likelihood for each capture history
  for(i in 1:Individuals){
    ch <- uch[i,1:(Time+1)] # the current capture history we are looking at
    # can rename ch to be "the current ch" because we no longer need ch, just uch
    transitions <- find.transitions2(ch)
    # print(transitions)
    # need to deal with all the capture histories that
    # only have one state in them first
    if(nrow(transitions)==1){
      current_state <- transitions[1,1]
      current_time <- transitions[1,3]
      if(transitions[1,4]==(Time+1)){ # if this occurs at the final time...do nothing
        next
      }else{ # else we are in a Pr_r0 situation
        current_state_index <- which(states==current_state) # what number is the current state
        current_age <- transitions[1,5]
        ll_i[i] <- ll_i[i] + Indicator(Pr_r02(current_state_index,current_time,current_age,phi,psi))
      }
    }else{ # therefore ch has multiple observations
      for(t in 1:nrow(transitions)){
        if(transitions[t,4]==(Time+1)){ # in the situation where the final time point is an observed state and not a zero
          next # then we just skip it
        }else if((transitions[t,4]-transitions[t,3])==1){ # aka are the states next to each other and not a state in the final time
          current_state_index <- which(states==transitions[t,1]) # what number is the current state
          current_time <- transitions[t,3]
          current_age <- transitions[t,5]
          # need to deal with if the transition is the final one, a state to zero
          if(transitions[t,2]=="0"){
            ll_i[i] <- ll_i[i] + Indicator(Pr_r02(current_state_index,current_time,current_age,phi,psi))
          }else{
            next_state_index <- which(states==transitions[t,2]) # what number is the next state
            ll_i[i] <- ll_i[i] + Indicator(Pr_rs2(current_state_index,next_state_index,current_time,current_age,phi,psi))
          } # end else
        } # all the states should now be next to each other so commented out the next bit
        # }else if((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]!="N"){ # we have a skipped a time and it was previously breeding
        #   current_state_index <- which(states==transitions[t,1]) # what number is the current state
        #   next_state_index <- which(states==transitions[t,2]) # what number is the next state
        #   current_time <- transitions[t,3]
        #   ll_i[i] <- ll_i[i] + Indicator(Pr_r0s(current_state_index,next_state_index,current_time,phi,psi))
        # }else if(((transitions[t,4]-transitions[t,3])==2 & transitions[t,1]=="N") | 
        #          ((transitions[t,4]-transitions[t,3])>2)){
        #   # we have a skipped a time and it was NOT previously breeding OR
        #   # we have skipped multiple times (and so theoretically should not be breeding) OR
        #   current_state_index <- which(states==transitions[t,1]) # what number is the current state
        #   current_time <- transitions[t,3]
        #   ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,current_time,phi,psi))
        # }
        # print(c(t,ll_i))
      } # end for t in transitions
    } # end else ch has multiple observations
  } # end i
  ll <- sum((uch$freq)*(ll_i)) # multiply the log-likelihood for each unique capture history
  # by how many times that capture history occurs and sum it up
  # to make the multinomial loglikelihood
  return(ll)
}

ll.il2 <- function(theta,phi.ind,delt.ind,kap.ind,rho.ind,gam.ind,ups.ind,struc,ch){
  states <- c("N","N_","B1","LB","L_B","LB_","L_B_","S")
  
  phi <- untrans2(logistic(theta[phi.ind]),struc$phi$age,struc$phi$time,struc$phi$state)
  delt <- untrans2(logistic(theta[delt.ind]),struc$delt$age,struc$delt$time,struc$delt$state)
  kap <- untrans2(logistic(theta[kap.ind]),struc$kap$age,struc$kap$time,struc$kap$state)
  rho <- untrans2(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
  gam <- untrans2(logistic(theta[gam.ind]),struc$gam$age,struc$gam$time,struc$gam$state)
  ups <- untrans2(logistic(theta[ups.ind]),struc$ups$age,struc$ups$time,struc$ups$state)
  
  psi <- make.psi2(delt,kap,rho,gam,ups)
  # print(psi)
  
  ll <- il2(ch,phi,psi)
  
  return(ll)
}