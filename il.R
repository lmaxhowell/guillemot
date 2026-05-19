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
  prob <- phi[r]*psi[r,skip]*(1-phi[skip]) + (1-phi[skip])
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
    # need to deal with all the capture histories that
    # only have one state in them first
    if(sum(ch!="0")==1){
      current_state <- ch[ch!="0"]
      if(which(ch!="0")==Time){ # if this occurs at the final time...do nothing
        next
      }else{ # else we are in a Pr_r0 situation
        current_state_index <- which(states==current_state) # what number is the current state
        ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,phi,psi))
      }
    }else{ # therefore ch has multiple observations
      for(t in 1:(Time-1)){# for each time, compare with the next time point 
        current_state <- unlist(ch[t]) # what state are we in right now?
        next_state <- unlist(ch[t+1]) # what is the next state?
        current_state_index <- which(states==current_state) # what number is the current state
        next_state_index <- which(states==next_state) # what number is the next state
        # print(c(current_state,next_state,current_state_index,next_state_index))
        if(current_state!="0"){ # only want to evaluate if we are in an observed state
          if(next_state=="0"){ # if next state is a zero then check
            # if we are in the final time as this is the only time
            # that Pr_r0 can happen, so we shall deal with that first
            if((t+1)==Time){
              ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,phi,psi))
            }else{ # otherwise its a Pr_r0s
              print(ch)
              plus <- 2
              next_state <- ch[t+plus] # what is the NEXT next state?
              while(next_state=="0"){
                plus <- plus + 1
                print(plus)
                next_state <- ch[t+plus] # what is the NEXT next state?
              }
              next_state_index <- which(states==next_state) # what number is the next state
              # next_state_index <- ifelse(next_state=="0",stop(print(ch)),which(states==next_state)) # what number is the next state
              ll_i[i] <- ll_i[i] + Indicator(Pr_r0s(current_state_index,next_state_index,phi,psi))
              next # might need a double next
            } # end else t+1==Time
          }else if(next_state!="0"){ # if it isnt a zero then we have an r->s transition
            ll_i[i] <- ll_i[i] + Indicator(Pr_rs(current_state_index,next_state_index,phi,psi))
          } # end else next state not equal zero
        } # end if current state not zero
      } # end t
    } # end i
  } # end else ch has multiple observations
  ll <- sum(uch$freq*ll_i) # multiply the log-likelihood for each unique capture history
  # by how many times that capture history occurs and sum it up
  # to make the multinomial loglikelihood
  return(ll)
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
    # need to deal with all the capture histories that
    # only have one state in them first
    if(sum(ch!="0")==1){
      current_state <- ch[ch!="0"]
      if(which(ch!="0")==Time){ # if this occurs at the final time...do nothing
        next
      }else{ # else we are in a Pr_r0 situation
        current_state_index <- which(states==current_state) # what number is the current state
        ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,phi,psi))
      }
    }else{ # therefore ch has multiple observations
      # want to trim the trailing zeros
      # find where the first and last state happen and what those states are
      # first_state_index <- head(which(ch!=0),1)
      last_state_index <- tail(which(ch!=0),1)
      # first_state <- ch[first_state_index]
      last_state <- ch[last_state_index]
      if(last_state_index!=(Time)){ # if there is a zero after the last observed state
        ch <- ch[1:(last_state_index+1)] # want only ONE trailing zero
      }
      for(t in 1:last_state_index){# for each time, compare with the next time point 
        current_state <- unlist(ch[t]) # what state are we in right now?
        next_state <- unlist(ch[t+1]) # what is the next state?
        current_state_index <- which(states==current_state) # what number is the current state
        next_state_index <- which(states==next_state) # what number is the next state
        # print(c(current_state,next_state,current_state_index,next_state_index))
        if(current_state!="0"){ # only want to evaluate if we are in an observed state
          if(next_state=="0"){ # if next state is a zero then check
            # if we are in the final time/final observed state as this is the only time
            # that Pr_r0 can happen, so we shall deal with that first
            if(t==last_state_index){
              ll_i[i] <- ll_i[i] + Indicator(Pr_r0(current_state_index,phi,psi))
            }else{ # otherwise its a Pr_r0s
              print(ch)
              plus <- 2
              next_state <- ch[t+plus] # what is the NEXT next state?
              while(next_state=="0"){
                plus <- plus + 1
                print(plus)
                next_state <- ch[t+plus] # what is the NEXT next state?
              }
              next_state_index <- which(states==next_state) # what number is the next state
              # next_state_index <- ifelse(next_state=="0",stop(print(ch)),which(states==next_state)) # what number is the next state
              ll_i[i] <- ll_i[i] + Indicator(Pr_r0s(current_state_index,next_state_index,phi,psi))
              next # might need a double next
            } # end else t+1==Time
          }else if(next_state!="0"){ # if it isnt a zero then we have an r->s transition
            ll_i[i] <- ll_i[i] + Indicator(Pr_rs(current_state_index,next_state_index,phi,psi))
          } # end else next state not equal zero
        } # end if current state not zero
      } # end t
    } # end i
  } # end else ch has multiple observations
  ll <- sum(uch$freq*ll_i) # multiply the log-likelihood for each unique capture history
  # by how many times that capture history occurs and sum it up
  # to make the multinomial loglikelihood
  return(ll)
}

il.ll <- function(theta,ch){
  
}

