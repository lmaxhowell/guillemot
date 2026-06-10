dat.simulate <- function(phi,delt,kap,rho,gam,ups,ni,seed=0){
  psi <- make.psi(delt,kap,rho,gam,ups)
  states <- c("N","N_","B1","LB","L_B","LB_","L_B_","S")
  
  Time <- length(ni)+1
  # ni is the number of individuals ringed at each occasion
  # an empty dataframe with a row for each individual
  # and no need to input Time as a parameter as can get this
  # from the length of ni. Plus one so we dont have to put
  # a zero at the end of the vector

  df <- as.data.frame(array(0,dim=c(sum(ni),Time)))
  
  if(seed!=0){
    set.seed(seed)
  }
  for(tn in 1:length(ni)){ # for each time that we release a new cohort
    # for all the individuals marked at this time, put an "N" in their first time
    # rows <- ifelse((tn-1)>0,(sum(ni[1:(tn-1)])+1):(sum(ni[1:(tn-1)])+ni[tn]),1:ni[1])
    if((tn-1)>0){
      rows <- (sum(ni[1:(tn-1)])+1):(sum(ni[1:(tn-1)])+ni[tn])
      addtoi <- sum(ni[1:(tn-1)]) # so that it gets the right line in the return df
    }else{
      rows <- 1:ni[1]
      addtoi <- 0 # so that it gets the right line in the return df
    }
    df[rows,tn] <- states[1]
    for(i in 1:ni[tn]){ # for every individual we release in that cohort
      current.state <- 1 # they start in state "N"
      for(t in tn:(Time-1)){ # for every time point from now until the end of the study
        survive <- rbinom(1,1,phi[current.state,t]) # do they survive
        if(survive==1){ # if they do survive
          theygo <- which(rmultinom(1,1,psi[current.state,,t])==1) # the state index they transition to
          df[i+addtoi,t+1] <- states[theygo]
          current.state <- theygo
        }else{ # they die and we skip to the next individual
          # df[i+addtoi,(t+1):Time]
          # previous.state <- df[i+addtoi,t-1]
          # need to remove any "skips" or "Nbars" in the previous time point
          # because its not observed again to infer this state
          # if(current.state %in% states[c(2,8)]){ # is the current state (aka the last state before it died)
          #   df[i+addtoi,t] <- "0" # replace the current "S" or "N_" with "0"
          # }
          # if(previous.state %in% states[c(2,8)]){
          #   df[i+addtoi,t-1] <- "0" # replace the current "S" or "N_" with "0"
          # }
          wdf <- tail(which(df[i+addtoi,] %in% states[-c(2,8)]),1) # what time did the last state occur that wasnt "S" or "N_"
          if(wdf<t){
            df[i+addtoi,(wdf+1):t] <- "0" # replace those "S" or "N_" with "0"
          }
          break # this is the skip to the next individual
        }
      } # end for every t
    } # end for every individuals
  } # end for every cohort release rn
  df[,Time][df[,Time]=="S"] <- "0" # remove any final skips that dont get written over by the bird dying
  df[,Time][df[,Time]=="N_"] <- "0" # remove any final skips that dont get written over by the bird dying
  return(df)
}

dat.simulate2 <- function(phi,delt,kap,rho,gam,eps,ni,seed=0){
  # ni should be a list of vectors of starting ages, the length of each vector
  # being the number of individuals released that time
  psi <- make.psi(delt,kap,rho,gam,eps)
  states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
  
  Time <- length(ni)+1
  # ni is the number of individuals ringed at each occasion
  # an empty dataframe with a row for each individual
  # and no need to input Time as a parameter as can get this
  # from the length of ni. Plus one so we dont have to put
  # a zero at the end of the vector
  
  # initialise return df. ncol is Time+1 so starting age can go at the end
  ni.s <- sapply(1:length(ni),function(x) length(ni[[x]]))
  df <- as.data.frame(array(0,dim=c(sum(ni.s),Time+1)))
  df[,Time+1] <- unlist(ni)
  
  
  if(seed!=0){
    set.seed(seed)
  }
  for(tn in 1:length(ni)){ # for each time that we release a new cohort
    # for all the individuals marked at this time, put an "N" in their first time
    # rows <- ifelse((tn-1)>0,(sum(ni[1:(tn-1)])+1):(sum(ni[1:(tn-1)])+ni[tn]),1:ni[1])
    if((tn-1)>0){
      rows <- (sum(ni.s[1:(tn-1)])+1):(sum(ni.s[1:(tn-1)])+ni.s[tn])
      addtoi <- sum(ni.s[1:(tn-1)]) # so that it gets the right line in the return df
    }else{
      rows <- 1:ni.s[1]
      addtoi <- 0 # so that it gets the right line in the return df
    }
    df[rows,tn] <- states[1]
    for(i in 1:ni.s[tn]){ # for every individual we release in that cohort
      current.state <- 1 # they start in state "N"
      current.age <- ni[[tn]][i] # what age do they start at
      for(t in tn:(Time-1)){ # for every time point from now until the end of the study
        # print(c(tn,t,i,current.age,current.state,phi[t,current.age,current.state]))
        survive <- rbinom(1,1,phi[t,current.age,current.state]) # do they survive
        if(survive==1){ # if they do survive
          theygo <- which(rmultinom(1,1,psi[current.state,,t,current.age])==1) # the state index they transition to
          df[i+addtoi,t+1] <- states[theygo]
          current.state <- theygo
          current.age <- current.age + 1
        }else{ # they die and we skip to the next individual
          # df[i+addtoi,(t+1):Time]
          # previous.state <- df[i+addtoi,t-1]
          # need to remove any "skips" or "Nbars" in the previous time point
          # because its not observed again to infer this state
          # if(current.state %in% states[c(2,8)]){ # is the current state (aka the last state before it died)
          #   df[i+addtoi,t] <- "0" # replace the current "S" or "N_" with "0"
          # }
          # if(previous.state %in% states[c(2,8)]){
          #   df[i+addtoi,t-1] <- "0" # replace the current "S" or "N_" with "0"
          # }
          break # this is the skip to the next individual
        }
      } # end for every t
      wdf <- tail(which(df[i+addtoi,1:Time] %in% states[-c(2,8)]),1) # what time did the last state occur that wasnt "S" or "N_"
      wdf2 <- tail(which(df[i+addtoi,1:Time] %in% states[c(2,8)]),1) # what time did the last state occur that was "S" or "N_"
      if(length(wdf2)>0){
        if(wdf2>wdf){
          df[i+addtoi,(wdf+1):Time] <- "0" # replace those "S" or "N_" with "0"
        } # end if wdf2>wdf
      } # end if wdf2>0
    } # end for every individuals
  } # end for every cohort release rn
  return(df)
}

dat.sim.wrap <- function(theta,phi.ind,delt.ind,kap.ind,rho.ind,gam.ind,eps.ind,struc,ni,seed=0){

  phi <- untrans(logistic(theta[phi.ind]),struc$phi$age,struc$phi$time,struc$phi$state)
  delt <- untrans(logistic(theta[delt.ind]),struc$delt$age,struc$delt$time,struc$delt$state)
  kap <- untrans(logistic(theta[kap.ind]),struc$kap$age,struc$kap$time,struc$kap$state)
  rho <- untrans(logistic(theta[rho.ind]),struc$rho$age,struc$rho$time,struc$rho$state)
  gam <- untrans(logistic(theta[gam.ind]),struc$gam$age,struc$gam$time,struc$gam$state)
  eps <- untrans(logistic(theta[eps.ind]),struc$eps$age,struc$eps$time,struc$eps$state)
  
  # psi <- make.psi(delt,kap,rho,gam,eps)
  # for(i in 1:16){
  #   print(sum(psi[,,1,1]==psi[,,i,1])==(prod(dim(psi[,,1,1]))))
  # }
  
  
  dat <- dat.simulate2(phi,delt,kap,rho,gam,eps,ni,seed)
  return(dat)
}