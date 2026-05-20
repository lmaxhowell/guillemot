dat.simulate <- function(phi,delt,kap,rho,gam,ni){
  psi <- make.psi(delt,kap,rho,gam)
  states <- c("N", "B1","LB","L_B","LB_","L_B_","S")
  
  Time <- length(ni)+1
  # an empty dataframe with a row for each individual
  # and no need to input Time as a parameter as can get this
  # from the length of ni. Plus one so we dont have to put
  # a zero at the end of the vector
  df <- as.data.frame(array(0,dim=c(sum(ni),Time)))
  
  for(tn in 1:length(ni)){ # for each time that we release a new cohort
    # for all the individuals marked at this time, put an "N" in their first time
    rows <- ifelse((tn-1)>0,(sum(ni[1:(tn-1)])+1):(sum(ni[1:(tn-1)])+ni[tn]),1:ni[1])
    df[rows,tn] <- states[1]
    for(i in 1:ni[tn]){ # for every individual we release in that cohort
      for(t in tn:(Time)){ # for every time point from now until the end of the study
        
      } # end for every t
    } # end for every individuals
  } # end for every cohort release rn
  return(df)
}