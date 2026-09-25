load("gf.storm.RData")
source("make.il.alt.R")

Time <- ncol(ch3)-1
states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")
struc <- list("phi"=list("age"=list(1,2,3,4:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "eps"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))))

op <- optim( rep(0,6),
             ll.il,
             phi.ind=1,
             delt.ind=2,
             kap.ind=3,
             rho.ind=4,
             gam.ind=5,
             eps.ind=6,
             struc=struc,
             ch=ch3,
             control=list(maxit=10000,fnscale=-1))

op <- optim( rep(0,6),
             ll.il.alt,
             phi.ind=1,
             delt.ind=2,
             kap.ind=3,
             rho.ind=4,
             gam.ind=5,
             eps.ind=6,
             struc=struc,
             mv=mv3,
             control=list(maxit=10000,fnscale=-1))
c("phi","delta","kappa","rho","gamma","eps")
logistic(op$par)

tb3 <- table(unlist(ch3[,1:16]))
gam <- (tb3["LB"]+tb3["L_B"])/(tb3["LB"]+tb3["L_B"] + tb3["LB_"]+tb3["L_B_"])
del <- 1-((tb3["LB"]+tb3["LB_"])/(tb3["LB"]+tb3["L_B"] + tb3["LB_"]+tb3["L_B_"]))


timer(op2 <- optim( rep(0,9),
                    ll.il,
                    phi.ind=1:4,
                    delt.ind=5,
                    kap.ind=6,
                    rho.ind=7,
                    gam.ind=8,
                    eps.ind=9,
                    struc=struc,
                    ch=ch3,
                    control=list(maxit=10000,fnscale=-1)))
# Time difference of 22.31095 mins

mv3 <- suff.stat(ch3)
timer(op3 <- optim( rep(0,9),
                    ll.il.alt,
                    phi.ind=1:4,
                    delt.ind=5,
                    kap.ind=6,
                    rho.ind=7,
                    gam.ind=8,
                    eps.ind=9,
                    struc=struc,
                    mv=mv3,
                    control=list(maxit=10000,fnscale=-1)))
# Time difference of 4.289034 mins (vs Time difference of 19.05396 mins
# when mv is not calculated prior to the function and you just give it ch3
# i.e. it has to calculate mv within the likelihood every time)
logistic(op3$par)

timer(op4 <- optim( rep(0,12),
                    ll.il.alt2,
                    ageclasses=4,
                    timeclasses=2,
                    beta.struc=list(1:11,12:15),
                    mv=mv3,
                    control=list(maxit=10000,fnscale=-1)))
data.frame("par"=c(paste0("alpha",1:3),paste0("beta",11:12),paste0("rho",c(23,4,5)),"delta","kappa","gamma","epsilon"),
           "mle"=c(op4$par[1:5],logistic(op4$par[6:12])))

phi.trans.s <- function(alpha,beta,beta.struc){
  # inputs should be vectors
  alpha <- c(0,alpha) # to constrain alpha_0=0
  Ages <- length(alpha)
  Time_1 <- length(unlist(beta.struc)) # cause this is Time-1
  phi <- array(0,dim=c(Time_1,Time_1+1,length(states)))
  for(a in 1:(Time_1+1)){
    for(t in 1:Time_1){
      # print(c(t,a,min(a,Ages)))
      for(i in 1:length(beta.struc)){
        if(t %in% beta.struc[[i]]){
          break
        }
      }
      phi[t,a,] <- alpha[min(a,Ages)]+beta[i]
    }
  }
  return(logistic(phi)[11:12,1:4,1])
}

phi.trans.s(op4$par[1:3],op4$par[4:5],list(1:11,12:15))



# get all the ones with no breeding states (i.e. chi is involved)
chitime <- which(sapply(1:nrow(ch3),function(x) sum(ch3[x,] %in% states[3:8])==0 ))
chil <- c()
chi <- c()
for(i in chitime){
  ws <- which(ch3[i,] %in% c("N","E"))
  t <- ws[length(ws)]
  a <- length(ws)
  chil <- c(chil,log(Chi(1,t,a,phi,psi)))
  chi <- c(chi,Chi(1,t,a,phi,psi))
  
}

View(ch3[chitime[which(chi==1)],]) # all these are correct and should be 1


ch3[7,]
1-phi[14,1,1] + phi[14,1,1]*psi[1,2,14,1]*(1-phi[15,2,2] + phi[15,2,2]*psi[2,2,15,2])
chi[7]

# a row corresponding to the minimum chi value calculated
ch3[17,]
Chi(1,3,3,phi,psi)
1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(
  1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(
    1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(
      1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(
        1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(
          1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]*(
            1-phi[9,9,2] + phi[9,9,2]*psi[2,2,9,9]*(
              1-phi[10,10,2] + phi[10,10,2]*psi[2,2,10,10]*(
                1-phi[11,11,2] + phi[11,11,2]*psi[2,2,11,11]*(
                  1-phi[12,12,2] + phi[12,12,2]*psi[2,2,12,12]*(
                    1-phi[13,13,2] + phi[13,13,2]*psi[2,2,13,13]*(
                      1-phi[14,14,2] + phi[14,14,2]*psi[2,2,14,14]*(
                        1-phi[15,15,2] + phi[15,15,2]*psi[2,2,15,15]))))))))))))

1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]

1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4])

1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]))

1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6])))

1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]))))

1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8])))))

1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]*1-phi[9,9,2] + phi[9,9,2]*psi[2,2,9,9])))))

Chi(1,3,3,phi,psi)
c(0.75,0.6875,0.671875,0.6679688,0.6669922,0.666748,0.666687,0.6666718,0.6666679,0.666667,0.6666667,0.6666667,0.6666667)
c(1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3],
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5])),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]))),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7])))),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]))))),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]*(1-phi[9,9,2] + phi[9,9,2]*psi[2,2,9,9])))))),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]*(1-phi[9,9,2] + phi[9,9,2]*psi[2,2,9,9]*(1-phi[10,10,2] + phi[10,10,2]*psi[2,2,10,10]))))))),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]*(1-phi[9,9,2] + phi[9,9,2]*psi[2,2,9,9]*(1-phi[10,10,2] + phi[10,10,2]*psi[2,2,10,10]*(1-phi[11,11,2] + phi[11,11,2]*psi[2,2,11,11])))))))),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]*(1-phi[9,9,2] + phi[9,9,2]*psi[2,2,9,9]*(1-phi[10,10,2] + phi[10,10,2]*psi[2,2,10,10]*(1-phi[11,11,2] + phi[11,11,2]*psi[2,2,11,11]*(1-phi[12,12,2] + phi[12,12,2]*psi[2,2,12,12]))))))))),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]*(1-phi[9,9,2] + phi[9,9,2]*psi[2,2,9,9]*(1-phi[10,10,2] + phi[10,10,2]*psi[2,2,10,10]*(1-phi[11,11,2] + phi[11,11,2]*psi[2,2,11,11]*(1-phi[12,12,2] + phi[12,12,2]*psi[2,2,12,12]*(1-phi[13,13,2] + phi[13,13,2]*psi[2,2,13,13])))))))))),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]*(1-phi[9,9,2] + phi[9,9,2]*psi[2,2,9,9]*(1-phi[10,10,2] + phi[10,10,2]*psi[2,2,10,10]*(1-phi[11,11,2] + phi[11,11,2]*psi[2,2,11,11]*(1-phi[12,12,2] + phi[12,12,2]*psi[2,2,12,12]*(1-phi[13,13,2] + phi[13,13,2]*psi[2,2,13,13]*(1-phi[14,14,2] + phi[14,14,2]*psi[2,2,14,14]))))))))))),
  1-phi[3,3,1] + phi[3,3,1]*psi[1,2,3,3]*(1-phi[4,4,2] + phi[4,4,2]*psi[2,2,4,4]*(1-phi[5,5,2] + phi[5,5,2]*psi[2,2,5,5]*(1-phi[6,6,2] + phi[6,6,2]*psi[2,2,6,6]*(1-phi[7,7,2] + phi[7,7,2]*psi[2,2,7,7]*(1-phi[8,8,2] + phi[8,8,2]*psi[2,2,8,8]*(1-phi[9,9,2] + phi[9,9,2]*psi[2,2,9,9]*(1-phi[10,10,2] + phi[10,10,2]*psi[2,2,10,10]*(1-phi[11,11,2] + phi[11,11,2]*psi[2,2,11,11]*(1-phi[12,12,2] + phi[12,12,2]*psi[2,2,12,12]*(1-phi[13,13,2] + phi[13,13,2]*psi[2,2,13,13]*(1-phi[14,14,2] + phi[14,14,2]*psi[2,2,14,14]*(1-phi[15,15,2] + phi[15,15,2]*psi[2,2,15,15]))))))))))))
)


struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "eps"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))))
phi <- untrans(logistic(0),struc$phi$age,struc$phi$time,struc$phi$state)
dat <- dat.simulate2(phi,phi,phi,phi,phi,phi,ni4[1:15],897423184)
mv.dat <- suff.stat(dat)
timer(op.sim <- optim( rep(0,6),
                    ll.il.alt,
                    phi.ind=1,
                    delt.ind=2,
                    kap.ind=3,
                    rho.ind=4,
                    gam.ind=5,
                    eps.ind=6,
                    struc=struc,
                    mv=mv.dat,
                    hessian=TRUE,
                    control=list(maxit=10000,fnscale=-1)))
op.sim$convergence
logistic(op.sim$par)
ci <- function(mle,hess,q=0.975){
  plusminus <- qnorm(q)*sqrt(diag(solve(-hess)))
  return(cbind(mle-plusminus,mle,mle+plusminus))
}
logistic(ci(op.sim$par,op.sim$hessian))

tb <- table(unlist(dat[,1:16]))
gam <- (tb["LB"]+tb["L_B"])/(tb["LB"]+tb["L_B"] + tb["LB_"]+tb["L_B_"])

mv0 <- suff.stat(cbind(ch3[,1:16],rep(0,nrow(ch3))))
timer(op0 <- optim( rep(0,9),
                    ll.il.alt,
                    phi.ind=1:4,
                    delt.ind=5,
                    kap.ind=6,
                    rho.ind=7,
                    gam.ind=8,
                    eps.ind=9,
                    struc=struc,
                    mv=mv0,
                    control=list(maxit=10000,fnscale=-1)))
op0$convergence
logistic(op0$par)
