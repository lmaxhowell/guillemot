load("gf.storm.RData")
source("make.il.alt.R")
source("il.R")

Time <- ncol(ch3)-1
states <- c("N","E","B1","LB","L_B","LB_","L_B_","S")

struc <- list("phi"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "eps"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "p"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))))

mv3 <- suff.stat(ch3)

ll.il.alt(rep(0,7),1,2,3,4,5,6,7,struc,mv3)

op <- optim( rep(0,7),
             ll.il.alt,
             phi.ind=1,
             delt.ind=2,
             kap.ind=3,
             rho.ind=4,
             gam.ind=5,
             eps.ind=6,
             p.ind=7,
             struc=struc,
             mv=mv3,
             control=list(maxit=10000,fnscale=-1))

op$convergence
logistic(op$par)

struc2 <- list("phi"=list("age"=list(1,2,3,4:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "delt"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "kap"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "rho"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "gam"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "eps"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))),
              "p"=list("age"=list(1:Time),"time"=list(1:(Time-1)),"state"=list(1:length(states))))

op2 <- optim( rep(0,10),
             ll.il.alt,
             phi.ind=1:4,
             delt.ind=5,
             kap.ind=6,
             rho.ind=7,
             gam.ind=8,
             eps.ind=9,
             p.ind=10,
             struc=struc2,
             mv=mv3,
             control=list(maxit=10000,fnscale=-1))
op2$convergence
logistic(op2$par)

op3 <- optim( rep(0,13),
              ll.il.alt2,
              ageclasses=4,
              timeclasses=2,
              beta.struc=list(1:12,13:15),
              mv=mv3,
              control=list(maxit=10000,fnscale=-1))
op3$convergence
data.frame("par"=c(paste0("alpha",1:3),paste0("beta",12:13),paste0("rho",c(23,4,5)),"delta","kappa","gamma","epsilon","p"),
           "mle"=c(op3$par[1:5],logistic(op3$par[6:13])))
logistic(phi.trans(op3$par[1:3],op3$par[4:5],list(1:12,13:15))[12:13,1:4,1])
