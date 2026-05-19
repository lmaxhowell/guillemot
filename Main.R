file.names <- list.files(path = getwd(), pattern = "\\.R$") # list vector of file names
file.names <- file.names[-grep("^make",file.names)] # any non-function scripts should start with the word make
file.names <- file.names[-which(file.names %in% c("Main.R"))] # remove current file and Main (if different)
lapply(file.names, source) # source all functions needed

load("matrix_AUK.RData")
load("age_mat.RData")

states <- c("N", "B1","LB","L_B","LB_","L_B_","S")
states_code <- c(1,2,3,5,4,6) # no skipped code needed so this vector is one shorter length than states
ch <- matrix_AUK
for(i in 1:length(states_code)){
  ch[which(ch==states_code[i],arr.ind=TRUE)] <- states[i]
}


# parameters are phi,delta,kappa,rho,gamma

# Time <- 7
# Ages <- Time
# Sexs <- 2

# phi <- array(dim=c(Time,Ages,Sexs))
# delta <- array(dim=c(Time,Ages,Sexs))
# kap <- array(dim=c(Time,Ages,Sexs))
# rho <- array(dim=c(Time,Ages,Sexs))
# gam <- array(dim=c(Time,Ages,Sexs))

phi <- array(rep(c(0.3,0.7),c(1,6)),dim=c(length(states)))
delt <- array(0.3,dim=c(length(states)))
kap <- array(0.2,dim=c(length(states)))
rho <- array(0.5,dim=c(length(states)))
gam <- array(0.8,dim=c(length(states)))

c("N","N","N","B1","LB","LB","0","L_B_","0")



psi <- make.psi(delt,kap,rho,gam)

il(ch,phi,psi)

theta <- logit(c(0.3,0.7,0.3,0.2,0.5,0.8))
struc <- list("phi"=list(1,2:length(states)),
              "delt"=list(1:length(states)),
              "kap"=list(1:length(states)),
              "rho"=list(1:length(states)),
              "gam"=list(1:length(states)))
ll.il(theta,1:2,3,4,5,6,struc,ch)



# library(plyr, include.only = c("count"))
# uch <- count(matrix_AUK)
# dch <- cbind(uch,"problem"=rep(FALSE,nrow(uch)))
# multiple_zeros <- rep(0,15)
# for(i in 1:nrow(uch)){
#   detected_states <- which(uch[i,1:16]!=0)
#   if(length(detected_states)>1){
#     for(t in length(detected_states):2){
#       counter <- detected_states[t]-detected_states[t-1]-1
#       if(counter>1){
#         multiple_zeros[counter] <- multiple_zeros[counter] + 1
#         dch[i,18] <- TRUE
#       }
#     } 
#   }
# }
# sum(dch$problem)/nrow(dch)

