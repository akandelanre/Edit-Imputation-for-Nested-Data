

###################################################################################
###################################################################################
###################### Error Location Model: Exploring Options ####################
###################################################################################
###################################################################################

###################################### START ######################################
###################### Phase One: Load Packages and Functions #####################
rm(list = ls())
Rcpp::sourceCpp('prZpost.cpp')
source("OtherFunctions.R")

##################### Phase Two: Prepare Data - Error Structure ###################
n <- 100
n_h <- round(c(1/2,1/3,1/6)*n); names(n_h) <- c("2","3","4")
N <- sum(as.numeric(as.character(names(n_h)))*(n_h))
q <- 2; p <- 3
epsilon_house_truth <- c(0.3,0.1)
E_house <- matrix(0,nrow=n,ncol=q)
for(k in 1:q){ E_house[,k] <- rbinom(n,1,epsilon_house_truth[k]) }
epsilon_indiv_truth <- matrix(c(0.4,0.05,0.2,0.001,0.5,0.45,0.15,0.35,0.01),ncol=3)
E_indiv <- matrix(0,nrow=N,ncol=p)
for(k in 1:p){
  E_indiv[(1:cumsum(n_h)[1]),k] <- rbinom(n_h[1],1,epsilon_indiv_truth[k,1]) 
  E_indiv[((cumsum(n_h)[1]+1):cumsum(n_h)[2]),k] <- rbinom(n_h[2],1,epsilon_indiv_truth[k,2]) 
  E_indiv[((cumsum(n_h)[2]+1):cumsum(n_h)[3]),k] <- rbinom(n_h[3],1,epsilon_indiv_truth[k,3]) 
  }

E_house <- data.frame(E_house); E_indiv <- data.frame(E_indiv)
for(i in 1:ncol(E_house)){ E_house[,i] = factor(E_house[,i],levels=c(0:1)) }
for(i in 1:ncol(E_indiv)){ E_indiv[,i] = factor(E_indiv[,i],levels=c(0:1)) }

############################### Phase Three: Fit Model #############################
TT <- 5
n_iter <- 10000
burn_in <- n_iter/2
MM <- 1
mc_thin <- 10
ModelFit <- Run_Model(E_house,E_indiv,TT,n_iter,burn_in,MM,mc_thin)

mean(ModelFit$Z_CLUST)
ModelFit$EPSILON_SUM[c(2,4,6),]/burn_in
epsilon_indiv_truth
ModelFit$EPSILON_SUM[c(8,10),]/burn_in
epsilon_house_truth

ModelFit$Imputations_house
E_house
ModelFit$Imputations_indiv
E_indiv

cbind(E_house,ModelFit$Imputations_house)
cbind(E_indiv,ModelFit$Imputations_indiv)





