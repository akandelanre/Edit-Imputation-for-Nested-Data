

###########################################################################
###########################################################################
################## Edit-Imputation Using the NDPMPM Model: ################
######################## Basic Rejection Sampler  #########################
###########################################################################
###########################################################################


###########################################################################
###########################################################################
################################## START ##################################
###########################################################################
###########################################################################


################### Phase One: One Time Data Preparation ##################
rm(list = ls())
###### 1: Import Data
House <- read.csv("Data/House.csv",header=T)
Indiv <- read.csv("Data/Indiv.csv",header=T)


###### 2: Remove Households with size < 2 and > 4
House <- House[which(House$NP >= 2 & House$NP <= 4),]


###### 3: Keep only Households with TEN == 1,2,or 3 and recode 1,2 as 1 and 3 as 2
House <- House[which(House$TEN == 1 | House$TEN == 2 | House$TEN == 3),]
House$TEN[which(House$TEN == 2)] <- 1
House$TEN[which(House$TEN == 3)] <- 2


###### 4: Take a sample of size 1,000 Households
set.seed(4321)
sample_size <- 1000
samp_index <- sort(sample(1:nrow(House),sample_size,replace=F))
House <- House[samp_index,]


###### 5: Pick the same households in the indiv data
pick_index <- is.element(Indiv$SERIALNO,House$SERIALNO)
Indiv <- Indiv[pick_index,]


###### 6: Recode within-household variables
###### 6a: First, the relationship variable
Indiv$RELP[which(Indiv$RELP == 11 | Indiv$RELP == 12 | Indiv$RELP == 13)] <- 12 #Boarder, roommate or partner
Indiv$RELP[which(Indiv$RELP == 14 | Indiv$RELP == 15)] <- 13 #Other non-relative or foster child
Indiv$RELP[which(Indiv$RELP == 10)] <- 11 #Other relative
Indiv$RELP[which(Indiv$RELP == 9)] <- 10 #Child-in-law
Indiv$RELP[which(Indiv$RELP == 8)] <- 9 #Parent-in-law
Indiv$RELP[which(Indiv$RELP == 7)] <- 8 #Grandchild
Indiv$RELP[which(Indiv$RELP == 6)] <- 7 #Parent
Indiv$RELP[which(Indiv$RELP == 5)] <- 6 #Sibling
Indiv$RELP[which(Indiv$RELP == 4)] <- 5 #Stepchild
Indiv$RELP[which(Indiv$RELP == 3)] <- 4 #Adopted child
Indiv$RELP[which(Indiv$RELP == 2)] <- 3 #Biological child
Indiv$RELP[which(Indiv$RELP == 1)] <- 2 #Spouse
Indiv$RELP[which(Indiv$RELP == 0)] <- 1 #Household head
###### 6b: Next, the race variable
Indiv$RAC3P[which(Indiv$RAC3P == 4 | Indiv$RAC3P == 8| Indiv$RAC3P == 9 | Indiv$RAC3P == 10)] <- 6
Indiv$RAC3P[which(Indiv$RAC3P == 5)] <- 4
Indiv$RAC3P[which(Indiv$RAC3P == 7)] <- 5
Indiv$RAC3P[which(Indiv$RAC3P >= 11 & Indiv$RAC3P <= 15)] <- 7
Indiv$RAC3P[which(Indiv$RAC3P >= 16 & Indiv$RAC3P <= 59)] <- 8
Indiv$RAC3P[which(Indiv$RAC3P >= 60 & Indiv$RAC3P <= 100)] <- 9
###### 6c: Next, the hisp variable
Indiv$HISP[which(Indiv$HISP >= 5 & Indiv$HISP <= 24)] <- 5
###### 6d: Lastly, age
Indiv$AGEP <- Indiv$AGEP + 1L


###### 7: Make household head into household level data
HHhead_data <- Indiv[which(Indiv$SPORDER==1),]
Indiv_minHH <- Indiv[-which(Indiv$SPORDER==1),]


###### 8: Combine Household and within-household data using the following ordering:
origdata <- data.frame(HHIndex = rep(c(1:sample_size),(House$NP-1L)),
                       WithinHHIndex = Indiv_minHH$SPORDER,
                       Gender = Indiv_minHH$SEX,Race = Indiv_minHH$RAC3P,Hisp = Indiv_minHH$HISP,
                       Age = Indiv_minHH$AGEP,Relate = Indiv_minHH$RELP,
                       Owner = rep(House$TEN,(House$NP-1L)),
                       HHGender = rep(HHhead_data$SEX,(House$NP-1L)),
                       HHRace = rep(HHhead_data$RAC3P,(House$NP-1L)),
                       HHHisp = rep(HHhead_data$HISP,(House$NP-1L)),
                       HHAge = rep(HHhead_data$AGEP,(House$NP-1L)),
                       HHRelate = rep(HHhead_data$RELP,(House$NP-1L)))


###### 9: Save and load back
#write.table(origdata,"Data/origdata.txt",row.names = FALSE)
#origdata <- read.table("Data/origdata.txt",header=T)


###### 10: Separate household and individual data
n_all <- length(unique(origdata$HHIndex))
X_indiv <- NULL
X_house <- NULL
for(i in 1:n_all){
  which_indiv <- which(origdata$HHIndex==i)
  X_indiv <- rbind(X_indiv,origdata[which_indiv,c("Gender","Race","Hisp","Age","Relate")])
  X_house <- rbind(X_house,cbind(length(which_indiv),origdata[
      which_indiv[length(which_indiv)],c("Owner","HHGender","HHRace","HHHisp","HHAge","HHRelate")]))
}
colnames(X_house) <- c("HHSize","Owner","HHGender","HHRace","HHHisp","HHAge","HHRelate")
X_house <- as.data.frame(X_house)


###### 11: Now create erroneous data using the measurement error model
set.seed(4321)
N <- nrow(X_indiv)
n <- nrow(X_house)
n_i <- as.numeric(as.character(X_house[,1]))
p <- ncol(X_indiv)
q <- ncol(X_house)
level_indiv = list(c(1:2),c(1:9),c(1:5),c(1:96),c(2:13))
level_house = list(c(1:3),c(1:2),c(1:2),c(1:9),c(1:5),c(16:96),c(1))
Y_house <- X_house; Y_indiv <- X_indiv
#a_epsilon_house <- c(0,1,1,3,3,2,0)
#b_epsilon_house <- c(1,10,15,20,15,10,1)
#a_epsilon_indiv <- c(1,1,3,3,2)
#b_epsilon_indiv <- c(5,10,15,10,20)
#epsilon_house <- rbeta(q,t(t(a_epsilon_house)),t(t(b_epsilon_house)))
#epsilon_indiv <- rbeta(p,t(t(a_epsilon_indiv)),t(t(b_epsilon_indiv)))
E_house <- E_indiv <- NULL
epsilon_indiv <- c(0.1,0.1,0.05,0.45,0.3)
epsilon_house <- c(0.0,0.2,0.1,0.05,0.35,0.5,0.0)
for(kq in 1:q){
  E_house <- cbind(E_house, rbinom(n,1,epsilon_house[kq]))
  for(i in 1:n){
    if(E_house[i,kq]==1){
      Y_house[i,kq] <- sample((level_house[[kq]])[-(1+X_house[i,kq]-min(level_house[[kq]]))],1)
    }
  }
}
for(kp in 1:p){
  E_indiv <- cbind(E_indiv, rbinom(N,1,epsilon_indiv[kp]))
  for(ii in 1:N){
    if(E_indiv[ii,kp]==1){
      Y_indiv[ii,kp] <- sample((level_indiv[[kp]])[-(1+X_indiv[ii,kp]-min(level_indiv[[kp]]))],1)
    }
  }
}
#head(E_house)
#head(X_house)
#head(Y_house)
#head(E_indiv)
#head(X_indiv)
#head(Y_indiv)


###### 12: Save!!!
write.table(Y_house, file = "Data/Y_house.txt",row.names = FALSE)
write.table(Y_indiv, file = "Data/Y_indiv.txt",row.names = FALSE)
write.table(X_house, file = "Results/Data_house_truth.txt",row.names = FALSE)
write.table(X_indiv, file = "Results/Data_indiv_truth.txt",row.names = FALSE)
write.table(epsilon_house, file = "Results/epsilon_house_truth.txt",row.names = FALSE)
write.table(epsilon_indiv, file = "Results/epsilon_indiv_truth.txt",row.names = FALSE)
#write.table(a_epsilon_house, file = "Results/a_epsilon_house_truth.txt",row.names = FALSE)
#write.table(a_epsilon_indiv, file = "Results/a_epsilon_indiv_truth.txt",row.names = FALSE)
#write.table(b_epsilon_house, file = "Results/b_epsilon_house_truth.txt",row.names = FALSE)
#write.table(b_epsilon_indiv, file = "Results/b_epsilon_indiv_truth.txt",row.names = FALSE)
############################ End of Phase One #############################


###########################################################################
###########################################################################
#############################   BLANK SPACE   #############################
###########################################################################
###########################################################################
###########################################################################


######################## Phase Two: Model Fitting #########################
rm(list = ls())
###### 1a: Load functions and packages
library(DirichletReg)
library(matrixStats)
library(coda)
Rcpp::sourceCpp('CppFunctions/prGpost.cpp')
Rcpp::sourceCpp('CppFunctions/prMpost.cpp')
Rcpp::sourceCpp('CppFunctions/checkSZ.cpp')
source("OtherFunctions/OtherFunctions.R")
###### 1b: Load prepared data; make sure data is in the right format
Y_house = read.table("Data/Y_house.txt",header=TRUE)
Y_indiv = read.table("Data/Y_indiv.txt",header=TRUE)
level_indiv = list(c(1:2),c(1:9),c(1:5),c(1:96),c(2:13))
level_house = list(c(1:3),c(1:2),c(1:2),c(1:9),c(1:5),c(16:96),c(1))
Y_house <- data.frame(Y_house)
Y_house_nf <- Y_house
for(i in 1:ncol(Y_house)){
  Y_house[,i] = factor(Y_house[,i],levels=level_house[[i]]) }
Y_indiv <- data.frame(Y_indiv)
Y_indiv_nf <- Y_indiv
for(i in 1:ncol(Y_indiv)){
  Y_indiv[,i] = factor(Y_indiv[,i],levels=level_indiv[[i]]) }


###### 2: Set global parameters and initialize chain
source("InitializeChain.R")

###### 3: Run MCMC
proc_total <- proc.time() 
source("GibbsSampler.R")
total_time <- (proc.time() - proc_total)[["elapsed"]]

#epsilon_indiv:
#c(0.1,0.1,0.05,0.45,0.3)
#epsilon_house:
#c(0.0,0.2,0.1,0.05,0.35,0.5,0.0)


###### 4: Save Results
if(hybrid_option){
  if(weight_option){
    MCMC_Results <- list(total_time_weighted_hybrid=total_time,
                         dp_imput_indiv_weighted_hybrid=dp_imput_indiv,
                         dp_imput_house_weighted_hybrid=dp_imput_house,
                         ALPHA_weighted_hybrid=ALPHA,BETA_weighted_hybrid=BETA,
                         N_ZERO_weighted_hybrid=N_ZERO,
                         M_CLUST_weighted_hybrid=M_CLUST,G_CLUST_weighted_hybrid=G_CLUST)
  } else {
    MCMC_Results <- list(total_time_hybrid=total_time,
                         dp_imput_indiv_hybrid=dp_imput_indiv,
                         dp_imput_house_hybrid=dp_imput_house,
                         ALPHA_hybrid=ALPHA,BETA_hybrid=BETA,N_ZERO_hybrid=N_ZERO,
                         M_CLUST_hybrid=M_CLUST,G_CLUST_hybrid=G_CLUST)
  }
} else {
  if(weight_option){
    MCMC_Results <- list(total_time_weighted=total_time,
                         dp_imput_indiv_weighted=dp_imput_indiv,
                         dp_imput_house_weighted=dp_imput_house,
                         ALPHA_weighted=ALPHA,BETA_weighted=BETA,N_ZERO_weighted=N_ZERO,
                         M_CLUST_weighted=M_CLUST,G_CLUST_weighted=G_CLUST)
  } else {
    MCMC_Results <- list(total_time=total_time,dp_imput_indiv=dp_imput_indiv,
                         dp_imput_house=dp_imput_house,
                         ALPHA=ALPHA,BETA=BETA,N_ZERO=N_ZERO,
                         M_CLUST=M_CLUST,G_CLUST=G_CLUST)
  }
}
writeFun <- function(LL){names.ll <- names(LL);for(i in names.ll){
    write.table(LL[[i]],paste0("Results/",i,".txt"),row.names = FALSE)}}
writeFun(MCMC_Results)
############################ End of Phase Two #############################

###########################################################################
###########################################################################
################################### END ###################################
###########################################################################
###########################################################################


