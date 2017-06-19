

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
set.seed(54321)
Rcpp::sourceCpp('CppFunctions/checkSZ.cpp')
###### 1: Import Data
House <- read.csv("Data/House.csv",header=T)
Indiv <- read.csv("Data/Indiv.csv",header=T)


###### 2: Remove Households with size < 2 and > 6
House <- House[which(House$NP >= 2 & House$NP <= 6),]


###### 3: Keep only Households with TEN == 1,2,or 3 and recode 1,2 as 1 and 3 as 2
House <- House[which(House$TEN == 1 | House$TEN == 2 | House$TEN == 3),]
House$TEN[which(House$TEN == 2)] <- 1
House$TEN[which(House$TEN == 3)] <- 2


###### 4: Take a sample of size 3000 Households
sample_size <- 3000
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
N <- nrow(X_indiv)
n <- nrow(X_house)
n_i <- as.numeric(as.character(X_house[,1]))
house_index <- rep(c(1:n),n_i)
p <- ncol(X_indiv)
q <- ncol(X_house)
level_indiv <- c(1:p)
level_indiv <- lapply(level_indiv, function(x) c(min(X_indiv[,x]):max(X_indiv[,x])) )
level_house <- c(1:q)
level_house <- lapply(level_house, function(x) c(min(X_house[,x]):max(X_house[,x])) )
Y_house <- X_house; Y_indiv <- X_indiv
struc_zero_variables_house <- which(is.element(colnames(X_house),c("HHGender","HHAge"))) ##gender is still included because I am still using 2012 data
struc_zero_variables_indiv <- which(is.element(colnames(X_indiv),c("Gender","Age","Relate"))) ##gender is still included because I am still using 2012 data
epsilon_indiv <- c(0.25,0.85,0.50)
epsilon_house <- c(0.35,0.60)
gamma <- 0.40
z_i <- rbinom(n,1,gamma)
Error_index_house <- which(z_i == 1)
E_house <- matrix(0,ncol=q,nrow=n)
E_indiv <- matrix(0,ncol=p,nrow=N)
for(i in Error_index_house){
  Error_index_indiv_i <- which(is.element(house_index,i)==TRUE)
  check_counter <- 1
  while(check_counter == 1){
    for(kq in struc_zero_variables_house){
      E_house[i,kq] <- rbinom(1,1,epsilon_house[which(struc_zero_variables_house==kq)])
      if(E_house[i,kq]==1){
        Y_house[i,kq] <- sample((level_house[[kq]])[-(1+X_house[i,kq]-min(level_house[[kq]]))],1)
      }
    }
    for(kp in struc_zero_variables_indiv){
      E_indiv[Error_index_indiv_i,kp] <- rbinom(length(Error_index_indiv_i),1,epsilon_indiv[which(struc_zero_variables_indiv==kp)])
      for(ii in 1:length(Error_index_indiv_i)){
        if(E_indiv[Error_index_indiv_i[ii],kp]==1){
          Y_indiv[Error_index_indiv_i[ii],kp] <- 
            sample((level_indiv[[kp]])[-(1+X_indiv[Error_index_indiv_i[ii],kp]-min(level_indiv[[kp]]))],1)
        }
      }
    }
    comb_to_check <- matrix(t(Y_indiv[Error_index_indiv_i,]),byrow=T,nrow=1)
    comb_to_check <- cbind(as.matrix(Y_house[i,(q-p+1):q]),comb_to_check) #add the household head before check
    check_counter <- checkSZ(comb_to_check,(length(Error_index_indiv_i)+1))
  }
}
E_house <- data.matrix(X_house)- data.matrix(Y_house)
E_house[E_house!=0] <- 1
E_indiv <- data.matrix(X_indiv)- data.matrix(Y_indiv)
E_indiv[E_indiv!=0] <- 1
#colSums(E_house)/length(Error_index_house)
#0.0000000 0.0000000 0.3826861 0.0000000 0.0000000 0.7265372 0.0000000 
#colSums(E_indiv)/length(which(is.element(house_index,Error_index_house)==TRUE))
#0.2830898 0.0000000 0.0000000 0.8960334 0.5991649 


###### 12: Add missing data
O_house <- matrix(1,ncol=q,nrow=n)
colnames(O_house) <- colnames(Y_house)
nonstruc_zero_variables_house <- c(1:ncol(Y_house))[-c(struc_zero_variables_house,which(colnames(X_house)=="HHSize"))]
O_house[,nonstruc_zero_variables_house] <- rbinom((n*length(nonstruc_zero_variables_house)),1,0.70)
O_indiv <- matrix(1,ncol=p,nrow=N)
colnames(O_indiv) <- colnames(Y_indiv)
nonstruc_zero_variables_indiv <- c(1:ncol(Y_indiv))[-struc_zero_variables_indiv]
O_indiv[,nonstruc_zero_variables_indiv] <- rbinom((N*length(nonstruc_zero_variables_indiv)),1,0.70)
Y_house[O_house==0] <- NA; Y_house$HHRelate <- 1;
Y_indiv[O_indiv==0] <- NA


###### 13: Save!!!
write.table(Y_house, file = "Data/Y_house.txt",row.names = FALSE)
write.table(Y_indiv, file = "Data/Y_indiv.txt",row.names = FALSE)
write.table(X_house, file = "Results/Data_house_truth.txt",row.names = FALSE)
write.table(X_indiv, file = "Results/Data_indiv_truth.txt",row.names = FALSE)
write.table(E_house, file = "Results/E_house_truth.txt",row.names = FALSE)
write.table(E_indiv, file = "Results/E_indiv_truth.txt",row.names = FALSE)
############################ End of Phase One #############################


###########################################################################
###########################################################################
#############################   BLANK SPACE   #############################
###########################################################################
###########################################################################
###########################################################################


######################## Phase Two: Model Fitting #########################
rm(list = ls())
###### 1: Load functions and packages
library(DirichletReg)
library(matrixStats)
library(coda)
Rcpp::sourceCpp('CppFunctions/prGpost.cpp')
Rcpp::sourceCpp('CppFunctions/prMpost.cpp')
Rcpp::sourceCpp('CppFunctions/checkSZ.cpp')
source("RFunctions/Functions.R")


###### 2: Load Data
Data_house = read.table("Data/Y_house.txt",header=TRUE)
Data_indiv = read.table("Data/Y_indiv.txt",header=TRUE)


###### 3: Set some parameters
HHSize_Index <- which(colnames(Data_house)=="HHSize")
HHRelate_Index <- which(colnames(Data_house)=="HHRelate")
struc_zero_house <- c("HHGender","HHAge")
struc_zero_indiv <- c("Gender","Age","Relate")


###### 4: Initialize chain
proc_total <- proc.time() 
source("RFunctions/InitializeChain.R")
(proc.time() - proc_total)[["elapsed"]]


###### 5: Run MCMC
proc_total <- proc.time() 
source("RFunctions/GibbsSampler.R")
total_time <- (proc.time() - proc_total)[["elapsed"]]
total_time


#colSums(E_house)/length(Error_index_house)
#0.0000000 0.0000000 0.3826861 0.0000000 0.0000000 0.7265372 0.0000000 
#colSums(E_indiv)/length(which(is.element(house_index,Error_index_house)==TRUE))
#0.2830898 0.0000000 0.0000000 0.8960334 0.5991649 


###### 5: Save Results
if(Weights$weight_option){
  MCMC_Results <- list(total_time_weighted=total_time,
                       dp_imput_indiv_weighted=PostSamples$dp_imput_indiv,
                       dp_imput_house_weighted=PostSamples$dp_imput_house,
                       ALPHA_weighted=PostSamples$ALPHA,BETA_weighted=PostSamples$BETA,
                       N_ZERO_weighted=PostSamples$N_ZERO,
                       M_CLUST_weighted=PostSamples$M_CLUST,G_CLUST_weighted=PostSamples$G_CLUST,
                       EPSILON_INDIV_weighted=PostSamples$EPSILON_INDIV,EPSILON_HOUSE_weighted=PostSamples$EPSILON_HOUSE)
} else {
  MCMC_Results <- list(total_time=total_time,dp_imput_indiv=PostSamples$dp_imput_indiv,
                       dp_imput_house=PostSamples$dp_imput_house,
                       ALPHA=PostSamples$ALPHA,BETA=PostSamples$BETA,N_ZERO=PostSamples$N_ZERO,
                       M_CLUST=PostSamples$M_CLUST,G_CLUST=PostSamples$G_CLUST,
                       EPSILON_INDIV=PostSamples$EPSILON_INDIV,EPSILON_HOUSE=PostSamples$EPSILON_HOUSE)
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


