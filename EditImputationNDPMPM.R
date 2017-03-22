

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
Rcpp::sourceCpp('CppFunctions/checkSZ.cpp')
###### 1: Import Data
House <- read.csv("Data/House.csv",header=T)
Indiv <- read.csv("Data/Indiv.csv",header=T)


###### 2: Remove Households with size < 2 and > 4
House <- House[which(House$NP >= 2 & House$NP <= 4),]


###### 3: Keep only Households with TEN == 1,2,or 3 and recode 1,2 as 1 and 3 as 2
House <- House[which(House$TEN == 1 | House$TEN == 2 | House$TEN == 3),]
House$TEN[which(House$TEN == 2)] <- 1
House$TEN[which(House$TEN == 3)] <- 2


###### 4: Take a sample of size 5000 Households
set.seed(4321)
sample_size <- 5000
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
house_index <- rep(c(1:n),n_i)
p <- ncol(X_indiv)
q <- ncol(X_house)
level_indiv = list(c(1:2),c(1:9),c(1:5),c(1:96),c(2:13))
level_house = list(c(1:3),c(1:2),c(1:2),c(1:9),c(1:5),c(16:96),c(1))
Y_house <- X_house; Y_indiv <- X_indiv
struc_zero_variables_house <- c(1,4) + 2 ##gender is still included because I am still using 2012 data
struc_zero_variables_indiv <- c(1,4,5) ##gender is still included because I am still using 2012 data
epsilon_indiv <- c(0.4,0.65,0.85)
epsilon_house <- c(0.3,0.75)
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
#0.0000000 0.0000000 0.3590361 0.0000000 0.0000000 0.7734940 0.0000000
#colSums(E_indiv)/length(which(is.element(house_index,Error_index_house)==TRUE))
#0.4235474 0.0000000 0.0000000 0.6574924 0.8853211

#colSums(E_house)/n
#0.000 0.000 0.149 0.000 0.000 0.321 0.000
#colSums(E_indiv)/N
#0.1736677 0.0000000 0.0000000 0.2695925 0.3630094

#gamma is 0.40
#sum(z_i) is 415

###### 12: Separate complete data
Data_indiv_cc <- Y_indiv[-z_i_index_indiv,]
Data_house_cc <- Y_house[-z_i_index_house,]


###### 13: Save!!!

###### 12: Save!!!
write.table(Y_house, file = "Data/Y_house.txt",row.names = FALSE)
write.table(Y_indiv, file = "Data/Y_indiv.txt",row.names = FALSE)
write.table(X_house, file = "Results/Data_house_truth.txt",row.names = FALSE)
write.table(X_indiv, file = "Results/Data_indiv_truth.txt",row.names = FALSE)
write.table(epsilon_house, file = "Results/epsilon_house_truth.txt",row.names = FALSE)
write.table(epsilon_indiv, file = "Results/epsilon_indiv_truth.txt",row.names = FALSE)
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

#colSums(E_indiv)/length(which(is.element(house_index,Error_index_house)==TRUE))
#0.4235474 0.6574924 0.8853211
#colSums(E_house)/length(Error_index_house)
#0.3590361 0.7734940

#colSums(E_indiv)/N
#0.1736677 0.2695925 0.3630094
#colSums(E_house)/n
#0.149 0.321 

#epsilon_indiv <- c(0.4,0.65,0.85)
#epsilon_house <- c(0.3,0.75)


###### 4: Save Results
if(weight_option){
  MCMC_Results <- list(total_time_weighted=total_time,
                       dp_imput_indiv_weighted=dp_imput_indiv,
                       dp_imput_house_weighted=dp_imput_house,
                       ALPHA_weighted=ALPHA,BETA_weighted=BETA,N_ZERO_weighted=N_ZERO,
                       M_CLUST_weighted=M_CLUST,G_CLUST_weighted=G_CLUST,
                       EPSILON_INDIV_weighted=EPSILON_INDIV,EPSILON_HOUSE_weighted=EPSILON_HOUSE)
} else {
  MCMC_Results <- list(total_time=total_time,dp_imput_indiv=dp_imput_indiv,
                       dp_imput_house=dp_imput_house,
                       ALPHA=ALPHA,BETA=BETA,N_ZERO=N_ZERO,
                       M_CLUST=M_CLUST,G_CLUST=G_CLUST,
                       EPSILON_INDIV=EPSILON_INDIV,EPSILON_HOUSE=EPSILON_HOUSE)
}
writeFun <- function(LL){names.ll <- names(LL);for(i in names.ll){
    write.table(LL[[i]],paste0("Results/",i,".txt"),row.names = FALSE)}}
writeFun(MCMC_Results)

write.table(Data_house_cc, file = "Results/Data_house_cc.txt",row.names = FALSE)
write.table(Data_indiv_cc, file = "Results/Data_indiv_cc.txt",row.names = FALSE)
############################ End of Phase Two #############################

###########################################################################
###########################################################################
################################### END ###################################
###########################################################################
###########################################################################


