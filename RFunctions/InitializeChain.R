

###########################################################################
###########################################################################
################## Edit-Imputation Using the NDPMPM Model: ################
########################## Initializing Chain  ############################
###########################################################################
###########################################################################


###### 1: Set global parameters for data, set data to the right format and fill missing values with starting values
GlobalPara <- SetGlobalPara(Data_house,Data_indiv,HHSize_Index)
GlobalPara$level_house <- list(c(min(Data_house[,1]):max(Data_house[,1])),c(1:2),c(1:2),c(1:9),c(1:5),c(16:96),c(1))
GlobalPara$level_indiv <- list(c(1:2),c(1:9),c(1:5),c(1:96),c(2:13))
AllData <- FormatData(Data_house,Data_indiv,GlobalPara,HHSize_Index)
remove(Data_house)
remove(Data_indiv)

###### 2: Set parameters for structural zeros; sample impossibles in batches before checking constraints
StrucZerosPara <- SetStrucZerosPara(AllData$Y_house,AllData$Y_indiv,struc_zero_house,struc_zero_indiv,
                                    GlobalPara$H,GlobalPara$n,imp_batch=1000,miss_batch=10,prop_batch=1.2)


###### 3: Weighting; set weight_option to true for weighting/capping option
Weights <- SetWeights(weight_option=FALSE,GlobalPara$H,rep(1/2,length(GlobalPara$H)))


###### 4: Check for erronous households and set error model parameters
ErrorIndicators <- SetErrorIndicators(AllData,GlobalPara,StrucZerosPara)
#E_house_truth <- read.table("Results/E_house_truth.txt",header=TRUE)
#E_indiv_truth <- read.table("Results/E_indiv_truth.txt",header=TRUE)
#ErrorIndicators$NA_house[E_indiv_truth==1] <- 1
#ErrorIndicators$NA_indiv[E_indiv_truth==1] <- 1
#ErrorIndicators$NA_house[E_indiv_truth==0] <- 0
#ErrorIndicators$NA_indiv[E_indiv_truth==0] <- 0


###### 7: Initialize chain
ModelPara <- SetModelPara(AllData,FF=20,SS=15,GlobalPara$n,GlobalPara$n_i,StrucZerosPara)

###### 9: Set MCMC parameters
MCMCPara <- list(n_iter=10000,burn_in=5000,MM=50,mc_thin=10)
MCMCPara$M_to_use_mc <- round(seq((MCMCPara$burn_in +1),MCMCPara$n_iter,length.out=MCMCPara$MM))


###### 10: Create empty matrices to save results
PostSamples <- list(dp_imput_house=NULL,dp_imput_indiv=NULL,ALPHA=NULL,BETA=NULL,PII=NULL,G_CLUST=NULL,M_CLUST=NULL,
                    N_ZERO=NULL,EPSILON_INDIV=NULL,EPSILON_HOUSE=NULL,conv_check=NULL)
#PostSamples$LAMBDA <- matrix(0,ncol=(ncol(ModelPara$lambda)*nrow(ModelPara$lambda)),
#                             nrow=(MCMCPara$n_iter-MCMCPara$burn_in))
#PostSamples$OMEGA <- matrix(0,ncol=(ncol(ModelPara$omega)*nrow(ModelPara$omega)),
#                            nrow=(MCMCPara$n_iter-MCMCPara$burn_in))

