


###########################################################################
###########################################################################
################################# Model Assessment ########################
###########################################################################
###########################################################################

################################# START ###################################
rm(list = ls())
source("RFunctions/ModelAssessmentFunctions.R")

###### 1: Get population data
Population <- GetData(min_HH_size=2,max_HH_size=6)
GlobalPara <- list()
GlobalPara$weight_option <- TRUE


###### 2: Load saved results
TrueSample <- list();
TrueSample$House <- read.table("Results/Data_house_truth.txt",header=TRUE)
TrueSample$Indiv <- read.table("Results/Data_indiv_truth.txt",header=TRUE)
Model <- list();
Model$House <- read.table("Results/dp_imput_house.txt",header=TRUE)
Model$Indiv <- read.table("Results/dp_imput_indiv.txt",header=TRUE)
if(GlobalPara$weight_option){
  Model_Weighted <- list();
  Model_Weighted$House <- read.table("Results/dp_imput_house_weighted.txt",header=TRUE)
  Model_Weighted$Indiv <- read.table("Results/dp_imput_indiv_weighted.txt",header=TRUE)
}


###### 3: Define other parameters
GlobalPara$mm <- 50;
GlobalPara$N <- nrow(TrueSample$Indiv)
GlobalPara$n <- nrow(TrueSample$House)
GlobalPara$n_i <- as.numeric(as.character((TrueSample$House)$HHSize))
GlobalPara$p <- ncol(TrueSample$Indiv)
GlobalPara$q <- ncol(TrueSample$House)
GlobalPara$house_index <- rep(c(1:GlobalPara$n),GlobalPara$n_i)
GlobalPara$H <- sort(unique(GlobalPara$n_i))
GlobalPara$N_Pop <- nrow(Population$Indiv)
GlobalPara$n_Pop <- nrow(Population$House)
GlobalPara$n_i_Pop <- as.numeric(as.character(Population$House[,1]))
GlobalPara$house_index_Pop <- rep(c(1:GlobalPara$n_Pop),GlobalPara$n_i_Pop)
GlobalPara$level_house <- list(c(min((TrueSample$House)$HHSize):max((TrueSample$House)$HHSize)),
                               c(1:2),c(1:2),c(1:9),c(1:5),c(16:96),c(1))
GlobalPara$level_indiv <- list(c(1:2),c(1:9),c(1:5),c(1:96),c(2:13))
GlobalPara$Estimands <- 
  c(apply(as.matrix(GlobalPara$H),1,function(x) paste0("Same race household, $n_i = ",as.character(x+1),"$")),
    "Spouse present","Black HH, own",
    "Spouse present, HH is White","Spouse present, HH is Black",
    "White couple", "Non-White couple, own","Same race couple", 
    "White-nonwhite couple", "Only one parent",
    "At least one biological child present","One grandchild present",
    "At least three generations present",
    "Couples with age difference less than 5",
    "HH older than spouse, White HH","White HH with Hispanic origin",
    "Only one adult female in house with at least one child less than 5", 
    "Only one adult Hispanic male in house with at least one child less than 10",
    "Only one adult Black female in house with at least one child less than 18",
    "HH over 35, no children present","Male HH, own",
    "Black HH younger than 40, own","White HH younger than 25, own",
    "Hispanic HH older than 50, own",
    "2 generations present, Black HH",
    "3 generations present, White couple",
    "At least 2 generations present, Hispanic couple",
    "At least one stepchild",
    "At least one adopted child, White couple",
    "At least one biological child, Hispanic couple",
    "At least two biological children, Black couple")
GlobalPara$samp_size_Pop <- matrix(0,nrow=length(GlobalPara$Estimands))
GlobalPara$samp_size_Pop[1:length(GlobalPara$H)] <- table(GlobalPara$n_i_Pop)
GlobalPara$samp_size_Pop[GlobalPara$samp_size_Pop==0] <- GlobalPara$n_Pop
GlobalPara$samp_size <- matrix(0,nrow=length(GlobalPara$Estimands))
GlobalPara$samp_size[1:length(GlobalPara$H)] <- table(GlobalPara$n_i)
GlobalPara$samp_size[GlobalPara$samp_size==0] <- GlobalPara$n
                          

###### 4: Calculate estimands of interest for population
PopulationResults <- GetAllProbs(Population$House,Population$Indiv,GlobalPara$level_house,GlobalPara$level_indiv)
PopulationResults$OtherProb <- GetOtherProbs(Population$House,Population$Indiv,GlobalPara$n_Pop,GlobalPara$p,GlobalPara$q,
                                             GlobalPara$house_index_Pop,GlobalPara$H,length(GlobalPara$Estimands))
PopulationResults$OtherProb <- PopulationResults$OtherProb/GlobalPara$samp_size_Pop


###### 5: Calculate estimands of interest for true sample data
TrueSampleResults <- GetAllProbs(TrueSample$House,TrueSample$Indiv,GlobalPara$level_house,GlobalPara$level_indiv)
TrueSampleResults$OtherProb <- GetOtherProbs(TrueSample$House,TrueSample$Indiv,GlobalPara$n,GlobalPara$p,GlobalPara$q,
                                               GlobalPara$house_index,GlobalPara$H,length(GlobalPara$Estimands))
TrueSampleResults$OtherProb <- TrueSampleResults$OtherProb/GlobalPara$samp_size
TrueSampleResults$OtherVar <- (TrueSampleResults$OtherProb*(1-TrueSampleResults$OtherProb))/GlobalPara$samp_size
TrueSampleResults$OtherCINT <- CalculateCI(TrueSampleResults$OtherProb,TrueSampleResults$OtherVar,imp_ind=F)


###### 6: Calculate estimands of interest for model
###### 6a: Regular model
ModelResults <- GetOtherProbsMI(Model$House,Model$Indiv,GlobalPara,GlobalPara$samp_size)
###### 6b: Weighted model
if(GlobalPara$weight_option){
  ModelResults_Weighted <- GetOtherProbsMI(Model_Weighted$House,Model_Weighted$Indiv,GlobalPara,GlobalPara$samp_size)
}


###### 7: Combine results, make plots and save!
###### 7a: Regular model
png("Results/AllProbabilities_ModelResults.png",pointsize=11,width = 13, height = 7,bg="white",units="in",res=150)
par(mfrow=c(1,3),oma=c(2,0,0,0))
plot(TrueSampleResults$MarginalProb,ModelResults$MarginalProb,pch = 2,col="red", las=2, main ="Marginal",
     xlab = "Sample Estimate", ylab = "Average From 50 Imputed Datasets",xlim = c(0,1.0),ylim = c(0,1.0))
abline(a=0, b=1,lty = 1, lwd = 1);
plot(TrueSampleResults$BivariateProb,ModelResults$BivariateProb,pch = 2,col="darkblue", las=2, main ="Bivariate",
     xlab = "Sample Estimate", ylab = "Average From 50 Imputed Datasets",xlim = c(0,1.0),ylim = c(0,1.0))
abline(a=0, b=1,lty = 1, lwd = 1);
plot(TrueSampleResults$TrivariateProb,ModelResults$TrivariateProb,pch = 2,col="darkgreen", las=2, main ="Trivariate",
     xlab = "Sample Estimate", ylab = "Average From 50 Imputed Datasets",xlim = c(0,1.0),ylim = c(0,1.0))
abline(a=0, b=1,lty = 1, lwd = 1);
dev.off()
###### 7b: Weighted model
if(GlobalPara$weight_option){
  png("Results/AllProbabilities_ModelResults_Weighted.png",pointsize=11,width = 13, height = 7,bg="white",units="in",res=150)
  par(mfrow=c(1,3),oma=c(2,0,0,0))
  plot(TrueSampleResults$MarginalProb,ModelResults_Weighted$MarginalProb,pch = 2,col="red", las=2, main ="Marginal",
       xlab = "Sample Estimate", ylab = "Average From 50 Imputed Datasets",xlim = c(0,1.0),ylim = c(0,1.0))
  abline(a=0, b=1,lty = 1, lwd = 1);
  plot(TrueSampleResults$BivariateProb,ModelResults_Weighted$BivariateProb,pch = 2,col="darkblue", las=2, main ="Bivariate",
       xlab = "Sample Estimate", ylab = "Average From 50 Imputed Datasets",xlim = c(0,1.0),ylim = c(0,1.0))
  abline(a=0, b=1,lty = 1, lwd = 1);
  plot(TrueSampleResults$TrivariateProb,ModelResults_Weighted$TrivariateProb,pch = 2,col="darkgreen", las=2, main ="Trivariate",
       xlab = "Sample Estimate", ylab = "Average From 50 Imputed Datasets",xlim = c(0,1.0),ylim = c(0,1.0))
  abline(a=0, b=1,lty = 1, lwd = 1);
  dev.off()
}

if(GlobalPara$weight_option){
  CompareProbs <- cbind(PopulationResults$OtherProb,
                        TrueSampleResults$OtherProb,ModelResults$OtherProb, ModelResults_Weighted$OtherProb,
                        TrueSampleResults$OtherCINT,ModelResults$OtherCINT, ModelResults_Weighted$OtherCINT,
                        ModelResults$OtherLambda, ModelResults_Weighted$OtherLambda,
                        ModelResults$OtherRE, ModelResults_Weighted$OtherRE)
  
  colnames(CompareProbs) <- 
    c("Pop. Truth","Orig. Data Q","Model Q","Weighted Sampler Q",
      "Orig. Data L","Orig. Data U","Model L","Model U","Weighted Sampler L","Weighted Sampler U",
      "Model M.Frac","Weighted Sampler M.Frac","Model RE","Weighted Sampler RE")
} else {
  CompareProbs <- cbind(PopulationResults$OtherProb,TrueSampleResults$OtherProb,ModelResults$OtherProb,
                        TrueSampleResults$OtherCINT,ModelResults$OtherCINT,
                        ModelResults$OtherLambda,ModelResults$OtherRE)
  colnames(CompareProbs) <- c("Pop. Truth","Orig. Data Q","Model Q","Orig. Data L","Orig. Data U","Model L","Model U",
                              "Model M.Frac","Model RE")
}
rownames(CompareProbs) <- GlobalPara$Estimands

write.table(CompareProbs,"Results/CompareProbs.txt",row.names = TRUE)
CompareProbs <- read.table("Results/CompareProbs.txt",header=TRUE)
round(CompareProbs,3)
#round(CompareProbs[order(CompareProbs$Pop..Truth,decreasing = T),c(7:10)],3)

#library(xtable)
#xtable(round(CompareProbs[,c(1:3)],3),digits=3)
#xtable(round(CompareProbs[,c(1,5:10)],3),digits=3)

################################## END ####################################

###########################################################################
###########################################################################
#############################   BLANK SPACE   #############################
###########################################################################
###########################################################################
###########################################################################

