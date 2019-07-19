

###########################################################################
###########################################################################
################## Edit-Imputation Using the NDPMPM Model: ################
########################### The Gibbs Sampler  ############################
###########################################################################
###########################################################################


###########################################################################
###########################################################################
################################## START ##################################
###########################################################################
###########################################################################


for(mc in 1:MCMCPara$n_iter){
  proc_t <- proc.time()


  ## Print iteration number
  cat(paste("Iteration ", mc,"\n", sep = ""))


  ## Sample structural zeros data
  StrucZerosPara$n_batch_sum <- StrucZerosPara$n_batch_sum + ceiling(StrucZerosPara$n_0*StrucZerosPara$prop_batch)
  StrucZerosPara$n_batch <- ceiling(StrucZerosPara$n_batch_sum/mc) #no. of batches of imputations to sample
  StrucZerosPara$n_0[] <- 0
  AllStrucZerosData <- SampleImpossibles(AllData,GlobalPara,StrucZerosPara,ModelPara,Weights,HHSize_Index,
                                           HHRelate_Index,GenSyn=F)
  StrucZerosPara$n_0 <- AllStrucZerosData$n_0


  ## Sample X, the true response
  AllData$X_house <- AllData$Y_house
  AllData$X_indiv <- AllData$Y_indiv
  StrucZerosPara$n_batch_imp_sum <- StrucZerosPara$n_batch_imp_sum +
    ceiling(StrucZerosPara$n_0_reject*StrucZerosPara$prop_batch)
  StrucZerosPara$n_batch_imp <- ceiling(StrucZerosPara$n_batch_imp_sum/mc) + 1 #no. of batches of imputations to sample
  StrucZerosPara$n_0_reject[] <- 0
  AllData <- SampleTrueResponse(AllData,GlobalPara,StrucZerosPara,ModelPara,ErrorIndicators,HHSize_Index,HHRelate_Index)
  StrucZerosPara$n_0_reject <- AllData$n_0_reject


  ## Sample E, the error indicators
  AllData$E_house <- data.matrix(AllData$X_house)- data.matrix(AllData$Y_house)
  AllData$E_house[AllData$E_house!=0] <- 1
  AllData$E_indiv <- data.matrix(AllData$X_indiv)- data.matrix(AllData$Y_indiv)
  AllData$E_indiv[AllData$E_indiv!=0] <- 1


  ## Sample epsilon
  ModelPara$epsilon_house <- SampleEpsilonHouse(AllData,ErrorIndicators,ModelPara,StrucZerosPara)
  ModelPara$epsilon_indiv <- SampleEpsilonIndiv(AllData,ErrorIndicators,ModelPara,StrucZerosPara)


  ## Sample G
  ModelPara$G <- SampleG(AllData,ModelPara,GlobalPara)
  ModelPara$rep_G <- rep(ModelPara$G,GlobalPara$n_i)


  ## Sample M
  ModelPara$M <- SampleM(AllData,ModelPara)


  ## Sample lambda
  ModelPara$lambda <- SampleLambda(AllData,ModelPara,GlobalPara,AllStrucZerosData,Weights)


  ## Sample phi
  ModelPara$phi <- SamplePhi(AllData,ModelPara,GlobalPara,AllStrucZerosData,Weights)


  ## Sample pii
  ModelPara$pii <- SamplePi(ModelPara,AllStrucZerosData,Weights)


  ## Sample omega
  ModelPara$omega <- SampleOmega(ModelPara,AllStrucZerosData,Weights)


  ## Sample alpha
  ModelPara$alpha <- rgamma(1,shape=(ModelPara$a_alpha+ModelPara$FF-1),
                            rate=(ModelPara$b_alpha-log(ModelPara$pii[ModelPara$FF])))


  ## Sample beta
  ModelPara$beta <- rgamma(1,shape=(ModelPara$a_beta+(ModelPara$FF*(ModelPara$SS-1))),
                           rate=(ModelPara$b_beta-sum(log(ModelPara$omega[,ModelPara$SS]))))


  ## Check number of occupied clusters
  ModelPara$G_all <- c(ModelPara$G,AllStrucZerosData$G_0)
  ModelPara$S_occup <- CheckOccupied(ModelPara,AllStrucZerosData)


  ## Print some summaries
  cat(paste("Number of Occupied Household Classes is ", length(unique(ModelPara$G_all)), "\n", sep = ''))
  cat(paste("Max Number of Occupied Individual Classes is ", max(ModelPara$S_occup), "\n", sep = ''))
  cat(paste("Number of Sampled Augmented Households is ", sum(StrucZerosPara$n_0), "\n", sep = ''))


  ## Save posterior sample and imputations
  if(mc > MCMCPara$burn_in){
    if((mc %% MCMCPara$mc_thin)==0){
      #PII <- rbind(PII,c(pii))
      PostSamples$ALPHA <- rbind(PostSamples$ALPHA,ModelPara$alpha)
      PostSamples$G_CLUST <- rbind(PostSamples$G_CLUST,length(unique(ModelPara$G)))
      PostSamples$M_CLUST <- rbind(PostSamples$M_CLUST,max(ModelPara$S_occup))
      PostSamples$BETA <- rbind(PostSamples$BETA,ModelPara$beta)
      #PostSamples$LAMBDA[(mc-burn_in),] <- c(ModelPara$lambda)
      #PostSamples$OMEGA[(mc-burn_in),] <- c(ModelPara$omega)
      PostSamples$N_ZERO <- rbind(PostSamples$N_ZERO,sum(StrucZerosPara$n_0))
      PostSamples$EPSILON_INDIV <- rbind(PostSamples$EPSILON_INDIV, ModelPara$epsilon_indiv)
      PostSamples$EPSILON_HOUSE <- rbind(PostSamples$EPSILON_HOUSE, ModelPara$epsilon_house)
    }
    if(sum(mc==MCMCPara$M_to_use_mc)==1){
      PostSamples$dp_imput_indiv <- rbind(PostSamples$dp_imput_indiv,AllData$X_indiv)
      PostSamples$dp_imput_house <- rbind(PostSamples$dp_imput_house,AllData$X_house)
    }
  }


  ## Print more summaries
  cat(paste("Number of Sampled Rejections for Missing Data is ", sum(StrucZerosPara$n_0_reject), "\n", sep = ''))
  cat(paste("Total (True) Number of Sampled Augmented Households is ",
            (sum(StrucZerosPara$n_0_reject)+sum(StrucZerosPara$n_0/Weights$struc_weight)),"\n", sep = ''))
  cat(paste("Epsilon_indiv:", round(ModelPara$epsilon_indiv,3), "\n", sep = ''))
  cat(paste("Epsilon_house:", round(ModelPara$epsilon_house,3), "\n", sep = ''))
  elapsed_time <- (proc.time() - proc_t)[["elapsed"]]
  cat(paste("Elapsed Time = ", elapsed_time, "\n\n", sep = ' '))


  ## Randomly pick one summary to monitor and plot it
  PostSamples$conv_check <- rbind(PostSamples$conv_check,t(ModelPara$lambda%*%ModelPara$pii))
  if(nrow(PostSamples$conv_check) > 1){
    plot(mcmc(PostSamples$conv_check[,sample(ncol(PostSamples$conv_check),1,replace=F)]),col="blue")
    #plot(1:length(PostSamples$conv_check[,sample(ncol(PostSamples$conv_check),1,replace=F)]),
    #     PostSamples$conv_check[,sample(ncol(PostSamples$conv_check),1,replace=F)],ylab="",xlab="Interations",
    #     col=rainbow(length(PostSamples$conv_check[,sample(ncol(PostSamples$conv_check),1,replace=F)])),type="b")
  }

}



