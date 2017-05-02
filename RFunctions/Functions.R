
###########################################################################
###########################################################################
################## Edit-Imputation Using the NDPMPM Model: ################
############################ All R Functions  #############################
###########################################################################
###########################################################################


## Function to sort a HH vector of within household data by relate variable
QuickSortHH <- function(X=x,P=p,R=r){# X =vector, p =no of variables, r =col index of relate var 
  XX = matrix(X,ncol=P,byrow = TRUE)
  XX = XX[order(XX[,R]),]
  return(matrix(t(XX),nrow=1))
}


## Function to set global parameters for data
SetGlobalPara <- function(Data_house,Data_indiv,HHSize_Index){
  N <- nrow(Data_indiv)
  n <- nrow(Data_house)
  n_i <- as.numeric(as.character(Data_house[,HHSize_Index]))
  p <- ncol(Data_indiv)
  q <- ncol(Data_house)
  house_index <- rep(c(1:n),n_i)
  n_i_index <- rep(n_i,n_i)
  H <- sort(unique(n_i))
  #level_indiv <- c(1:p)
  #level_indiv <- lapply(level_indiv, function(x) c(min(Data_indiv[,x],na.rm=T):max(Data_indiv[,x],na.rm=T)) )
  #level_house <- c(1:q)
  #level_house <- lapply(level_house, function(x) c(min(Data_house[,x],na.rm=T):max(Data_house[,x],na.rm=T)) )
  
  return(list(N=N,n=n,n_i=n_i,p=p,q=q,house_index=house_index,n_i_index=n_i_index,H=H))
              #level_indiv=level_indiv,level_house=level_house
}


## Function to make data into data frame, keep a nonfactor copy and fill missing values with starting values
FormatData <- function(Data_house,Data_indiv,GlobalPara,HHSize_Index){
Data_house <- data.frame(Data_house)
NA_miss_house <- Data_house
for(i in 1:GlobalPara$q){
  Data_house[,i] = factor(Data_house[,i],levels=GlobalPara$level_house[[i]]) }
Data_indiv <- data.frame(Data_indiv)
NA_miss_indiv <- Data_indiv
for(i in 1:GlobalPara$p){
  Data_indiv[,i] = factor(Data_indiv[,i],levels=GlobalPara$level_indiv[[i]]) }
#if(sum(is.na(NA_miss_indiv)) > 0){
#  for (ii in 1:GlobalPara$p){
#    Data_indiv[is.na(Data_indiv[,ii]),ii] <- 
#      sample(GlobalPara$level_indiv[[ii]],length(Data_indiv[is.na(Data_indiv[,ii]),ii]),
#             replace=T,prob=summary(na.omit(Data_indiv[,ii])))
#  }
#}
#if(sum(is.na(NA_miss_house)) > 0){
#  for (jj in c(1:GlobalPara$q)[-HHSize_Index]){
#    Data_house[is.na(Data_house[,jj]),jj] <- 
#      sample(GlobalPara$level_house[[jj]],length(Data_house[is.na(Data_house[,jj]),jj]),
#             replace=T,prob=summary(na.omit(Data_house[,jj])))
#  }
#}

d_k_house <- d_k_indiv <- ini_marg_house <- ini_marg_indiv <- NULL
for(k in 1:GlobalPara$q){
  d_k_house <- cbind(d_k_house,nlevels(Data_house[,k]))
  ini_marg_k <- as.data.frame(table(Data_house[,k]))$Freq/sum(table(Data_house[,k]))
  ini_marg_k <- matrix(ini_marg_k,ncol=1)
  ini_marg_house <- rbind(ini_marg_house,ini_marg_k)  }
for(k in 1:GlobalPara$p){
  d_k_indiv <- cbind(d_k_indiv,nlevels(Data_indiv[,k]))
  ini_marg_k <- as.data.frame(table(Data_indiv[,k]))$Freq/sum(table(Data_indiv[,k]))
  ini_marg_k <- matrix(ini_marg_k,ncol=1)
  ini_marg_indiv <- rbind(ini_marg_indiv,ini_marg_k)  }
d_k_indiv_cum <- 1+cumsum(c(0,d_k_indiv[,-GlobalPara$p]))
d_k_house_cum <- 1+cumsum(c(0,d_k_house[,-GlobalPara$q]))

FFF_house <- matrix(rep(cumsum(c(0,d_k_house[,-GlobalPara$q])),each=GlobalPara$n),ncol=GlobalPara$q)
FFF_indiv <- matrix(rep(cumsum(c(0,d_k_indiv[,-GlobalPara$p])),each=GlobalPara$N),ncol=GlobalPara$p)

return(list(Y_house=Data_house,Y_indiv=Data_indiv,NA_miss_house=NA_miss_house,NA_miss_indiv=NA_miss_indiv,
            ini_marg_house=ini_marg_house,ini_marg_indiv=ini_marg_indiv,
            d_k_house_cum=d_k_house_cum,d_k_indiv_cum=d_k_indiv_cum,d_k_house=d_k_house,d_k_indiv=d_k_indiv,
            FFF_house=FFF_house,FFF_indiv=FFF_indiv))
}


## Function to set structural zeros parameters
SetStrucZerosPara <- function(Y_house,Y_indiv,struc_zero_house,struc_zero_indiv,H,n,imp_batch,miss_batch,prop_batch){
  struc_zero_variables_house <- which(is.element(colnames(Y_house),struc_zero_house)) 
  struc_zero_variables_indiv <- which(is.element(colnames(Y_indiv),struc_zero_indiv))
  nonstruc_zero_variables_indiv <- c(1:ncol(Y_indiv))[-struc_zero_variables_indiv]
  nonstruc_zero_variables_house <- c(1:ncol(Y_house))[-struc_zero_variables_house]
  n_batch_sum <- rep(imp_batch,length(H)) 
  n_0 <- rep(1,length(H))
  n_batch_imp_sum <- rep(miss_batch,n)
  n_0_reject <- rep(1,n)
  
  return(list(struc_zero_variables_house=struc_zero_variables_house,
              struc_zero_variables_indiv=struc_zero_variables_indiv,
              nonstruc_zero_variables_house=nonstruc_zero_variables_house,
              nonstruc_zero_variables_indiv=nonstruc_zero_variables_indiv,n_batch_sum=n_batch_sum,n_0=n_0,
              n_batch_imp_sum=n_batch_imp_sum,n_0_reject=n_0_reject,prop_batch=prop_batch) )
}


## Function to make weights for capping
SetWeights <- function(weight_option,H,struc_weight){
  if(weight_option){
    struc_weight <- struc_weight
  } else {
    struc_weight <- rep(1,length(H)) #set weights: must be ordered & no household size must be excluded
  }
  struc_weight <- as.matrix(struc_weight)
  rownames(struc_weight) <- as.character(H)
  
  return(list(struc_weight=struc_weight))
}


## Function to check errorneous households
SetErrorIndicators <- function(AllData,GlobalPara,StrucZerosPara){
  z_i <- matrix(0,ncol=1,nrow=GlobalPara$n)
  NA_error_house <- matrix(0,ncol=GlobalPara$q,nrow=GlobalPara$n,byrow=T)
  NA_error_indiv <- matrix(0,ncol=GlobalPara$p,nrow=GlobalPara$N,byrow=T)
  for(hh_size in GlobalPara$H){
    house_index_hh <- which(GlobalPara$n_i == hh_size)
    comb_to_check <- AllData$Y_indiv[which(is.element(GlobalPara$house_index,house_index_hh)==TRUE),]
    comb_to_check <- apply(comb_to_check,2,function(x) as.numeric(as.character(x)))
    comb_to_check <- matrix(t(comb_to_check),byrow=T,ncol=(GlobalPara$p*hh_size))
    #add the household head before check
    comb_to_check <- cbind(AllData$Y_house[house_index_hh,(GlobalPara$q-GlobalPara$p+1):GlobalPara$q],comb_to_check)
    comb_to_check <- apply(comb_to_check,2,function(x) as.numeric(as.character(x)))
    z_i[house_index_hh] <- ifelse(checkSZ(comb_to_check,(hh_size+1))==1,0,1)
    NA_error_house[house_index_hh,StrucZerosPara$struc_zero_variables_house] <- 
      rep(z_i[house_index_hh],length(StrucZerosPara$struc_zero_variables_house))
    NA_error_indiv[which(is.element(GlobalPara$house_index,house_index_hh)==TRUE),StrucZerosPara$struc_zero_variables_indiv] <- 
      rep(rep(z_i[house_index_hh],GlobalPara$n_i[house_index_hh]),length(StrucZerosPara$struc_zero_variables_indiv))
  }
  z_i_index_house <- which(z_i == 1)
  z_i_index_indiv <- which(is.element(GlobalPara$house_index,z_i_index_house)==TRUE)
  NA_error_house[NA_error_house==1] <- NA
  NA_error_indiv[NA_error_indiv==1] <- NA
  NA_error_house[is.na(AllData$NA_miss_house)] <- NA #For missing data
  NA_error_indiv[is.na(AllData$NA_miss_indiv)] <- NA #For missing data
  
  return(list(z_i_index_house=z_i_index_house,z_i_index_indiv=z_i_index_indiv,
              NA_error_house=NA_error_house,NA_error_indiv=NA_error_indiv))
}


## Function to set model parameters for data
SetModelPara <- function(AllData,FF,SS,n,n_i,StrucZerosPara){
  alpha <- beta <- 1
  a_kdk <- 1
  a_alpha <- b_alpha <- a_beta <- b_beta <- 0.25
  lambda <- matrix(rep(AllData$ini_marg_house,FF),ncol=FF)
  phi <- matrix(0,nrow=length(AllData$ini_marg_indiv),ncol=FF*SS) #make phi matrix and not array for c++
  for(gm in 1:(FF*SS)){
    phi[,gm] <- AllData$ini_marg_indiv   }
  U <- matrix(rbeta(FF,1,alpha),nrow=FF)
  V <- matrix(rbeta((FF*SS),1,beta),nrow=FF,ncol=SS)
  U[FF]<-1
  V[,SS]<-1
  one_min_U <- 1L-U
  one_min_U <- c(1,cumprod(one_min_U[1:(FF-1)]))
  one_min_V <- 1L-V
  one_min_V <- cbind(1,t(apply(one_min_V[,-SS],1,cumprod)))
  pii <- U*one_min_U
  omega <- V*one_min_V
  a_epsilon_house <- rep(1,length(StrucZerosPara$struc_zero_variables_house))
  b_epsilon_house <- rep(1,length(StrucZerosPara$struc_zero_variables_house))
  a_epsilon_indiv <- rep(1,length(StrucZerosPara$struc_zero_variables_indiv))
  b_epsilon_indiv <- rep(1,length(StrucZerosPara$struc_zero_variables_indiv))
  epsilon_house <- rbeta(length(StrucZerosPara$struc_zero_variables_house),t(t(a_epsilon_house)),t(t(b_epsilon_house)))
  epsilon_indiv <- rbeta(length(StrucZerosPara$struc_zero_variables_indiv),t(t(a_epsilon_indiv)),t(t(b_epsilon_indiv)))
  pr_G <- matrix(pii,byrow=T,ncol=FF,nrow=n)
  Ran_unif_G <- runif(nrow(pr_G))
  cumul_G <- pr_G%*%upper.tri(diag(ncol(pr_G)),diag=TRUE)
  G <- rowSums(Ran_unif_G>cumul_G) + 1L
  rep_G <- rep(G,n_i)
  pr_M <- omega[rep_G,]
  Ran_unif_M <- runif(nrow(pr_M))
  cumul_M <- pr_M%*%upper.tri(diag(ncol(pr_M)),diag=TRUE)
  M <- rowSums(Ran_unif_M>cumul_M) + 1L
  
  return(list(FF=FF,SS=SS,alpha=alpha,beta=beta,a_kdk=a_kdk,a_alpha=a_alpha,b_alpha=b_alpha,a_beta=a_beta,b_beta=b_beta,
              lambda=lambda,phi=phi,pii=pii,omega=omega,epsilon_house=epsilon_house,epsilon_indiv=epsilon_indiv,
              a_epsilon_house=a_epsilon_house,b_epsilon_house=b_epsilon_house,
              a_epsilon_indiv=a_epsilon_indiv,b_epsilon_indiv=b_epsilon_indiv,G=G,rep_G=rep_G,M=M) )
}


## Function to sample impossible household data for a given household size
SampleImpossibles_h <- 
  function(AllData,GlobalPara,StrucZerosPara,ModelPara,Weights,HHSize_Index,HHRelate_Index,hh_size,GenSyn){
  H <- GlobalPara$H; n_i <- GlobalPara$n_i; struc_weight <- Weights$struc_weight;
  lambda <- ModelPara$lambda; level_house <- GlobalPara$level_house; level_indiv <- GlobalPara$level_indiv;
  pii <- ModelPara$pii; phi <- ModelPara$phi; FF <- ModelPara$FF; omega <- ModelPara$omega;
  n_batch <- StrucZerosPara$n_batch; p <- GlobalPara$p; q <- GlobalPara$q;
  d_k_house_cum <- AllData$d_k_house_cum; d_k_indiv_cum <- AllData$d_k_indiv_cum;
  d_k_house <- AllData$d_k_house; d_k_indiv <- AllData$d_k_indiv; n_0 <- StrucZerosPara$n_0
  
  G_0 <- M_0 <- X_house_struc <- X_indiv_struc <- X_house_valid <- X_indiv_valid <- NULL
  ss <- which(H==hh_size); t_0 <- 0; t_1 <- 0
  n_possibles <- ceiling(length(which(n_i == hh_size))*struc_weight[as.character(hh_size),])
  while(t_1 < n_possibles){
    quick_d_k_house_index <- d_k_house_cum[HHSize_Index]:cumsum(d_k_house)[HHSize_Index]
    pr_G_t <- lambda[(quick_d_k_house_index[which(level_house[[HHSize_Index]]==hh_size)]),]*pii
    G_t <- sample(FF,n_batch[ss],prob=pr_G_t,replace=T)
    rep_G_t <- rep(G_t,each=hh_size)
    pr_M_post_t <- omega[rep_G_t,]
    Ran_unif_M_t <- runif(nrow(pr_M_post_t))
    cumul_M_t <- pr_M_post_t%*%upper.tri(diag(ncol(pr_M_post_t)),diag=TRUE)
    M_t <- rowSums(Ran_unif_M_t>cumul_M_t) + 1L
    X_house_t <- matrix(0,nrow=length(G_t),ncol=q)
    X_house_t[,HHSize_Index] <- hh_size
    lambda_g_t <- t(lambda[,G_t])
    for(tt in c(1:q)[-c(HHSize_Index,HHRelate_Index)]){
      pr_house_t <- lambda_g_t[,d_k_house_cum[tt]:cumsum(d_k_house)[tt]]
      Ran_unif_t <- runif(nrow(pr_house_t))
      cumul_t <- pr_house_t%*%upper.tri(diag(ncol(pr_house_t)),diag=TRUE)
      level_house_t <- level_house[[tt]]
      X_house_t[,tt] <- level_house_t[rowSums(Ran_unif_t > cumul_t) + 1L] 
    }
    X_house_t[,HHRelate_Index] <- 1 #relate is always 1 for household head
    X_indiv_t <- matrix(0,nrow=length(M_t),ncol=p)
    phi_m_g_t <- t(phi[,(rep_G_t+((M_t-1)*FF))])
    for(ttt in 1:p){
      pr_indiv_t <- phi_m_g_t[,d_k_indiv_cum[ttt]:cumsum(d_k_indiv)[ttt]]
      Ran_unif_t <- runif(nrow(pr_indiv_t))
      cumul_t <- pr_indiv_t%*%upper.tri(diag(ncol(pr_indiv_t)),diag=TRUE)
      level_indiv_t <- level_indiv[[ttt]]
      X_indiv_t[,ttt] <- level_indiv_t[rowSums(Ran_unif_t > cumul_t) + 1L] 
    }
    comb_to_check <- X_indiv_t
    comb_to_check <- matrix(t(comb_to_check),byrow=T,nrow=n_batch[ss])
    comb_to_check <- cbind(X_house_t[,(q-p+1):q],comb_to_check) #add the household head before check; arrangement key!
    check_counter <- checkSZ(comb_to_check,(hh_size+1))
    t_1 <- t_1 + sum(check_counter);
    
    if(t_1 <= n_possibles){
      t_0 <- t_0 + (n_batch[ss] - sum(check_counter))
      X_indiv_struc <- rbind(X_indiv_struc,X_indiv_t[which(rep(check_counter,each=hh_size)==0),])
      X_house_struc <- rbind(X_house_struc,X_house_t[which(check_counter==0),])
      if(GenSyn){
        X_indiv_valid <- rbind(X_indiv_valid,X_indiv_t[which(rep(check_counter,each=hh_size)==1),])
        X_house_valid <- rbind(X_house_valid,X_house_t[which(check_counter==1),])
      }
      M_0 <- c(M_0,M_t[which(rep(check_counter,each=hh_size)==0)])
      G_0 <- c(G_0,G_t[which(check_counter==0)])
    } else {
      t_needed <- sum(check_counter) - (t_1 - n_possibles)
      index_needed <- which(cumsum(check_counter)==t_needed)[1]
      check_counter_needed <- check_counter[1:index_needed]
      t_0 <- t_0 + (length(check_counter_needed) - sum(check_counter_needed))
      X_indiv_struc <- rbind(X_indiv_struc,X_indiv_t[which(rep(check_counter_needed,each=hh_size)==0),])
      X_house_struc <- rbind(X_house_struc,X_house_t[which(check_counter_needed==0),])
      if(GenSyn){
        X_indiv_valid <- rbind(X_indiv_valid,X_indiv_t[which(rep(check_counter_needed,each=hh_size)==1),])
        X_house_valid <- rbind(X_house_valid,X_house_t[which(check_counter_needed==1),])
      }
      M_0 <- c(M_0,M_t[which(rep(check_counter_needed,each=hh_size)==0)])
      G_0 <- c(G_0,G_t[which(check_counter_needed==0)])
    }
  }
  n_0[ss] <- n_0[ss] + t_0
  
  if(GenSyn){
    return(list(X_house_struc=X_house_struc,X_indiv_struc=X_indiv_struc,
                X_house_valid=X_house_valid,X_indiv_valid=X_indiv_valid,G_0=G_0,M_0=M_0,n_0=n_0))
  } else {
    return(list(X_house_struc=X_house_struc,X_indiv_struc=X_indiv_struc,G_0=G_0,M_0=M_0,n_0=n_0))
  }
}


## Function to sample all impossible household data
SampleImpossibles <- 
  function(AllData,GlobalPara,StrucZerosPara,ModelPara,Weights,HHSize_Index,HHRelate_Index,GenSyn){
  
    AllStrucZerosData <- list()
    for(hh_size in GlobalPara$H){
      ImpossiblesData <- SampleImpossibles_h(AllData,GlobalPara,StrucZerosPara,ModelPara,Weights,HHSize_Index,
                                             HHRelate_Index,hh_size,GenSyn)
      AllStrucZerosData$X_house_struc <- rbind(AllStrucZerosData$X_house_struc,ImpossiblesData$X_house_struc)
      AllStrucZerosData$X_indiv_struc <- rbind(AllStrucZerosData$X_indiv_struc,ImpossiblesData$X_indiv_struc)
      #AllStrucZerosData$X_house_valid <- rbind(AllStrucZerosData$X_house_struc,ImpossiblesData$X_house_valid)
      #AllStrucZerosData$X_indiv_valid <- rbind(AllStrucZerosData$X_indiv_struc,ImpossiblesData$X_indiv_valid)
      AllStrucZerosData$G_0 <- c(AllStrucZerosData$G_0,ImpossiblesData$G_0)
      AllStrucZerosData$M_0 <- c(AllStrucZerosData$M_0,ImpossiblesData$M_0)
      AllStrucZerosData$n_0[which(GlobalPara$H==hh_size)] <- ImpossiblesData$n_0[which(GlobalPara$H==hh_size)]
    }
    AllStrucZerosData$n_i_0 <- as.numeric(as.character(AllStrucZerosData$X_house_struc[,HHSize_Index]))
    AllStrucZerosData$house_index_0 <- rep(c(1:sum(AllStrucZerosData$n_0)),AllStrucZerosData$n_i_0)
    AllStrucZerosData$n_i_index_0 <- rep(AllStrucZerosData$n_i_0,AllStrucZerosData$n_i_0)
    AllStrucZerosData$rep_G_0 <- rep(AllStrucZerosData$G_0,AllStrucZerosData$n_i_0)
    row.names(AllStrucZerosData$X_house_struc) <- NULL; 
    #row.names(AllStrucZerosData$X_house_valid) <- NULL
    AllStrucZerosData$X_house_struc <- as.data.frame(AllStrucZerosData$X_house_struc)
    #AllStrucZerosData$X_house_valid <- as.data.frame(AllStrucZerosData$X_house_valid)
    for(ii in 1:GlobalPara$q){
      AllStrucZerosData$X_house_struc[,ii] <- 
        factor(AllStrucZerosData$X_house_struc[,ii],levels=GlobalPara$level_house[[ii]])
      #AllStrucZerosData$X_house_valid[,ii] <- 
      factor(AllStrucZerosData$X_house_valid[,ii],levels=GlobalPara$level_house[[ii]])
    }
    AllStrucZerosData$X_indiv_struc <- as.data.frame(AllStrucZerosData$X_indiv_struc)
    #AllStrucZerosData$X_indiv_valid <- as.data.frame(AllStrucZerosData$X_indiv_valid)
    for(iii in 1:GlobalPara$p){
      AllStrucZerosData$X_indiv_struc[,iii] <- 
        factor(AllStrucZerosData$X_indiv_struc[,iii],levels=GlobalPara$level_indiv[[iii]])
      #AllStrucZerosData$X_indiv_valid[,iii] <- 
      factor(AllStrucZerosData$X_indiv_valid[,iii],levels=GlobalPara$level_indiv[[iii]])
    }
    colnames(AllStrucZerosData$X_house_struc) <- colnames(AllData$Y_house)
    colnames(AllStrucZerosData$X_indiv_struc) <- colnames(AllData$Y_indiv)
    #colnames(AllStrucZerosData$X_house_valid) <- colnames(AllData$Y_house)
    #colnames(AllStrucZerosData$X_indiv_valid) <- colnames(AllData$Y_indiv)
    
    return(AllStrucZerosData)
}


## Function to sample true responses for a given household
SampleTrueResponse_h <- function(AllData,GlobalPara,StrucZerosPara,ModelPara,ErrorIndicators,sss){
  lambda <- ModelPara$lambda; phi <- ModelPara$phi; FF <- ModelPara$FF; G <- ModelPara$G; M <- ModelPara$M; 
  level_house <- GlobalPara$level_house; level_indiv <- GlobalPara$level_indiv; p <- GlobalPara$p; q <- GlobalPara$q;
  house_index <- GlobalPara$house_index; 
  d_k_house_cum <- AllData$d_k_house_cum; d_k_indiv_cum <- AllData$d_k_indiv_cum;
  d_k_house <- AllData$d_k_house; d_k_indiv <- AllData$d_k_indiv;
  NA_miss_house <- AllData$NA_miss_house; NA_miss_indiv <- AllData$NA_miss_indiv;
  n_batch_imp <- StrucZerosPara$n_batch_imp; n_0_reject <- StrucZerosPara$n_0_reject;
  struc_zero_variables_house <- StrucZerosPara$struc_zero_variables_house;
  struc_zero_variables_indiv <- StrucZerosPara$struc_zero_variables_indiv;
  NA_error_house <- ErrorIndicators$NA_error_house; NA_error_indiv <- ErrorIndicators$NA_error_indiv
  epsilon_house <- ModelPara$epsilon_house; epsilon_indiv <- ModelPara$epsilon_indiv
  
  
  another_index <- which(house_index==sss); 
  X_house_sss_prop <- NA_miss_house[sss,]
  X_house_sss_prop <- matrix(rep(t(X_house_sss_prop),n_batch_imp[sss]),byrow=T,ncol=q)
  X_house_sss_prop_nf <- X_house_sss_prop
  X_indiv_sss_prop <- NA_miss_indiv[another_index,]
  X_indiv_sss_prop <- matrix(rep(t(X_indiv_sss_prop),n_batch_imp[sss]),byrow=T,ncol=p)
  X_indiv_sss_prop_nf <- X_indiv_sss_prop
  NA_error_house_sss <- NA_error_house[sss,]
  NA_error_house_sss <- matrix(rep(t(NA_error_house_sss),n_batch_imp[sss]),byrow=T,ncol=q)
  NA_error_indiv_sss <- NA_error_indiv[another_index,]
  NA_error_indiv_sss <- matrix(rep(t(NA_error_indiv_sss),n_batch_imp[sss]),byrow=T,ncol=p)
  G_prop <- rep(G[sss],n_batch_imp[sss])
  M_prop <- rep(M[another_index],n_batch_imp[sss])
  lambda_g <- t(lambda[,G_prop])
  phi_m_g <- t(phi[,(G[sss]+((M_prop-1)*FF))])
  check_counter_sss <- 0;
  while(check_counter_sss < 1){
    #First sample household data
    for(kkk in struc_zero_variables_house){
      level_house_k <- level_house[[kkk]]
      q_house_k <- 1/(d_k_house[kkk] - 1)
      if(length(which(is.na(NA_error_house_sss[,kkk])==TRUE))>0){
        corr_factor_Y_house_k <- matrix(0,ncol=d_k_house[kkk],nrow=length(G_prop))
        corr_factor_Y_house_k <- cbind(corr_factor_Y_house_k,X_house_sss_prop_nf[,kkk])
        corr_factor_house_k <- t(apply(corr_factor_Y_house_k,1,function(x) 
          replace(x,(1+x[(d_k_house[kkk]+1)]-min(level_house_k)),1-epsilon_house[which(struc_zero_variables_house==kkk)])))
        corr_factor_house_k <- corr_factor_house_k[,-ncol(corr_factor_house_k)]
        corr_factor_house_k[corr_factor_house_k==0] <- epsilon_house[which(struc_zero_variables_house==kkk)]*q_house_k
        pr_X_house_k <- lambda_g[,d_k_house_cum[kkk]:cumsum(d_k_house)[kkk]]*corr_factor_house_k
        pr_X_house_k <- t(apply(pr_X_house_k,1,function(x) x/sum(x)))
        #pr_X_house_k <- pr_X_house_k[which(is.na(NA_error_house_sss[,kkk])==TRUE),]
        Ran_unif_X_house_k <- runif(nrow(pr_X_house_k))
        cumul_X_house_k <- pr_X_house_k%*%upper.tri(diag(ncol(pr_X_house_k)),diag=TRUE)
        X_house_sss_prop[is.na(NA_error_house_sss[,kkk]),kkk] <- 
          level_house_k[rowSums(Ran_unif_X_house_k>cumul_X_house_k) + 1L]   
      } else if(length(which(NA_error_house_sss[,kkk]==1))>0){
        corr_factor_house_k <- matrix(1,ncol=d_k_house[kkk],nrow=length(G_prop))
        corr_factor_Y_house_k <- cbind(corr_factor_house_k,X_house_sss_prop_nf[sss,kkk])
        corr_factor_house_k <- t(apply(corr_factor_Y_house_k,1,function(x) replace(x,(1+x[(d_k_house[kkk]+1)]-min(level_house_k)),0)))
        corr_factor_house_k <- corr_factor_house_k[,-ncol(corr_factor_house_k)]
        corr_factor_house_k[corr_factor_house_k==1] <- epsilon_house[which(struc_zero_variables_house==kkk)]*q_house_k
        pr_X_house_k <- lambda_g[,d_k_house_cum[kkk]:cumsum(d_k_house)[kkk]]*corr_factor_house_k
        pr_X_house_k <- t(apply(pr_X_house_k,1,function(x) x/sum(x)))
        #pr_X_house_k <- pr_X_house_k[which(NA_error_house_sss[,kkk]==1),]
        Ran_unif_X_house_k <- runif(nrow(pr_X_house_k))
        cumul_X_house_k <- pr_X_house_k%*%upper.tri(diag(ncol(pr_X_house_k)),diag=TRUE)
        X_house_sss_prop[which(NA_error_house_sss[,kkk]==1),kkk] <- 
          level_house_k[rowSums(Ran_unif_X_house_k>cumul_X_house_k) + 1L]   
      }
    }
    #Then sample individual data
    for(kkkk in struc_zero_variables_indiv){
      level_indiv_k <- level_indiv[[kkkk]]
      q_indiv_k <- 1/(d_k_indiv[kkkk] - 1)
      if(length(which(is.na(NA_error_indiv_sss[,kkkk])==TRUE))>0){
        corr_factor_indiv_k <- matrix(0,ncol=d_k_indiv[kkkk],nrow=length(M_prop))
        corr_factor_Y_indiv_k <- cbind(corr_factor_indiv_k,X_indiv_sss_prop_nf[,kkkk])
        corr_factor_indiv_k <- t(apply(corr_factor_Y_indiv_k,1,function(x) 
          replace(x,(1+x[(d_k_indiv[kkkk]+1)]-min(level_indiv_k)),1-epsilon_indiv[which(struc_zero_variables_indiv==kkkk)])))
        corr_factor_indiv_k <- corr_factor_indiv_k[,-ncol(corr_factor_indiv_k)]
        corr_factor_indiv_k[corr_factor_indiv_k==0] <- epsilon_indiv[which(struc_zero_variables_indiv==kkkk)]*q_indiv_k
        pr_X_indiv_k <- phi_m_g[,d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]*corr_factor_indiv_k
        pr_X_indiv_k <- t(apply(pr_X_indiv_k,1,function(x) x/sum(x)))
        pr_X_indiv_k <- pr_X_indiv_k[which(is.na(NA_error_indiv_sss[,kkkk])==TRUE),]
        Ran_unif_X_indiv_k <- runif(nrow(pr_X_indiv_k))
        cumul_X_indiv_k <- pr_X_indiv_k%*%upper.tri(diag(ncol(pr_X_indiv_k)),diag=TRUE)
        X_indiv_sss_prop[is.na(NA_error_indiv_sss[,kkkk]),kkkk] <- level_indiv_k[rowSums(Ran_unif_X_indiv_k>cumul_X_indiv_k) + 1L]
      } else if(length(which(NA_error_indiv_sss[,kkkk]==1))>0){
        corr_factor_indiv_k <- matrix(1,ncol=d_k_indiv[kkkk],nrow=length(M_prop))
        corr_factor_Y_indiv_k <- cbind(corr_factor_indiv_k,X_indiv_sss_prop_nf[,kkkk])
        corr_factor_indiv_k <- t(apply(corr_factor_Y_indiv_k,1,function(x) replace(x,(1+x[(d_k_indiv[kkkk]+1)]-min(level_indiv_k)),0)))
        corr_factor_indiv_k <- corr_factor_indiv_k[,-ncol(corr_factor_indiv_k)]
        corr_factor_indiv_k[corr_factor_indiv_k==1] <- epsilon_indiv[which(struc_zero_variables_indiv==kkkk)]*q_indiv_k
        pr_X_indiv_k <- phi_m_g[,d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]*corr_factor_indiv_k
        pr_X_indiv_k <- t(apply(pr_X_indiv_k,1,function(x) x/sum(x)))
        pr_X_indiv_k <- pr_X_indiv_k[which(NA_error_indiv_sss[,kkkk]==1),]
        Ran_unif_X_indiv_k <- runif(nrow(pr_X_indiv_k))
        cumul_X_indiv_k <- pr_X_indiv_k%*%upper.tri(diag(ncol(pr_X_indiv_k)),diag=TRUE)
        X_indiv_sss_prop[which(NA_error_indiv_sss[,kkkk]==1),kkkk] <- level_indiv_k[rowSums(Ran_unif_X_indiv_k>cumul_X_indiv_k) + 1L]
      }
    }
    #Check edit rules
    comb_to_check <- matrix(t(X_indiv_sss_prop),nrow=n_batch_imp[sss],byrow=TRUE)
    comb_to_check <- cbind(X_house_sss_prop[,(q-p+1):q],comb_to_check)
    check_counter <- checkSZ(comb_to_check,(length(another_index) + 1))
    check_counter_sss <- check_counter_sss + sum(check_counter)
    if(length(which(check_counter==1))>0){
      n_0_reject[sss] <- n_0_reject[sss] + length(which(check_counter[1:which(check_counter==1)[1]]==0))
    } else{
      n_0_reject[sss] <- n_0_reject[sss] + n_batch_imp[sss]
    }
  }
  X_house <- X_house_sss_prop[which(check_counter==1)[1],]
  X_indiv <- matrix(comb_to_check[which(check_counter==1)[1],-c(1:p)],byrow=TRUE,ncol=p) #remove household head
  
  return(list(X_house=X_house, X_indiv=X_indiv, n_0_reject=n_0_reject))
}


## Function to sample true responses
SampleTrueResponse <- function(AllData,GlobalPara,StrucZerosPara,ModelPara,ErrorIndicators,HHSize_Index,HHRelate_Index){
  lambda <- ModelPara$lambda; phi <- ModelPara$phi; FF <- ModelPara$FF; 
  G <- ModelPara$G; M <- ModelPara$M; rep_G <- ModelPara$rep_G;
  level_house <- GlobalPara$level_house; level_indiv <- GlobalPara$level_indiv; p <- GlobalPara$p; q <- GlobalPara$q;
  d_k_house_cum <- AllData$d_k_house_cum; d_k_indiv_cum <- AllData$d_k_indiv_cum;
  d_k_house <- AllData$d_k_house; d_k_indiv <- AllData$d_k_indiv;
  nonstruc_zero_variables_house <- StrucZerosPara$nonstruc_zero_variables_house;
  nonstruc_zero_variables_indiv <- StrucZerosPara$nonstruc_zero_variables_indiv;
  NA_miss_house <- AllData$NA_miss_house; NA_miss_indiv <- AllData$NA_miss_indiv;
  
  #Sample errorneous data (struc zeros variables only!!!)
  AllData$n_0_reject <- StrucZerosPara$n_0_reject;
  for(sss in ErrorIndicators$z_i_index_house){
    ErrorData <- SampleTrueResponse_h(AllData,GlobalPara,StrucZerosPara,ModelPara,ErrorIndicators,sss)
    AllData$X_house[sss,] <- ErrorData$X_house
    AllData$X_indiv[which(GlobalPara$house_index==sss),] <- ErrorData$X_indiv
    AllData$n_0_reject[sss] <- ErrorData$n_0_reject[sss]
  }
  
  #Sample missing data (non-struc zeros variables only!!!)
  #Household data
  if(sum(is.na(NA_miss_house[,nonstruc_zero_variables_house])) > 0){
    lambda_g <- t(lambda[,G])
    for(kkk in nonstruc_zero_variables_house[!is.element(nonstruc_zero_variables_house,c(HHSize_Index,HHRelate_Index))]){
      if(length(which(is.na(NA_miss_house[,kkk])==TRUE))>0){
        pr_X_miss_q <- lambda_g[which(is.na(NA_miss_house[,kkk])==TRUE),d_k_house_cum[kkk]:cumsum(d_k_house)[kkk]]
        Ran_unif_miss_q <- runif(nrow(pr_X_miss_q))
        cumul_miss_q <- pr_X_miss_q%*%upper.tri(diag(ncol(pr_X_miss_q)),diag=TRUE)
        level_house_q <- level_house[[kkk]]
        AllData$X_house[is.na(NA_miss_house[,kkk]),kkk] <- level_house_q[rowSums(Ran_unif_miss_q>cumul_miss_q) + 1L]    
      }
    }
  }
  #Individual data
  if(sum(is.na(NA_miss_indiv[,nonstruc_zero_variables_indiv])) > 0){
    phi_m_g <- t(phi[,(rep_G+((M-1)*FF))])
    for(kkkk in nonstruc_zero_variables_indiv){
      if(length(which(is.na(NA_miss_indiv[,kkkk])==TRUE))>0){
        pr_X_miss_p <- phi_m_g[which(is.na(NA_miss_indiv[,kkkk])==TRUE),d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]
        Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
        cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
        level_indiv_p <- level_indiv[[kkkk]]
        AllData$X_indiv[is.na(NA_miss_indiv[,kkkk]),kkkk] <- level_indiv_p[rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L]
      }
    }
  }
  
  for(kq in 1:q){
    AllData$X_house[,kq] <- factor(AllData$X_house[,kq],levels=level_house[[kq]]) }
  for(kp in 1:p){
    AllData$X_indiv[,kp] <- factor(AllData$X_indiv[,kp],levels=level_indiv[[kp]]) }
  
  return(AllData)
}


## Function to sample household epsilon
SampleEpsilonHouse <- function(AllData,ErrorIndicators,ModelPara,StrucZerosPara){
  a_epsilon_house_star <- ModelPara$a_epsilon_house + 
    colSums(AllData$E_house[ErrorIndicators$z_i_index_house,StrucZerosPara$struc_zero_variables_house])
  b_epsilon_house_star <- ModelPara$b_epsilon_house + length(ErrorIndicators$z_i_index_house) - 
    colSums(AllData$E_house[ErrorIndicators$z_i_index_house,StrucZerosPara$struc_zero_variables_house])
  epsilon_house <- rbeta(length(StrucZerosPara$struc_zero_variables_house),t(t(a_epsilon_house_star)),t(t(b_epsilon_house_star)))
  
  return(epsilon_house)
}


## Function to sample individual epsilon
SampleEpsilonIndiv <- function(AllData,ErrorIndicators,ModelPara,StrucZerosPara){
  a_epsilon_indiv_star <- ModelPara$a_epsilon_indiv + 
    colSums(AllData$E_indiv[ErrorIndicators$z_i_index_indiv,StrucZerosPara$struc_zero_variables_indiv])
  b_epsilon_indiv_star <- ModelPara$b_epsilon_indiv + length(ErrorIndicators$z_i_index_indiv) - 
    colSums(AllData$E_indiv[ErrorIndicators$z_i_index_indiv,StrucZerosPara$struc_zero_variables_indiv])
  epsilon_indiv <- rbeta(length(StrucZerosPara$struc_zero_variables_indiv),t(t(a_epsilon_indiv_star)),t(t(b_epsilon_indiv_star)))
  
  return(epsilon_indiv)
}


## Function to sample household-level class
SampleG <- function(AllData,ModelPara,GlobalPara){
  ## First create indexes for phi and lambda; data.matrix function won't mess anything up as long as columns are coded as factors
  phi_index <- data.matrix(AllData$X_indiv)+AllData$FFF_indiv #has to be within loop for MI/EI
  lambda_index <- data.matrix(AllData$X_house)+AllData$FFF_house #has to be within loop  for MI/EI
  pr_G_post <- prGpost(phi_index,lambda_index,ModelPara$phi,ModelPara$lambda,ModelPara$omega,c(ModelPara$pii),ModelPara$FF,ModelPara$SS,GlobalPara$n_i)
  Ran_unif_G <- runif(nrow(pr_G_post))
  cumul_G <- pr_G_post%*%upper.tri(diag(ncol(pr_G_post)),diag=TRUE)
  G <- rowSums(Ran_unif_G>cumul_G) + 1L
  
  return(G)
}


## Function to sample individual-level class
SampleM <- function(AllData,ModelPara){
  phi_index <- data.matrix(AllData$X_indiv)+AllData$FFF_indiv
  pr_M_post <- prMpost(phi_index,ModelPara$phi,ModelPara$omega,ModelPara$rep_G,ModelPara$FF,ModelPara$SS)
  Ran_unif_M <- runif(nrow(pr_M_post))
  cumul_M <- pr_M_post%*%upper.tri(diag(ncol(pr_M_post)),diag=TRUE)
  M <- rowSums(Ran_unif_M>cumul_M) + 1L
  
  return(M)
}


## Function to sample household-level multinomial probabilities
SampleLambda <- function(AllData,ModelPara,GlobalPara,AllStrucZerosData,Weights){
  lambda <- ModelPara$lambda
  
  for(kk in 1:GlobalPara$q){
    lambda_count_table <- table(factor(ModelPara$G,levels=c(1:ModelPara$FF)),AllData$X_house[,kk])
    for(w_i in 1:length(Weights$struc_weight)){
      hh_size <- as.numeric(rownames(Weights$struc_weight)[w_i])
      w_i_index <-  AllStrucZerosData$n_i_0==hh_size
      G_0_w_i <- AllStrucZerosData$G_0[w_i_index]
      X_house_struc_w_i <- AllStrucZerosData$X_house_struc[w_i_index,]
      lambda_count_table <- 
        lambda_count_table + (table(factor(G_0_w_i,levels=c(1:ModelPara$FF)),X_house_struc_w_i[,kk])/Weights$struc_weight[w_i])
    }
    lambda[AllData$d_k_house_cum[kk]:cumsum(AllData$d_k_house)[kk],] <- 
      t(rdirichlet(ModelPara$FF,matrix(ModelPara$a_kdk + lambda_count_table,nrow=ModelPara$FF)))
  }
  
  return(lambda)
}


## Function to sample individual-level multinomial probabilities
SamplePhi <- function(AllData,ModelPara,GlobalPara,AllStrucZerosData,Weights){
  phi <- ModelPara$phi
  
  for(gg in 1:ModelPara$SS){
    for(ggg in 1:GlobalPara$p){
      phi_count_table <- table(factor(ModelPara$rep_G[which(ModelPara$M==gg)],levels=c(1:ModelPara$FF)),AllData$X_indiv[which(ModelPara$M==gg),ggg])
      for(w_i in 1:length(Weights$struc_weight)){
        hh_size <- as.numeric(rownames(Weights$struc_weight)[w_i])
        w_i_index <- AllStrucZerosData$n_i_index_0==hh_size
        rep_G_0_w_i <- AllStrucZerosData$rep_G_0[w_i_index]
        M_0_w_i <- AllStrucZerosData$M_0[w_i_index]
        X_indiv_struc_w_i <- AllStrucZerosData$X_indiv_struc[w_i_index,]
        phi_count_table <- phi_count_table + 
          (table(factor(rep_G_0_w_i[which(M_0_w_i==gg)],levels=c(1:ModelPara$FF)),
                 X_indiv_struc_w_i[which(M_0_w_i==gg),ggg])/Weights$struc_weight[w_i])
      }
      phi[AllData$d_k_indiv_cum[ggg]:cumsum(AllData$d_k_indiv)[ggg],(c(1:ModelPara$FF)+((gg-1)*ModelPara$FF))] <- 
        t(rdirichlet(ModelPara$FF,matrix(ModelPara$a_kdk + phi_count_table,nrow=ModelPara$FF)))
    }
  }
  
  return(phi)
}


## Function to sample household-level latent class probabilities
SamplePi <- function(ModelPara,AllStrucZerosData,Weights){
  pii <- ModelPara$pii; U <- pii;
  
  n_f <- matrix(summary(factor(ModelPara$G,levels=c(1:ModelPara$FF))),ncol=1)
  for(w_i in 1:length(Weights$struc_weight)){
    hh_size <- as.numeric(rownames(Weights$struc_weight)[w_i])
    w_i_index <- AllStrucZerosData$n_i_0==hh_size
    G_0_w_i <- AllStrucZerosData$G_0[w_i_index]
    n_f <- n_f + (matrix(summary(factor(G_0_w_i,levels=c(1:ModelPara$FF))),ncol=1)/Weights$struc_weight[w_i])
  }
  U[]<-1
  U[1:(ModelPara$FF-1),1] <- rbeta((ModelPara$FF-1),(1L+n_f[1:(ModelPara$FF-1)]),(ModelPara$alpha+(sum(n_f)-cumsum(n_f[-ModelPara$FF]))))
  if(length(which(U[-ModelPara$FF]==1))>0){
    U[which(U[-ModelPara$FF]==1)] <- 0.99999
  }
  one_min_U <- 1L-U
  one_min_U_prod <- c(1,cumprod(one_min_U[1:(ModelPara$FF-1)]))
  pii <- U*one_min_U_prod
  
  return(pii)
}


## Function to sample individual-level latent class probabilities
SampleOmega <- function(ModelPara,AllStrucZerosData,Weights){
  omega <- ModelPara$omega; V <- omega;
  
  M_G <- table(factor(ModelPara$rep_G,levels=c(1:ModelPara$FF)),factor(ModelPara$M,levels=c(1:ModelPara$SS)))
  for(w_i in 1:length(Weights$struc_weight)){
    hh_size <- as.numeric(rownames(Weights$struc_weight)[w_i])
    w_i_index <- AllStrucZerosData$n_i_index_0==hh_size
    rep_G_0_w_i <- AllStrucZerosData$rep_G_0[w_i_index]
    M_0_w_i <- AllStrucZerosData$M_0[w_i_index]
    M_G <- M_G + (table(factor(rep_G_0_w_i,levels=c(1:ModelPara$FF)),factor(M_0_w_i,levels=c(1:ModelPara$SS)))/Weights$struc_weight[w_i])
  }
  n_gm <- as.data.frame(M_G)$Freq
  V[]<-1
  no_V_to_sim <- (ModelPara$FF*(ModelPara$SS-1))
  V[,1:(ModelPara$SS-1)] <- rbeta(no_V_to_sim,(1L+n_gm[1:no_V_to_sim]),
                        (ModelPara$beta + c( matrix(rowSums(M_G),ncol=ModelPara$SS-1,nrow=ModelPara$FF)-
                                     t(apply(M_G,1,cumsum))[,1:ModelPara$SS-1])))
  if(length(which(V[,-ModelPara$SS]==1))>0){
    V[which(V[,-ModelPara$SS]==1)] <- 0.99999
  }
  one_min_V <- 1L-V
  one_min_V_prod <- cbind(1,t(apply(one_min_V[,-ModelPara$SS],1,cumprod)))
  omega <- V*one_min_V_prod
  
  return(omega)
}


## Function to check number of occupied individual-level clusters
CheckOccupied <- function(ModelPara,AllStrucZerosData){
  rep_G_all <- c(ModelPara$rep_G,AllStrucZerosData$rep_G_0); 
  G_all <- ModelPara$G_all; M_all <- c(ModelPara$M,AllStrucZerosData$M_0)
  S_occup <- NULL
  for(occ in sort(unique(G_all))){
    S_occup <- rbind(S_occup,dim(table(rep_G_all[which(rep_G_all==occ)],M_all[which(rep_G_all==occ)]))[2])
  }
  
  return(S_occup)
}







