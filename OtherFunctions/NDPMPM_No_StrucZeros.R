
###################################### START ######################################
############### Function for fitting NDPMPM without structural zeros and with missing data
############### Returns a list containing all parameters and MM imputed datasets of same size:: CR
############### The MM imputed datasets are stacked by rows:: CR
############### Also returns n_prop posterior draws for households with individuals that have missing entries/data
############### The posterior draws for the households are stacked by rows
library(DirichletReg)
library(matrixStats)
library(coda)
Rcpp::sourceCpp('CppFunctions/prGpost.cpp')
Rcpp::sourceCpp('CppFunctions/prMpost.cpp')
Rcpp::sourceCpp('CppFunctions/checkSZ.cpp')

fit_NDPMPM <- function(Data_house,Data_indiv,FF,SS,n_iter,burn_in,MM,n_prop,struc_zero,valid_prop,mc_thin,
                       save_imp,save_prop){
  ####### struc_zero = T to fit the model without augmentation but ensure sampled values are valid at every iteration
  ####### valid_prop = T to fit the model without augmentation but ensure sampled values are valid only for proposals/imputations
  ####### Data_house is the household data
  ####### Data_indiv is the individual data
  ####### FF and SS are the number of clusters for household and individuals respectively
  ####### n_iter is the number of MCMC iterations, burn_in is the number of burn-in
  ####### MM is the number of imputed datasets to generate -- spaced mc_thin times apart in the chain
  ####### N_prop is the number of proposed completions to generate -- spaced mc_thin times apart in the chain
  ####### Hyper-priors are set as in the paper
  ####### Parameters are initialized at observed estimates
  
  ###### 1: Set values for n, N, p, q, and so on
  N <- nrow(Data_indiv)
  n <- nrow(Data_house)
  n_i <- as.numeric(as.character(Data_house[,1]))
  p <- ncol(Data_indiv)
  q <- ncol(Data_house)
  house_index <- rep(c(1:n),n_i) 
  level_indiv <- vector("list",p)
  for(i in 1:p){
      level_indiv[[i]] <- c(as.numeric(levels(Data_indiv[,i])))
  }
  level_house <- vector("list",q)
  for(i in 1:q){
      level_house[[i]] <- c(as.numeric(levels(Data_house[,i])))
  }
  
  ###### 2: Fill missing values with plausible entries
  NA_indiv <- Data_indiv;
  if(sum(is.na(NA_indiv)) > 0){
    for (ii in 1:ncol(Data_indiv)){
      Data_indiv[is.na(Data_indiv[,ii]),ii] <- 
        sample(level_indiv[[ii]],length(Data_indiv[is.na(Data_indiv[,ii]),ii]),replace=T,
               prob=summary(na.omit(Data_indiv[,ii])))
    }
  }
  NA_house <- Data_house;
  if(sum(is.na(NA_house)) > 0){
    for (jj in 2:ncol(Data_house)){
      Data_house[is.na(Data_house[,jj]),jj] <- 
        sample(level_house[[jj]],length(Data_house[is.na(Data_house[,jj]),jj]),replace=T,
               prob=summary(na.omit(Data_house[,jj])))
    }
  }
  struc_zero_variables <- c(1,4,5)
  nonstruc_zero_variables <- c(2,3)
  Indiv_miss_index_HH <- sort(unique(house_index[which(complete.cases(NA_indiv[,struc_zero_variables])==FALSE)]))
  n_miss <- length(Indiv_miss_index_HH)
  Indiv_miss_index <- which(is.element(house_index,Indiv_miss_index_HH)==TRUE)
  
  
  ###### 3: Calculate observed proportions and number of categories for each variable
  d_k_house <- d_k_indiv <- ini_marg_house <- ini_marg_indiv <- NULL
  for(k in 1:q){
    d_k_house <- cbind(d_k_house,nlevels(Data_house[,k]))
    ini_marg_k <- as.data.frame(table(Data_house[,k]))$Freq/n
    ini_marg_k <- matrix(ini_marg_k,ncol=1)
    ini_marg_house <- rbind(ini_marg_house,ini_marg_k)
  }
  for(k in 1:p){
    d_k_indiv <- cbind(d_k_indiv,nlevels(Data_indiv[,k]))
    ini_marg_k <- as.data.frame(table(Data_indiv[,k]))$Freq/N
    ini_marg_k <- matrix(ini_marg_k,ncol=1)
    ini_marg_indiv <- rbind(ini_marg_indiv,ini_marg_k)
  }
  
  ###### 4: Initialize chain
  alpha = beta = 1
  a_kdk = 1
  a_alpha <- b_alpha <- a_beta <- b_beta <- 0.25
  pr_G_post <- matrix(0,ncol=FF,nrow=n)
  pii = matrix(0,nrow=FF)
  lambda = matrix(rep(ini_marg_house,FF),ncol=FF)
  pr_M_post = matrix(0,ncol=SS,nrow=N)
  omega = matrix(0,nrow=FF,ncol=SS)
  phi = matrix(0,nrow=length(ini_marg_indiv),ncol=FF*SS) #make phi matrix and not array for c++
  for(gm in 1:(FF*SS)){
    phi[,gm] = ini_marg_indiv
  }
  U = matrix(rbeta(FF,1,alpha),nrow=FF)
  V = matrix(rbeta((FF*SS),1,beta),nrow=FF,ncol=SS)
  U[FF]=1
  V[,SS]=1
  one_min_U = 1L-U
  one_min_U = c(1,cumprod(one_min_U[1:(FF-1)]))
  one_min_V = 1L-V
  one_min_V = cbind(1,t(apply(one_min_V[,-SS],1,cumprod)))
  pii = U*one_min_U
  omega = V*one_min_V
  d_k_house_cum = 1+cumsum(c(0,d_k_house[,-q]))
  d_k_indiv_cum = 1+cumsum(c(0,d_k_indiv[,-p]))
  dp_imput_house_nz <- dp_imput_indiv_nz <- NULL
  M_to_use_mc <- sort(sample(seq((burn_in +1),n_iter,by=mc_thin),MM,replace=F))
  n_prop_to_use_mc <- sort(sample(seq((burn_in +1),n_iter,by=mc_thin),n_prop,replace=F))
  FFF_indiv = matrix(rep(cumsum(c(0,d_k_indiv[,-p])),each=N),ncol=p)
  FFF_house = matrix(rep(cumsum(c(0,d_k_house[,-q])),each=n),ncol=q)
  ALPHA = BETA = G_CLUST = M_CLUST = NULL
  #PII = NULL
  #LAMBDA = matrix(0,ncol=ncol(lambda),nrow=nrow(lambda))
  #OMEGA = matrix(0,ncol=ncol(omega),nrow=nrow(omega))
  #PHI <- matrix(0,ncol=ncol(phi),nrow=nrow(phi))
  DATA_INDIV_MISS <- NULL
  
  n_batch_imp_init <- rep(50,n_miss) #sample imputations in batches before checking constraints
  n_0_reject <- rep(0,n_miss)
  prop_batch <- 1.2
  
  ###### 5: Run MCMC
  for(mc in 1:n_iter){
    cat(paste("Iteration ", mc,"\n", sep = ""))
    proc_t <- proc.time()
    #sample G
    phi_index = data.matrix(Data_indiv)+FFF_indiv
    lambda_index = data.matrix(Data_house)+FFF_house
    pr_G_post = prGpost(phi_index,lambda_index,phi,lambda,omega,c(pii),FF,SS,n_i)
    Ran_unif_G = runif(nrow(pr_G_post))
    cumul_G = pr_G_post%*%upper.tri(diag(ncol(pr_G_post)),diag=TRUE)
    G = rowSums(Ran_unif_G>cumul_G) + 1L
    
    #sample M
    rep_G = rep(G,n_i)
    pr_M_post = prMpost(phi_index,phi,omega,rep_G,FF,SS)
    Ran_unif_M = runif(nrow(pr_M_post))
    cumul_M = pr_M_post%*%upper.tri(diag(ncol(pr_M_post)),diag=TRUE)
    M = rowSums(Ran_unif_M>cumul_M) + 1L
    
    #sample phi
    for(gg in 1:SS){
      for(ggg in 1:p){
        phi[d_k_indiv_cum[ggg]:cumsum(d_k_indiv)[ggg],(c(1:FF)+((gg-1)*FF))] = 
          t(rdirichlet(FF,matrix(a_kdk + table(factor(rep_G[which(M==gg)],levels=c(1:FF)),
                                               Data_indiv[which(M==gg),ggg]),nrow=FF)))
      }
    }
    
    #sample lambda
    for(kk in 1:q){
      lambda[d_k_house_cum[kk]:cumsum(d_k_house)[kk],] = 
        t(rdirichlet(FF,matrix(a_kdk + table(factor(G,levels=c(1:FF)),Data_house[,kk]),nrow=FF)))
    }
    
    #sample U and pii
    n_f = matrix(summary(factor(G,levels=c(1:FF))),ncol=1)
    U[FF]=1
    U[1:(FF-1),1] = rbeta((FF-1),(1L+n_f[1:(FF-1)]),(alpha+(n-cumsum(n_f[-FF]))))
    if(length(which(U[-FF]==1))>0){
      U[which(U[-FF]==1)] = 0.99999
    }
    one_min_U = 1L-U
    one_min_U_prod = c(1,cumprod(one_min_U[1:(FF-1)]))
    pii = U*one_min_U_prod
    
    #sample V and omega
    M_G = table(factor(rep_G,levels=c(1:FF)),factor(M,levels=c(1:SS)))
    n_gm = as.data.frame(M_G)$Freq
    V[,SS]=1
    no_V_to_sim = (FF*(SS-1))
    V[,1:(SS-1)] = rbeta(no_V_to_sim,(1L+n_gm[1:no_V_to_sim]),
                         (beta + c( matrix(rowSums(M_G),ncol=SS-1,nrow=FF)-
                                      t(apply(M_G,1,cumsum))[,1:SS-1])))
    if(length(which(V[,-SS]==1))>0){
      V[which(V[,-SS]==1)] = 0.99999
    }
    one_min_V = 1L-V
    one_min_V_prod = cbind(1,t(apply(one_min_V[,-SS],1,cumprod)))
    omega = V*one_min_V_prod
    
    #sample alpha
    alpha = rgamma(1,shape=(a_alpha+FF-1),rate=(b_alpha-log(pii[FF])))
    
    #sample beta
    beta = rgamma(1,shape=(a_beta+(FF*(SS-1))),rate=(b_beta-sum(log(omega[,SS]))))
    
    #check number of occupied clusters
    S_occup = NULL
    for(occ in sort(unique(G))){
      S_occup = rbind(S_occup,dim(table(rep_G[which(rep_G==occ)],M[which(rep_G==occ)]))[2])
    }
    cat(paste("Number of Occupied Household Classes is ", length(unique(G)), "\n", sep = ''))
    cat(paste("Max Number of Occupied Individual Classes is ", max(S_occup), "\n", sep = ''))
    
    
    #sample missing data
    #first household
    if(sum(is.na(NA_house)) > 0){
      lambda_g <- t(lambda[,G])
      for(kkk in 2:q){
        if(length(which(is.na(NA_house[,kkk])==TRUE))>0){
          pr_X_miss_q <- lambda_g[which(is.na(NA_house[,kkk])==TRUE), d_k_house_cum[kkk]:cumsum(d_k_house)[kkk]]
          Ran_unif_miss_q <- runif(nrow(pr_X_miss_q))
          cumul_miss_q <- pr_X_miss_q%*%upper.tri(diag(ncol(pr_X_miss_q)),diag=TRUE)
          level_house_q <- level_house[[kkk]]
          Data_house[which(is.na(NA_house[,kkk])==TRUE),kkk] <- level_house_q[rowSums(Ran_unif_miss_q>cumul_miss_q) + 1L]    
        }
      }
    }
    #now individuals
    #first sample non structural zero variables
    if(sum(is.na(NA_indiv)) > 0){
      phi_m_g <- matrix(0,nrow=N,ncol=dim(phi)[1])
      for(jjj in 1:N){
        phi_m_g[jjj,] <- phi[,(rep_G[jjj]+((M[jjj]-1)*FF))]
      }
      for(kkkk in nonstruc_zero_variables){
        if(length(which(is.na(NA_indiv[,kkkk])==TRUE))>0){
          pr_X_miss_p <- phi_m_g[which(is.na(NA_indiv[,kkkk])==TRUE),d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]
          Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
          cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
          level_indiv_p <- level_indiv[[kkkk]]
          Data_indiv[which(is.na(NA_indiv[,kkkk])==TRUE),kkkk] <- level_indiv_p[rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L]
        }
      }
      if(!struc_zero){
        for(kkkk in struc_zero_variables){
          if(length(which(is.na(NA_indiv[,kkkk])==TRUE))>0){
            pr_X_miss_p <- phi_m_g[which(is.na(NA_indiv[,kkkk])==TRUE),d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]
            Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
            cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
            level_indiv_p <- level_indiv[[kkkk]]
            Data_indiv[which(is.na(NA_indiv[,kkkk])==TRUE),kkkk] <-
              level_indiv_p[rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L]
          }
        }
      } else{
        for(sss in 1:n_miss){
          n_batch_imp <- n_batch_imp_init[sss] + ceiling(n_0_reject[sss]*prop_batch) #no. of batches of imp.s to sample
          n_0_reject[sss] <- 0
          another_index <- which(is.element(house_index,Indiv_miss_index_HH[sss])==TRUE)
          n_another_index <- length(another_index) + 1
          NA_indiv_prop <- Data_indiv[another_index,]
          NA_indiv_prop[,struc_zero_variables] <- NA_indiv[another_index,struc_zero_variables]
          NA_indiv_prop <- apply(NA_indiv_prop,2,function(x) as.numeric(as.character(x)))
          NA_indiv_prop <- matrix(rep(t(NA_indiv_prop),n_batch_imp),byrow=TRUE,ncol=p)
          rep_G_prop <- rep(rep_G[another_index],n_batch_imp)
          rep_M_prop <- rep(M[another_index],n_batch_imp)
          phi_m_g <- t(phi[,(rep_G_prop+((rep_M_prop-1)*FF))])
          check_counter_sss <- 0;
          while(check_counter_sss < 1){
            Data_indiv_prop <- NA_indiv_prop
            for(kkkk in struc_zero_variables){
              if(length(which(is.na(NA_indiv_prop[,kkkk])==TRUE))>0){
                pr_X_miss_p <- matrix(t(phi_m_g[which(is.na(NA_indiv_prop[,kkkk])==TRUE),
                                                d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]),
                                      nrow=length(which(is.na(NA_indiv_prop[,kkkk])==TRUE)),byrow=T)
                Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
                cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
                level_indiv_p <- level_indiv[[kkkk]]
                Data_indiv_prop[is.na(NA_indiv_prop[,kkkk]),kkkk] <- 
                  level_indiv_p[rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L]
              }
            }
            #Check edit rules
            comb_to_check <- matrix(t(Data_indiv_prop),nrow=n_batch_imp,byrow=TRUE)
            NA_house_prop <- Data_house[Indiv_miss_index_HH[sss],(q-p+1):q] 
            comb_to_check_HH <- matrix(rep(apply(NA_house_prop,2,function(x) as.numeric(as.character(x))),
                                          n_batch_imp),nrow=n_batch_imp,byrow = T)
            comb_to_check <- cbind(comb_to_check_HH,comb_to_check)
            check_counter <- checkSZ(comb_to_check,n_another_index)
            check_counter_sss <- check_counter_sss + sum(check_counter)
            if(length(which(check_counter==1))>0){
              n_0_reject[sss] <- n_0_reject[sss] + length(which(check_counter[1:which(check_counter==1)[1]]==0))
            } else{
              n_0_reject[sss] <- n_0_reject[sss] + n_batch_imp
            }
          }
          Data_indiv[another_index,] <-
            matrix(comb_to_check[which(check_counter==1)[1],-c(1:p)],byrow=TRUE,ncol=p) #remove household head
        }
      }
    }
    
    #save and sample synthetic data
    if(mc > burn_in){
      #PII = rbind(PII,c(pii))
      ALPHA = rbind(ALPHA,alpha)
      BETA = rbind(BETA,beta)
      G_CLUST <- rbind(G_CLUST,length(unique(G)))
      M_CLUST <- rbind(M_CLUST,max(S_occup))
      #LAMBDA = LAMBDA + lambda
      #OMEGA = OMEGA + omega
      #PHI <- PHI + phi
      
      if(sum(mc==M_to_use_mc)==1 && save_imp){
        dp_imput_indiv_nz <- rbind(dp_imput_indiv_nz,Data_indiv)  
        dp_imput_house_nz <- rbind(dp_imput_house_nz,Data_house)
      }
      
      if(sum(mc==n_prop_to_use_mc)==1 && save_prop){
        if(!struc_zero && valid_prop){
          if(sum(is.na(NA_indiv[,struc_zero_variables])) > 0){
            n_miss_struc <- length(sort(unique(house_index[which(complete.cases(NA_indiv[,struc_zero_variables])==FALSE)])))
            for(sss in 1:n_miss_struc){
              n_batch_imp <- n_batch_imp_init[sss] + ceiling(n_0_reject[sss]*prop_batch) #no. of batches of imp.s to sample
              n_0_reject[sss] <- 0
              another_index <- which(is.element(house_index,Indiv_miss_index_HH[sss])==TRUE)
              n_another_index <- length(another_index) + 1
              NA_indiv_prop <- Data_indiv[another_index,]
              NA_indiv_prop[,struc_zero_variables] <- NA_indiv[another_index,struc_zero_variables]
              NA_indiv_prop <- apply(NA_indiv_prop,2,function(x) as.numeric(as.character(x)))
              NA_indiv_prop <- matrix(rep(t(NA_indiv_prop),n_batch_imp),byrow=TRUE,ncol=p)
              rep_G_prop <- rep(rep_G[another_index],n_batch_imp)
              rep_M_prop <- rep(M[another_index],n_batch_imp)
              phi_m_g <- t(phi[,(rep_G_prop+((rep_M_prop-1)*FF))])
              check_counter_sss <- 0;
              while(check_counter_sss < 1){
                Data_indiv_prop <- NA_indiv_prop
                for(kkkk in struc_zero_variables){
                  if(length(which(is.na(NA_indiv_prop[,kkkk])==TRUE))>0){
                    pr_X_miss_p <- matrix(t(phi_m_g[which(is.na(NA_indiv_prop[,kkkk])==TRUE),
                                                    d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]),
                                          nrow=length(which(is.na(NA_indiv_prop[,kkkk])==TRUE)),byrow=T)
                    Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
                    cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
                    level_indiv_p <- level_indiv[[kkkk]]
                    Data_indiv_prop[is.na(NA_indiv_prop[,kkkk]),kkkk] <- 
                      level_indiv_p[rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L]
                  }
                }
                #Check edit rules
                comb_to_check <- matrix(t(Data_indiv_prop),nrow=n_batch_imp,byrow=TRUE)
                NA_house_prop <- Data_house[Indiv_miss_index_HH[sss],(q-p+1):q] 
                comb_to_check_HH <- matrix(rep(apply(NA_house_prop,2,function(x) as.numeric(as.character(x))),
                                               n_batch_imp),nrow=n_batch_imp,byrow = T)
                comb_to_check <- cbind(comb_to_check_HH,comb_to_check)
                check_counter <- checkSZ(comb_to_check,n_another_index)
                check_counter_sss <- check_counter_sss + sum(check_counter)
                if(length(which(check_counter==1))>0){
                  n_0_reject[sss] <- n_0_reject[sss] + length(which(check_counter[1:which(check_counter==1)[1]]==0))
                } else{
                  n_0_reject[sss] <- n_0_reject[sss] + n_batch_imp
                }
              }
              Data_indiv[another_index,] <-
                matrix(comb_to_check[which(check_counter==1)[1],-c(1:p)],byrow=TRUE,ncol=p) #remove household head
            }
          }
        }
        DATA_INDIV_MISS <- rbind(DATA_INDIV_MISS,Data_indiv[Indiv_miss_index,])
      }
    }
    
    #print
    elapsed_time <- (proc.time() - proc_t)[["elapsed"]]
    cat(paste("Number of Sampled Rejections for Missing Data is ",
              ifelse(sum(mc==n_prop_to_use_mc)==1 |sum(mc==M_to_use_mc)==1 |struc_zero,sum(n_0_reject),0),"\n",sep =''))
    cat(paste("Elapsed Time = ", elapsed_time, "\n\n", sep = ' '))
  }
  #list(dp_imput_indiv_nz=dp_imput_indiv_nz,dp_imput_house_nz=dp_imput_house_nz,DATA_INDIV_MISS=DATA_INDIV_MISS,
  #     DATA_HOUSE_MISS=DATA_HOUSE_MISS,PII=PII,ALPHA=ALPHA,BETA=BETA,LAMBDA=LAMBDA,OMEGA=OMEGA,PHI=PHI)
  if(save_imp && !save_prop){
    list(dp_imput_indiv_nz=dp_imput_indiv_nz,dp_imput_house_nz=dp_imput_house_nz,
         ALPHA_nz=ALPHA,BETA_nz=BETA,M_CLUST_nz=M_CLUST,G_CLUST_nz=G_CLUST)
  } else if(save_prop && !save_imp){
    list(DATA_INDIV_MISS=DATA_INDIV_MISS,
         ALPHA_nz=ALPHA,BETA_nz=BETA,M_CLUST_nz=M_CLUST,G_CLUST_nz=G_CLUST)
  } else if(save_prop && save_imp){
    list(dp_imput_indiv_nz=dp_imput_indiv_nz,dp_imput_house_nz=dp_imput_house_nz,
         DATA_INDIV_MISS=DATA_INDIV_MISS,
         ALPHA_nz=ALPHA,BETA_nz=BETA,M_CLUST_nz=M_CLUST,G_CLUST_nz=G_CLUST)
  }
}




