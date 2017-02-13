

###################################################################################
###################################################################################
###################### Error Location Model: Exploring Options ####################
###################################################################################
###################################################################################

###################################### START ######################################
###################### Phase One: Load Packages and Functions #####################
rm(list = ls())
Rcpp::sourceCpp('prZpost.cpp')
#library(DirichletReg)

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

################### Phase Three: Make function for Gibbs Sampler ##################
Run_Model <- function(E_house,E_indiv,TT,n_iter,burn_in){
  
  ###### 1: Set global parameters for data
  N <- nrow(E_indiv)
  n <- nrow(E_house)
  n_i <- as.numeric(as.character(rep(names(n_h),n_h)))
  p <- ncol(E_indiv)
  q <- ncol(E_house)
  #house_index <- rep(c(1:n),n_i)
  #n_i_index <- rep(n_i,n_i)
  
  ###### 2: Initialize parameters
  a_epsilon_j <- 1; b_epsilon_j <- 1
  gamma <- 1
  a_gamma <- b_gamma <- 0.25
  epsilon_house <- matrix(c(0.9,0.1),ncol=TT,nrow=(q*2))
  epsilon_indiv <- matrix(c(0.9,0.1),ncol=TT,nrow=(p*2))
  r <- matrix(rbeta(TT,1,gamma),nrow=TT)
  r[TT] <- 1
  one_min_r <- 1L-r
  one_min_r <- c(1,cumprod(one_min_r[1:(TT-1)]))
  zeta <- r*one_min_r
  epsilon_indiv_index <- data.matrix(E_indiv)+ matrix(rep(cumsum(c(0,rep(2,(p-1)))),each=N),ncol=p)
  epsilon_house_index <- data.matrix(E_house)+ matrix(rep(cumsum(c(0,rep(2,(q-1)))),each=n),ncol=q)
  d_k_indiv_cum <- 1+cumsum(c(0,rep(2,(p-1))))
  d_k_house_cum <- 1+cumsum(c(0,rep(2,(q-1))))
  ZETA <- GAMMA <- Z_CLUST <- NULL
  EPSILON_SUM <- 0
  EPSILON <- matrix(0,ncol=(TT*((p+q)*2)),nrow=(n_iter-burn_in))
  
  ###### 3: The MCMC
  ## Begin MCMC
  for(mc in 1:n_iter){
    proc_t <- proc.time()
    
    ## Sample the z_i's
    post_pr_z <- prZpost(epsilon_indiv_index,epsilon_house_index,epsilon_indiv,epsilon_house,zeta,TT,n_i)
    Ran_unif_z <- runif(nrow(post_pr_z))
    cumul_z <- post_pr_z%*%upper.tri(diag(ncol(post_pr_z)),diag=TRUE)
    z <- rowSums(Ran_unif_z>cumul_z) + 1L
    remove(post_pr_z); remove(Ran_unif_z); remove(cumul_z);
    
    ## Sample Epsilon_house
    count_table_house <- matrix(apply(E_house,2,function(x) t(table(factor(z,levels=c(1:TT)),x))),byrow=T,ncol=2)
    beta_sample_house <- 
      rbeta((TT*q),matrix(a_epsilon_j + count_table_house[,1],nrow=(TT*q)),matrix(b_epsilon_j + count_table_house[,2],nrow=(TT*q)))
    for(kq in 1:q){
      beta_sample_house_kq <- beta_sample_house[1:TT+ TT*(kq-1)]
      epsilon_house[d_k_house_cum[kq]:cumsum(rep(2,q))[kq],] <- rbind(beta_sample_house_kq,1-beta_sample_house_kq)
    }
    remove(count_table_house)
    
    ## Sample Epsilon_indiv
    rep_z <- rep(z,n_i)
    count_table_indiv <- matrix(apply(E_indiv,2,function(x) t(table(factor(rep_z,levels=c(1:TT)),x))),byrow=T,ncol=2)
    beta_sample_indiv <- 
      rbeta((TT*p),matrix(a_epsilon_j + count_table_indiv[,1],nrow=(TT*p)),matrix(b_epsilon_j + count_table_indiv[,2],nrow=(TT*p)))
    for(kp in 1:p){
      beta_sample_indiv_kp <- beta_sample_indiv[1:TT+ TT*(kp-1)]
      epsilon_indiv[d_k_indiv_cum[kp]:cumsum(rep(2,p))[kp],] <- rbind(beta_sample_indiv_kp,1-beta_sample_indiv_kp)
    }
    remove(count_table_indiv)
    
    ## Sample r and zeta
    n_t <- matrix(summary(factor(z,levels=c(1:TT))),ncol=1)
    r[TT]<-1
    r[1:(TT-1),1] <- rbeta((TT-1),(1L+n_t[1:(TT-1)]),(gamma+(sum(n_t)-cumsum(n_t[-TT]))))
    if(length(which(r[-TT]==1))>0){
      r[which(r[-TT]==1)] <- 0.999999999
    }
    one_min_r <- 1L-r
    one_min_r_prod <- c(1,cumprod(one_min_r[1:(TT-1)]))
    zeta <- r*one_min_r_prod
    remove(n_t)
    
    ## Sample gamma
    gamma <- rgamma(1,shape=(a_gamma+TT-1),rate=(b_gamma-log(zeta[TT])))
    
    ## Sample synthetic data
    E_house_new <- E_house
    epsilon_house_t <- t(epsilon_house_t[,z])
    for(kkq in 1:q){
      pr_new_E_q <- epsilon_house_t
    }
    
    
    for(kkk in 2:q){
      if(length(which(is.na(NA_house[,kkk])==TRUE))>0){
        pr_X_miss_q <- lambda_g[which(is.na(NA_house[,kkk])==TRUE),d_k_house_cum[kkk]:cumsum(d_k_house)[kkk]]
        Ran_unif_miss_q <- runif(nrow(pr_X_miss_q))
        cumul_miss_q <- pr_X_miss_q%*%upper.tri(diag(ncol(pr_X_miss_q)),diag=TRUE)
        level_house_q <- level_house[[kkk]]
        Data_house[which(is.na(NA_house[,kkk])==TRUE),kkk] <- level_house_q[rowSums(Ran_unif_miss_q>cumul_miss_q) + 1L]    
      }
    }
    
    ## Print a few summaries
    if(sum(mc == seq(1,n_iter,by=50))==1){
      cat(paste("Iteration ", mc,"\n", sep = ""))
      cat(paste("Number of Occupied Latent Classes is ", length(unique(z)), "\n", sep = ''))
      elapsed_time <- (proc.time() - proc_t)[["elapsed"]]
      cat(paste("Elapsed Time = ", elapsed_time, "\n\n", sep = ' '))
    } else if(mc==n_iter){
      cat(paste("Iteration ", mc,"\n", sep = ""))
      cat(paste("Number of Occupied Latent Classes is ", length(unique(z)), "\n", sep = ''))
      elapsed_time <- (proc.time() - proc_t)[["elapsed"]]
      cat(paste("Elapsed Time = ", elapsed_time, "\n\n", sep = ' '))
    }
    
    ## Save posterior samples
    if(mc > burn_in){
      ZETA <- rbind(ZETA,c(zeta))
      GAMMA <- rbind(GAMMA,gamma)
      Z_CLUST <- rbind(Z_CLUST,length(unique(z)))
      EPSILON[(mc-burn_in),] <- c(rbind(epsilon_indiv,epsilon_house))
      EPSILON_SUM <- EPSILON_SUM + rbind(epsilon_indiv,epsilon_house)
      
      if(sum(mc==M_to_use_mc)==1){
        dp_imput_indiv <- rbind(dp_imput_indiv,Data_indiv)  
        dp_imput_house <- rbind(dp_imput_house,Data_house)
      }
    }
  }
  ## End MCMC
  ## Return results
  list(ZETA=ZETA, GAMMA=GAMMA, Z_CLUST=Z_CLUST, EPSILON=EPSILON, EPSILON_SUM=EPSILON_SUM)
}

############################### Phase Four: Fit Model #############################
TT <- 5
n_iter <- 10000
burn_in <- n_iter/2
ModelFit <- Run_Model(E_house,E_indiv,TT,n_iter,burn_in)

mean(ModelFit$Z_CLUST)
ModelFit$EPSILON_SUM[c(2,4,6),]/burn_in
epsilon_indiv_truth
ModelFit$EPSILON_SUM[c(8,10),]/burn_in
epsilon_house_truth




