
###################################################################################
###################################################################################
####################### Error Location Model: Other Functions #####################
###################################################################################
###################################################################################

################### One: Function for Gibbs Sampler ##################
Run_Model <- function(E_house,E_indiv,TT,n_iter,burn_in,MM,mc_thin){
  
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
  M_to_use_mc <- sort(sample(seq((burn_in +1),n_iter,by=mc_thin),MM,replace=F))
  ZETA <- GAMMA <- Z_CLUST <- NULL
  EPSILON_SUM <- 0
  EPSILON <- matrix(0,ncol=(TT*((p+q)*2)),nrow=(n_iter-burn_in))
  Imputations_indiv <- Imputations_house <- vector("list", MM)
  imp_index <- 0
  
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
    
    ## Sample Epsilon_house and synthetic household data
    if(sum(mc==M_to_use_mc)==1){
      E_house_new <- E_house 
    }
    count_table_house <- matrix(apply(E_house,2,function(x) t(table(factor(z,levels=c(1:TT)),x))),byrow=T,ncol=2)
    beta_sample_house <- 
      rbeta((TT*q),matrix(a_epsilon_j + count_table_house[,1],nrow=(TT*q)),matrix(b_epsilon_j + count_table_house[,2],nrow=(TT*q)))
    for(kq in 1:q){
      beta_sample_house_kq <- beta_sample_house[1:TT+ TT*(kq-1)]
      epsilon_house[d_k_house_cum[kq]:cumsum(rep(2,q))[kq],] <- rbind(beta_sample_house_kq,1-beta_sample_house_kq)
      if(sum(mc==M_to_use_mc)==1){
        pr_new_E_house_q <- t(epsilon_house[d_k_house_cum[kq]:cumsum(rep(2,q))[kq],z])
        Ran_unif_new_house_q <- runif(nrow(pr_new_E_house_q))
        cumul_new_house_q <- pr_new_E_house_q%*%upper.tri(diag(ncol(pr_new_E_house_q)),diag=TRUE)
        E_house_new[,kq] <- rowSums(Ran_unif_new_house_q > cumul_new_house_q) 
      }
    }
    remove(count_table_house)
    
    ## Sample Epsilon_indiv and synthetic individual data
    if(sum(mc==M_to_use_mc)==1){
      E_indiv_new <- E_indiv
    }
    rep_z <- rep(z,n_i)
    count_table_indiv <- matrix(apply(E_indiv,2,function(x) t(table(factor(rep_z,levels=c(1:TT)),x))),byrow=T,ncol=2)
    beta_sample_indiv <- 
      rbeta((TT*p),matrix(a_epsilon_j + count_table_indiv[,1],nrow=(TT*p)),matrix(b_epsilon_j + count_table_indiv[,2],nrow=(TT*p)))
    for(kp in 1:p){
      beta_sample_indiv_kp <- beta_sample_indiv[1:TT+ TT*(kp-1)]
      epsilon_indiv[d_k_indiv_cum[kp]:cumsum(rep(2,p))[kp],] <- rbind(beta_sample_indiv_kp,1-beta_sample_indiv_kp)
      if(sum(mc==M_to_use_mc)==1){
        pr_new_E_indiv_p <- t(epsilon_indiv[d_k_indiv_cum[kp]:cumsum(rep(2,p))[kp],rep_z])
        Ran_unif_new_indiv_p <- runif(nrow(pr_new_E_indiv_p))
        cumul_new_indiv_p <- pr_new_E_indiv_p%*%upper.tri(diag(ncol(pr_new_E_indiv_p)),diag=TRUE)
        E_indiv_new[,kp] <- rowSums(Ran_unif_new_indiv_p > cumul_new_indiv_p) 
      }
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
        imp_index <- imp_index + 1
        Imputations_indiv[[imp_index]] <- E_indiv_new
        Imputations_house[[imp_index]] <- E_house_new
      }
    }
  }
  ## End MCMC
  ## Return results
  list(ZETA=ZETA, GAMMA=GAMMA, Z_CLUST=Z_CLUST, EPSILON=EPSILON, EPSILON_SUM=EPSILON_SUM,
       Imputations_indiv=Imputations_indiv,Imputations_house=Imputations_house)
}


