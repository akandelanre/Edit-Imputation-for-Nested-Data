

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


for(mc in 1:n_iter){
  proc_t <- proc.time()
  
  
  ## Print iteration number
  cat(paste("Iteration ", mc,"\n", sep = ""))
  
  
  ## Sample structural zeros data
  G_0 <- NULL; M_0 <- NULL
  X_indiv_struc <- NULL; X_house_struc <- NULL
  X_indiv_valid <- NULL; X_house_valid <- NULL
  n_batch_sum <- n_batch_sum + ceiling(n_0*prop_batch)
  n_batch <- ceiling(n_batch_sum/mc) #no. of batches of imputations to sample
  n_0[] <- 0
  for(hh_size in H){
    ss <- which(H==hh_size)
    t_0 <- 0; t_1 <- 0
    n_possibles <- ceiling(length(which(n_i == hh_size))*struc_weight[as.character(hh_size),])
    while(t_1 < n_possibles){
      pr_G_t <- lambda[which(level_house[[1]]==hh_size),]*pii #1 is the location of HHsize in level_house
      G_t <- sample(FF,n_batch[ss],prob=pr_G_t,replace=T)
      rep_G_t <- rep(G_t,each=hh_size)
      pr_M_post_t <- omega[rep_G_t,]
      Ran_unif_M_t <- runif(nrow(pr_M_post_t))
      cumul_M_t <- pr_M_post_t%*%upper.tri(diag(ncol(pr_M_post_t)),diag=TRUE)
      M_t <- rowSums(Ran_unif_M_t>cumul_M_t) + 1L
      X_house_t <- hh_size
      lambda_g_t <- t(lambda[,G_t])
      for(tt in 2:(q-1)){ #q is the HH relate variable, can only take one value
        pr_house_t <- lambda_g_t[,d_k_house_cum[tt]:cumsum(d_k_house)[tt]]
        Ran_unif_t <- runif(nrow(pr_house_t))
        cumul_t <- pr_house_t%*%upper.tri(diag(ncol(pr_house_t)),diag=TRUE)
        level_house_t <- level_house[[tt]]
        X_house_t <- cbind(X_house_t,level_house_t[rowSums(Ran_unif_t > cumul_t) + 1L])    
      }
      X_house_t <- cbind(X_house_t,1) #relate is always 1 for household head
      X_indiv_t <- NULL
      phi_m_g_t <- t(phi[,(rep_G_t+((M_t-1)*FF))])
      for(ttt in 1:p){
        pr_indiv_t <- phi_m_g_t[,d_k_indiv_cum[ttt]:cumsum(d_k_indiv)[ttt]]
        Ran_unif_t <- runif(nrow(pr_indiv_t))
        cumul_t <- pr_indiv_t%*%upper.tri(diag(ncol(pr_indiv_t)),diag=TRUE)
        level_indiv_t <- level_indiv[[ttt]]
        X_indiv_t <- cbind(X_indiv_t,level_indiv_t[rowSums(Ran_unif_t > cumul_t) + 1L])  
      }
      comb_to_check <- X_indiv_t
      comb_to_check <- matrix(t(comb_to_check),byrow=T,nrow=n_batch[ss])
      comb_to_check <- cbind(X_house_t[,(q-p+1):q],comb_to_check) #add the household head before check
      check_counter <- checkSZ(comb_to_check,(hh_size+1)) 
      X_indiv_t <- matrix(t(comb_to_check[,-c(1:p)]),byrow=T,ncol=p)
      
      t_1 <- t_1 + sum(check_counter);
      if(t_1 <= n_possibles){
        t_0 <- t_0 + (n_batch[ss] - sum(check_counter))
        X_indiv_struc <- rbind(X_indiv_struc,X_indiv_t[which(rep(check_counter,each=hh_size)==0),])
        X_house_struc <- rbind(X_house_struc,X_house_t[which(check_counter==0),])
        X_indiv_valid <- rbind(X_indiv_valid,X_indiv_t[which(rep(check_counter,each=hh_size)==1),])
        X_house_valid <- rbind(X_house_valid,X_house_t[which(check_counter==1),])
        M_0 <- c(M_0,M_t[which(rep(check_counter,each=hh_size)==0)])
        G_0 <- c(G_0,G_t[which(check_counter==0)])
      } else {
        t_needed <- sum(check_counter) - (t_1 - n_possibles)
        index_needed <- which(cumsum(check_counter)==t_needed)[1]
        check_counter_needed <- check_counter[1:index_needed]
        t_0 <- t_0 + (length(check_counter_needed) - sum(check_counter_needed))
        X_indiv_struc <- rbind(X_indiv_struc,X_indiv_t[which(rep(check_counter_needed,each=hh_size)==0),])
        X_house_struc <- rbind(X_house_struc,X_house_t[which(check_counter_needed==0),])
        X_indiv_valid <- rbind(X_indiv_valid,X_indiv_t[which(rep(check_counter_needed,each=hh_size)==1),])
        X_house_valid <- rbind(X_house_valid,X_house_t[which(check_counter_needed==1),])
        M_0 <- c(M_0,M_t[which(rep(check_counter_needed,each=hh_size)==0)])
        G_0 <- c(G_0,G_t[which(check_counter_needed==0)])
      }
    }
    n_0[ss] <- n_0[ss] + t_0
  }
  rep_G_0 <- rep(G_0,X_house_struc[,1])
  row.names(X_house_struc) <- NULL; row.names(X_house_valid) <- NULL
  X_house_struc <- as.data.frame(X_house_struc)
  X_house_valid <- as.data.frame(X_house_valid)
  for(ii in 1:q){
    X_house_struc[,ii] <- factor(X_house_struc[,ii],levels=level_house[[ii]])
    X_house_valid[,ii] <- factor(X_house_valid[,ii],levels=level_house[[ii]])
  }
  X_indiv_struc <- as.data.frame(X_indiv_struc)
  X_indiv_valid <- as.data.frame(X_indiv_valid)
  for(iii in 1:p){
    X_indiv_struc[,iii] <- factor(X_indiv_struc[,iii],levels=level_indiv[[iii]])
    X_indiv_valid[,iii] <- factor(X_indiv_valid[,iii],levels=level_indiv[[iii]])
  }
  colnames(X_house_struc) <- colnames(Y_house)
  colnames(X_indiv_struc) <- colnames(Y_indiv)
  colnames(X_house_valid) <- colnames(Y_house)
  colnames(X_indiv_valid) <- colnames(Y_indiv)
  n_i_0 <- as.numeric(as.character(X_house_struc[,1]))
  house_index_0 <- rep(c(1:sum(n_0)),n_i_0)
  n_i_index_0 <- rep(n_i_0,n_i_0)
  
  
  ## Free up some memory
  remove(pr_M_post_t); remove(Ran_unif_M_t); remove(cumul_M_t)
  remove(X_house_t); remove(lambda_g_t); remove(X_indiv_t); remove(phi_m_g_t)
  remove(pr_house_t); remove(pr_indiv_t); remove(Ran_unif_t); remove(cumul_t)
  remove(comb_to_check); remove(check_counter); remove(check_counter_needed)
  
  
  ## Sample X, the true response
  X_house <- Y_house
  X_indiv <- Y_indiv
  n_batch_imp_sum <- n_batch_imp_sum + ceiling(n_0_reject*prop_batch)
  n_batch_imp <- ceiling(n_batch_imp_sum/mc) + 1 #no. of batches of imputations to sample
  n_0_reject[] <- 0
  #First sample erroneous data (struc zeros variables only!!!)
  for(sss in z_i_index_house){
    another_index <- which(house_index==sss)
    n_another_index <- length(another_index) + 1
    NA_house_prop <- X_house[sss,]
    #NA_house_prop[is.na(NA_house[sss,])] <- NA
    NA_house_prop <- apply(NA_house_prop,2,function(x) as.numeric(as.character(x)))
    NA_house_prop <- matrix(rep(t(NA_house_prop),n_batch_imp[sss]),byrow=TRUE,ncol=q)
    NA_indiv_prop <- X_indiv[another_index,]
    #NA_indiv_prop[is.na(NA_indiv[another_index,])] <- NA
    NA_indiv_prop <- apply(NA_indiv_prop,2,function(x) as.numeric(as.character(x)))
    NA_indiv_prop <- matrix(rep(t(NA_indiv_prop),n_batch_imp[sss]),byrow=TRUE,ncol=p)
    G_prop <- rep(M[sss],n_batch_imp[sss])
    rep_G_prop <- rep(rep_G[another_index],n_batch_imp[sss])
    M_prop <- rep(M[another_index],n_batch_imp[sss])
    lambda_g <- t(lambda[,G_prop])
    phi_m_g <- t(phi[,(rep_G_prop+((M_prop-1)*FF))])
    check_counter_sss <- 0;
    NA_house_sss <- NA_house[sss,]
    NA_house_sss <- matrix(rep(t(NA_house_sss),n_batch_imp[sss]),byrow=TRUE,ncol=q)
    NA_indiv_sss <- NA_indiv[another_index,]
    NA_indiv_sss <- matrix(rep(t(NA_indiv_sss),n_batch_imp[sss]),byrow=TRUE,ncol=p)
    while(check_counter_sss < 1){
      #First sample household data
      X_house_prop <- NA_house_prop
      for(kkk in struc_zero_variables_house){
        level_house_k <- level_house[[kkk]]
        q_house_k <- 1/(d_k_house[kkk] - 1)
        if(length(which(is.na(NA_house_sss[,kkk])==TRUE))>0){
          corr_factor_house_k <- matrix(0,ncol=d_k_house[kkk],nrow=length(G_prop))
          corr_factor_Y_house_k <- cbind(corr_factor_house_k,Y_house_nf[sss,kkk])
          corr_factor_house_k <- t(apply(corr_factor_Y_house_k,1,function(x)
            replace(x,(1+x[(d_k_house[kkk]+1)]-min(level_house_k)),1-epsilon_house[which(struc_zero_variables_house==kkk)])))
          corr_factor_house_k <- corr_factor_house_k[,-ncol(corr_factor_house_k)]
          corr_factor_house_k[corr_factor_house_k==0] <- epsilon_house[which(struc_zero_variables_house==kkk)]*q_house_k
          pr_X_house_k <- lambda_g[,d_k_house_cum[kkk]:cumsum(d_k_house)[kkk]]*corr_factor_house_k
          pr_X_house_k <- t(apply(pr_X_house_k,1,function(x) x/sum(x)))
          pr_X_house_k <- pr_X_house_k[which(is.na(NA_house_sss[,kkk])==TRUE),]
          Ran_unif_X_house_k <- runif(nrow(pr_X_house_k))
          cumul_X_house_k <- pr_X_house_k%*%upper.tri(diag(ncol(pr_X_house_k)),diag=TRUE)
          X_house_prop[is.na(NA_house_sss[,kkk]),kkk] <- level_house_k[rowSums(Ran_unif_X_house_k>cumul_X_house_k) + 1L]   
        } else if(length(which(NA_house_sss[,kkk]==1))>0){
          corr_factor_house_k <- matrix(1,ncol=d_k_house[kkk],nrow=length(G_prop))
          corr_factor_Y_house_k <- cbind(corr_factor_house_k,Y_house_nf[sss,kkk])
          corr_factor_house_k <- t(apply(corr_factor_Y_house_k,1,function(x) replace(x,(1+x[(d_k_house[kkk]+1)]-min(level_house_k)),0)))
          corr_factor_house_k <- corr_factor_house_k[,-ncol(corr_factor_house_k)]
          corr_factor_house_k[corr_factor_house_k==1] <- epsilon_house[which(struc_zero_variables_house==kkk)]*q_house_k
          pr_X_house_k <- lambda_g[,d_k_house_cum[kkk]:cumsum(d_k_house)[kkk]]*corr_factor_house_k
          pr_X_house_k <- t(apply(pr_X_house_k,1,function(x) x/sum(x)))
          pr_X_house_k <- pr_X_house_k[which(NA_house_sss[,kkk]==1),]
          Ran_unif_X_house_k <- runif(nrow(pr_X_house_k))
          cumul_X_house_k <- pr_X_house_k%*%upper.tri(diag(ncol(pr_X_house_k)),diag=TRUE)
          X_house_prop[which(NA_house_sss[,kkk]==1),kkk] <- level_house_k[rowSums(Ran_unif_X_house_k>cumul_X_house_k) + 1L]   
        }
      }
      #Then sample individual data
      X_indiv_prop <- NA_indiv_prop
      for(kkkk in struc_zero_variables_indiv){
        level_indiv_k <- level_indiv[[kkkk]]
        q_indiv_k <- 1/(d_k_indiv[kkkk] - 1)
        if(length(which(is.na(NA_indiv_sss[,kkkk])==TRUE))>0){
          corr_factor_indiv_k <- matrix(0,ncol=d_k_indiv[kkkk],nrow=length(M_prop))
          corr_factor_Y_indiv_k <- cbind(corr_factor_indiv_k,Y_indiv_nf[another_index,kkkk])
          corr_factor_indiv_k <- t(apply(corr_factor_Y_indiv_k,1,function(x) 
            replace(x,(1+x[(d_k_indiv[kkkk]+1)]-min(level_indiv_k)),1-epsilon_indiv[which(struc_zero_variables_indiv==kkkk)])))
          corr_factor_indiv_k <- corr_factor_indiv_k[,-ncol(corr_factor_indiv_k)]
          corr_factor_indiv_k[corr_factor_indiv_k==0] <- epsilon_indiv[which(struc_zero_variables_indiv==kkkk)]*q_indiv_k
          pr_X_indiv_k <- phi_m_g[,d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]*corr_factor_indiv_k
          pr_X_indiv_k <- t(apply(pr_X_indiv_k,1,function(x) x/sum(x)))
          pr_X_indiv_k <- pr_X_indiv_k[which(is.na(NA_indiv_sss[,kkkk])==TRUE),]
          Ran_unif_X_indiv_k <- runif(nrow(pr_X_indiv_k))
          cumul_X_indiv_k <- pr_X_indiv_k%*%upper.tri(diag(ncol(pr_X_indiv_k)),diag=TRUE)
          X_indiv_prop[is.na(NA_indiv_sss[,kkkk]),kkkk] <- level_indiv_k[rowSums(Ran_unif_X_indiv_k>cumul_X_indiv_k) + 1L]
        } else if(length(which(NA_indiv_sss[,kkkk]==1))>0){
          corr_factor_indiv_k <- matrix(1,ncol=d_k_indiv[kkkk],nrow=length(M_prop))
          corr_factor_Y_indiv_k <- cbind(corr_factor_indiv_k,Y_indiv_nf[another_index,kkkk])
          corr_factor_indiv_k <- t(apply(corr_factor_Y_indiv_k,1,function(x) replace(x,(1+x[(d_k_indiv[kkkk]+1)]-min(level_indiv_k)),0)))
          corr_factor_indiv_k <- corr_factor_indiv_k[,-ncol(corr_factor_indiv_k)]
          corr_factor_indiv_k[corr_factor_indiv_k==1] <- epsilon_indiv[which(struc_zero_variables_indiv==kkkk)]*q_indiv_k
          pr_X_indiv_k <- phi_m_g[,d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]*corr_factor_indiv_k
          pr_X_indiv_k <- t(apply(pr_X_indiv_k,1,function(x) x/sum(x)))
          pr_X_indiv_k <- pr_X_indiv_k[which(NA_indiv_sss[,kkkk]==1),]
          Ran_unif_X_indiv_k <- runif(nrow(pr_X_indiv_k))
          cumul_X_indiv_k <- pr_X_indiv_k%*%upper.tri(diag(ncol(pr_X_indiv_k)),diag=TRUE)
          X_indiv_prop[which(NA_indiv_sss[,kkkk]==1),kkkk] <- level_indiv_k[rowSums(Ran_unif_X_indiv_k>cumul_X_indiv_k) + 1L]
        }
      }
      #Check edit rules
      comb_to_check <- matrix(t(X_indiv_prop),nrow=n_batch_imp[sss],byrow=TRUE)
      comb_to_check <- cbind(X_house_prop[,(q-p+1):q],comb_to_check)
      check_counter <- checkSZ(comb_to_check,n_another_index)
      check_counter_sss <- check_counter_sss + sum(check_counter)
      if(length(which(check_counter==1))>0){
        n_0_reject[sss] <- n_0_reject[sss] + length(which(check_counter[1:which(check_counter==1)[1]]==0))
      } else{
        n_0_reject[sss] <- n_0_reject[sss] + n_batch_imp[sss]
      }
    }
    X_house[sss,] <- X_house_prop[which(check_counter==1)[1],]
    X_indiv[another_index,] <- matrix(comb_to_check[which(check_counter==1)[1],-c(1:p)],byrow=TRUE,ncol=p) #remove household head
  }
  #Now sample missing data (non-struc zeros variables only!!!)
  #Household data
  if(sum(is.na(NA_house[,nonstruc_zero_variables_house])) > 0){
    lambda_g <- t(lambda[,G])
    for(kkk in nonstruc_zero_variables_house){
      if(length(which(is.na(NA_house[,kkk])==TRUE))>0){
        pr_X_miss_q <- lambda_g[which(is.na(NA_house[,kkk])==TRUE),d_k_house_cum[kkk]:cumsum(d_k_house)[kkk]]
        Ran_unif_miss_q <- runif(nrow(pr_X_miss_q))
        cumul_miss_q <- pr_X_miss_q%*%upper.tri(diag(ncol(pr_X_miss_q)),diag=TRUE)
        level_house_q <- level_house[[kkk]]
        X_house[is.na(NA_house[,kkk]),kkk] <- level_house_q[rowSums(Ran_unif_miss_q>cumul_miss_q) + 1L]    
      }
    }
  }
  #Individual data
  if(sum(is.na(NA_indiv[,nonstruc_zero_variables_indiv])) > 0){
    phi_m_g <- t(phi[,(rep_G+((M-1)*FF))])
    for(kkkk in nonstruc_zero_variables_indiv){
      if(length(which(is.na(NA_indiv[,kkkk])==TRUE))>0){
        pr_X_miss_p <- phi_m_g[which(is.na(NA_indiv[,kkkk])==TRUE),d_k_indiv_cum[kkkk]:cumsum(d_k_indiv)[kkkk]]
        Ran_unif_miss_p <- runif(nrow(pr_X_miss_p))
        cumul_miss_p <- pr_X_miss_p%*%upper.tri(diag(ncol(pr_X_miss_p)),diag=TRUE)
        level_indiv_p <- level_indiv[[kkkk]]
        X_indiv[is.na(NA_indiv[,kkkk]),kkkk] <- level_indiv_p[rowSums(Ran_unif_miss_p>cumul_miss_p) + 1L]
      }
    }
  }
  for(kq in 1:q){
    X_house[,kq] <- factor(X_house[,kq],levels=level_house[[kq]]) }
  for(kp in 1:p){
    X_indiv[,kp] <- factor(X_indiv[,kp],levels=level_indiv[[kp]]) }
  
  
  ## Sample E, the error indicators
  E_house <- data.matrix(X_house)- data.matrix(Y_house)
  E_house[E_house!=0] <- 1
  E_indiv <- data.matrix(X_indiv)- data.matrix(Y_indiv)
  E_indiv[E_indiv!=0] <- 1
  
  
  ## Sample epsilon
  a_epsilon_house_star <- a_epsilon_house + colSums(E_house[z_i_index_house,struc_zero_variables_house])
  b_epsilon_house_star <- b_epsilon_house + length(z_i_index_house) - a_epsilon_house_star
  a_epsilon_indiv_star <- a_epsilon_indiv + colSums(E_indiv[z_i_index_indiv,struc_zero_variables_indiv])
  b_epsilon_indiv_star <- b_epsilon_indiv + length(z_i_index_indiv) - a_epsilon_indiv_star
  epsilon_house <- rbeta(length(struc_zero_variables_house),t(t(a_epsilon_house_star)),t(t(b_epsilon_house_star)))
  epsilon_indiv <- rbeta(length(struc_zero_variables_indiv),t(t(a_epsilon_indiv_star)),t(t(b_epsilon_indiv_star)))
  
  
  ## Sample G
  ## First create indexes for phi and lambda; data.matrix function won't mess anything up as long as columns are coded as factors
  phi_index <- data.matrix(X_indiv)+FFF_indiv #has to be within loop for MI/EI
  lambda_index <- data.matrix(X_house)+FFF_house #has to be within loop  for MI/EI
  pr_G_post <- prGpost(phi_index,lambda_index,phi,lambda,omega,c(pii),FF,SS,n_i)
  Ran_unif_G <- runif(nrow(pr_G_post))
  cumul_G <- pr_G_post%*%upper.tri(diag(ncol(pr_G_post)),diag=TRUE)
  G <- rowSums(Ran_unif_G>cumul_G) + 1L
  remove(pr_G_post); remove(Ran_unif_G); remove(cumul_G); remove(lambda_index)
  
  
  ## Sample M
  rep_G <- rep(G,n_i)
  pr_M_post <- prMpost(phi_index,phi,omega,rep_G,FF,SS)
  Ran_unif_M <- runif(nrow(pr_M_post))
  cumul_M <- pr_M_post%*%upper.tri(diag(ncol(pr_M_post)),diag=TRUE)
  M <- rowSums(Ran_unif_M>cumul_M) + 1L
  remove(pr_M_post); remove(Ran_unif_M); remove(cumul_M); remove(phi_index)
  
  
  ## Sample phi
  rep_G_all <- c(rep_G,rep_G_0)
  M_all <- c(M,M_0)
  for(gg in 1:SS){
    for(ggg in 1:p){
      phi_count_table <- table(factor(rep_G[which(M==gg)],levels=c(1:FF)),X_indiv[which(M==gg),ggg])
      for(w_i in 1:length(struc_weight)){
        hh_size <- as.numeric(rownames(struc_weight)[w_i])
        w_i_index <- n_i_index_0==hh_size
        rep_G_0_w_i <- rep_G_0[w_i_index]
        M_0_w_i <- M_0[w_i_index]
        X_indiv_struc_w_i <- X_indiv_struc[w_i_index,]
        phi_count_table <- phi_count_table + 
          (table(factor(rep_G_0_w_i[which(M_0_w_i==gg)],levels=c(1:FF)),
                 X_indiv_struc_w_i[which(M_0_w_i==gg),ggg])/struc_weight[w_i])
      }
      phi[d_k_indiv_cum[ggg]:cumsum(d_k_indiv)[ggg],(c(1:FF)+((gg-1)*FF))] <- 
        t(rdirichlet(FF,matrix(a_kdk + phi_count_table,nrow=FF)))
    }
  }
  
  
  ## Sample lambda
  G_all <- c(G,G_0)
  for(kk in 1:q){
    lambda_count_table <- table(factor(G,levels=c(1:FF)),X_house[,kk])
    for(w_i in 1:length(struc_weight)){
      hh_size <- as.numeric(rownames(struc_weight)[w_i])
      w_i_index <- n_i_0==hh_size
      G_0_w_i <- G_0[w_i_index]
      X_house_struc_w_i <- X_house_struc[w_i_index,]
      lambda_count_table <- 
        lambda_count_table + (table(factor(G_0_w_i,levels=c(1:FF)),X_house_struc_w_i[,kk])/struc_weight[w_i])
    }
    lambda[d_k_house_cum[kk]:cumsum(d_k_house)[kk],] <- t(rdirichlet(FF,matrix(a_kdk + lambda_count_table,nrow=FF)))
  }
  
  
  ## Sample U and pii
  n_f <- matrix(summary(factor(G,levels=c(1:FF))),ncol=1)
  for(w_i in 1:length(struc_weight)){
    hh_size <- as.numeric(rownames(struc_weight)[w_i])
    w_i_index <- n_i_0==hh_size
    G_0_w_i <- G_0[w_i_index]
    n_f <- n_f + (matrix(summary(factor(G_0_w_i,levels=c(1:FF))),ncol=1)/struc_weight[w_i])
  }
  U[FF]<-1
  U[1:(FF-1),1] <- rbeta((FF-1),(1L+n_f[1:(FF-1)]),(alpha+(sum(n_f)-cumsum(n_f[-FF]))))
  if(length(which(U[-FF]==1))>0){
    U[which(U[-FF]==1)] <- 0.99999
  }
  one_min_U <- 1L-U
  one_min_U_prod <- c(1,cumprod(one_min_U[1:(FF-1)]))
  pii <- U*one_min_U_prod
  remove(n_f);
  
  
  ## Sample V and omega
  M_G <- table(factor(rep_G,levels=c(1:FF)),factor(M,levels=c(1:SS)))
  for(w_i in 1:length(struc_weight)){
    hh_size <- as.numeric(rownames(struc_weight)[w_i])
    w_i_index <- n_i_index_0==hh_size
    rep_G_0_w_i <- rep_G_0[w_i_index]
    M_0_w_i <- M_0[w_i_index]
    M_G <- M_G + (table(factor(rep_G_0_w_i,levels=c(1:FF)),factor(M_0_w_i,levels=c(1:SS)))/struc_weight[w_i])
  }
  n_gm <- as.data.frame(M_G)$Freq
  V[,SS]<-1
  no_V_to_sim <- (FF*(SS-1))
  V[,1:(SS-1)] <- rbeta(no_V_to_sim,(1L+n_gm[1:no_V_to_sim]),
                        (beta + c( matrix(rowSums(M_G),ncol=SS-1,nrow=FF)-
                                     t(apply(M_G,1,cumsum))[,1:SS-1])))
  if(length(which(V[,-SS]==1))>0){
    V[which(V[,-SS]==1)] <- 0.99999
  }
  one_min_V <- 1L-V
  one_min_V_prod <- cbind(1,t(apply(one_min_V[,-SS],1,cumprod)))
  omega <- V*one_min_V_prod
  remove(M_G); remove(n_gm)
  
  
  ## Sample alpha
  alpha <- rgamma(1,shape=(a_alpha+FF-1),rate=(b_alpha-log(pii[FF])))
  
  
  ## Sample beta
  beta <- rgamma(1,shape=(a_beta+(FF*(SS-1))),rate=(b_beta-sum(log(omega[,SS]))))
  
  
  ## Check number of occupied clusters
  S_occup <- NULL
  for(occ in sort(unique(G_all))){
    S_occup <- rbind(S_occup,dim(table(rep_G_all[which(rep_G_all==occ)],M_all[which(rep_G_all==occ)]))[2])
  }
  cat(paste("Number of Occupied Household Classes is ", length(unique(G_all)), "\n", sep = ''))
  cat(paste("Max Number of Occupied Individual Classes is ", max(S_occup), "\n", sep = ''))
  cat(paste("Number of Sampled Augmented Households is ", sum(n_0), "\n", sep = ''))
  remove(rep_G_all); remove(M_all); remove(G_all)
  
  
  ## Save posterior sample and imputations
  if(mc > burn_in){
    if((mc %% mc_thin)==0){
      #PII <- rbind(PII,c(pii))
      ALPHA <- rbind(ALPHA,alpha)
      G_CLUST <- rbind(G_CLUST,length(unique(G)))
      M_CLUST <- rbind(M_CLUST,max(S_occup))
      BETA <- rbind(BETA,beta)
      #LAMBDA[(mc-burn_in),] <- c(lambda)
      #OMEGA[(mc-burn_in),] <- c(omega)
      N_ZERO <- rbind(N_ZERO,sum(n_0))
      EPSILON_INDIV <- rbind(EPSILON_INDIV, epsilon_indiv)
      EPSILON_HOUSE <- rbind(EPSILON_HOUSE, epsilon_house)
    }
    if(sum(mc==M_to_use_mc)==1){
      dp_imput_indiv <- rbind(dp_imput_indiv,X_indiv)  
      dp_imput_house <- rbind(dp_imput_house,X_house)
    }
  }
  
  
  ## Print some summaries
  cat(paste("Number of Sampled Rejections for Missing Data is ", sum(n_0_reject), "\n", sep = ''))
  cat(paste("Total (True) Number of Sampled Augmented Households is ",
            (sum(n_0_reject)+sum(n_0/struc_weight)),"\n", sep = ''))
  cat(paste("Epsilon_indiv:", round(epsilon_indiv,3), "\n", sep = ''))
  cat(paste("Epsilon_house:", round(epsilon_house,3), "\n", sep = ''))
  elapsed_time <- (proc.time() - proc_t)[["elapsed"]]
  cat(paste("Elapsed Time = ", elapsed_time, "\n\n", sep = ' '))
  
  
  ## Randomly pick one summary to monitor and plot it
  conv_check <- rbind(conv_check,t(lambda%*%pii))
  if(nrow(conv_check) > 1){
    plot(mcmc(conv_check[,sample(ncol(conv_check),1,replace=F)]),
         col="blue")
    #plot(1:length(conv_check[,sample(ncol(conv_check),1,replace=F)]),
    #     conv_check[,sample(ncol(conv_check),1,replace=F)],ylab="",xlab="Interations",
    #     col=rainbow(length(conv_check[,sample(ncol(conv_check),1,replace=F)])),type="b")
  }
  
}



