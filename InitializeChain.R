

###########################################################################
###########################################################################
################## Edit-Imputation Using the NDPMPM Model: ################
########################## Initializing Chain  ############################
###########################################################################
###########################################################################


###########################################################################
###########################################################################
################################## START ##################################
###########################################################################
###########################################################################


###### 1: Set global parameters for data
N <- nrow(Y_indiv)
n <- nrow(Y_house)
n_i <- as.numeric(as.character(Y_house[,1]))
p <- ncol(Y_indiv)
q <- ncol(Y_house)
house_index <- rep(c(1:n),n_i)
n_i_index <- rep(n_i,n_i)
struc_zero_variables_house <- c(1,4) + (q-p) ##gender is still included because I am still using 2012 data
struc_zero_variables_indiv <- c(1,4,5) ##gender is still included because I am still using 2012 data
nonstruc_zero_variables_indiv <- c(1:ncol(Y_indiv))[-struc_zero_variables_indiv]
nonstruc_zero_variables_house <- c(1:ncol(Y_house))[-struc_zero_variables_house]
H <- sort(unique(n_i))
NA_miss_indiv <- Y_indiv; NA_miss_house <- Y_house; # For missing data


###### 2: Fill missing values with starting values
if(sum(is.na(NA_miss_indiv)) > 0){
  for (ii in 1:ncol(Y_indiv)){
    Y_indiv[is.na(Y_indiv[,ii]),ii] <- 
      sample(level_indiv[[ii]],length(Y_indiv[is.na(Y_indiv[,ii]),ii]),replace=T,prob=summary(na.omit(Y_indiv[,ii])))
  }
}
if(sum(is.na(NA_miss_house)) > 0){
  for (jj in 2:ncol(Y_house)){
    Y_house[is.na(Y_house[,jj]),jj] <- 
      sample(level_house[[jj]],length(Y_house[is.na(Y_house[,jj]),jj]),replace=T,prob=summary(na.omit(Y_house[,jj])))
  }
}


###### 3: Calculate observed proportions and number of categories for each variable
d_k_house <- d_k_indiv <- ini_marg_house <- ini_marg_indiv <- NULL
for(k in 1:q){
  d_k_house <- cbind(d_k_house,nlevels(Y_house[,k]))
  ini_marg_k <- as.data.frame(table(Y_house[,k]))$Freq/sum(table(Y_house[,k]))
  ini_marg_k <- matrix(ini_marg_k,ncol=1)
  ini_marg_house <- rbind(ini_marg_house,ini_marg_k)  }
for(k in 1:p){
  d_k_indiv <- cbind(d_k_indiv,nlevels(Y_indiv[,k]))
  ini_marg_k <- as.data.frame(table(Y_indiv[,k]))$Freq/sum(table(Y_indiv[,k]))
  ini_marg_k <- matrix(ini_marg_k,ncol=1)
  ini_marg_indiv <- rbind(ini_marg_indiv,ini_marg_k)  }
d_k_indiv_cum <- 1+cumsum(c(0,d_k_indiv[,-p]))
d_k_house_cum <- 1+cumsum(c(0,d_k_house[,-q]))
FFF_indiv <- matrix(rep(cumsum(c(0,d_k_indiv[,-p])),each=N),ncol=p)
FFF_house <- matrix(rep(cumsum(c(0,d_k_house[,-q])),each=n),ncol=q)


###### 4: Set parameters for structural zeros
n_batch_sum <- rep(1000,length(H)) #sample impossibles in batches before checking constraints
n_0 <- rep(0,length(H))
n_batch_imp_sum <- rep(10,n)
n_0_reject <- rep(0,n)
prop_batch <- 1.2


###### 5: Weighting
weight_option <- FALSE #set to true for weighting/capping option
if(weight_option){
  struc_weight <- c(1/2,1/2,1/3) #set weights: must be ordered & no household size must be excluded
} else {
  struc_weight <- rep(1,length(level_house[[1]])) #set weights: must be ordered & no household size must be excluded
}
struc_weight <- as.matrix(struc_weight)
rownames(struc_weight) <- as.character(unique(sort(n_i)))


###### 6: Check for erronous households...
z_i <- matrix(0,ncol=1,nrow=n)
NA_house <- matrix(0,ncol=q,nrow=n,byrow=T)
NA_indiv <- matrix(0,ncol=p,nrow=N,byrow=T)
for(hh_size in H){
  house_index_hh <- which(n_i == hh_size)
  comb_to_check <- Y_indiv[which(is.element(house_index,house_index_hh)==TRUE),]
  comb_to_check <- apply(comb_to_check,2,function(x) as.numeric(as.character(x)))
  comb_to_check <- matrix(t(comb_to_check),byrow=T,ncol=(p*hh_size))
  comb_to_check <- cbind(Y_house[house_index_hh,(q-p+1):q],comb_to_check) #add the household head before check
  comb_to_check <- apply(comb_to_check,2,function(x) as.numeric(as.character(x)))
  z_i[house_index_hh] <- ifelse(checkSZ(comb_to_check,(hh_size+1))==1,0,1)
  NA_house[house_index_hh,struc_zero_variables_house] <- rep(z_i[house_index_hh],length(struc_zero_variables_house))
  NA_indiv[which(is.element(house_index,house_index_hh)==TRUE),struc_zero_variables_indiv] <- 
    rep(rep(z_i[house_index_hh],n_i[house_index_hh]),length(struc_zero_variables_indiv))
}
z_i_index_house <- which(z_i == 1)
z_i_index_indiv <- which(is.element(house_index,z_i_index_house)==TRUE)
#NA_house <- read.table("Results/E_house_truth.txt",header=TRUE)
#NA_indiv <- read.table("Results/E_indiv_truth.txt",header=TRUE)
NA_house[NA_house==1] <- NA
NA_indiv[NA_indiv==1] <- NA
NA_house[is.na(NA_miss_house)] <- NA #For missing data
NA_indiv[is.na(NA_miss_indiv)] <- NA #For missing data


###### 7: Initialize chain
FF <- 20
SS <- 10
alpha <- beta <- 1
a_kdk <- 1
a_alpha <- b_alpha <- a_beta <- b_beta <- 0.25
lambda <- matrix(rep(ini_marg_house,FF),ncol=FF)
phi <- matrix(0,nrow=length(ini_marg_indiv),ncol=FF*SS) #make phi matrix and not array for c++
for(gm in 1:(FF*SS)){
  phi[,gm] <- ini_marg_indiv   }
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
a_epsilon_house <- b_epsilon_house <- rep(1,length(struc_zero_variables_house))
a_epsilon_indiv <- b_epsilon_indiv <- rep(1,length(struc_zero_variables_indiv))
epsilon_house <- rbeta(length(struc_zero_variables_house),t(t(a_epsilon_house)),t(t(b_epsilon_house)))
epsilon_indiv <- rbeta(length(struc_zero_variables_indiv),t(t(a_epsilon_indiv)),t(t(b_epsilon_indiv)))
pr_G <- matrix(pii,byrow=T,ncol=FF,nrow=n)
Ran_unif_G <- runif(nrow(pr_G))
cumul_G <- pr_G%*%upper.tri(diag(ncol(pr_G)),diag=TRUE)
G <- rowSums(Ran_unif_G>cumul_G) + 1L
rep_G <- rep(G,n_i)
pr_M <- omega[rep_G,]
Ran_unif_M <- runif(nrow(pr_M))
cumul_M <- pr_M%*%upper.tri(diag(ncol(pr_M)),diag=TRUE)
M <- rowSums(Ran_unif_M>cumul_M) + 1L


###### 8: Free up some memory
remove(pr_G); remove(Ran_unif_G); remove(cumul_G);
remove(pr_M); remove(Ran_unif_M); remove(cumul_M);


###### 9: Set MCMC parameters
n_iter <- 10000
burn_in <- 0.5*n_iter
MM <- 50
mc_thin <- 10
#M_to_use_mc <- sort(sample(seq((burn_in +1),n_iter,by=mc_thin),MM,replace=F))
M_to_use_mc <- round(seq((burn_in +1),n_iter,length.out=MM))


###### 10: Create empty matrices to save results
dp_imput_house <- dp_imput_indiv <- NULL
ALPHA <- BETA <- PII <- G_CLUST <- M_CLUST <- N_ZERO <- NULL
EPSILON_INDIV <- EPSILON_HOUSE <- NULL
#LAMBDA <- matrix(0,ncol=(ncol(lambda)*nrow(lambda)),nrow=(n_iter-burn_in))
#OMEGA <- matrix(0,ncol=(ncol(omega)*nrow(omega)),nrow=(n_iter-burn_in))
conv_check <- NULL






