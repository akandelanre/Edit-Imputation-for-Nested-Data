

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
struc_zero_variables <- c(1,4,5)
nonstruc_zero_variables <- c(2,3)
H <- sort(unique(n_i))

###### 2: Calculate observed proportions and number of categories for each variable
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


###### 3: Set parameters for structural zeros
n_batch_init <- 1000 #sample impossibles in batches before checking constraints
n_0 <- rep(0,length(level_house[[1]]))
n_batch_imp_init <- 50 #sample imputations in batches before checking constraints
n_0_reject <- rep(0,n)
prop_batch <- 1.5


###### 4: Weighting
weight_option <- FALSE #set to true for weighting/capping option
if(weight_option){
  struc_weight <- c(1/2,1/2,1/3) #set weights: must be ordered & no household size must be excluded
} else {
  struc_weight <- rep(1,length(level_house[[1]])) #set weights: must be ordered & no household size must be excluded
}
struc_weight <- as.matrix(struc_weight)
rownames(struc_weight) <- as.character(unique(sort(n_i)))


###### 5: Initialize chain
FF <- 30
SS <- 15
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
#a_epsilon_house <- c(0,1,1,1,1,1,0); b_epsilon_house <- rep(1,q)
#a_epsilon_indiv <- c(1,1,1,1,1); b_epsilon_indiv <- rep(1,p)
#epsilon_house <- rbeta(q,t(t(a_epsilon_house)),t(t(b_epsilon_house)))
#epsilon_indiv <- rbeta(p,t(t(a_epsilon_indiv)),t(t(b_epsilon_indiv)))
epsilon_indiv <- c(0.1,0.1,0.05,0.45,0.3)
epsilon_house <- c(0.0,0.2,0.1,0.05,0.35,0.5,0.0)
#E_house <- matrix(rbinom((n*q),1,epsilon_house),ncol=q,byrow=T)
#E_indiv <- matrix(rbinom((N*p),1,epsilon_indiv),ncol=p,byrow=T)
pr_G <- matrix(pii,byrow=T,ncol=FF,nrow=n)
Ran_unif_G <- runif(nrow(pr_G))
cumul_G <- pr_G%*%upper.tri(diag(ncol(pr_G)),diag=TRUE)
G <- rowSums(Ran_unif_G>cumul_G) + 1L
rep_G <- rep(G,n_i)
pr_M <- omega[rep_G,]
Ran_unif_M <- runif(nrow(pr_M))
cumul_M <- pr_M%*%upper.tri(diag(ncol(pr_M)),diag=TRUE)
M <- rowSums(Ran_unif_M>cumul_M) + 1L


###### 7: Free up some memory
remove(pr_G); remove(Ran_unif_G); remove(cumul_G);
remove(pr_M); remove(Ran_unif_M); remove(cumul_M);


###### 8: Set MCMC parameters
n_iter <- 10000
burn_in <- 0.5*n_iter
MM <- 5
mc_thin <- 10
#M_to_use_mc <- sort(sample(seq((burn_in +1),n_iter,by=mc_thin),MM,replace=F))
M_to_use_mc <- round(seq((burn_in +1),n_iter,length.out=MM))

###### 7: Create empty matrices to save results
dp_imput_house <- dp_imput_indiv <- NULL
ALPHA <- BETA <- PII <- G_CLUST <- M_CLUST <- N_ZERO <- NULL
EPSILON_INDIV <- EPSILON_HOUSE <- NULL
#LAMBDA <- matrix(0,ncol=(ncol(lambda)*nrow(lambda)),nrow=(n_iter-burn_in))
#OMEGA <- matrix(0,ncol=(ncol(omega)*nrow(omega)),nrow=(n_iter-burn_in))
conv_check <- NULL






