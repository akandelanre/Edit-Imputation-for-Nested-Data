


###########################################################################
###########################################################################
####### MI using NDPMPM: Use rejection sampler to sample missing data #####
####### Add weighting option for sampling impossibles #####################
####### Add the Rejection sampler/sampling proposals hybrid option ########
####### Also make household head into a household level variable ##########
###########################################################################
###########################################################################

###################################### START ######################################
########################## Step 1: Model Assessment ##########################
rm(list = ls())


###### 1a: Load saved results
Data_house_truth <- read.table("Results/Data_house_truth.txt",header=TRUE)
Data_indiv_truth <- read.table("Results/Data_indiv_truth.txt",header=TRUE)
Data_house_cc <- read.table("Results/Data_house_cc.txt",header=TRUE)
Data_indiv_cc <- read.table("Results/Data_indiv_cc.txt",header=TRUE)
dp_imput_house <- read.table("Results/dp_imput_house.txt",header=TRUE)
dp_imput_indiv <- read.table("Results/dp_imput_indiv.txt",header=TRUE)
dp_imput_house_nz <- read.table("Results/dp_imput_house_nz.txt",header=TRUE)
dp_imput_indiv_nz <- read.table("Results/dp_imput_indiv_nz.txt",header=TRUE)
dp_imput_house_hybrid <- read.table("Results/dp_imput_house_hybrid.txt",header=TRUE)
dp_imput_indiv_hybrid <- read.table("Results/dp_imput_indiv_hybrid.txt",header=TRUE)
dp_imput_house_weighted_hybrid <- read.table("Results/dp_imput_house_weighted_hybrid.txt",header=TRUE)
dp_imput_indiv_weighted_hybrid <- read.table("Results/dp_imput_indiv_weighted_hybrid.txt",header=TRUE)
dp_imput_house_weighted <- read.table("Results/dp_imput_house_weighted.txt",header=TRUE)
dp_imput_indiv_weighted <- read.table("Results/dp_imput_indiv_weighted.txt",header=TRUE)


###### 1b: Define parameters
mm <- 50;
N <- nrow(Data_indiv_truth)
n <- nrow(Data_house_truth)
n_i <- as.numeric(as.character(Data_house_truth[,1]))
p <- ncol(Data_indiv_truth)
q <- ncol(Data_house_truth)
house_index <- rep(c(1:n),n_i)
N_cc <- nrow(Data_indiv_cc)
n_cc <- nrow(Data_house_cc)
n_i_cc <- as.numeric(as.character(Data_house_cc[,1]))
house_index_cc <- rep(c(1:n_cc),n_i_cc)


###### 1c: Define categories for RELATE, GENDER & RACE variable
HEAD <- 1
SPOUSE <- 2
BIOLOGICALCHILD <- 3
ADOPTEDCHILD <- 4
STEPCHILD <- 5
SIBLING <- 6
PARENT <- 7
GRANDCHILD <- 8
PARENTINLAW <- 9
CHILDINLAW <- 10
WHITE <- 1; BLACK <- 2
MALE <- 1; FEMALE <- 2
NONHISP <- 1; OWN <- 1


###### 2: Calculate probabilities that depend on relationship variable
###### 2a: Original Data
Probs <- matrix(0,nrow=21)
n_row_2 <- length(which(n_i==1))
n_row_3 <- length(which(n_i==2))
n_row_4 <- length(which(n_i==3))
samp_size <- matrix(c(n_row_2,n_row_3,n_row_4,rep(n,18)),21,1)
for(kk in 1:n){
  hh_check <- Data_house_truth[kk,(q-p+1):q]
  colnames(hh_check) <- colnames(Data_indiv_truth)
  hh_check <- rbind(hh_check,Data_indiv_truth[which(house_index==kk),])
  hh_check <- data.frame(Owner=Data_house_truth[kk,"Owner"],hh_check)
  if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
    Probs[1] <- Probs[1] + 1 #All same race, n_i = 2
  }
  if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
     hh_check[2,"Race"]==hh_check[3,"Race"]){
    Probs[2] <- Probs[2] + 1 #All same race, n_i = 3
  }
  if(nrow(hh_check)==4 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
     hh_check[2,"Race"]==hh_check[3,"Race"]&&hh_check[3,"Race"]==hh_check[4,"Race"]){
    Probs[3] <- Probs[3] + 1 #All same race, n_i = 4
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1){
    Probs[4] <- Probs[4] + 1 #Spouse present
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs[5] <- Probs[5] + 1 #Spouse with white HH
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
    Probs[6] <- Probs[6] + 1 #Spouse with black HH
  }
  if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs[7] <- Probs[7] + 1 #White Couple
  }
  if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1 && hh_check[1,"Owner"] == OWN){
    Probs[8] <- Probs[8] + 1 #White Couple, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
    Probs[9] <- Probs[9] + 1 #Same race couple
  }
  if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
     sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1){
    Probs[10] <- Probs[10] + 1 #White-nonwhite couple
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs[11] <- Probs[11] + 1 #Non-white couple, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==PARENT,"Gender"]==FEMALE)==1 && length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
    Probs[12] <- Probs[12] + 1 #Only one parent, mother
  }
  if(length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
    Probs[13] <- Probs[13] + 1 #Only one parent   
  }
  if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 | sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1){
    Probs[14] <- Probs[14] + 1 #Children present  
  }
  if(sum(hh_check[,"Relate"]==SIBLING)>=1){
    Probs[15] <- Probs[15] + 1 #Siblings present  
  }
  if(sum(hh_check[,"Relate"]==GRANDCHILD)>=1){
    Probs[16] <- Probs[16] + 1 #Grandchild present  
  }
  if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0)+ 
     ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 |sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0)+
     ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
    Probs[17] <- Probs[17] + 1 #Three generations present  
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1){
    Probs[18] <- Probs[18] + 1 #Same age couple
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
     sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
    Probs[19] <- Probs[19] + 1 #HH older than spouse, white HH
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]==NONHISP)==1){
    Probs[20] <- Probs[20] + 1 #Non Hisp HH
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
     hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
    Probs[21] <- Probs[21] + 1 #White HH with hisp origin
  }
}
Probs <- Probs/samp_size
V <- (Probs*(1-Probs))/samp_size
CIntLower <- Probs + (qnorm(0.025)*sqrt(V))
CIntUpper <- Probs - (qnorm(0.025)*sqrt(V))
CInt <- cbind(CIntLower,CIntUpper)


###### 2b: Complete-Case Data
Probs_cc <- matrix(0,nrow=21)
n_row_2_cc <- length(which(n_i_cc==1))
n_row_3_cc <- length(which(n_i_cc==2))
n_row_4_cc <- length(which(n_i_cc==3))
samp_size_cc <- matrix(c(n_row_2_cc,n_row_3_cc,n_row_4_cc,rep(n_cc,18)),21,1)
for(kk in 1:n_cc){
  hh_check <- Data_house_cc[kk,(q-p+1):q]
  colnames(hh_check) <- colnames(Data_indiv_cc)
  hh_check <- rbind(hh_check,Data_indiv_cc[which(house_index_cc==kk),])
  hh_check <- data.frame(Owner=Data_house_cc[kk,"Owner"],hh_check)
  if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
    Probs_cc[1] <- Probs_cc[1] + 1 #All same race, n_i = 2
  }
  if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
     hh_check[2,"Race"]==hh_check[3,"Race"]){
    Probs_cc[2] <- Probs_cc[2] + 1 #All same race, n_i = 3
  }
  if(nrow(hh_check)==4 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
     hh_check[2,"Race"]==hh_check[3,"Race"]&&hh_check[3,"Race"]==hh_check[4,"Race"]){
    Probs_cc[3] <- Probs_cc[3] + 1 #All same race, n_i = 4
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1){
    Probs_cc[4] <- Probs_cc[4] + 1 #Spouse present
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs_cc[5] <- Probs_cc[5] + 1 #Spouse with white HH
  }
  if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
    Probs_cc[6] <- Probs_cc[6] + 1 #Spouse with black HH
  }
  if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
    Probs_cc[7] <- Probs_cc[7] + 1 #White Couple
  }
  if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
     sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1 && hh_check[1,"Owner"] == OWN){
    Probs_cc[8] <- Probs_cc[8] + 1 #White Couple, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
    Probs_cc[9] <- Probs_cc[9] + 1 #Same race couple
  }
  if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
     sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1){
    Probs_cc[10] <- Probs_cc[10] + 1 #White-nonwhite couple
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1 &&
     hh_check[1,"Owner"]==OWN){
    Probs_cc[11] <- Probs_cc[11] + 1 #Non-white couple, own
  }
  if(sum(hh_check[hh_check[,"Relate"]==PARENT,"Gender"]==FEMALE)==1 && length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
    Probs_cc[12] <- Probs_cc[12] + 1 #Only one parent, mother
  }
  if(length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
    Probs_cc[13] <- Probs_cc[13] + 1 #Only one parent   
  }
  if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 | sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1){
    Probs_cc[14] <- Probs_cc[14] + 1 #Children present  
  }
  if(sum(hh_check[,"Relate"]==SIBLING)>=1){
    Probs_cc[15] <- Probs_cc[15] + 1 #Siblings present  
  }
  if(sum(hh_check[,"Relate"]==GRANDCHILD)>=1){
    Probs_cc[16] <- Probs_cc[16] + 1 #Grandchild present  
  }
  if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0)+ 
     ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 |sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0)+
     ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
    Probs_cc[17] <- Probs_cc[17] + 1 #Three generations present  
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1){
    Probs_cc[18] <- Probs_cc[18] + 1 #Same age couple
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
     sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
    Probs_cc[19] <- Probs_cc[19] + 1 #HH older than spouse, white HH
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]==NONHISP)==1){
    Probs_cc[20] <- Probs_cc[20] + 1 #Non Hisp HH
  }
  if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
         hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
    Probs_cc[21] <- Probs_cc[21] + 1 #White HH with hisp origin
  }
}
Probs_cc <- Probs_cc/samp_size_cc
V_cc <- (Probs_cc*(1-Probs_cc))/samp_size_cc
CIntLower_cc <- Probs_cc + (qnorm(0.025)*sqrt(V_cc))
CIntUpper_cc <- Probs_cc - (qnorm(0.025)*sqrt(V_cc))
CInt_cc <- cbind(CIntLower_cc,CIntUpper_cc)


###### 2c: Regular Rejection Sampler
Probs_syn <- matrix(0,nrow=21,ncol=mm)
V_syn <- matrix(0,nrow=21,ncol=mm)
for(k in 1:mm){
  k_imp_house <- dp_imput_house[((n*(k-1))+1):(n*k),]
  k_imp_indiv <- dp_imput_indiv[((N*(k-1))+1):(N*k),]
  new_n_i <- as.numeric(as.character(k_imp_house[,"HHSize"]))
  new_house_index <- rep(c(1:n),new_n_i)
  for(kk in 1:n){
    hh_check <- k_imp_house[kk,(q-p+1):q]
    colnames(hh_check) <- colnames(k_imp_indiv)
    hh_check <- rbind(hh_check,k_imp_indiv[which(new_house_index==kk),])
    hh_check <- data.frame(Owner=k_imp_house[kk,"Owner"],hh_check)
    if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
      Probs_syn[1,k] <- Probs_syn[1,k] + 1 #All same race, n_i = 2
    }
    if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]){
      Probs_syn[2,k] <- Probs_syn[2,k] + 1 #All same race, n_i = 3
    }
    if(nrow(hh_check)==4 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]&&hh_check[3,"Race"]==hh_check[4,"Race"]){
      Probs_syn[3,k] <- Probs_syn[3,k] + 1 #All same race, n_i = 4
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1){
      Probs_syn[4,k] <- Probs_syn[4,k] + 1 #Spouse present
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn[5,k] <- Probs_syn[5,k] + 1 #Spouse with white HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs_syn[6,k] <- Probs_syn[6,k] + 1 #Spouse with black HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn[7,k] <- Probs_syn[7,k] + 1 #White Couple
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1 && hh_check[1,"Owner"] == OWN){
      Probs_syn[8,k] <- Probs_syn[8,k] + 1 #White Couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
      Probs_syn[9,k] <- Probs_syn[9,k] + 1 #Same race couple
    }
    if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
       sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1){
      Probs_syn[10,k] <- Probs_syn[10,k] + 1 #White-nonwhite couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn[11,k] <- Probs_syn[11,k] + 1 #Non-white couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==PARENT,"Gender"]==FEMALE)==1 && length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn[12,k] <- Probs_syn[12,k] + 1 #Only one parent, mother
    }
    if(length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn[13,k] <- Probs_syn[13,k] + 1 #Only one parent   
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 | sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1){
      Probs_syn[14,k] <- Probs_syn[14,k] + 1 #Children present  
    }
    if(sum(hh_check[,"Relate"]==SIBLING)>=1){
      Probs_syn[15,k] <- Probs_syn[15,k] + 1 #Siblings present  
    }
    if(sum(hh_check[,"Relate"]==GRANDCHILD)>=1){
      Probs_syn[16,k] <- Probs_syn[16,k] + 1 #Grandchild present  
    }
    if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0)+ 
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 |sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0)+
       ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
      Probs_syn[17,k] <- Probs_syn[17,k] + 1 #Three generations present  
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1){
      Probs_syn[18,k] <- Probs_syn[18,k] + 1 #Same age couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn[19,k] <- Probs_syn[19,k] + 1 #HH older than spouse, white HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]==NONHISP)==1){
      Probs_syn[20,k] <- Probs_syn[20,k] + 1 #Non Hisp HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
           hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn[21,k] <- Probs_syn[21,k] + 1 #White HH with hisp origin
    }
  }
  Probs_syn[,k] <- Probs_syn[,k]/samp_size
  V_syn[,k] <- (Probs_syn[,k]*(1-Probs_syn[,k]))/samp_size
}
dp_qbar_syn <- rowMeans(Probs_syn)
dp_b_syn <- apply(Probs_syn,1,var)
dp_ubar_syn <- rowMeans(V_syn)
dp_t_syn <- dp_ubar_syn + (dp_b_syn*(mm+1)/mm) #dp_t_syn <- dp_ubar_syn + (dp_b/mm) for synthetic data
dp_r_syn <- dp_ubar_syn/dp_b_syn
dp_v_syn <- (mm-1)*((1+((mm/(mm+1))*dp_r_syn))^2) #dp_v_syn <- (mm-1)*((1+((mm)*dp_r_syn))^2) for synthetic data
CIntLower_dp_syn <- dp_qbar_syn + (qt(0.025,dp_v_syn)*sqrt(dp_t_syn))
CIntUpper_dp_syn <- dp_qbar_syn - (qt(0.025,dp_v_syn)*sqrt(dp_t_syn))
CInt_syn <- cbind(CIntLower_dp_syn,CIntUpper_dp_syn)


###### 2d: No structural zeros model
Probs_syn_nz <- matrix(0,nrow=21,ncol=mm)
V_syn_nz <- matrix(0,nrow=21,ncol=mm)
for(k in 1:mm){
  k_imp_house <- dp_imput_house_nz[((n*(k-1))+1):(n*k),]
  k_imp_indiv <- dp_imput_indiv_nz[((N*(k-1))+1):(N*k),]
  new_n_i <- as.numeric(as.character(k_imp_house[,"HHSize"]))
  new_house_index <- rep(c(1:n),new_n_i)
  for(kk in 1:n){
    hh_check <- k_imp_house[kk,(q-p+1):q]
    colnames(hh_check) <- colnames(k_imp_indiv)
    hh_check <- rbind(hh_check,k_imp_indiv[which(new_house_index==kk),])
    hh_check <- data.frame(Owner=k_imp_house[kk,"Owner"],hh_check)
    if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
      Probs_syn_nz[1,k] <- Probs_syn_nz[1,k] + 1 #All same race, n_i = 2
    }
    if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]){
      Probs_syn_nz[2,k] <- Probs_syn_nz[2,k] + 1 #All same race, n_i = 3
    }
    if(nrow(hh_check)==4 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]&&hh_check[3,"Race"]==hh_check[4,"Race"]){
      Probs_syn_nz[3,k] <- Probs_syn_nz[3,k] + 1 #All same race, n_i = 4
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1){
      Probs_syn_nz[4,k] <- Probs_syn_nz[4,k] + 1 #Spouse present
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_nz[5,k] <- Probs_syn_nz[5,k] + 1 #Spouse with white HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs_syn_nz[6,k] <- Probs_syn_nz[6,k] + 1 #Spouse with black HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_nz[7,k] <- Probs_syn_nz[7,k] + 1 #White Couple
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1 && hh_check[1,"Owner"] == OWN){
      Probs_syn_nz[8,k] <- Probs_syn_nz[8,k] + 1 #White Couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
      Probs_syn_nz[9,k] <- Probs_syn_nz[9,k] + 1 #Same race couple
    }
    if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
       sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1){
      Probs_syn_nz[10,k] <- Probs_syn_nz[10,k] + 1 #White-nonwhite couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn_nz[11,k] <- Probs_syn_nz[11,k] + 1 #Non-white couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==PARENT,"Gender"]==FEMALE)==1 && length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn_nz[12,k] <- Probs_syn_nz[12,k] + 1 #Only one parent, mother
    }
    if(length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn_nz[13,k] <- Probs_syn_nz[13,k] + 1 #Only one parent   
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 | sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1){
      Probs_syn_nz[14,k] <- Probs_syn_nz[14,k] + 1 #Children present  
    }
    if(sum(hh_check[,"Relate"]==SIBLING)>=1){
      Probs_syn_nz[15,k] <- Probs_syn_nz[15,k] + 1 #Siblings present  
    }
    if(sum(hh_check[,"Relate"]==GRANDCHILD)>=1){
      Probs_syn_nz[16,k] <- Probs_syn_nz[16,k] + 1 #Grandchild present  
    }
    if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0)+ 
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 |sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0)+
       ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
      Probs_syn_nz[17,k] <- Probs_syn_nz[17,k] + 1 #Three generations present  
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1){
      Probs_syn_nz[18,k] <- Probs_syn_nz[18,k] + 1 #Same age couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_nz[19,k] <- Probs_syn_nz[19,k] + 1 #HH older than spouse, white HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]==NONHISP)==1){
      Probs_syn_nz[20,k] <- Probs_syn_nz[20,k] + 1 #Non Hisp HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
           hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_nz[21,k] <- Probs_syn_nz[21,k] + 1 #White HH with hisp origin
    }
  }
  Probs_syn_nz[,k] <- Probs_syn_nz[,k]/samp_size
  V_syn_nz[,k] <- (Probs_syn_nz[,k]*(1-Probs_syn_nz[,k]))/samp_size
}
dp_qbar_syn_nz <- rowMeans(Probs_syn_nz)
dp_b_syn_nz <- apply(Probs_syn_nz,1,var)
dp_ubar_syn_nz <- rowMeans(V_syn_nz)
dp_t_syn_nz <- dp_ubar_syn_nz + (dp_b_syn_nz*(mm+1)/mm) #dp_t_syn_nz <- dp_ubar_syn_nz + (dp_b/mm) for synthetic data
dp_r_syn_nz <- dp_ubar_syn_nz/dp_b_syn_nz
dp_v_syn_nz <- (mm-1)*((1+((mm/(mm+1))*dp_r_syn_nz))^2) #dp_v_syn_nz <- (mm-1)*((1+((mm)*dp_r_syn_nz))^2) for synthetic data
CIntLower_dp_syn_nz <- dp_qbar_syn_nz + (qt(0.025,dp_v_syn_nz)*sqrt(dp_t_syn_nz))
CIntUpper_dp_syn_nz <- dp_qbar_syn_nz - (qt(0.025,dp_v_syn_nz)*sqrt(dp_t_syn_nz))
CInt_syn_nz <- cbind(CIntLower_dp_syn_nz,CIntUpper_dp_syn_nz)


###### 2e: Hybrid Sampler
Probs_syn_hybrid <- matrix(0,nrow=21,ncol=mm)
V_syn_hybrid <- matrix(0,nrow=21,ncol=mm)
for(k in 1:mm){
  k_imp_house <- dp_imput_house_hybrid[((n*(k-1))+1):(n*k),]
  k_imp_indiv <- dp_imput_indiv_hybrid[((N*(k-1))+1):(N*k),]
  new_n_i <- as.numeric(as.character(k_imp_house[,"HHSize"]))
  new_house_index <- rep(c(1:n),new_n_i)
  for(kk in 1:n){
    hh_check <- k_imp_house[kk,(q-p+1):q]
    colnames(hh_check) <- colnames(k_imp_indiv)
    hh_check <- rbind(hh_check,k_imp_indiv[which(new_house_index==kk),])
    hh_check <- data.frame(Owner=k_imp_house[kk,"Owner"],hh_check)
    if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
      Probs_syn_hybrid[1,k] <- Probs_syn_hybrid[1,k] + 1 #All same race, n_i = 2
    }
    if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]){
      Probs_syn_hybrid[2,k] <- Probs_syn_hybrid[2,k] + 1 #All same race, n_i = 3
    }
    if(nrow(hh_check)==4 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]&&hh_check[3,"Race"]==hh_check[4,"Race"]){
      Probs_syn_hybrid[3,k] <- Probs_syn_hybrid[3,k] + 1 #All same race, n_i = 4
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1){
      Probs_syn_hybrid[4,k] <- Probs_syn_hybrid[4,k] + 1 #Spouse present
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_hybrid[5,k] <- Probs_syn_hybrid[5,k] + 1 #Spouse with white HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs_syn_hybrid[6,k] <- Probs_syn_hybrid[6,k] + 1 #Spouse with black HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_hybrid[7,k] <- Probs_syn_hybrid[7,k] + 1 #White Couple
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1 && hh_check[1,"Owner"] == OWN){
      Probs_syn_hybrid[8,k] <- Probs_syn_hybrid[8,k] + 1 #White Couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
      Probs_syn_hybrid[9,k] <- Probs_syn_hybrid[9,k] + 1 #Same race couple
    }
    if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
       sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1){
      Probs_syn_hybrid[10,k] <- Probs_syn_hybrid[10,k] + 1 #White-nonwhite couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn_hybrid[11,k] <- Probs_syn_hybrid[11,k] + 1 #Non-white couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==PARENT,"Gender"]==FEMALE)==1 && length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn_hybrid[12,k] <- Probs_syn_hybrid[12,k] + 1 #Only one parent, mother
    }
    if(length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn_hybrid[13,k] <- Probs_syn_hybrid[13,k] + 1 #Only one parent   
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 | sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1){
      Probs_syn_hybrid[14,k] <- Probs_syn_hybrid[14,k] + 1 #Children present  
    }
    if(sum(hh_check[,"Relate"]==SIBLING)>=1){
      Probs_syn_hybrid[15,k] <- Probs_syn_hybrid[15,k] + 1 #Siblings present  
    }
    if(sum(hh_check[,"Relate"]==GRANDCHILD)>=1){
      Probs_syn_hybrid[16,k] <- Probs_syn_hybrid[16,k] + 1 #Grandchild present  
    }
    if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0)+ 
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 |sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0)+
       ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
      Probs_syn_hybrid[17,k] <- Probs_syn_hybrid[17,k] + 1 #Three generations present  
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1){
      Probs_syn_hybrid[18,k] <- Probs_syn_hybrid[18,k] + 1 #Same age couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_hybrid[19,k] <- Probs_syn_hybrid[19,k] + 1 #HH older than spouse, white HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]==NONHISP)==1){
      Probs_syn_hybrid[20,k] <- Probs_syn_hybrid[20,k] + 1 #Non Hisp HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
           hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_hybrid[21,k] <- Probs_syn_hybrid[21,k] + 1 #White HH with hisp origin
    }
  }
  Probs_syn_hybrid[,k] <- Probs_syn_hybrid[,k]/samp_size
  V_syn_hybrid[,k] <- (Probs_syn_hybrid[,k]*(1-Probs_syn_hybrid[,k]))/samp_size
}
dp_qbar_syn_hybrid <- rowMeans(Probs_syn_hybrid)
dp_b_syn_hybrid <- apply(Probs_syn_hybrid,1,var)
dp_ubar_syn_hybrid <- rowMeans(V_syn_hybrid)
dp_t_syn_hybrid <- dp_ubar_syn_hybrid + (dp_b_syn_hybrid*(mm+1)/mm) #dp_t_syn_hybrid <- dp_ubar_syn_hybrid + (dp_b/mm) for synthetic data
dp_r_syn_hybrid <- dp_ubar_syn_hybrid/dp_b_syn_hybrid
dp_v_syn_hybrid <- (mm-1)*((1+((mm/(mm+1))*dp_r_syn_hybrid))^2) #dp_v_syn_hybrid <- (mm-1)*((1+((mm)*dp_r_syn_hybrid))^2) for synthetic data
CIntLower_dp_syn_hybrid <- dp_qbar_syn_hybrid + (qt(0.025,dp_v_syn_hybrid)*sqrt(dp_t_syn_hybrid))
CIntUpper_dp_syn_hybrid <- dp_qbar_syn_hybrid - (qt(0.025,dp_v_syn_hybrid)*sqrt(dp_t_syn_hybrid))
CInt_syn_hybrid <- cbind(CIntLower_dp_syn_hybrid,CIntUpper_dp_syn_hybrid)


###### 2f: Weighted Hybrid Sampler
Probs_syn_weighted_hybrid <- matrix(0,nrow=21,ncol=mm)
V_syn_weighted_hybrid <- matrix(0,nrow=21,ncol=mm)
for(k in 1:mm){
  k_imp_house <- dp_imput_house_weighted_hybrid[((n*(k-1))+1):(n*k),]
  k_imp_indiv <- dp_imput_indiv_weighted_hybrid[((N*(k-1))+1):(N*k),]
  new_n_i <- as.numeric(as.character(k_imp_house[,"HHSize"]))
  new_house_index <- rep(c(1:n),new_n_i)
  for(kk in 1:n){
    hh_check <- k_imp_house[kk,(q-p+1):q]
    colnames(hh_check) <- colnames(k_imp_indiv)
    hh_check <- rbind(hh_check,k_imp_indiv[which(new_house_index==kk),])
    hh_check <- data.frame(Owner=k_imp_house[kk,"Owner"],hh_check)
    if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
      Probs_syn_weighted_hybrid[1,k] <- Probs_syn_weighted_hybrid[1,k] + 1 #All same race, n_i = 2
    }
    if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]){
      Probs_syn_weighted_hybrid[2,k] <- Probs_syn_weighted_hybrid[2,k] + 1 #All same race, n_i = 3
    }
    if(nrow(hh_check)==4 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]&&hh_check[3,"Race"]==hh_check[4,"Race"]){
      Probs_syn_weighted_hybrid[3,k] <- Probs_syn_weighted_hybrid[3,k] + 1 #All same race, n_i = 4
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1){
      Probs_syn_weighted_hybrid[4,k] <- Probs_syn_weighted_hybrid[4,k] + 1 #Spouse present
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_weighted_hybrid[5,k] <- Probs_syn_weighted_hybrid[5,k] + 1 #Spouse with white HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs_syn_weighted_hybrid[6,k] <- Probs_syn_weighted_hybrid[6,k] + 1 #Spouse with black HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_weighted_hybrid[7,k] <- Probs_syn_weighted_hybrid[7,k] + 1 #White Couple
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1 && hh_check[1,"Owner"] == OWN){
      Probs_syn_weighted_hybrid[8,k] <- Probs_syn_weighted_hybrid[8,k] + 1 #White Couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
      Probs_syn_weighted_hybrid[9,k] <- Probs_syn_weighted_hybrid[9,k] + 1 #Same race couple
    }
    if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
       sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1){
      Probs_syn_weighted_hybrid[10,k] <- Probs_syn_weighted_hybrid[10,k] + 1 #White-nonwhite couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn_weighted_hybrid[11,k] <- Probs_syn_weighted_hybrid[11,k] + 1 #Non-white couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==PARENT,"Gender"]==FEMALE)==1 && length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn_weighted_hybrid[12,k] <- Probs_syn_weighted_hybrid[12,k] + 1 #Only one parent, mother
    }
    if(length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn_weighted_hybrid[13,k] <- Probs_syn_weighted_hybrid[13,k] + 1 #Only one parent   
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 | sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1){
      Probs_syn_weighted_hybrid[14,k] <- Probs_syn_weighted_hybrid[14,k] + 1 #Children present  
    }
    if(sum(hh_check[,"Relate"]==SIBLING)>=1){
      Probs_syn_weighted_hybrid[15,k] <- Probs_syn_weighted_hybrid[15,k] + 1 #Siblings present  
    }
    if(sum(hh_check[,"Relate"]==GRANDCHILD)>=1){
      Probs_syn_weighted_hybrid[16,k] <- Probs_syn_weighted_hybrid[16,k] + 1 #Grandchild present  
    }
    if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0)+ 
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 |sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0)+
       ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
      Probs_syn_weighted_hybrid[17,k] <- Probs_syn_weighted_hybrid[17,k] + 1 #Three generations present  
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1){
      Probs_syn_weighted_hybrid[18,k] <- Probs_syn_weighted_hybrid[18,k] + 1 #Same age couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_weighted_hybrid[19,k] <- Probs_syn_weighted_hybrid[19,k] + 1 #HH older than spouse, white HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]==NONHISP)==1){
      Probs_syn_weighted_hybrid[20,k] <- Probs_syn_weighted_hybrid[20,k] + 1 #Non Hisp HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
           hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_weighted_hybrid[21,k] <- Probs_syn_weighted_hybrid[21,k] + 1 #White HH with hisp origin
    }
  }
  Probs_syn_weighted_hybrid[,k] <- Probs_syn_weighted_hybrid[,k]/samp_size
  V_syn_weighted_hybrid[,k] <- (Probs_syn_weighted_hybrid[,k]*(1-Probs_syn_weighted_hybrid[,k]))/samp_size
}
dp_qbar_syn_weighted_hybrid <- rowMeans(Probs_syn_weighted_hybrid)
dp_b_syn_weighted_hybrid <- apply(Probs_syn_weighted_hybrid,1,var)
dp_ubar_syn_weighted_hybrid <- rowMeans(V_syn_weighted_hybrid)
dp_t_syn_weighted_hybrid <- dp_ubar_syn_weighted_hybrid + (dp_b_syn_weighted_hybrid*(mm+1)/mm)
#dp_t_syn_weighted_hybrid <- dp_ubar_syn_weighted_hybrid + (dp_b/mm) for synthetic data
dp_r_syn_weighted_hybrid <- dp_ubar_syn_weighted_hybrid/dp_b_syn_weighted_hybrid
dp_v_syn_weighted_hybrid <- (mm-1)*((1+((mm/(mm+1))*dp_r_syn_weighted_hybrid))^2)
#dp_v_syn_weighted_hybrid <- (mm-1)*((1+((mm)*dp_r_syn_weighted_hybrid))^2) for synthetic data
CIntLower_dp_syn_weighted_hybrid <- dp_qbar_syn_weighted_hybrid + (qt(0.025,dp_v_syn_weighted_hybrid)*sqrt(dp_t_syn_weighted_hybrid))
CIntUpper_dp_syn_weighted_hybrid <- dp_qbar_syn_weighted_hybrid - (qt(0.025,dp_v_syn_weighted_hybrid)*sqrt(dp_t_syn_weighted_hybrid))
CInt_syn_weighted_hybrid <- cbind(CIntLower_dp_syn_weighted_hybrid,CIntUpper_dp_syn_weighted_hybrid)


###### 2g: Weighted Sampler
Probs_syn_weighted <- matrix(0,nrow=21,ncol=mm)
V_syn_weighted <- matrix(0,nrow=21,ncol=mm)
for(k in 1:mm){
  k_imp_house <- dp_imput_house_weighted[((n*(k-1))+1):(n*k),]
  k_imp_indiv <- dp_imput_indiv_weighted[((N*(k-1))+1):(N*k),]
  new_n_i <- as.numeric(as.character(k_imp_house[,"HHSize"]))
  new_house_index <- rep(c(1:n),new_n_i)
  for(kk in 1:n){
    hh_check <- k_imp_house[kk,(q-p+1):q]
    colnames(hh_check) <- colnames(k_imp_indiv)
    hh_check <- rbind(hh_check,k_imp_indiv[which(new_house_index==kk),])
    hh_check <- data.frame(Owner=k_imp_house[kk,"Owner"],hh_check)
    if(nrow(hh_check)==2 && hh_check[1,"Race"]==hh_check[2,"Race"]){
      Probs_syn_weighted[1,k] <- Probs_syn_weighted[1,k] + 1 #All same race, n_i = 2
    }
    if(nrow(hh_check)==3 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]){
      Probs_syn_weighted[2,k] <- Probs_syn_weighted[2,k] + 1 #All same race, n_i = 3
    }
    if(nrow(hh_check)==4 && hh_check[1,"Race"]==hh_check[2,"Race"]&&
       hh_check[2,"Race"]==hh_check[3,"Race"]&&hh_check[3,"Race"]==hh_check[4,"Race"]){
      Probs_syn_weighted[3,k] <- Probs_syn_weighted[3,k] + 1 #All same race, n_i = 4
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1){
      Probs_syn_weighted[4,k] <- Probs_syn_weighted[4,k] + 1 #Spouse present
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_weighted[5,k] <- Probs_syn_weighted[5,k] + 1 #Spouse with white HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs_syn_weighted[6,k] <- Probs_syn_weighted[6,k] + 1 #Spouse with black HH
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs_syn_weighted[7,k] <- Probs_syn_weighted[7,k] + 1 #White Couple
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1 && hh_check[1,"Owner"] == OWN){
      Probs_syn_weighted[8,k] <- Probs_syn_weighted[8,k] + 1 #White Couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
      Probs_syn_weighted[9,k] <- Probs_syn_weighted[9,k] + 1 #Same race couple
    }
    if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
       sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1){
      Probs_syn_weighted[10,k] <- Probs_syn_weighted[10,k] + 1 #White-nonwhite couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1 && sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs_syn_weighted[11,k] <- Probs_syn_weighted[11,k] + 1 #Non-white couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==PARENT,"Gender"]==FEMALE)==1 && length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn_weighted[12,k] <- Probs_syn_weighted[12,k] + 1 #Only one parent, mother
    }
    if(length(hh_check[hh_check[,"Relate"]==PARENT,"Gender"])==1){
      Probs_syn_weighted[13,k] <- Probs_syn_weighted[13,k] + 1 #Only one parent   
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 | sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1){
      Probs_syn_weighted[14,k] <- Probs_syn_weighted[14,k] + 1 #Children present  
    }
    if(sum(hh_check[,"Relate"]==SIBLING)>=1){
      Probs_syn_weighted[15,k] <- Probs_syn_weighted[15,k] + 1 #Siblings present  
    }
    if(sum(hh_check[,"Relate"]==GRANDCHILD)>=1){
      Probs_syn_weighted[16,k] <- Probs_syn_weighted[16,k] + 1 #Grandchild present  
    }
    if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0)+ 
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 |sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0)+
       ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
      Probs_syn_weighted[17,k] <- Probs_syn_weighted[17,k] + 1 #Three generations present  
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1){
      Probs_syn_weighted[18,k] <- Probs_syn_weighted[18,k] + 1 #Same age couple
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_weighted[19,k] <- Probs_syn_weighted[19,k] + 1 #HH older than spouse, white HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]==NONHISP)==1){
      Probs_syn_weighted[20,k] <- Probs_syn_weighted[20,k] + 1 #Non Hisp HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
           hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs_syn_weighted[21,k] <- Probs_syn_weighted[21,k] + 1 #White HH with hisp origin
    }
  }
  Probs_syn_weighted[,k] <- Probs_syn_weighted[,k]/samp_size
  V_syn_weighted[,k] <- (Probs_syn_weighted[,k]*(1-Probs_syn_weighted[,k]))/samp_size
}
dp_qbar_syn_weighted <- rowMeans(Probs_syn_weighted)
dp_b_syn_weighted <- apply(Probs_syn_weighted,1,var)
dp_ubar_syn_weighted <- rowMeans(V_syn_weighted)
dp_t_syn_weighted <- dp_ubar_syn_weighted + (dp_b_syn_weighted*(mm+1)/mm)
#dp_t_syn_weighted <- dp_ubar_syn_weighted + (dp_b/mm) for synthetic data
dp_r_syn_weighted <- dp_ubar_syn_weighted/dp_b_syn_weighted
dp_v_syn_weighted <- (mm-1)*((1+((mm/(mm+1))*dp_r_syn_weighted))^2)
#dp_v_syn_weighted <- (mm-1)*((1+((mm)*dp_r_syn_weighted))^2) for synthetic data
CIntLower_dp_syn_weighted <- dp_qbar_syn_weighted + (qt(0.025,dp_v_syn_weighted)*sqrt(dp_t_syn_weighted))
CIntUpper_dp_syn_weighted <- dp_qbar_syn_weighted - (qt(0.025,dp_v_syn_weighted)*sqrt(dp_t_syn_weighted))
CInt_syn_weighted <- cbind(CIntLower_dp_syn_weighted,CIntUpper_dp_syn_weighted)


###### 3: Combine and save!!!
#CompareProbs <- cbind(Probs,Probs_cc,dp_qbar_syn,dp_qbar_syn_weighted,dp_qbar_syn_hybrid,dp_qbar_syn_weighted_hybrid,dp_qbar_syn_nz,
#                      CInt,CInt_cc,CInt_syn,CInt_syn_weighted,CInt_syn_hybrid,CInt_syn_weighted_hybrid,CInt_syn_nz)
CompareProbs <- cbind(Probs,dp_qbar_syn,CInt,CInt_syn)
#colnames(CompareProbs) <- c("Orig. Data Q","CC Data Q","Rej. Sampler Q","Weighted Sampler Q","Hybrid Sampler Q",
#                            "Weighted Hybrid Sampler Q","No Zeros Q",
#                            "Orig. Data L","Orig. Data U","CC Data L","CC Data U","Rej. Sampler L","Rej. Sampler U",
#                            "Weighted Sampler L","Weighted Sampler U",
#                            "Hybrid Sampler L","Hybrid Sampler U","Weighted Hybrid Sampler L","Weighted Hybrid Sampler U",
#                            "No Zeros L","No Zeros U")
colnames(CompareProbs) <- c("Orig. Data Q","Full Model Q","Orig. Data L","Orig. Data U",
                            "Full Model L","Full Model U")
rownames(CompareProbs) <- c("$n_i = 2$","$n_i = 3$","$n_i = 4$","SP present","SP with white HH","SP with black HH","White couple",
                            "White couple, own","Same race couple","White-nonwhite couple","Nonwhite couple, own","Only mother present",
                            "Only one parent present","Children present","Siblings present","Grandchild present","Three generations present",
                            "Same age couple","White HH, older than SP","Nonhisp HH","White, Hisp HH")
write.table(CompareProbs,"Results/CompareProbs.txt",row.names = TRUE)
CompareProbs <- read.table("Results/CompareProbs.txt",header=TRUE)
round(CompareProbs,3)

round(CompareProbs[,c(1:3)],3)
round(CompareProbs[,c(4:6)],3)
round(CompareProbs[,c(7:9)],3)

library(xtable)
xtable(round(CompareProbs[,c(1,3:7)],3),digits=3)
xtable(round(CompareProbs[,c(8,seq(10,ncol(CompareProbs),by=2))],3),digits=3)
xtable(round(CompareProbs[,c(9,seq(11,ncol(CompareProbs),by=2))],3),digits=3)

xtable(round(CompareProbs[,c(8,9,12,13,16,17,20,21)],3),digits=3)
xtable(round(CompareProbs[,c(12:19)],3),digits=3)
########################## End of Step 1 ########################## 

###########################################################################
###########################################################################
#############################***Blank Space***#############################
###########################################################################
###########################################################################
###########################################################################

library(coda)
plot(mcmc(read.table("Results/N_ZERO.txt",header = T)))
plot(mcmc(read.table("Results/ALPHA.txt",header = T)))
plot(mcmc(read.table("Results/BETA.txt",header = T)))
plot(mcmc(read.table("Results/M_CLUST.txt",header = T)))
plot(mcmc(read.table("Results/G_CLUST.txt",header = T)))

plot(mcmc(read.table("Results/N_ZERO_weighted.txt",header = T)))
plot(mcmc(read.table("Results/ALPHA_weighted.txt",header = T)))
plot(mcmc(read.table("Results/BETA_weighted.txt",header = T)))
plot(mcmc(read.table("Results/M_CLUST_weighted.txt",header = T)))
plot(mcmc(read.table("Results/G_CLUST_weighted.txt",header = T)))

summary(rbind(read.table("Results/M_CLUST_nz.txt",header = T),read.table("Results/M_CLUST_nz.txt",header = T),
          read.table("Results/M_CLUST_weighted.txt",header = T),read.table("Results/M_CLUST_weighted_hybrid.txt",header = T),
          read.table("Results/M_CLUST.txt",header = T)))
        
summary(rbind(read.table("Results/G_CLUST_nz.txt",header = T),read.table("Results/G_CLUST_nz.txt",header = T),
              read.table("Results/G_CLUST_weighted.txt",header = T),read.table("Results/G_CLUST_weighted_hybrid.txt",header = T),
              read.table("Results/G_CLUST.txt",header = T)))



abs(CompareProbs$Rej..Sampler.L - CompareProbs$Hybrid.Sampler.L) /CompareProbs$Rej..Sampler.L











