
###########################################################################
###########################################################################
################## Edit-Imputation Using the NDPMPM Model: ################
####################### Model Assessment Functions  #######################
###########################################################################
###########################################################################


## Function to process true population data
GetData <- function(min_HH_size,max_HH_size){
  House <- read.csv("Data/House.csv",header=T)
  Indiv <- read.csv("Data/Indiv.csv",header=T)
  
  ###### 1: Remove Households with size < min_HH_size and > max_HH_size
  House <- House[which(House$NP >= min_HH_size & House$NP <= max_HH_size),]
  
  
  ###### 2: Keep only Households with TEN == 1,2,or 3 and recode 1,2 as 1 and 3 as 2
  House <- House[which(House$TEN == 1 | House$TEN == 2 | House$TEN == 3),]
  House$TEN[which(House$TEN == 2)] <- 1
  House$TEN[which(House$TEN == 3)] <- 2
  
  
  ###### 3: Pick the same households in the indiv data
  pick_index <- is.element(Indiv$SERIALNO,House$SERIALNO)
  Indiv <- Indiv[pick_index,]
  
  
  ###### 4: Recode within-household variables
  ###### 4a: First, the relationship variable
  Indiv$RELP[which(Indiv$RELP == 11 | Indiv$RELP == 12 | Indiv$RELP == 13)] <- 12 #Boarder, roommate or partner
  Indiv$RELP[which(Indiv$RELP == 14 | Indiv$RELP == 15)] <- 13 #Other non-relative or foster child
  Indiv$RELP[which(Indiv$RELP == 10)] <- 11 #Other relative
  Indiv$RELP[which(Indiv$RELP == 9)] <- 10 #Child-in-law
  Indiv$RELP[which(Indiv$RELP == 8)] <- 9 #Parent-in-law
  Indiv$RELP[which(Indiv$RELP == 7)] <- 8 #Grandchild
  Indiv$RELP[which(Indiv$RELP == 6)] <- 7 #Parent
  Indiv$RELP[which(Indiv$RELP == 5)] <- 6 #Sibling
  Indiv$RELP[which(Indiv$RELP == 4)] <- 5 #Stepchild
  Indiv$RELP[which(Indiv$RELP == 3)] <- 4 #Adopted child
  Indiv$RELP[which(Indiv$RELP == 2)] <- 3 #Biological child
  Indiv$RELP[which(Indiv$RELP == 1)] <- 2 #Spouse
  Indiv$RELP[which(Indiv$RELP == 0)] <- 1 #Household head
  ###### 4b: Next, the race variable
  Indiv$RAC3P[which(Indiv$RAC3P == 4 | Indiv$RAC3P == 8| Indiv$RAC3P == 9 | Indiv$RAC3P == 10)] <- 6
  Indiv$RAC3P[which(Indiv$RAC3P == 5)] <- 4
  Indiv$RAC3P[which(Indiv$RAC3P == 7)] <- 5
  Indiv$RAC3P[which(Indiv$RAC3P >= 11 & Indiv$RAC3P <= 15)] <- 7
  Indiv$RAC3P[which(Indiv$RAC3P >= 16 & Indiv$RAC3P <= 59)] <- 8
  Indiv$RAC3P[which(Indiv$RAC3P >= 60 & Indiv$RAC3P <= 100)] <- 9
  ###### 4c: Next, the hisp variable
  Indiv$HISP[which(Indiv$HISP >= 5 & Indiv$HISP <= 24)] <- 5
  ###### 4d: Lastly, age
  Indiv$AGEP <- Indiv$AGEP + 1L
  
  
  ###### 5: Make household head into household level data
  HHhead_data <- Indiv[which(Indiv$SPORDER==1),]
  Indiv_minHH <- Indiv[-which(Indiv$SPORDER==1),]
  
  
  ###### 6: Combine household head with household level data, rename all columns
  Indiv_variables <- c("SEX","RAC3P","HISP","AGEP","RELP")
  X_house <- cbind((House$NP-1L),House$TEN,HHhead_data[,Indiv_variables])
  X_indiv <- Indiv_minHH[,Indiv_variables]
  colnames(X_house) <- c("HHSize","Owner","HHGender","HHRace","HHHisp","HHAge","HHRelate")
  colnames(X_indiv) <- c("Gender","Race","Hisp","Age","Relate")
  
  return(list(House=X_house,Indiv=X_indiv))
}


## Function to calculate all probabilities
GetAllProbs <- function(Data_house,Data_indiv,level_house,level_indiv){
  n_i <- as.numeric(as.character(Data_house$HHSize))
  FullData <- apply(Data_house,2,function(x) rep(x,n_i))
  row.names(FullData) <- NULL
  FullData <- cbind(FullData,Data_indiv)
  AllLevels <- c(level_house,level_indiv)
  n_col <- ncol(FullData)
  FullData <- data.frame(FullData)
  
  MarginalProb <- NULL; MarginalVar <- NULL; 
  for(j in 1:n_col){
    FullData[,j] <- factor(FullData[,j],levels=AllLevels[[j]])
    MarginalProb_j <- matrix(table(FullData[,j])/sum(table(FullData[,j])),ncol=1)
    MarginalProb <- rbind(MarginalProb,MarginalProb_j)
    MarginalVar_j <- (MarginalProb_j*(1-MarginalProb_j))/sum(table(FullData[,j]))
    MarginalVar <- rbind(MarginalVar,MarginalVar_j)
  }
  
  BivariateProb <- NULL; BivariateVar <- NULL; 
  AllPossBivCombs = combn(n_col,2)
  for(jj in 1:ncol(AllPossBivCombs)){
    bivcomb_j <- AllPossBivCombs[,jj]
    BivariateProb_j <- matrix(table(FullData[,bivcomb_j])/sum(table(FullData[,bivcomb_j])),ncol=1)
    BivariateProb <- rbind(BivariateProb,BivariateProb_j)
    BivariateVar_j <- (BivariateProb_j*(1-BivariateProb_j))/sum(table(FullData[,bivcomb_j]))
    BivariateVar <- rbind(BivariateVar,BivariateVar_j)
  }
  
  TrivariateProb <- NULL; TrivariateVar <- NULL; 
  AllPossTrivCombs = combn(n_col,3)
  for(jj in 1:ncol(AllPossTrivCombs)){
    trivcomb_j <- AllPossTrivCombs[,jj]
    TrivariateProb_j <- matrix(table(FullData[,trivcomb_j])/sum(table(FullData[,trivcomb_j])),ncol=1)
    TrivariateProb <- rbind(TrivariateProb,TrivariateProb_j)
    TrivariateVar_j <- (TrivariateProb_j*(1-TrivariateProb_j))/sum(table(FullData[,trivcomb_j]))
    TrivariateVar <- rbind(TrivariateVar,TrivariateVar_j)
  }
  
  return(list(MarginalProb=MarginalProb,MarginalVar=MarginalVar,
              BivariateProb=BivariateProb,BivariateVar=BivariateVar,
              TrivariateProb=TrivariateProb,TrivariateVar=TrivariateVar))
}


## Function to calculate other quantities/estimands of interest
GetOtherProbs <- function(Data_house,Data_indiv,n,p,q,house_index,H,n_estimates){
  
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
  
  Probs <- matrix(0,nrow=n_estimates)
  for(kk in 1:n){
    hh_check <- Data_house[kk,(q-p+1):q]
    colnames(hh_check) <- colnames(Data_indiv)
    hh_check <- rbind(hh_check,Data_indiv[which(house_index==kk),])
    hh_check <- data.frame(Owner=Data_house[kk,"Owner"],hh_check)
    
    for(hh in H){
      if(nrow(hh_check)== (hh+1) && length(unique(hh_check$Race))==1){
        Probs[which(H==hh)] <- Probs[which(H==hh)] + 1 #Same race household
      }
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1){
      Probs[1 + length(H)] <- Probs[1 + length(H)] + 1 #Spouse present
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs[2 + length(H)] <- Probs[2 + length(H)] + 1 #Black HH, own
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 && 
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs[3 + length(H)] <- Probs[3 + length(H)] + 1 #Spouse present, HH is white
    }
    if(sum(hh_check[,"Relate"]==SPOUSE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs[4 + length(H)] <- Probs[4 + length(H)] + 1 #Spouse present, HH is black
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs[5 + length(H)] <- Probs[5 + length(H)] + 1 #White Couple
    }
    if(sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]!=WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]!=WHITE)==1 &&
       hh_check[1,"Owner"] == OWN){
      Probs[6 + length(H)] <- Probs[6 + length(H)] + 1 #Non-White Couple, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==hh_check[hh_check[,"Relate"]==SPOUSE,"Race"])==1){
      Probs[7 + length(H)] <- Probs[7 + length(H)] + 1 #Same race couple
    }
    if((sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1 &&
        sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]!=WHITE)==1) |
       (sum(hh_check[hh_check[,"Relate"]==SPOUSE,"Race"]==WHITE)==1 &&
        sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]!=WHITE)==1)){
      Probs[8 + length(H)] <- Probs[8 + length(H)] + 1 #White-nonwhite couple
    }
    if(sum(hh_check[,"Relate"]==PARENT)==1){
      Probs[9 + length(H)] <- Probs[9 + length(H)] + 1 #Only one parent   
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1){
      Probs[10 + length(H)] <- Probs[10 + length(H)] + 1 #At least one biological child present
    }
    if(sum(hh_check[,"Relate"]==GRANDCHILD)==1){
      Probs[11 + length(H)] <- Probs[11 + length(H)] + 1 #One grandchild present
    }
    if(ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
               sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
       ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=2){
      Probs[12 + length(H)] <- Probs[12 + length(H)] + 1 #At least three generations present
    }
    if(sum(abs(hh_check[hh_check[,"Relate"]==HEAD,"Age"]-hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])<5)==1){
      Probs[13 + length(H)] <- Probs[13 + length(H)] + 1 #Couples with age difference less than 5
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>hh_check[hh_check[,"Relate"]==SPOUSE,"Age"])==1 &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs[14 + length(H)] <- Probs[14 + length(H)] + 1 #HH older than spouse, white HH
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Hisp"]!=NONHISP &&
           hh_check[hh_check[,"Relate"]==HEAD,"Race"]==WHITE)==1){
      Probs[15 + length(H)] <- Probs[15 + length(H)] + 1 #White HH with hisp origin
    }
    if(sum(hh_check[hh_check[,"Gender"]==FEMALE,"Age"]>=21)==1 &&
       ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<5)>=1 |
               sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<5)>=1 |
               sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<5)>=1),1,0)==1){
      Probs[16 + length(H)] <- Probs[16 + length(H)] + 1
      #Only one adult female (>=21) in house with at least one child less than 5
    }
    if(sum(hh_check[hh_check[,"Gender"]==MALE & hh_check[,"Hisp"]!=NONHISP,"Age"]>=21)==1 &&
       ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<10)>=1 |
               sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<10)>=1 |
               sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<10)>=1),1,0)==1){
      Probs[17 + length(H)] <- Probs[17 + length(H)] + 1
      #Only one adult male hispanic (>=21) in house with at least one child less than 10
    }
    if(sum(hh_check[hh_check[,"Gender"]==FEMALE & hh_check[,"Race"]==BLACK,"Age"]>=21)==1 &&
       ifelse((sum(hh_check[hh_check[,"Relate"]==BIOLOGICALCHILD,"Age"]<18)>=1 |
               sum(hh_check[hh_check[,"Relate"]==ADOPTEDCHILD,"Age"]<18)>=1 |
               sum(hh_check[hh_check[,"Relate"]==STEPCHILD,"Age"]<18)>=1),1,0)==1){
      Probs[18 + length(H)] <- Probs[18 + length(H)] + 1
      #Only one adult female black (>=21) in house with at least one child less than 18
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Age"]>35)==1 &&
       ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)==0 &
               sum(hh_check[,"Relate"]==ADOPTEDCHILD)==0 & sum(hh_check[,"Relate"]==STEPCHILD)==0),1,0)==1){
      Probs[19 + length(H)] <- Probs[19 + length(H)] + 1 #HH over 35, no children present
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD,"Gender"]==MALE)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs[20 + length(H)] <- Probs[20 + length(H)] + 1 #Male HH, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK,"Age"]<40)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs[21 + length(H)] <- Probs[21 + length(H)] + 1 #Black HH younger than 40, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE,"Age"]<25)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs[22 + length(H)] <- Probs[22 + length(H)] + 1 #White HH younger than 25, own
    }
    if(sum(hh_check[hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP,"Age"]>50)==1 &&
       hh_check[1,"Owner"]==OWN){
      Probs[23 + length(H)] <- Probs[23 + length(H)] + 1 #Hispanic HH older than 50, own
    }
    if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
        ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
                sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
        ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==1) &&
       sum(hh_check[hh_check[,"Relate"]==HEAD,"Race"]==BLACK)==1){
      Probs[24 + length(H)] <- Probs[24 + length(H)] + 1 #2 generations present, black HH
    }
    if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
        ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
                sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
        ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)==2) &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs[25 + length(H)] <- Probs[25 + length(H)] + 1 #3 generations present, white couple
    }
    if((ifelse(sum(hh_check[,"Relate"]==PARENT)>=1,1,0) + 
        ifelse((sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 |
                sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 | sum(hh_check[,"Relate"]==STEPCHILD)>=1),1,0) +
        ifelse(sum(hh_check[,"Relate"]==GRANDCHILD)>=1,1,0)>=1) &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
      Probs[26 + length(H)] <- Probs[26 + length(H)] + 1 #At least 2 generations present, hispanic couple
    }
    if(sum(hh_check[,"Relate"]==STEPCHILD)>=1){
      Probs[27 + length(H)] <- Probs[27 + length(H)] + 1 #At least one stepchild
    }
    if(sum(hh_check[,"Relate"]==ADOPTEDCHILD)>=1 &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==WHITE)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==WHITE)==1){
      Probs[28 + length(H)] <- Probs[28 + length(H)] + 1 #At least one adopted child, white couple
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=1 &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Hisp"]!=NONHISP)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Hisp"]!=NONHISP)==1){
      Probs[29 + length(H)] <- Probs[29 + length(H)] + 1 #At least one biological child, hispanic couple
    }
    if(sum(hh_check[,"Relate"]==BIOLOGICALCHILD)>=2 &&
       sum(hh_check[,"Relate"]==SPOUSE & hh_check[,"Race"]==BLACK)==1 &&
       sum(hh_check[,"Relate"]==HEAD & hh_check[,"Race"]==BLACK)==1){
      Probs[30 + length(H)] <- Probs[30 + length(H)] + 1 #At least two biological children, Black couple
    }
  }
  
  return(Probs)
}


## Function to calculate CI
CalculateCI <- function(Probs,V,imp_ind){
  
  if(!imp_ind){
    CIntLower <- Probs + (qnorm(0.025)*sqrt(V))
    CIntUpper <- Probs - (qnorm(0.025)*sqrt(V))
    CInt <- cbind(CIntLower,CIntUpper)
    
    return(CINT=CInt)
  } else {
    mm <- ncol(Probs)
    qbar <- rowMeans(Probs)
    b <- apply(Probs,1,var)
    ubar <- rowMeans(V)
    t <- ubar + (b*(mm+1)/mm) #t <- ubar + (b/mm) for synthetic data
    r <- ubar/b
    v <- (mm-1)*((1+((mm/(mm+1))*r))^2) #v <- (mm-1)*((1+((mm)*r))^2) for synthetic data
    CIntLower <- qbar + (qt(0.025,v)*sqrt(t))
    CIntUpper <- qbar - (qt(0.025,v)*sqrt(t))
    CInt <- cbind(CIntLower,CIntUpper)
    r <- (1+(1/mm))/r
    lambda <- (r+(2/(v + 3)))/(r + 1)
    RE <- 1/(1 + (lambda/mm))
    
    return(list(PE=qbar,VAR=t,CINT=CInt,Lambda=lambda,RE=RE))
  }
}


## Function to calculate the other probabilities for multiply-imputed datasets
GetOtherProbsMI <- function(House,Indiv,GlobalPara,samp_size){
  Results <- list()
  MarginalProb <- MarginalVar <- BivariateProb <- BivariateVar <- TrivariateProb <- TrivariateVar <- NULL;
  OtherProb <- matrix(0,nrow=length(GlobalPara$Estimands),ncol=GlobalPara$mm)
  OtherVar <- matrix(0,nrow=length(GlobalPara$Estimands),ncol=GlobalPara$mm)
  
  for(k in 1:GlobalPara$mm){
    Data_house <- House[((GlobalPara$n*(k-1))+1):(GlobalPara$n*k),]
    Data_indiv <- Indiv[((GlobalPara$N*(k-1))+1):(GlobalPara$N*k),]
    new_n_i <- as.numeric(as.character(Data_house$HHSize))
    new_house_index <- rep(c(1:GlobalPara$n),new_n_i)
    
    Results_k <- GetAllProbs(Data_house,Data_indiv,GlobalPara$level_house,GlobalPara$level_indiv)
    MarginalProb <- cbind(MarginalProb,Results_k$MarginalProb)
    MarginalVar <- cbind(MarginalVar,Results_k$MarginalVar)
    BivariateProb <- cbind(BivariateProb,Results_k$BivariateProb)
    BivariateVar <- cbind(BivariateVar,Results_k$BivariateVar)
    TrivariateProb <- cbind(TrivariateProb,Results_k$TrivariateProb)
    TrivariateVar <- cbind(TrivariateVar,Results_k$TrivariateVar)
    
    OtherProb_k <- GetOtherProbs(Data_house,Data_indiv,GlobalPara$n,GlobalPara$p,GlobalPara$q,
                                 new_house_index,GlobalPara$H,length(GlobalPara$Estimands))
    OtherProb[,k] <- OtherProb_k/samp_size
    OtherVar[,k] <- (OtherProb[,k]*(1-OtherProb[,k]))/samp_size
  }
  
  Marginal <- CalculateCI(MarginalProb,MarginalVar,imp_ind=T)
  Results$MarginalProb <- Marginal$PE; Results$MarginalVar <- Marginal$VAR
  Bivariate <- CalculateCI(BivariateProb,BivariateVar,imp_ind=T)
  Results$BivariateProb <- Bivariate$PE; Results$BivariateVar <- Bivariate$VAR
  Trivariate <- CalculateCI(TrivariateProb,TrivariateVar,imp_ind=T)
  Results$TrivariateProb <- Trivariate$PE; Results$TrivariateVar <- Trivariate$VAR
  
  Others <- CalculateCI(OtherProb,OtherVar,imp_ind=T)
  Results$OtherProb <- Others$PE; Results$OtherVar <- Others$VAR; Results$OtherCINT <- Others$CINT
  Results$OtherLambda <- Others$Lambda; Results$OtherRE <- Others$RE
  
  
  return(Results)
}





